/*  bam_consensus.c -- consensus subcommand.

    Copyright (C) 1993 Medical Research Council (Gap4/5 source)
    Copyright (C) 2016-2020 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/*
 * The origins of this algorithm come from the Staden Package, initially
 * the Gap4 consensus algorithm, revised for Gap5, and substantially
 * rewritten again for the qual-lossy Crumble compression tool which
 * is where this algorithm has forked from.
 */

#include <config.h>

#include <stdio.h>
#include <stdlib.h>

#include <htslib/sam.h>

#include "samtools.h"
#include "sam_opts.h"


typedef struct {
    samFile *fp;
    sam_hdr_t *h;
    hts_idx_t *idx;
    hts_itr_t *iter;
    char *reg;
    int use_qual;
    int min_depth;
    double call_fract;
    double het_fract;
} consensus_opts;

static int readaln(void *data, bam1_t *b) {
    consensus_opts *dat = (consensus_opts *)data;
    if (dat->iter)
	return sam_itr_next(dat->fp, dat->iter, b);
    else
	return sam_read1(dat->fp, dat->h, b);
}


// call1 / score1 / depth1 is the highest scoring allele.
// call2 / score2 / depth2 is the second highest scoring allele.
//
// Het_fract:  score2/score1
// Call_fract: score1 or score1+score2 over total score
// Min_depth:  minimum total depth of utilised bases (depth1+depth2)
// Min_score:  minimum total score of utilised bases (score1+score2)
//
// Eg het_fract 0.66, call_fract 0.75 and min_depth 10.
// 11A, 2C, 2G (14 total depth) is A.
// 9A, 2C, 2G  (12 total depth) is N as depth(A) < 10.
// 11A, 5C, 5G (21 total depth) is N as 11/21 < 0.75 (call_fract)
//
// 
// 6A, 5G, 1C  (12 total depth) is AG het as depth(A)+depth(G) >= 10
//                              and 5/6 >= 0.66 and 11/12 >= 0.75.
//
// 6A, 5G, 4C  (15 total depth) is N as (6+5)/15 < 0.75 (call_fract).
//
//
// Note for the purpose of deletions, a base/del has an ambiguity
// code of lower-case base (otherwise it is uppercase).
static int consensus(const bam_pileup1_t *plp, int nplp,
		     consensus_opts *opts, int *qual) {
    int i, min_qual = 0;

    // Ignore ambiguous bases in seq for now, so we don't treat R, Y, etc
    // as part of one base and part another.  Based on BAM seqi values.
    // We also use freq[16] as "*" for gap.
    int freq[17] = {0};  // base frequency, aka depth
    int score[17] = {0}; // summation of base qualities

    // Accumulate
    for (i = 0; i < nplp; i++) {
        const bam_pileup1_t *p = plp+i;
        int q = bam_get_qual(p->b)[p->qpos];
        if (q < min_qual)
            // Should we still record these in freq[] somewhere so
            // we can use them in the fracts?
            // Difference between >= X% of high-qual bases calling Y
            // and >= X% of all bases are high-quality Y calls.
            continue;

        int b = p->is_del ? 16 : bam_seqi(bam_get_seq(p->b), p->qpos);
        freq[b]++;
        score[b] += opts->use_qual ? q : 1;
    }

    // Total usable depth
    int tdepth = 0, tscore = 0;
    for (i = 0; i < 5; i++) {
        tdepth += freq[1<<i];
        tscore += score[1<<i];
    }

    // Best and second best potential calls
    int call1  = 15, call2 = 15;
    int depth1 = 0,  depth2 = 0;
    int score1 = 0,  score2 = 0;
    for (i = 0; i < 5; i++) {
        int c = 1<<i; // A C G T *
        if (score1 < score[c]) {
            depth2 = depth1;
            score2 = score1;
            call2  = call1;
            depth1 = freq[c];
            score1 = score[c];
            call1  = c;
        } else if (score2 < score[c]) {
            depth2 = freq[c];
            score2 = score[c];
            call2  = c;
        }
    }

    // Work out which best and second best are usable as a call
    int used_score = score1;
    int used_depth = depth1;
    int used_base  = call1;
    if (score2 >= opts->het_fract * score1) {
        used_base  |= call2;
        used_score += score1;
        used_depth += depth2;
    }

    // N is too shallow, or insufficient proportion of total
    if (used_depth < opts->min_depth ||
        used_score < opts->call_fract * tscore) {
        used_depth = 0;
        // But note shallow gaps are still called gaps, not N, as
        // we're still more confident there is no base than it is
        // A, C, G or T.
        used_base = call1 == 16 /*&& depth1 >= call_fract * depth*/
            ? 16 : 0; // * or N
    }

    // Our final call.  "?" shouldn't be possible to generate
    const char *het =
        "NACMGRSVTWYHKDB?"
        "*ac?g???t???????";

    //printf("%c %d\n", het[used_base], used_depth);
    if (qual)
        *qual = 100.0 * used_score / tscore;
    return het[used_base];
}

// FIXME: move to header file if we intend to keep this interaction.
extern int pileup_seq(FILE *fp, const bam_pileup1_t *p, hts_pos_t pos,
		      hts_pos_t ref_len, const char *ref, kstring_t *ks,
		      int rev_del, int no_ins, int no_ins_mods,
		      int no_del, int no_ends);

void consensus_pileup(consensus_opts *opts, const bam_pileup1_t *p,
		      int np, int tid, int pos) {
    kstring_t ks = {0,0};
    int cq, cb = consensus(p, np, opts, &cq);
    printf("%s\t%d\t%c\t%d\t",
	   sam_hdr_tid2name(opts->h, tid), pos, cb, cq);

    int j;
    for (j = 0; j < np; j++)
	pileup_seq(stdout, p+j, pos, 0, NULL, &ks, 0, 2, 1, 2, 1);
    putchar('\n');
}

// Iterate over the pileup
static int consensus_loop(consensus_opts *opts) {
    bam_plp_t iter;
    const bam_pileup1_t *p;
    int tid, pos, n;

    iter = bam_plp_init(readaln, (void *)opts);
    while ((p = bam_plp_auto(iter, &tid, &pos, &n)) != 0) {
	if (opts->iter && (pos <= opts->iter->beg || pos > opts->iter->end))
	    continue;
	consensus_pileup(opts, p, n, tid, pos);
    }
    bam_plp_destroy(iter);

    return 0;
}

static void usage_exit(FILE *fp, int exit_status) {
    fprintf(fp, "Usage: samtools consensus [options] <in.bam>\n");
    fprintf(fp, "\nOptions:\n");
    fprintf(fp, "   -r, --region REG    Limit query to REG. Requires an index\n");
    fprintf(fp, "   -q, --use-qual      Use quality values in calculation\n");
    fprintf(fp, "   -m, --min-depth D   Minimum depth of D [20]\n");
    fprintf(fp, "   -m, --call-fract C  At least C portion of bases must agree [0.75]\n");
    fprintf(fp, "   -m, --het-fract C   Minimum fraction of 2nd-most to most common base [0.66]\n");
    sam_global_opt_help(fp, "-.---@-.");
    exit(exit_status);
}

int main_consensus(int argc, char **argv) {
    int c;

    consensus_opts opts = {
	.use_qual   = 0,
	.min_depth  = 20,
	.call_fract = 0.75,
	.het_fract  = 0.66,
    };

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', '-', '-', '@'),
	{"use-qual",   no_argument,       NULL, 'q'},
	{"min-depth",  required_argument, NULL, 'd'},
	{"call-fract", required_argument, NULL, 'c'},
	{"het-fract",  required_argument, NULL, 'H'},
	{"region",     required_argument, NULL, 'r'},
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "@:qd:c:H:r:", lopts, NULL)) >= 0) {
	switch (c) {
	case 'q': opts.use_qual=1; break;
	case 'd': opts.min_depth = atoi(optarg); break;
	case 'c': opts.call_fract = atof(optarg); break;
	case 'H': opts.het_fract = atof(optarg); break;
	case 'r': opts.reg = optarg; break;

        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
        case '?':
            usage_exit(stderr, EXIT_FAILURE);
        }
    }

    if (argc != optind+1) {
        if (argc == optind) usage_exit(stdout, EXIT_SUCCESS);
        else usage_exit(stderr, EXIT_FAILURE);
    }
    opts.fp = sam_open_format(argv[optind], "r", &ga.in);
    if (opts.fp == NULL) {
        print_error_errno("consensus", "Cannot open input file \"%s\"",
			  argv[optind]);
        return 1;
    }
    if (ga.nthreads > 0)
        hts_set_threads(opts.fp, ga.nthreads);

    if (hts_set_opt(opts.fp, CRAM_OPT_DECODE_MD, 0)) {
        fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
        return 1;
    }

    if (!(opts.h = sam_hdr_read(opts.fp))) {
        fprintf(stderr, "Failed to read header for \"%s\"\n", argv[optind]);
        return 1;
    }

    if (opts.reg) {
	opts.idx = sam_index_load(opts.fp, argv[optind]);
	if (!opts.idx) {
	    print_error("consensus", "Cannot load index for input file \"%s\"",
			argv[optind]);
	    return 1;
	}
	opts.iter = sam_itr_querys(opts.idx, opts.h, opts.reg);
	if (!opts.iter) {
	    print_error("consensus", "Failed to parse region \"%s\"",
			opts.reg);
	    return 1;
	}
    }

    consensus_loop(&opts);

    if (opts.iter)
	hts_itr_destroy(opts.iter);
    if (opts.idx)
	hts_idx_destroy(opts.idx);

    if (sam_close(opts.fp) < 0) {
	print_error_errno("consensus", "Closing input file \"%s\"",
			  argv[optind]);
	return 1;
    }

    sam_hdr_destroy(opts.h);
    sam_global_args_free(&ga);

    return 0;
}

