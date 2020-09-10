/*  bam_consensus.c -- consensus subcommand.

    Copyright (C) 1998-2001,2003 Medical Research Council (Gap4/5 source)
    Copyright (C) 2003-2005,2007-2020 Genome Research Ltd.

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
#include <math.h>
#include <limits.h>
#include <float.h>

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
    int gap5;
} consensus_opts;

static int readaln(void *data, bam1_t *b) {
    consensus_opts *dat = (consensus_opts *)data;
    if (dat->iter)
        return sam_itr_next(dat->fp, dat->iter, b);
    else
        return sam_read1(dat->fp, dat->h, b);
}


/* --------------------------------------------------------------------------
 * A bayesian consensus algorithm that analyses the data to work out
 * which hypothesis of pure A/C/G/T/absent and all combinations of two
 * such bases meets the observations.
 *
 * This has its origins in Gap4 (homozygous) -> Gap5 (heterozygous)
 * -> Crumble (tidied up to use htslib's pileup) -> here.
 * 
 */

#define CONS_DISCREP    4
#define CONS_ALL        15

#define CONS_MQUAL      16

typedef struct {
    /* the most likely base call - we never call N here */
    /* A=0, C=1, G=2, T=3, *=4 */
    int call;

    /* The most likely heterozygous base call */
    /* Use "ACGT*"[het / 5] vs "ACGT*"[het % 5] for the combination */
    int het_call;

    /* Log-odds for het_call */
    int het_phred;

    /* Single phred style call */
    unsigned char phred;

    /* Sequence depth */
    int depth;

    /* Discrepancy search score */
    float discrep;
} consensus_t;

#define P_HET 1e-6

#define LOG10            2.30258509299404568401
#define TENOVERLOG10     4.34294481903251827652
#define TENLOG2OVERLOG10 3.0103


/* Sequencing technologies for seq_t.seq_tech; 5 bits, so max=31 */
#define STECH_UNKNOWN    0
#define STECH_SANGER     1
#define STECH_SOLEXA     2
#define STECH_SOLID      3
#define STECH_454        4
#define STECH_HELICOS    5
#define STECH_IONTORRENT 6
#define STECH_PACBIO     7
#define STECH_ONT        8
#define STECH_LAST       8 // highest value

double tech_undercall[] = {
    1.00, // unknown
    1.00, // sanger
    1.00, // solexa/illumina
    1.00, // solid
    1.00, // 454
    1.00, // helicos
    1.00, // iontorrent
    1.00, // pacbio
    1.63, // ont
};

#ifdef __GNUC__
#define ALIGNED(x) __attribute((aligned(x)))
#else
#define ALIGNED(x)
#endif

static double prior[25]    ALIGNED(16);  /* Sum to 1.0 */
static double lprior15[15] ALIGNED(16);  /* 15 combinations of {ACGT*} */

/* Precomputed matrices for the consensus algorithm */
static double pMM[9][101] ALIGNED(16);
static double p__[9][101] ALIGNED(16);
static double p_M[9][101] ALIGNED(16);
static double po_[9][101] ALIGNED(16);
static double poM[9][101] ALIGNED(16);
static double poo[9][101] ALIGNED(16);
static double puu[9][101] ALIGNED(16);
static double pum[9][101] ALIGNED(16);
static double pmm[9][101] ALIGNED(16);

static double e_tab_a[1002]  ALIGNED(16);
static double *e_tab = &e_tab_a[500];
static double e_tab2_a[1002] ALIGNED(16);
static double *e_tab2 = &e_tab2_a[500];
static double e_log[501]     ALIGNED(16);

/*
 * Lots of confusing matrix terms here, so some definitions will help.
 *
 * M = match base
 * m = match pad
 * _ = mismatch
 * o = overcall
 * u = undercall
 *
 * We need to distinguish between homozygous columns and heterozygous columns,
 * done using a flat prior.  This is implemented by treating every observation
 * as coming from one of two alleles, giving us a 2D matrix of possibilities
 * (the hypotheses) for each and every call (the observation).
 *
 * So pMM[] is the chance that given a call 'x' that it came from the
 * x/x allele combination.  Similarly p_o[] is the chance that call
 * 'x' came from a mismatch (non-x) / overcall (consensus=*) combination.
 *
 * Examples with observation (call) C and * follows
 *
 *  C | A  C  G  T  *          * | A  C  G  T  * 
 *  -----------------          ----------------- 
 *  A | __ _M __ __ o_         A | uu uu uu uu um
 *  C | _M MM _M _M oM         C | uu uu uu uu um
 *  G | __ _M __ __ o_         G | uu uu uu uu um
 *  T | __ _M __ __ o_         T | uu uu uu uu um
 *  * | o_ oM o_ o_ oo         * | um um um um mm
 *
 * In calculation terms, the _M is half __ and half MM, similarly o_ and um.
 *
 * Relative weights of substitution vs overcall vs undercall are governed on a
 * per base basis using the P_OVER and P_UNDER scores (subst is 1-P_OVER-P_UNDER).
 *
 * The heterozygosity weight though is a per column calculation as we're
 * trying to model whether the column is pure or mixed. Hence this is done
 * once via a prior and has no affect on the individual matrix cells.
 */

static void consensus_init(double p_het) {
    int i, t;

    for (i = -500; i <= 500; i++)
        e_tab[i] = exp(i);
    for (i = -500; i <= 500; i++)
        e_tab2[i] = exp(i/10.);
    for (i = 0; i <= 500; i++)
        e_log[i] = log(i);

    // Heterozygous locations
    for (i = 0; i < 25; i++)
        prior[i] = p_het / 20;
    prior[0] = prior[6] = prior[12] = prior[18] = prior[24] = (1-p_het)/5;

    lprior15[0]  = log(prior[0]);
    lprior15[1]  = log(prior[1]*2);
    lprior15[2]  = log(prior[2]*2);
    lprior15[3]  = log(prior[3]*2);
    lprior15[4]  = log(prior[4]*2);
    lprior15[5]  = log(prior[6]);
    lprior15[6]  = log(prior[7]*2);
    lprior15[7]  = log(prior[8]*2);
    lprior15[8]  = log(prior[9]*2);
    lprior15[9]  = log(prior[12]);
    lprior15[10] = log(prior[13]*2);
    lprior15[11] = log(prior[14]*2);
    lprior15[12] = log(prior[18]);
    lprior15[13] = log(prior[19]*2);
    lprior15[14] = log(prior[24]);


    // Rewrite as new form
    for (t = STECH_UNKNOWN; t <= STECH_LAST; t++) {
        for (i = 1; i < 101; i++) {
            double prob = 1 - pow(10, -i / 10.0);

            // May want to multiply all these by 5 so pMM[i] becomes close
            // to -0 for most data. This makes the sums increment very slowly,
            // keeping bit precision in the accumulator.
            pMM[t][i] = log(prob/5);
            p__[t][i] = log((1-prob)/20);
            p_M[t][i] = log((exp(pMM[t][i]) + exp(p__[t][i]))/2);

            puu[t][i] = p__[t][i];

            poM[t][i] = p_M[t][i] *= tech_undercall[t];
            po_[t][i] = p__[t][i] *= tech_undercall[t];
            poo[t][i] = p__[t][i] *= tech_undercall[t];
            pum[t][i] = p_M[t][i] *= tech_undercall[t];
            pmm[t][i] = pMM[t][i] *= tech_undercall[t];
        }

        pMM[t][0] = pMM[t][1];
        p__[t][0] = p__[t][1];
        p_M[t][0] = p_M[t][1];

        pmm[t][0] = pmm[t][1];
        poo[t][0] = poo[t][1];
        po_[t][0] = po_[t][1];
        poM[t][0] = poM[t][1];
        puu[t][0] = puu[t][1];
        pum[t][0] = pum[t][1];
    }
}

static inline double fast_exp(double y) {
    if (y >= -50 && y <= 50)
        return e_tab2[(int)(y*10)];

    if (y < -500)
        y = -500;
    if (y > 500)
        y = 500;

    return e_tab[(int)y];
}

/*Taylor (deg 3) implementation of the log: http://www.flipcode.com/cgi-bin/fcarticles.cgi?show=63828*/
static inline double fast_log2(double val)
{
   register int64_t *const     exp_ptr = ((int64_t*)&val);
   register int64_t            x = *exp_ptr;
   register const int      log_2 = ((x >> 52) & 2047) - 1024;
   x &= ~(2047LL << 52);
   x += 1023LL << 52;
   *exp_ptr = x;

   val = ((-1.0f/3) * val + 2) * val - 2.0f/3;

   return val + log_2;
}

static inline double fast_log (double val) {
    return fast_log2(val)*0.69314718;
}

//#define fast_exp exp
//#define fast_log log

#define ph_log(x) (-TENLOG2OVERLOG10*fast_log2((x)))


/*
 * As per calculate_consensus_bit_het but for a single pileup column.
 */
int calculate_consensus_pileup(int flags,
                               const bam_pileup1_t *p,
                               int np,
                               consensus_t *cons) {
    int i, j;
    static int init_done =0;
    static double q2p[101], mqual_pow[256];
    double min_e_exp = DBL_MIN_EXP * log(2) + 1;

    double S[15] ALIGNED(16) = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double sumsC[6] = {0,0,0,0,0,0}, sumsE = 0;
    int depth = 0;

    /* Map the 15 possible combinations to 1-base or 2-base encodings */
    static int map_sing[15] ALIGNED(16) =
        {0, 5, 5, 5, 5,
            1, 5, 5, 5,
               2, 5, 5,
                  3, 5,
                     4};
    static int map_het[15] ALIGNED(16) =
        {0,  1,  2,  3,  4,
             6,  7,  8,  9,
                12, 13, 14,
                    18, 19,
                        24};

    if (!init_done) {
        init_done = 1;
        consensus_init(P_HET);

        for (i = 0; i <= 100; i++) {
            q2p[i] = pow(10, -i/10.0);
        }

        for (i = 0; i < 255; i++) {
            //mqual_pow[i] = 1-pow(10, -(i+.01)/10.0);
            //mqual_pow[i] = 1-pow(10, -(i/3+.1)/10.0);
            mqual_pow[i] = 1-pow(10, -(i/2+.05)/10.0);
        }
        // unknown mqual
        mqual_pow[255] = mqual_pow[10];
    }

    /* Initialise */
    int counts[6] = {0};

    /* Accumulate */
    int n;
    //printf("-----\n");

    // FIXME: also seed with unknown alleles so low coverage data is
    // less confident.

    for (n = 0; n < np; n++) {
        if (p[n].is_refskip)
            continue;

        bam1_t *b = p[n].b;
        uint8_t base = bam_seqi(bam_get_seq(b), p[n].qpos);
        uint8_t qual = bam_get_qual(b)[p[n].qpos];
        const int stech = STECH_SOLEXA;

        // =ACM GRSV TWYH KDBN
        static int L[16] = {
            5,0,1,5, 2,5,5,5, 3,5,5,5, 5,5,5,5
        };

        // convert from sam base to acgt*n order.
        base = L[base];
        if (p[n].is_del) base = 4;

        double MM, __, _M, qe;

        // Correction for mapping quality.  Maybe speed up via lookups?
        // Cannot nullify mapping quality completely.  Lots of (true)
        // SNPs means low mapping quality.  (Ideally need to know
        // hamming distance to next best location.)

        if (flags & CONS_MQUAL) {
            double _p = mqual_pow[qual];
            double _m = mqual_pow[b->core.qual];

            //printf("%c %d -> %d, %f %f\n", "ACGT*N"[base], qual, (int)(-TENOVERLOG10 * log(1-(_m * _p + (1 - _m)/4))), _p, _m);
            qual = ph_log(1-(_m * _p + (1 - _m)/4));
        }

        /* Quality 0 should never be permitted as it breaks the math */
        if (qual < 1)
            qual = 1;

        __ = p__[stech][qual];
        MM = pMM[stech][qual] - __;
        _M = p_M[stech][qual] - __;

        if (flags & CONS_DISCREP) {
            qe = q2p[qual];
            sumsE += qe;
            sumsC[base] += 1 - qe;
        }

        counts[base]++;

        switch (base) {
        case 0:
            S[0] += MM; S[1 ]+= _M; S[2 ]+= _M; S[3 ]+= _M; S[4 ]+= _M;
            break;

        case 1:
            S[1 ]+= _M; S[5 ]+= MM; S[6 ]+= _M; S[7 ]+= _M; S[8 ]+= _M;
            break;

        case 2:
            S[2 ]+= _M; S[6 ]+= _M; S[9 ]+= MM; S[10]+= _M; S[11]+= _M; 
            break;

        case 3:
            S[3 ]+= _M; S[7 ]+= _M; S[10]+= _M; S[12]+= MM; S[13]+= _M; 
            break;

        case 4:
            S[4 ]+= _M; S[8 ]+= _M; S[11]+= _M; S[13]+= _M; S[14]+= MM;
            break;

        case 5: /* N => equal weight to all A,C,G,T but not a pad */
            S[0] += MM; S[1 ]+= MM; S[2 ]+= MM; S[3 ]+= MM; S[4 ]+= _M;
                        S[5 ]+= MM; S[6 ]+= MM; S[7 ]+= MM; S[8 ]+= _M;
                                    S[9 ]+= MM; S[10]+= MM; S[11]+= _M; 
                                                S[12]+= MM; S[13]+= _M; 
            break;
        }

        depth++;
    }


    /* and speculate */
    {
        double shift, max, max_het, norm[15];
        int call = 0, het_call = 0, ph;
        double tot1, tot2;

        /*
         * Scale numbers so the maximum score is 0. This shift is essentially 
         * a multiplication in non-log scale to both numerator and denominator,
         * so it cancels out. We do this to avoid calling exp(-large_num) and
         * ending up with norm == 0 and hence a 0/0 error.
         *
         * Can also generate the base-call here too.
         */
        shift = -DBL_MAX;
        max = -DBL_MAX;
        max_het = -DBL_MAX;

        for (j = 0; j < 15; j++) {
            S[j] += lprior15[j];
            if (shift < S[j])
                shift = S[j];

            /* Only call pure AA, CC, GG, TT, ** for now */
            if (j != 0 && j != 5 && j != 9 && j != 12 && j != 14) {
                if (max_het < S[j]) {
                    max_het = S[j];
                    het_call = j;
                }
                continue;
            }

    if (max < S[j]) {
                max = S[j];
                call = j;
            }
        }

        /*
         * Shift and normalise.
         * If call is, say, b we want p = b/(a+b+c+...+n), but then we do
         * p/(1-p) later on and this has exceptions when p is very close
         * to 1.
         *
         * Hence we compute b/(a+b+c+...+n - b) and
         * rearrange (p/norm) / (1 - (p/norm)) to be p/norm2.
         */
        for (j = 0; j < 15; j++) {
            S[j] -= shift;
            double e = fast_exp(S[j]);
            S[j] = (S[j] > min_e_exp) ? e : DBL_MIN;
            norm[j] = 0;
        }

        tot1 = tot2 = 0;
        for (j = 0; j < 15; j++) {
            norm[j]    += tot1;
            norm[14-j] += tot2;
            tot1 += S[j];
            tot2 += S[14-j];
        }

        /* And store result */
        if (depth && depth != counts[5] /* all N */) {
            double m;

            cons->depth = depth;

            cons->call     = map_sing[call];
            if (norm[call] == 0) norm[call] = DBL_MIN;
            ph = ph_log(norm[call]) + .5;
            cons->phred = ph > 255 ? 255 : (ph < 0 ? 0 : ph);
            //cons->call_prob1 = norm[call]; // p = 1 - call_prob1

            cons->het_call = map_het[het_call];
            if (norm[het_call] == 0) norm[het_call] = DBL_MIN;
            ph = TENLOG2OVERLOG10 * (fast_log2(S[het_call]) - fast_log2(norm[het_call])) + .5;

            cons->het_phred = ph;
            //cons->het_prob_n = S[het_call]; // p = prob_n / prob_d
            //cons->het_prob_d = norm[het_call];

            /* Compute discrepancy score */
            if (flags & CONS_DISCREP) {
                m = sumsC[0]+sumsC[1]+sumsC[2]+sumsC[3]+sumsC[4];
                double c;
                if (cons->het_phred > 0)
                    c = sumsC[cons->het_call%5] + sumsC[cons->het_call/5];
                else
                    c = sumsC[cons->call];;
                cons->discrep = (m-c)/sqrt(m);
//              printf("Discrep = %f,  %f %f %f %f %f\n", cons->discrep,
//                     sumsC[0], sumsC[1], sumsC[2], sumsC[3], sumsC[4]);
//              if (cons->discrep > 1)
//                  printf("XYZZY\n");
            }
        } else {
            cons->call = 5; /* N */
            cons->het_call = 0;
            cons->het_phred = 0;
            cons->phred = 0;
            cons->depth = 0;
            cons->discrep = 0;
        }
    }

    return 0;
}


/* --------------------------------------------------------------------------
 * A simple summing algorithm, either pure base frequency, or by
 * weighting them according to their quality values.
 *
 * This is crude, but easy to understand and fits with several
 * standard pileup criteria (eg COG-UK / CLIMB Covid-19 seq project).
 *
 *
 * call1 / score1 / depth1 is the highest scoring allele.
 * call2 / score2 / depth2 is the second highest scoring allele.
 *
 * Het_fract:  score2/score1
 * Call_fract: score1 or score1+score2 over total score
 * Min_depth:  minimum total depth of utilised bases (depth1+depth2)
 * Min_score:  minimum total score of utilised bases (score1+score2)
 *
 * Eg het_fract 0.66, call_fract 0.75 and min_depth 10.
 * 11A, 2C, 2G (14 total depth) is A.
 * 9A, 2C, 2G  (12 total depth) is N as depth(A) < 10.
 * 11A, 5C, 5G (21 total depth) is N as 11/21 < 0.75 (call_fract)
 *
 * 
 * 6A, 5G, 1C  (12 total depth) is AG het as depth(A)+depth(G) >= 10
 *                              and 5/6 >= 0.66 and 11/12 >= 0.75.
 *
 * 6A, 5G, 4C  (15 total depth) is N as (6+5)/15 < 0.75 (call_fract).
 *
 *
 * Note for the purpose of deletions, a base/del has an ambiguity
 * code of lower-case base (otherwise it is uppercase).
 */

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


/* --------------------------------------------------------------------------
 * Main processing logic
 */

// FIXME: move to header file if we intend to keep this interaction.
extern int pileup_seq(FILE *fp, const bam_pileup1_t *p, hts_pos_t pos,
                      hts_pos_t ref_len, const char *ref, kstring_t *ks,
                      int rev_del, int no_ins, int no_ins_mods,
                      int no_del, int no_ends);

void consensus_pileup(consensus_opts *opts, const bam_pileup1_t *p,
                      int np, int tid, int pos) {
    kstring_t ks = {0,0};

    int cq, cb;
    if (opts->gap5) {
        consensus_t cons;
        calculate_consensus_pileup(CONS_ALL, p, np, &cons);
        if (cons.het_phred > 0) {
            cb = "AMRWa" // 5x5 matrix with ACGT* per row / col
                 "MCSYc" 
                 "RSGKg"
                 "WYKTt"
                 "acgt*"[cons.het_call];
            cq = cons.het_phred > 255 ? 255 : cons.het_phred;
        } else {
            cb = "ACGT*"[cons.call];
            cq = cons.phred;
        }
    } else {
        cb = consensus(p, np, opts, &cq);
    }

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

    while ((c = getopt_long(argc, argv, "@:qd:c:H:r:5", lopts, NULL)) >= 0) {
        switch (c) {
        case 'q': opts.use_qual=1; break;
        case 'd': opts.min_depth = atoi(optarg); break;
        case 'c': opts.call_fract = atof(optarg); break;
        case 'H': opts.het_fract = atof(optarg); break;
        case 'r': opts.reg = optarg; break;
        case '5': opts.gap5 = 1; break;

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

