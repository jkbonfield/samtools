/*  bam_sort2.c -- experimental pipelined, cache-aware sort.

    Copyright (C) 2026 Genome Research Ltd.

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
 * Pipeline
 * --------
 * Two record batches, each with half of the -m budget.  The main thread
 * reads records into one batch while a sorter thread processes the other:
 *
 *   reader:  sam_read1 -> copy into batch arena -> hand off when full
 *   sorter:  extract keys (parallel) -> sort keys (parallel)
 *            -> write one sorted run (parallel serialise + compress)
 *
 * The final batch stays in memory and joins the on-disk runs in the merge.
 *
 * Temporary runs
 * --------------
 * A run is a valid BGZF file.  Records are framed as
 *     uint32 body_len, uint32 key_len,
 *     body: key (padded to 8) | s2_rec_t + pad | data
 * Workers each take a chunk of consecutive records (~128 KB), serialise and
 * compress it into BGZF blocks with bgzf_compress(), and append the chunks
 * in order.  Every chunk starts on a block boundary, so the first key and
 * file offset of each chunk form a sparse in-memory index of the run.
 *
 * Final merge
 * -----------
 * For BAM output, the key space is cut into ranges at quantiles of the run
 * indexes.  A splitter is a (key, source) pair, so even long runs of equal
 * keys, such as unmapped reads, split across ranges without changing tie
 * order.  Each worker thread seeks every run to the start of a range,
 * merges it with its own loser tree, encodes BAM and compresses BGZF blocks
 * itself; ranges are appended to the output in order.  Other output
 * formats, or records BAM cannot hold directly (more than 65535 CIGAR ops,
 * positions beyond 2^31), use a serial merge through sam_write1().
 *
 * Ordering matches "samtools sort" as follows (ties keep input order):
 *   coordinate:  tid (unmapped last), pos, reverse-strand flag
 *   -n:          strnum_cmp(qname), then READ1/READ2, primary,
 *                supplementary, secondary (as samtools develop)
 *   -N:          strcmp(qname),     then the same
 * The name tie-break is name_flag_byte() below.
 */

#if defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE            // syscall() for memfd_create
#endif
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <inttypes.h>
#include <limits.h>
#include <pthread.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <sys/mman.h>
#if defined(__linux__)
#include <sys/syscall.h>
#endif

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "htslib/hts.h"
#include "htslib/thread_pool.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "bam.h"   // unclipped_start() etc., as bam_sort.c uses
#include "sam_opts.h"
#include "samtools.h"
#include "sort2_engine.h"

#define S2_DEFAULT_MEM_PER_THREAD ((size_t)768 << 20)
#define S2_RUN_MAGIC "S2R1"
#define S2_CHUNK_BYTES (128 * 1024)   // run index granularity

#define ALIGN8(x) (((x) + 7) & ~(size_t)7)


/* Sort orders, as samtools sort: coordinate (default), -n natural name, -N
 * byte name, -t TAG then coordinate (TAG_COORD) or name (TAG_NAME, with -n
 * or -N), -M minimiser (MINHASH), --template-coordinate (TEMPLATE). */
enum s2_order { ORDER_COORD, ORDER_NAME_NAT, ORDER_NAME_BYTE,
                ORDER_TAG_COORD, ORDER_TAG_NAME, ORDER_MINHASH,
                ORDER_TEMPLATE };

/* Numeric orders have keys recomputed from the record; the others have
 * byte-comparable keys, which are stored alongside records in run files. */
static inline int order_numeric(enum s2_order o) {
    return o == ORDER_COORD || o == ORDER_MINHASH;
}
enum s2_merge { MERGE_SERIAL, MERGE_HTS };

/* A record as held in a batch arena or a run file: the core struct and data
 * length, padding to 8, then the variable-length data. */
typedef struct {
    bam1_core_t core;
    uint32_t l_data;
    uint32_t pad;
} s2_rec_t;
#define REC_DATA_OFF ALIGN8(sizeof(s2_rec_t))

/* A sort key.  Coordinate order: (k_hi, k64) = (tid', ((pos+1) << 1) | rev),
 * where tid' maps unmapped past every real tid; positions may use the full
 * 64-bit range SAM allows.  Name order: bytes p[0..len). */
typedef struct {
    uint64_t k64;
    const uint8_t *p;
    uint32_t len;
    uint64_t k_hi;
    uint64_t k_x;              // minimiser order: hash position and strand
} s2_key_t;

/* A merge position: records order by (key, source), source being the run
 * index (the in-memory batch is last).  Used for range splitters. */
typedef struct {
    s2_key_t k;
    int src;
} s2_bound_t;

typedef struct {
    uint64_t coff;             // compressed offset of the chunk's first block
    uint64_t ord;              // ordinal of the chunk's first record
    s2_key_t k;                // first key in the chunk (name bytes owned)
} s2_idx_t;

/* A sorted run file.  Records are read until nrec have been seen. */
typedef struct {
    char *fn;
    s2_idx_t *idx;
    size_t nidx;
    uint64_t nrec;
    uint64_t ubytes;           // uncompressed framed bytes
    int unlinked;              // already removed from the filesystem
} s2_run_t;


typedef struct {
    uint8_t *arena;
    size_t used, cap;
    size_t *off;               // arena offset of each record, input order
    size_t n, m_off;
    size_t charged;            // memory charged against the batch limit

    // name keys: kptr[i]/klen[i] for record i, stored in per-thread buffers
    const uint8_t **kptr;
    uint32_t *klen;
    uint8_t **tbuf;
    size_t *tbuf_cap;

    s2_kv_t *kv, *tmp, *sorted;
    size_t m_kv;
    int wide;                  // coordinate keys too large to pack
    int eof;                   // set on the last batch
} s2_batch_t;

typedef struct {
    enum s2_order order;
    int nthreads;              // compute threads (>= 1)
    size_t mem_total, batch_limit;
    int tmp_level;
    int max_runs;
    int stats;
    char *tmp_prefix;

    hts_tpool *pool;           // htslib I/O pool, may be NULL
    int error;                 // sticky; protected by lock

    // runs on disk, in input order
    s2_run_t *runs;
    int n_runs, m_runs, run_seq;

    int phase1_done;           // sorter finished the last batch

    // final merge
    int serial_merge;          // --serial-merge
    int merge_kind;            // MERGE_SERIAL or MERGE_HTS
    const char *out_fn, *out_mode;
    void *out_opts;            // format-specific output options
    size_t drain_recs;         // see pm_worker; SORT2_TEST_DRAIN overrides

    // reader <-> sorter hand-off
    s2_batch_t batch[2];
    pthread_mutex_t lock;
    pthread_cond_t cond;
    s2_batch_t *full;          // batch waiting for the sorter
    int free_mask;             // bit i set if batch[i] is free
    s2_batch_t *final;         // sorted last batch, left by the sorter

    // statistics
    size_t n_records, n_batches, n_intermediate, out_ranges;
    uint64_t tmp_cbytes, tmp_ubytes;   // temp data written (compressed, raw)
    uint64_t tmp_uread;                // raw temp bytes read back by merges
    uint64_t scratch_bytes;            // output bytes via scratch files
    char *idx_fn;                      // --write-index: index file name
    char sort_tag[3];                  // -t TAG
    int tag_natural;                   // -t with -n: natural name order
    int kmer;                          // -M -K: minimiser k-mer size
    int try_rev, no_squash, window;    // -R, -H, -w
    void *kmer_h;                      // -I: reference minimiser index
    int placed;                        // some read has tid >= 0 (for @HD)
    void *libs;                        // --template-coordinate: RG -> LB
    hts_idx_t *idx;                    // built by the parallel BAM merge
    double t_close, t_idx_save;
    int out_threads;
    const char *serial_why;    // why the final merge was not parallel
    double t_read_stall, t_sort_idle, t_keys, t_sort, t_runs, t_inter,
           t_final;
} s2_state_t;

static long fd_budget(long want);

/* CPU and I/O counters at a point in time, for per-phase --stats.  rchar
 * and wchar (Linux /proc/self/io) count bytes passed through read and
 * write system calls, page cache hits included; -1 if unavailable. */
typedef struct {
    double user, sys;
    int64_t rchar, wchar;
} s2_usage_t;

static void usage_now(s2_usage_t *u) {
    struct rusage ru;
    memset(u, 0, sizeof(*u));
    u->rchar = u->wchar = -1;
    if (getrusage(RUSAGE_SELF, &ru) == 0) {
        u->user = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6;
        u->sys  = ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
    }
#if defined(__linux__)
    FILE *f = fopen("/proc/self/io", "r");
    if (f) {
        char key[64];
        long long v;
        while (fscanf(f, "%63[^:]: %lld\n", key, &v) == 2) {
            if (strcmp(key, "rchar") == 0)
                u->rchar = v;
            else if (strcmp(key, "wchar") == 0)
                u->wchar = v;
        }
        fclose(f);
    }
#endif
}

static double now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static void set_error(s2_state_t *st) {
    pthread_mutex_lock(&st->lock);
    st->error = 1;
    pthread_cond_broadcast(&st->cond);
    pthread_mutex_unlock(&st->lock);
}

/* ---------------------------------------------------------------------- */
/* Keys */

// Must match the tie-break in bam_sort.c's name comparison.
/* As samtools' bam1_cmp_core(): ((f&0xc0)<<8)|((f&0x100)<<3)|((f&0x800)>>3),
 * whose bits all lie in its high byte. */
static inline uint8_t name_flag_byte(const bam1_core_t *c) {
    return (c->flag & 0xc0) | ((c->flag & 0x100) >> 5) | ((c->flag & 0x800) >> 11);
}

static inline uint32_t name_len(const s2_rec_t *r) {
    return r->core.l_qname - 1 - r->core.l_extranul;
}

/* (tid' << 33) | ((pos+1) << 1) | rev.  tid' maps unmapped (-1) past every
 * real tid.  Returns -1 if pos+1 does not fit in 32 bits. */
/* Full coordinate key.  pos <= HTS_POS_MAX < 2^63 - 2^31, so
 * ((pos+1) << 1) | rev fits in 64 bits. */
static inline void coord_key2(const bam1_core_t *c, uint64_t *hi,
                              uint64_t *lo) {
    *hi = c->tid < 0 ? 0x7fffffffu : (uint32_t) c->tid;
    uint64_t p1 = c->pos < 0 ? 0 : (uint64_t) c->pos + 1;
    *lo = (p1 << 1) | !!(c->flag & BAM_FREVERSE);
}

/* Numeric key (hi, lo, x) for COORD and MINHASH.  An unplaced read in
 * minimiser order carries its minimiser in pos/mpos and the minimiser's
 * position in isize (set by minimise_record(), as samtools sort does, so
 * they survive temporary files); like samtools it sorts after placed reads,
 * by minimiser, then position, then strand. */
static inline void num_key(enum s2_order o, const bam1_core_t *c,
                           uint64_t *hi, uint64_t *lo, uint64_t *x) {
    if (o == ORDER_MINHASH && c->tid < 0) {
        *hi = 0x7fffffffu;
        *lo = ((uint64_t) c->pos << 31) | (uint32_t) c->mpos;
        // Larger isize (hash further left) first, then strand.
        *x = ((uint64_t) (0xffffffffu - (uint32_t) c->isize) << 1)
             | !!(c->flag & BAM_FREVERSE);
        return;
    }
    coord_key2(c, hi, lo);
    *x = 0;
}

/* Packed single-word key, (tid' << 33) | lo, when pos + 1 < 2^32 (all BAM
 * data).  Returns -1 if the position is too large to pack. */
static inline int coord_key(const bam1_core_t *c, uint64_t *key) {
    uint64_t hi, lo;
    coord_key2(c, &hi, &lo);
    if (lo >> 33)
        return -1;
    *key = (hi << 33) | lo;
    return 0;
}

/* ---------------------------------------------------------------------- */
/* Minimiser order (-M), as samtools sort's minhash() and reverse_complement() */

/* The minimiser functions below are samtools develop's minhash(),
 * build_minhash_index(), minhash_with_idx() and minhash_with_idx_squash()
 * (the last two merged, with a flag), apart from names and a 64-bit mask
 * constant.  The index is read-only once built, so the sort threads share
 * it safely. */
#define S2_XOR 0xdead7878beef7878ULL
#define S2_UNIQ_BIT  60
#define S2_UNIQ_TEST(x) (((x) & (1ULL<<S2_UNIQ_BIT))==0)
#define S2_UNIQ_MASK ((1ULL<<S2_UNIQ_BIT)-1)
KHASH_MAP_INIT_INT64(s2kmer, int64_t)

static uint64_t s2_minhash(bam1_t *b, int kmer, int window, int *curr_pos,
                           int *end, int *is_rev, int try_fwd, int try_rev,
                           int no_squash) {
    uint64_t hashf = 0, minhashf = UINT64_MAX;
    int minhashpf = *curr_pos, i, j;
    uint64_t mask = ((uint64_t) 1 << (2 * kmer)) - 1;
    uint8_t *seq = bam_get_seq(b);
    int len = b->core.l_qseq;
    uint64_t xor = S2_XOR & mask;
    if (is_rev) *is_rev = 0;
    // bam_seqi() codes to 0123 forward / reverse hashes (=ACM GRSV TWYH KDBN)
    static const unsigned char L[16] = { 0,0,1,0, 2,0,0,0, 3,0,0,0, 0,0,0,0 };
    uint64_t R[16] = { 0,3,2,0, 1,0,0,0, 0,0,0,0, 0,0,0,0 };
    for (i = 0; i < 16; i++)
        R[i] <<= 2 * (kmer - 1);
    int i_start = *curr_pos;
    int i_end = i_start + window < len ? i_start + window : len;
    int last_base = -1;
    if (try_fwd) {
        for (i = i_start, j = 0; j < kmer - 1 && i < i_end; i++) {
            int base = bam_seqi(seq, i);
            if (no_squash || last_base != base) {   // collapse homopolymers
                last_base = base;
                hashf = (hashf << 2) | L[base];
                j++;
            }
        }
        if (no_squash) {
            for (; i < i_end; i++) {
                int base = bam_seqi(seq, i);
                hashf = (hashf << 2) | L[base];
                uint64_t hashfx = (hashf ^ S2_XOR) & mask;
                if (minhashf > hashfx)
                    minhashf = hashfx, minhashpf = i;
            }
        } else {
            for (; i < i_end; i++) {
                int base = bam_seqi(seq, i);
                if (last_base != base) {
                    last_base = base;
                    hashf = (hashf << 2) | L[base];
                    uint64_t hashfx = (hashf ^ S2_XOR) & mask;
                    if (minhashf > hashfx)
                        minhashf = hashfx, minhashpf = i;
                }
            }
        }
    }
    if (try_rev) {
        uint64_t hashr = 0, minhashr = UINT64_MAX;
        int minhashpr = *curr_pos;
        int last_base = -1;
        for (i = i_start, j = 0; j < kmer - 1 && i < len; i++) {
            int base = bam_seqi(seq, i);
            if (no_squash || last_base != base) {
                last_base = base;
                hashr = (hashr >> 2) | R[base];
                j++;
            }
        }
        if (no_squash) {
            for (; i < i_end; i++) {
                int base = bam_seqi(seq, i);
                hashr = (hashr >> 2) | R[base];
                if (minhashr > (hashr ^ xor))
                    minhashr = (hashr ^ xor), minhashpr = len - i + kmer - 2;
            }
        } else {
            for (; i < i_end; i++) {
                int base = bam_seqi(seq, i);
                if (last_base != base) {
                    last_base = base;
                    hashr = (hashr >> 2) | R[base];
                    if (minhashr > (hashr ^ xor))
                        minhashr = (hashr ^ xor), minhashpr = len - i + kmer - 2;
                }
            }
        }
        if (minhashr < minhashf) {
            minhashf = minhashr;
            minhashpf = minhashpr;
            if (is_rev) *is_rev = 1;
        }
    }
    *curr_pos = minhashpf - (kmer - 1);
    if (end) *end = (i_end == len);
    return minhashf;
}

// -I ref.fa: hash -> position in the concatenated reference (+ dup bit).
static int s2_build_minhash_index(s2_state_t *st, const char *fn) {
    int ret = 1, r;
    samFile *in = sam_open(fn, "r");
    sam_hdr_t *h = NULL;
    bam1_t *b = NULL;
    khash_t(s2kmer) *kh;
    if (!in) {
        perror(fn);
        return 1;
    }
    if (!(kh = kh_init(s2kmer)))
        goto err;
    st->kmer_h = kh;
    if (!(h = sam_hdr_read(in)) || !(b = bam_init1()))
        goto err;
    uint64_t tpos = 0;
    while ((r = sam_read1(in, h, b)) >= 0) {
        int pos = 0, end = 0, kret;
        if (b->core.l_qseq < st->window)
            continue;
        while (!end) {
            int last_pos = pos;
            uint64_t hashf = s2_minhash(b, st->kmer, st->window, &pos, &end,
                                        NULL, 1, 0, st->no_squash);
            khiter_t k = kh_put(s2kmer, kh, hashf, &kret);
            if (kret < 0)
                goto err;
            kh_value(kh, k) = tpos + pos + (((uint64_t) !kret) << S2_UNIQ_BIT);
            pos = last_pos + st->kmer > pos + 1 ? last_pos + st->kmer : pos + 1;
        }
        tpos += b->core.l_qseq;
    }
    if (r < -1)
        goto err;
    ret = 0;
 err:
    if (b) bam_destroy1(b);
    if (h) sam_hdr_destroy(h);
    if (sam_close(in) < 0)
        ret = 1;
    return ret;
}

static void s2_minhash_index_free(void *vp) {
    if (vp)
        kh_destroy(s2kmer, (khash_t(s2kmer) *) vp);
}

/* Minimiser against the -I index: prefer k-mers uniquely in the index, then
 * duplicated ones, then any; the key is then the position in the reference.
 * With squash set, homopolymers are collapsed (minhash_with_idx_squash()). */
static uint64_t s2_minhash_idx(const void *vkh, bam1_t *b, int kmer,
                               int *pos, int *rev, int try_rev, int squash) {
    const khash_t(s2kmer) *kh = vkh;
    uint64_t hashf = 0, minhashf = UINT64_MAX, minhashfi = UINT64_MAX;
    uint64_t minhashfd = UINT64_MAX;
    int minhashpf = 0, minhashpfi = 0, minhashpfd = 0, i, j;
    uint64_t mask = ((uint64_t) 1 << (2 * kmer)) - 1;
    unsigned char *seq = bam_get_seq(b);
    int len = b->core.l_qseq;
    const uint64_t xor = S2_XOR & mask;
    static const unsigned char L[16] = { 0,0,1,0, 2,0,0,0, 3,0,0,0, 0,0,0,0 };
    uint64_t R[16] = { 0,3,2,0, 1,0,0,0, 0,0,0,0, 0,0,0,0 };
    for (i = 0; i < 16; i++)
        R[i] <<= 2 * (kmer - 1);
    int last_base = -1;
    for (i = j = 0; j < kmer - 1 && i < len; i++) {
        int base = bam_seqi(seq, i);
        if (squash) {
            if (base == last_base)
                continue;
            last_base = base;
        }
        j++;
        hashf = (hashf << 2) | L[base];
    }
    int found_f = 0, found_r = 0;
    for (; i < len; i++) {
        int base = bam_seqi(seq, i);
        if (squash) {
            if (base == last_base)
                continue;
            last_base = base;
        }
        hashf = ((hashf << 2) | L[base]) & mask;
        const uint64_t hashfx = hashf ^ xor;
        // Priority: 1. unique key in the index, 2. duplicated key, 3. other
        int index = 0;
        if (minhashfi > hashfx || (found_f < 2 && minhashfd > hashfx)) {
            khiter_t k = kh_get(s2kmer, kh, hashfx);
            if (k != kh_end(kh))
                index = S2_UNIQ_TEST(kh_value(kh, k)) ? 2 : 1;
        }
        found_f |= index;
        switch (index) {
        case 2: minhashfi = hashfx, minhashpfi = i; break;
        case 1: minhashfd = hashfx, minhashpfd = i; break;
        default:
            if (minhashf > hashfx)
                minhashf = hashfx, minhashpf = i;
        }
    }
    if (minhashfi != UINT64_MAX)
        minhashf = minhashfi, minhashpf = minhashpfi;
    else if (minhashfd != UINT64_MAX)
        minhashf = minhashfd, minhashpf = minhashpfd;
    int dir = 0;
    if (try_rev) {
        uint64_t hashr = 0, minhashr = UINT64_MAX, minhashri = UINT64_MAX;
        uint64_t minhashrd = UINT64_MAX;
        int minhashpr = 0, minhashpri = 0, minhashprd = 0;
        last_base = -1;
        for (i = j = 0; j < kmer - 1 && i < len; i++) {
            int base = bam_seqi(seq, i);
            if (squash) {
                if (base == last_base)
                    continue;
                last_base = base;
            }
            j++;
            hashr = (hashr >> 2) | R[base];
        }
        for (; i < len; i++) {
            int base = bam_seqi(seq, i);
            if (squash) {
                if (base == last_base)
                    continue;
                last_base = base;
            }
            hashr = (hashr >> 2) | R[base];
            const uint64_t hashrx = hashr ^ xor;
            int index = 0;
            if (minhashri > hashrx || (found_r < 2 && minhashrd > hashrx)) {
                khiter_t k = kh_get(s2kmer, kh, hashrx);
                if (k != kh_end(kh))
                    index = S2_UNIQ_TEST(kh_value(kh, k)) ? 2 : 1;
            }
            found_r |= index;
            switch (index) {
            case 2: minhashri = hashrx, minhashpri = i; break;
            case 1: minhashrd = hashrx, minhashprd = i; break;
            default:
                if (minhashr > hashrx)
                    minhashr = hashrx, minhashpr = i;
            }
        }
        if (minhashri != UINT64_MAX)
            minhashr = minhashri, minhashpr = minhashpri;
        else if (minhashrd != UINT64_MAX)
            minhashr = minhashrd, minhashpr = minhashprd;
        // Pick reverse if it maps better.
        if ((minhashf > minhashr) || (!found_f && found_r)) {
            if (!found_f || found_r) {
                minhashf = minhashr;
                minhashpf = b->core.l_qseq - minhashpr + kmer - 2;
                dir = 1;
            }
        }
    }
    // For indexed k-mers the key is the k-mer's reference position.
    khiter_t k = kh_get(s2kmer, kh, minhashf);
    if (k != kh_end(kh))
        minhashf = kh_value(kh, k) & S2_UNIQ_MASK;
    *rev = dir;
    *pos = minhashpf;
    return minhashf != UINT64_MAX ? minhashf : 0;
}

/* samtools sort's complement table, verbatim (so that output is identical,
 * including its handling of the rarer IUPAC codes). */
static const char s2_comp[256] = {
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    '@','T','V','G', 'H','E','F','C', 'D','I','H','M', 'L','K','N','O',
    'P','Q','Y','S', 'A','A','B','W', 'X','Y','Z','[','\\','[','^','_',
    '`','t','v','g', 'h','e','f','c', 'd','i','j','m', 'l','k','n','o',
    'p','q','y','s', 'a','a','b','w', 'x','y','z','{', '|','}','~',127,
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
    'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
};

// Reverse complement sequence and qualities in place, flipping BAM_FREVERSE.
static void s2_reverse_complement(bam1_t *b) {
    uint8_t comp16[16];
    for (int c = 0; c < 16; c++)
        comp16[c] = seq_nt16_table[(unsigned char) s2_comp[(unsigned char) seq_nt16_str[c]]];
    uint8_t *seq = bam_get_seq(b), *qual = bam_get_qual(b);
    int i, j;
    for (i = 0, j = b->core.l_qseq - 1; i < j; i++, j--) {
        int a = bam_seqi(seq, i), z = bam_seqi(seq, j);
        bam_set_seqi(seq, i, comp16[z]);
        bam_set_seqi(seq, j, comp16[a]);
        uint8_t t = qual[i]; qual[i] = qual[j]; qual[j] = t;
    }
    if (i == j)
        bam_set_seqi(seq, i, comp16[bam_seqi(seq, i)]);
    b->core.flag ^= BAM_FREVERSE;
}

/* ---------------------------------------------------------------------- */
/* Byte-comparable keys for -t TAG and --template-coordinate */

static inline uint8_t *put_be32(uint8_t *p, uint32_t v) {
    p[0] = v >> 24; p[1] = v >> 16; p[2] = v >> 8; p[3] = v;
    return p + 4;
}
static inline uint8_t *put_be64(uint8_t *p, uint64_t v) {
    for (int i = 7; i >= 0; i--)
        *p++ = (uint8_t) (v >> (8 * i));
    return p;
}
// Signed integers and doubles in an order-preserving unsigned form.
static inline uint64_t ord_i64(int64_t v) { return (uint64_t) v ^ ((uint64_t) 1 << 63); }
static inline uint64_t ord_f64(double d) {
    uint64_t u;
    if (d == 0)
        d = 0.0;               // -0.0 compares equal to 0.0
    memcpy(&u, &d, 8);
    return (u >> 63) ? ~u : (u | ((uint64_t) 1 << 63));
}

/* -t TAG: as samtools bam1_cmp_by_tag().  Reads without the tag first; then
 * by normalised type when types differ ('A' < 'B' < 'H'/'Z' < numbers, with
 * integers and floats compared numerically); then the value (integers as
 * int64, floats as doubles, characters, strings by strcmp; B arrays equal);
 * then coordinate, or name with -n/-N. */
static size_t tag_key_bound(const s2_state_t *st, bam1_t *v) {
    const uint8_t *aux = bam_aux_get(v, st->sort_tag);
    size_t n = 1 + 16 + 13 + s2_name_key_bound(v->core.l_qname) + 2;
    if (aux && (*aux == 'Z' || *aux == 'H'))
        n += strlen((const char *) aux + 1) + 1;
    return n;
}

static size_t tag_key(const s2_state_t *st, bam1_t *v, uint8_t *o) {
    uint8_t *p = o;
    const uint8_t *aux = bam_aux_get(v, st->sort_tag);
    if (!aux) {
        *p++ = 0;
    } else {
        uint8_t t = *aux;
        if (strchr("cCsSiI", t)) {
            *p++ = 'c';        // numbers sort as one type
            int64_t iv = bam_aux2i(aux);
            p = put_be64(p, ord_f64((double) iv));
            p = put_be64(p, ord_i64(iv));
        } else if (t == 'f' || t == 'd') {
            *p++ = 'c';
            double d = bam_aux2f(aux);
            // An integral float ties with the equal integer, as there.
            int64_t iv = (d > -9.2e18 && d < 9.2e18) ? (int64_t) d : 0;
            p = put_be64(p, ord_f64(d));
            p = put_be64(p, ord_i64(iv));
        } else if (t == 'A') {
            *p++ = 'A';
            *p++ = (uint8_t) bam_aux2A(aux);
        } else if (t == 'Z' || t == 'H') {
            *p++ = 'H';
            const char *z = bam_aux2Z(aux);
            size_t l = strlen(z);
            memcpy(p, z, l);
            p += l;
            *p++ = 0;
        } else {
            *p++ = t;          // e.g. B arrays: compare equal
        }
    }
    if (st->order == ORDER_TAG_NAME) {
        p += s2_name_key(bam_get_qname(v), v->core.l_qname - 1
                         - v->core.l_extranul, st->tag_natural,
                         name_flag_byte(&v->core), p);
    } else {
        // As bam1_cmp_core(): tid as uint64 (unmapped last), pos+1, strand.
        p = put_be32(p, (uint32_t) v->core.tid);
        p = put_be64(p, (uint64_t) (v->core.pos + 1));
        *p++ = bam_is_rev(v);
    }
    return p - o;
}

/* --template-coordinate: as samtools template_coordinate_key() and
 * bam1_cmp_template_coordinate(): tid1, tid2, pos1, pos2 (unclipped 5'
 * ends, the mate's from its MC tag; soft clips only), neg1, neg2 (negative
 * first), library, CB, MI (ignoring a trailing "/x"), name, and whether
 * this is the upper end of the pair. */
// RG ID -> LB, sorted by ID (read groups are few).
typedef struct {
    size_t n;
    char **id, **lb;
} s2_libs_t;

static int str_qcmp(const void *a, const void *b) {
    return strcmp(*(char *const *) a, *(char *const *) b);
}

static const char *template_library(const s2_state_t *st, bam1_t *v) {
    const uint8_t *rg = bam_aux_get(v, "RG");
    const s2_libs_t *L = st->libs;
    if (rg && rg[0] == 'Z' && L && L->n) {
        const char *key = (const char *) rg + 1;
        char **hit = bsearch(&key, L->id, L->n, sizeof(char *), str_qcmp);
        if (hit)
            return L->lb[hit - L->id];
    }
    return "";
}

static size_t template_key_bound(const s2_state_t *st, bam1_t *v) {
    const uint8_t *mi = bam_aux_get(v, "MI"), *cb = bam_aux_get(v, "CB");
    size_t n = 4 + 4 + 8 + 8 + 2 + strlen(template_library(st, v)) + 1
        + v->core.l_qname + 1 + 1 + 1;
    if (cb && *cb == 'Z')
        n += strlen((const char *) cb + 1);
    if (mi && *mi == 'Z')
        n += strlen((const char *) mi + 1) + 1;
    else
        n += 1;
    return n;
}

// Returns the key length, or 0 (with a message) if the key can't be built.
static size_t template_key(const s2_state_t *st, bam1_t *b, uint8_t *o) {
    int tid1 = INT32_MAX, tid2 = INT32_MAX;
    hts_pos_t pos1 = HTS_POS_MAX, pos2 = HTS_POS_MAX;
    int neg1 = 0, neg2 = 0;
    const char *mid = "", *lib = template_library(st, b);
    const uint8_t *data;
    if (!(b->core.flag & BAM_FUNMAP)) {
        tid1 = b->core.tid;
        neg1 = bam_is_rev(b);
        pos1 = neg1 ? unclipped_end(b, 0) : unclipped_start(b, 0);
    }
    if ((b->core.flag & BAM_FPAIRED) && !(b->core.flag & BAM_FMUNMAP)) {
        char *cigar;
        if (!(data = bam_aux_get(b, "MC"))) {
            print_error("sort2", "no MC tag on \"%s\". Please run samtools "
                        "fixmate on the file first", bam_get_qname(b));
            return 0;
        }
        if (!(cigar = bam_aux2Z(data))) {
            print_error("sort2", "MC tag wrong type on \"%s\". Please use "
                        "the MC tag provided by samtools fixmate",
                        bam_get_qname(b));
            return 0;
        }
        tid2 = b->core.mtid;
        neg2 = bam_is_mrev(b);
        pos2 = neg2 ? unclipped_other_end(b->core.mpos, cigar, 0)
                    : unclipped_other_start(b->core.mpos, cigar, 0);
    }
    const char *cid = "";
    if ((data = bam_aux_get(b, "CB")) && !(cid = bam_aux2Z(data))) {
        print_error("sort2", "CB tag wrong type (not a string) on \"%s\"",
                    bam_get_qname(b));
        return 0;
    }
    if ((data = bam_aux_get(b, "MI")) && !(mid = bam_aux2Z(data))) {
        print_error("sort2", "MI tag wrong type (not a string) on \"%s\"",
                    bam_get_qname(b));
        return 0;
    }
    // Lower end first; is_upper records whether this read was swapped.
    int is_upper = !(tid1 < tid2 || (tid1 == tid2 && pos1 < pos2)
                     || (tid1 == tid2 && pos1 == pos2 && !neg1));
    if (is_upper) {
        int t = tid1; tid1 = tid2; tid2 = t;
        hts_pos_t q = pos1; pos1 = pos2; pos2 = q;
        int n = neg1; neg1 = neg2; neg2 = n;
    }
    uint8_t *p = o;
    p = put_be32(p, (uint32_t) tid1);       // tids are >= 0 or INT32_MAX
    p = put_be32(p, (uint32_t) tid2);
    p = put_be64(p, ord_i64(pos1));
    p = put_be64(p, ord_i64(pos2));
    *p++ = neg1 ? 0 : 1;                    // negative strand first
    *p++ = neg2 ? 0 : 1;
    size_t l = strlen(lib);
    memcpy(p, lib, l); p += l; *p++ = 0;
    l = strlen(cid);
    memcpy(p, cid, l); p += l; *p++ = 0;
    l = strlen(mid);
    if (l >= 2 && mid[l - 2] == '/')
        l -= 2;                             // trailing "/x" is ignored
    memcpy(p, mid, l); p += l; *p++ = 0;
    const char *qn = bam_get_qname(b);
    l = strlen(qn);
    memcpy(p, qn, l); p += l; *p++ = 0;
    *p++ = (uint8_t) is_upper;
    return p - o;
}

// RG ID -> LB, for --template-coordinate (as samtools lookup_libraries():
// the first read group with each ID wins).
static void template_libraries_free(void *vp) {
    s2_libs_t *L = vp;
    if (!L)
        return;
    for (size_t i = 0; i < L->n; i++) {
        free(L->id[i]);
        free(L->lb[i]);
    }
    free(L->id);
    free(L->lb);
    free(L);
}

static void *template_libraries(sam_hdr_t *h) {
    int n = sam_hdr_count_lines(h, "RG");
    s2_libs_t *L = calloc(1, sizeof(*L));
    if (n < 0 || !L)
        goto fail;
    // Pairs of (ID, LB) sorted together: sort indexes by ID, stably.
    char **id = calloc(n ? n : 1, sizeof(char *)),
         **lb = calloc(n ? n : 1, sizeof(char *));
    L->id = id;
    L->lb = lb;
    if (!id || !lb)
        goto fail;
    kstring_t ks = { 0, 0, NULL };
    for (int i = 0; i < n; i++) {
        const char *rid = sam_hdr_line_name(h, "RG", i);
        ks.l = 0;
        if (!rid || sam_hdr_find_tag_pos(h, "RG", i, "LB", &ks) < 0 || !ks.s)
            continue;
        int dup = 0;
        for (size_t j = 0; j < L->n; j++)
            if (strcmp(id[j], rid) == 0)
                dup = 1;
        if (dup)
            continue;
        if (!(id[L->n] = strdup(rid)) || !(lb[L->n] = strdup(ks.s))) {
            free(id[L->n]);
            free(ks.s);
            goto fail;
        }
        L->n++;
    }
    free(ks.s);
    // Insertion sort of the pairs by ID (few read groups).
    for (size_t i = 1; i < L->n; i++)
        for (size_t j = i; j > 0 && strcmp(id[j - 1], id[j]) > 0; j--) {
            char *t = id[j]; id[j] = id[j - 1]; id[j - 1] = t;
            t = lb[j]; lb[j] = lb[j - 1]; lb[j - 1] = t;
        }
    return L;
 fail:
    template_libraries_free(L);
    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Batches */

static void *big_alloc(size_t sz) {
#ifdef MADV_HUGEPAGE
    void *p = NULL;
    // Huge pages cut TLB misses in the radix scatter and gather passes.
    if (sz >= ((size_t)4 << 20) && posix_memalign(&p, (size_t)2 << 20, sz) == 0) {
        madvise(p, sz, MADV_HUGEPAGE);
        return p;
    }
#endif
    return malloc(sz);
}

static void batch_reset(s2_batch_t *b) {
    b->wide = 0;
    b->used = 0;
    b->n = 0;
    b->charged = 0;
    b->sorted = NULL;
    b->eof = 0;
}

static void batch_free(s2_batch_t *b, int nthreads) {
    free(b->arena);
    free(b->off);
    free(b->kptr);
    free(b->klen);
    if (b->tbuf)
        for (int t = 0; t < nthreads; t++)
            free(b->tbuf[t]);
    free(b->tbuf);
    free(b->tbuf_cap);
    free(b->kv);
    free(b->tmp);
    memset(b, 0, sizeof(*b));
}

// Bytes charged per record beyond the record itself.
static size_t record_overhead(const s2_state_t *st, const bam1_t *r) {
    size_t o = sizeof(size_t) + 2 * sizeof(s2_kv_t);
    if (!order_numeric(st->order))
        o += sizeof(uint8_t *) + sizeof(uint32_t)
            + s2_name_key_bound(r->core.l_qname) / 2
            + (st->order == ORDER_NAME_NAT || st->order == ORDER_NAME_BYTE
               ? 0 : 48);
    return o;
}

/* Append a record.  Returns 1 if added, 0 if the batch is full (and not
 * empty), -1 on error. */
static int batch_add(s2_state_t *st, s2_batch_t *b, const bam1_t *r) {
    size_t need = REC_DATA_OFF + ALIGN8((size_t) r->l_data);
    size_t charge = need + record_overhead(st, r);
    if (b->n > 0 && b->charged + charge > st->batch_limit)
        return 0;

    if (b->used + need > b->cap) {
        // Only an oversized first record (or a tiny -m) gets here.
        size_t ncap = b->used + need;
        if (ncap < st->batch_limit)
            ncap = st->batch_limit;
        uint8_t *a = realloc(b->arena, ncap);
        if (!a)
            return -1;
        b->arena = a;
        b->cap = ncap;
    }
    if (b->n == b->m_off) {
        size_t m = b->m_off ? b->m_off + b->m_off / 2 : 65536;
        size_t *o = realloc(b->off, m * sizeof(*o));
        if (!o)
            return -1;
        b->off = o;
        b->m_off = m;
    }

    s2_rec_t *h = (s2_rec_t *) (b->arena + b->used);
    memset(h, 0, REC_DATA_OFF);
    h->core = r->core;
    h->l_data = r->l_data;
    memcpy(b->arena + b->used + REC_DATA_OFF, r->data, r->l_data);
    b->off[b->n++] = b->used;
    b->used += need;
    b->charged += charge;
    return 1;
}

static inline s2_rec_t *batch_rec(const s2_batch_t *b, size_t i) {
    return (s2_rec_t *) (b->arena + b->off[i]);
}

typedef struct {
    s2_state_t *st;
    s2_batch_t *b;
    int fail;                  // a worker ran out of memory
    int wide;                  // a coordinate key needs more than 64 bits
    pthread_mutex_t m;         // guards fail and wide (rarely taken)
} kx_arg_t;

static void kx_flag(kx_arg_t *x, int *flag, int v) {
    pthread_mutex_lock(&x->m);
    *flag = v;
    pthread_mutex_unlock(&x->m);
}

/* Wide coordinate keys: kv[].key = ((pos+1) << 1) | rev here, then tid' in
 * a second pass; two stable radix sorts give (tid', pos, rev) order. */
static void kx_coord_lo(void *vp, int t, int nt) {
    kx_arg_t *x = vp;
    s2_batch_t *b = x->b;
    size_t lo, hi;
    s2_chunk(b->n, t, nt, &lo, &hi);
    for (size_t i = lo; i < hi; i++) {
        uint64_t khi;
        coord_key2(&batch_rec(b, i)->core, &khi, &b->kv[i].key);
        b->kv[i].ref = i;
    }
}

typedef struct {
    s2_batch_t *b;
    s2_kv_t *a;
} kx_tid_t;

static void kx_coord_tid(void *vp, int t, int nt) {
    kx_tid_t *x = vp;
    size_t lo, hi;
    s2_chunk(x->b->n, t, nt, &lo, &hi);
    for (size_t i = lo; i < hi; i++) {
        uint64_t khi, klo;
        coord_key2(&batch_rec(x->b, x->a[i].ref)->core, &khi, &klo);
        x->a[i].key = khi;
    }
}

static void kx_coord(void *vp, int t, int nt) {
    kx_arg_t *x = vp;
    s2_batch_t *b = x->b;
    size_t lo, hi;
    int wide = 0;
    s2_chunk(b->n, t, nt, &lo, &hi);
    for (size_t i = lo; i < hi; i++) {
        if (i + 8 < hi)
            S2_PREFETCH(b->arena + b->off[i + 8]);
        if (coord_key(&batch_rec(b, i)->core, &b->kv[i].key) < 0)
            wide = 1;
        b->kv[i].ref = i;
    }
    if (wide)
        kx_flag(x, &x->wide, 1);
}

static inline void rec_view(bam1_t *v, s2_rec_t *r);

/* Byte keys for -t and --template-coordinate.  Fails with x->fail = 2 if a
 * template key can't be built (the message has been printed). */
static void kx_bytes(void *vp, int t, int nt) {
    kx_arg_t *x = vp;
    s2_batch_t *b = x->b;
    const s2_state_t *st = x->st;
    int failed = 0, is_tag = st->order != ORDER_TEMPLATE;
    size_t lo, hi, pos = 0;
    s2_chunk(b->n, t, nt, &lo, &hi);
    uint8_t *buf = b->tbuf[t];
    size_t cap = b->tbuf_cap[t];
    bam1_t v;
    memset(&v, 0, sizeof(v));
    for (size_t i = lo; i < hi && !failed; i++) {
        rec_view(&v, batch_rec(b, i));
        size_t bound = is_tag ? tag_key_bound(st, &v)
                              : template_key_bound(st, &v);
        if (pos + bound > cap) {
            size_t ncap = (pos + bound) * 3 / 2 + 4096;
            uint8_t *nb = realloc(buf, ncap);
            if (!nb) {
                failed = 1;
                break;
            }
            buf = nb;
            cap = ncap;
        }
        size_t l = is_tag ? tag_key(st, &v, buf + pos)
                          : template_key(st, &v, buf + pos);
        if (!l) {
            failed = 2;
            break;
        }
        b->klen[i] = (uint32_t) l;
        b->kv[i].key = pos;    // offset for now; buffer may still move
        b->kv[i].ref = i;
        pos += l;
    }
    b->tbuf[t] = buf;
    b->tbuf_cap[t] = cap;
    if (failed) {
        kx_flag(x, &x->fail, failed);
        return;
    }
    for (size_t i = lo; i < hi; i++)
        b->kptr[i] = buf + b->kv[i].key;
}

/* -M: find each unplaced read's minimiser, reverse complementing the read
 * if it was on the reverse strand, and store it in pos/mpos/isize, exactly
 * as samtools sort's worker_minhash().  Then the first sort key (hash
 * position and strand; see num_key()). */
static void kx_minhash(void *vp, int t, int nt) {
    kx_arg_t *x = vp;
    s2_batch_t *b = x->b;
    size_t lo, hi;
    s2_chunk(b->n, t, nt, &lo, &hi);
    bam1_t v;
    memset(&v, 0, sizeof(v));
    for (size_t i = lo; i < hi; i++) {
        s2_rec_t *r = batch_rec(b, i);
        if (r->core.tid == -1) {
            const s2_state_t *st = x->st;
            rec_view(&v, r);
            int pos = 0, rev = 0;
            uint64_t mh = st->kmer_h
                ? s2_minhash_idx(st->kmer_h, &v, st->kmer, &pos, &rev,
                                 st->try_rev, !st->no_squash)
                : s2_minhash(&v, st->kmer, v.core.l_qseq, &pos, NULL, &rev,
                             1, st->try_rev, st->no_squash);
            if (rev)
                s2_reverse_complement(&v);   // seq/qual in place
            // As samtools worker_minhash().
            if (!st->kmer_h) {
                mh += 1ULL << 30;
                pos = 65535 - pos >= 0 ? 65535 - pos : 0;
            } else {
                mh -= pos;
                pos = 0;
            }
            r->core.flag = v.core.flag;
            r->core.pos = (mh >> 31) & 0x7fffffff;
            r->core.mpos = mh & 0x7fffffff;
            r->core.isize = pos;
        }
        uint64_t khi, klo, kx;
        num_key(ORDER_MINHASH, &r->core, &khi, &klo, &kx);
        b->kv[i].key = kx;
        b->kv[i].ref = i;
    }
}

// Replace each sorted entry's key with one component of its numeric key.
typedef struct {
    s2_batch_t *b;
    s2_kv_t *a;
    enum s2_order o;
    int which;                 // 0 = hi, 1 = lo
} kx_field_t;

static void kx_field(void *vp, int t, int nt) {
    kx_field_t *x = vp;
    size_t lo, hi;
    s2_chunk(x->b->n, t, nt, &lo, &hi);
    for (size_t i = lo; i < hi; i++) {
        uint64_t khi, klo, kx;
        num_key(x->o, &batch_rec(x->b, x->a[i].ref)->core, &khi, &klo, &kx);
        x->a[i].key = x->which ? klo : khi;
    }
}

static void kx_name(void *vp, int t, int nt) {
    kx_arg_t *x = vp;
    s2_batch_t *b = x->b;
    int natural = x->st->order == ORDER_NAME_NAT;
    int failed = 0;
    size_t lo, hi, pos = 0;
    s2_chunk(b->n, t, nt, &lo, &hi);
    uint8_t *buf = b->tbuf[t];
    size_t cap = b->tbuf_cap[t];
    for (size_t i = lo; i < hi; i++) {
        if (i + 8 < hi)
            S2_PREFETCH(b->arena + b->off[i + 8]);
        s2_rec_t *r = batch_rec(b, i);
        uint32_t nl = name_len(r);
        size_t bound = s2_name_key_bound(nl);
        if (pos + bound > cap) {
            size_t ncap = (pos + bound) * 3 / 2 + 4096;
            uint8_t *nb = realloc(buf, ncap);
            if (!nb) {
                failed = 1;
                break;
            }
            buf = nb;
            cap = ncap;
        }
        const char *qn = (const char *) b->arena + b->off[i] + REC_DATA_OFF;
        b->klen[i] = s2_name_key(qn, nl, natural, name_flag_byte(&r->core),
                                 buf + pos);
        b->kv[i].key = pos;    // offset for now; buffer may still move
        b->kv[i].ref = i;
        pos += b->klen[i];
    }
    b->tbuf[t] = buf;
    b->tbuf_cap[t] = cap;
    if (failed) {
        kx_flag(x, &x->fail, 1);
        return;
    }
    for (size_t i = lo; i < hi; i++)
        b->kptr[i] = buf + b->kv[i].key;
}

static int batch_sort(s2_state_t *st, s2_batch_t *b) {
    double t0 = now();
    if (b->n > b->m_kv) {
        free(b->kv);
        free(b->tmp);
        b->kv  = big_alloc(b->n * sizeof(s2_kv_t));
        b->tmp = big_alloc(b->n * sizeof(s2_kv_t));
        b->m_kv = b->n;
        if (!order_numeric(st->order)) {
            free(b->kptr);
            free(b->klen);
            b->kptr = malloc(b->n * sizeof(*b->kptr));
            b->klen = malloc(b->n * sizeof(*b->klen));
        }
        if (!b->kv || !b->tmp
            || (!order_numeric(st->order) && (!b->kptr || !b->klen))) {
            b->m_kv = 0;
            print_error("sort2", "couldn't allocate memory for sort keys");
            return -1;
        }
    }
    if (!order_numeric(st->order) && !b->tbuf) {
        b->tbuf = calloc(st->nthreads, sizeof(*b->tbuf));
        b->tbuf_cap = calloc(st->nthreads, sizeof(*b->tbuf_cap));
        if (!b->tbuf || !b->tbuf_cap) {
            print_error("sort2", "couldn't allocate memory for name keys");
            return -1;
        }
    }

    int nt = b->n >= 65536 ? st->nthreads : 1;
    kx_arg_t x;
    memset(&x, 0, sizeof(x));
    x.st = st;
    x.b = b;
    pthread_mutex_init(&x.m, NULL);
    void (*kx)(void *, int, int) =
        st->order == ORDER_COORD ? kx_coord
        : st->order == ORDER_MINHASH ? kx_minhash
        : st->order == ORDER_NAME_NAT || st->order == ORDER_NAME_BYTE ? kx_name
        : kx_bytes;
    s2_parallel(nt, kx, &x);
    if (x.fail) {
        pthread_mutex_destroy(&x.m);
        if (x.fail == 1)
            print_error("sort2", "couldn't allocate memory for sort keys");
        return -1;
    }
    double t1 = now();

    if (st->order == ORDER_COORD && x.wide) {
        // Positions beyond 32 bits (SAM allows them): sort by position and
        // strand, then stably by tid.
        b->wide = 1;
        s2_parallel(nt, kx_coord_lo, &x);
        s2_kv_t *r = s2_radix_sort(b->kv, b->tmp, b->n, st->nthreads);
        kx_tid_t xt = { b, r };
        s2_parallel(nt, kx_coord_tid, &xt);
        b->sorted = s2_radix_sort(r, r == b->kv ? b->tmp : b->kv, b->n,
                                  st->nthreads);
    } else if (st->order == ORDER_COORD) {
        b->sorted = s2_radix_sort(b->kv, b->tmp, b->n, st->nthreads);
    } else if (st->order == ORDER_MINHASH) {
        // Three stable passes: hash position and strand, then the hash (or
        // position), then the tid.
        b->wide = 1;
        s2_kv_t *r = s2_radix_sort(b->kv, b->tmp, b->n, st->nthreads);
        kx_field_t xf = { b, r, ORDER_MINHASH, 1 };
        s2_parallel(nt, kx_field, &xf);
        r = s2_radix_sort(r, r == b->kv ? b->tmp : b->kv, b->n, st->nthreads);
        xf.a = r;
        xf.which = 0;
        s2_parallel(nt, kx_field, &xf);
        b->sorted = s2_radix_sort(r, r == b->kv ? b->tmp : b->kv, b->n,
                                  st->nthreads);
    } else {
        s2_name_sort(b->kv, b->tmp, b->n, b->kptr, b->klen, st->nthreads);
        b->sorted = b->kv;
    }
    double t2 = now();
    pthread_mutex_destroy(&x.m);
    st->t_keys += t1 - t0;
    st->t_sort += t2 - t1;
    return 0;
}

/* ---------------------------------------------------------------------- */
/* Keys and bounds */

static inline int kcmp(enum s2_order o, const s2_key_t *a, const s2_key_t *b) {
    if (order_numeric(o)) {
        if (a->k_hi != b->k_hi)
            return a->k_hi < b->k_hi ? -1 : 1;
        if (a->k64 != b->k64)
            return a->k64 < b->k64 ? -1 : 1;
        return (a->k_x > b->k_x) - (a->k_x < b->k_x);
    }
    return s2_key_cmp(a->p, a->len, b->p, b->len);
}

// Sign of (k, src) relative to bound b.
static inline int bcmp_src(enum s2_order o, const s2_key_t *k, int src,
                           const s2_bound_t *b) {
    int c = kcmp(o, k, &b->k);
    if (c)
        return c;
    return (src > b->src) - (src < b->src);
}

// Key of the i-th record in a batch's sorted order.
static inline s2_key_t batch_key(const s2_state_t *st, const s2_batch_t *b,
                                 size_t i) {
    s2_key_t k = {0, NULL, 0, 0, 0};
    const s2_kv_t *kv = &b->sorted[i];
    if (order_numeric(st->order)) {
        if (!b->wide) {
            k.k_hi = kv->key >> 33;
            k.k64 = kv->key & (((uint64_t) 1 << 33) - 1);
        } else {
            num_key(st->order, &batch_rec(b, kv->ref)->core, &k.k_hi,
                    &k.k64, &k.k_x);
        }
    } else {
        k.p = b->kptr[kv->ref];
        k.len = b->klen[kv->ref];
    }
    return k;
}

/* ---------------------------------------------------------------------- */
/* BGZF block helpers */

static const uint8_t bgzf_eof_block[28] = {
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00,
    0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00
};

typedef struct {
    uint8_t *p;
    size_t len, cap;
} s2_buf_t;

static int buf_reserve(s2_buf_t *b, size_t extra) {
    if (b->len + extra <= b->cap)
        return 0;
    size_t ncap = (b->len + extra) + (b->len + extra) / 2 + 4096;
    uint8_t *np = realloc(b->p, ncap);
    if (!np)
        return -1;
    b->p = np;
    b->cap = ncap;
    return 0;
}

// Compress src (at most BGZF_BLOCK_SIZE bytes) as one block appended to out.
static int buf_add_block(s2_buf_t *out, const uint8_t *src, size_t slen,
                         int level) {
    if (buf_reserve(out, BGZF_MAX_BLOCK_SIZE) < 0)
        return -1;
    size_t dlen = BGZF_MAX_BLOCK_SIZE;
    if (bgzf_compress(out->p + out->len, &dlen, src, slen, level) != 0)
        return -1;
    out->len += dlen;
    return 0;
}

// Compress all of src as consecutive blocks appended to out.
static int buf_add_blocks(s2_buf_t *out, const uint8_t *src, size_t slen,
                          int level) {
    for (size_t o = 0; o < slen; o += BGZF_BLOCK_SIZE) {
        size_t l = slen - o < BGZF_BLOCK_SIZE ? slen - o : BGZF_BLOCK_SIZE;
        if (buf_add_block(out, src + o, l, level) < 0)
            return -1;
    }
    return 0;
}

static int write_all(int fd, const void *p, size_t len) {
    const uint8_t *c = p;
    while (len > 0) {
        ssize_t w = write(fd, c, len);
        if (w < 0) {
            if (errno == EINTR)
                continue;
            return -1;
        }
        c += w;
        len -= w;
    }
    return 0;
}

/* ---------------------------------------------------------------------- */
/* Run files */

static inline size_t frame_len(const s2_key_t *k, enum s2_order o,
                               const s2_rec_t *r) {
    uint32_t klen = order_numeric(o) ? 0 : k->len;
    return 8 + ALIGN8((size_t) klen) + REC_DATA_OFF + r->l_data;
}

// Append one framed record; buf must have frame_len() bytes reserved.
static inline void frame_put(uint8_t *o, const s2_key_t *k, enum s2_order ord,
                             const s2_rec_t *r) {
    // Numeric keys are recomputed from the record when read back.
    uint32_t klen = order_numeric(ord) ? 0 : k->len;
    size_t kpad = ALIGN8((size_t) klen), rlen = REC_DATA_OFF + r->l_data;
    uint32_t hdr[2] = { (uint32_t) (kpad + rlen), klen };
    memcpy(o, hdr, 8);
    if (klen)
        memcpy(o + 8, k->p, klen);
    memset(o + 8 + klen, 0, kpad - klen);
    memcpy(o + 8 + kpad, r, rlen);
}

/* Remove a run file from the filesystem.  Called as soon as nothing else
 * needs to open it by name; its space is released when the last reader
 * closes it.  The time taken is reported in --stats. */
static double g_t_unlink;      // only one thread unlinks at a time
static void run_unlink(s2_run_t *r) {
    if (r->fn && !r->unlinked) {
        double t0 = now();
        unlink(r->fn);
        r->unlinked = 1;
        g_t_unlink += now() - t0;
    }
}

static void run_free(s2_run_t *r, int do_unlink) {
    if (!r->fn)
        return;
    if (do_unlink)
        run_unlink(r);
    free(r->fn);
    for (size_t i = 0; i < r->nidx; i++)
        free((void *) r->idx[i].k.p);
    free(r->idx);
    memset(r, 0, sizeof(*r));
}

// Create a run file, writing the magic block.  Returns fd or -1.
static int run_create(s2_state_t *st, s2_run_t *run) {
    memset(run, 0, sizeof(*run));
    run->fn = malloc(strlen(st->tmp_prefix) + 32);
    if (!run->fn)
        return -1;
    sprintf(run->fn, "%s.%04d.s2r", st->tmp_prefix, st->run_seq++);
    int fd = open(run->fn, O_WRONLY | O_CREAT | O_TRUNC, 0600);
    if (fd < 0) {
        print_error_errno("sort2", "couldn't create temporary file \"%s\"",
                          run->fn);
        free(run->fn);
        run->fn = NULL;
        return -1;
    }
    s2_buf_t b = {0};
    if (buf_add_block(&b, (const uint8_t *) S2_RUN_MAGIC, 4, st->tmp_level) < 0
        || write_all(fd, b.p, b.len) < 0) {
        print_error_errno("sort2", "failed writing temporary file \"%s\"",
                          run->fn);
        free(b.p);
        close(fd);
        run_free(run, 1);
        return -1;
    }
    free(b.p);
    return fd;
}

static int run_add(s2_state_t *st, s2_run_t *run) {
    if (st->n_runs == st->m_runs) {
        int m = st->m_runs ? st->m_runs * 2 : 32;
        s2_run_t *r = realloc(st->runs, m * sizeof(*r));
        if (!r)
            return -1;
        st->runs = r;
        st->m_runs = m;
    }
    st->runs[st->n_runs++] = *run;
    return 0;
}

// Fill an index entry, copying name key bytes.
static int idx_set(const s2_state_t *st, s2_idx_t *e, uint64_t coff,
                   uint64_t ord, const s2_key_t *k) {
    e->coff = coff;
    e->ord = ord;
    e->k = *k;
    if (!order_numeric(st->order)) {
        uint8_t *p = malloc(k->len ? k->len : 1);
        if (!p)
            return -1;
        memcpy(p, k->p, k->len);
        e->k.p = p;
    }
    return 0;
}

/* Parallel run writer.  Workers claim groups of S2_GROUP_CHUNKS chunks in
 * order and serialise and compress them privately.  Under the lock a worker
 * only waits for the previous group's size to be known and reserves its own
 * file range; the pwrite() itself happens outside the lock, so several
 * large writes can be in flight at once. */
#define S2_GROUP_CHUNKS 16    // ~2 MB uncompressed per write

static int pwrite_all(int fd, const void *p, size_t len, off_t off) {
    const uint8_t *c = p;
    while (len > 0) {
        ssize_t w = pwrite(fd, c, len, off);
        if (w < 0) {
            if (errno == EINTR)
                continue;
            return -1;
        }
        c += w;
        off += w;
        len -= w;
    }
    return 0;
}

typedef struct {
    s2_state_t *st;
    s2_batch_t *b;
    int fd;
    size_t crec, nchunk, ngroup;
    pthread_mutex_t m;
    pthread_cond_t cv;
    size_t next_claim, next_write;
    uint64_t off;              // end of the space reserved so far
    uint64_t *chunk_off;
    uint64_t ubytes;
    int fail;
} rw_t;

static void rw_worker(void *vp, int t, int nt) {
    (void) t; (void) nt;
    rw_t *w = vp;
    s2_state_t *st = w->st;
    s2_batch_t *b = w->b;
    s2_buf_t u = {0}, c = {0};
    uint64_t ubytes = 0, rel[S2_GROUP_CHUNKS];
    for (;;) {
        pthread_mutex_lock(&w->m);
        size_t g = w->fail ? w->ngroup : w->next_claim++;
        pthread_mutex_unlock(&w->m);
        if (g >= w->ngroup)
            break;

        size_t k0 = g * S2_GROUP_CHUNKS, k1 = k0 + S2_GROUP_CHUNKS;
        if (k1 > w->nchunk)
            k1 = w->nchunk;
        c.len = 0;
        int err = 0;
        for (size_t k = k0; k < k1 && !err; k++) {
            size_t lo = k * w->crec, hi = lo + w->crec;
            if (hi > b->n)
                hi = b->n;
            u.len = 0;
            for (size_t i = lo; i < hi; i++) {
                if (i + 8 < hi)
                    S2_PREFETCH(b->arena + b->off[b->sorted[i + 8].ref]);
                s2_key_t key = batch_key(st, b, i);
                const s2_rec_t *r = batch_rec(b, b->sorted[i].ref);
                size_t fl = frame_len(&key, st->order, r);
                if (buf_reserve(&u, fl) < 0) {
                    err = 1;
                    break;
                }
                frame_put(u.p + u.len, &key, st->order, r);
                u.len += fl;
            }
            // Each chunk starts a new block, so its offset is an index point.
            rel[k - k0] = c.len;
            if (!err && buf_add_blocks(&c, u.p, u.len, st->tmp_level) < 0)
                err = 1;
            ubytes += u.len;
        }

        // Reserve file space in group order.
        uint64_t base = 0;
        pthread_mutex_lock(&w->m);
        while (!err && !w->fail && w->next_write != g)
            pthread_cond_wait(&w->cv, &w->m);
        if (!err && !w->fail) {
            base = w->off;
            w->off += c.len;
            w->next_write++;
        }
        if (err)
            w->fail = 1;
        int give_up = w->fail;
        pthread_cond_broadcast(&w->cv);
        pthread_mutex_unlock(&w->m);
        if (give_up)
            break;

        if (pwrite_all(w->fd, c.p, c.len, (off_t) base) < 0) {
            pthread_mutex_lock(&w->m);
            w->fail = 1;
            pthread_cond_broadcast(&w->cv);
            pthread_mutex_unlock(&w->m);
            break;
        }
        for (size_t k = k0; k < k1; k++)
            w->chunk_off[k] = base + rel[k - k0];
    }
    pthread_mutex_lock(&w->m);
    w->ubytes += ubytes;
    pthread_mutex_unlock(&w->m);
    free(u.p);
    free(c.p);
}

static size_t chunk_records(const s2_state_t *st, const s2_batch_t *b) {
    if (b->n == 0)
        return 1;
    size_t avg = b->used / b->n + 16 + (order_numeric(st->order) ? 8 : 40);
    size_t crec = S2_CHUNK_BYTES / avg;
    return crec < 16 ? 16 : crec;
}

static int write_run(s2_state_t *st, s2_batch_t *b) {
    double t0 = now();
    s2_run_t run;
    int fd = run_create(st, &run);
    if (fd < 0)
        return -1;

    rw_t w;
    memset(&w, 0, sizeof(w));
    w.st = st;
    w.b = b;
    w.fd = fd;
    w.crec = chunk_records(st, b);
    w.nchunk = (b->n + w.crec - 1) / w.crec;
    w.ngroup = (w.nchunk + S2_GROUP_CHUNKS - 1) / S2_GROUP_CHUNKS;
    w.chunk_off = malloc((w.nchunk ? w.nchunk : 1) * sizeof(uint64_t));
    run.idx = calloc(w.nchunk ? w.nchunk : 1, sizeof(*run.idx));
    off_t magic_end = lseek(fd, 0, SEEK_CUR);
    int ret = -1;
    if (!w.chunk_off || !run.idx || magic_end < 0)
        goto out;
    w.off = magic_end;
    pthread_mutex_init(&w.m, NULL);
    pthread_cond_init(&w.cv, NULL);
    int nt = (int) (w.ngroup < (size_t) st->nthreads
                    ? w.ngroup : (size_t) st->nthreads);
    s2_parallel(nt > 0 ? nt : 1, rw_worker, &w);
    pthread_mutex_destroy(&w.m);
    pthread_cond_destroy(&w.cv);
    if (w.fail || pwrite_all(fd, bgzf_eof_block, sizeof(bgzf_eof_block),
                             (off_t) w.off) < 0)
        goto out;

    for (size_t k = 0; k < w.nchunk; k++) {
        s2_key_t key = batch_key(st, b, k * w.crec);
        if (idx_set(st, &run.idx[k], w.chunk_off[k], k * w.crec, &key) < 0)
            goto out;
        run.nidx = k + 1;
    }
    run.ubytes = w.ubytes;
    st->tmp_cbytes += w.off + sizeof(bgzf_eof_block);
    st->tmp_ubytes += w.ubytes;
    run.nrec = b->n;
    ret = 0;

 out:
    if (close(fd) < 0)
        ret = -1;
    free(w.chunk_off);
    if (ret == 0 && run_add(st, &run) < 0)
        ret = -1;
    if (ret < 0) {
        print_error_errno("sort2", "failed writing temporary file \"%s\"",
                          run.fn ? run.fn : "");
        run_free(&run, 1);
    }
    st->t_runs += now() - t0;
    return ret;
}

/* Serial run writer, used by intermediate merges. */
typedef struct {
    s2_state_t *st;
    s2_run_t run;
    int fd;
    uint64_t off;
    s2_buf_t u, c;
    size_t m_idx;
} cw_t;

static int cw_open(s2_state_t *st, cw_t *w) {
    memset(w, 0, sizeof(*w));
    w->st = st;
    w->fd = run_create(st, &w->run);
    if (w->fd < 0)
        return -1;
    off_t o = lseek(w->fd, 0, SEEK_CUR);
    if (o < 0)
        return -1;
    w->off = o;
    return 0;
}

static int cw_flush(cw_t *w) {
    if (!w->u.len)
        return 0;
    w->c.len = 0;
    if (buf_add_blocks(&w->c, w->u.p, w->u.len, w->st->tmp_level) < 0
        || write_all(w->fd, w->c.p, w->c.len) < 0)
        return -1;
    w->off += w->c.len;
    w->run.ubytes += w->u.len;
    w->u.len = 0;
    return 0;
}

static int cw_put(cw_t *w, const s2_key_t *k, const s2_rec_t *r) {
    if (w->u.len == 0) {   // new chunk: index its first key
        if (w->run.nidx == w->m_idx) {
            size_t m = w->m_idx ? w->m_idx * 2 : 256;
            s2_idx_t *x = realloc(w->run.idx, m * sizeof(*x));
            if (!x)
                return -1;
            w->run.idx = x;
            w->m_idx = m;
        }
        if (idx_set(w->st, &w->run.idx[w->run.nidx], w->off, w->run.nrec,
                    k) < 0)
            return -1;
        w->run.nidx++;
    }
    size_t fl = frame_len(k, w->st->order, r);
    if (buf_reserve(&w->u, fl) < 0)
        return -1;
    frame_put(w->u.p + w->u.len, k, w->st->order, r);
    w->u.len += fl;
    w->run.nrec++;
    return w->u.len >= S2_CHUNK_BYTES ? cw_flush(w) : 0;
}

// Finish the run; on success it is moved into *out.
static int cw_close(cw_t *w, int ok, s2_run_t *out) {
    if (ok && (cw_flush(w) < 0
               || write_all(w->fd, bgzf_eof_block, sizeof(bgzf_eof_block)) < 0))
        ok = 0;
    if (w->fd >= 0 && close(w->fd) < 0)
        ok = 0;
    free(w->u.p);
    free(w->c.p);
    if (!ok) {
        run_free(&w->run, 1);
        return -1;
    }
    w->st->tmp_cbytes += w->off + sizeof(bgzf_eof_block);
    w->st->tmp_ubytes += w->run.ubytes;
    *out = w->run;
    return 0;
}

/* ---------------------------------------------------------------------- */
/* Merge sources */

typedef struct {
    int src;                   // source index: run number, or n_runs for mem
    BGZF *fp;                  // on-disk run, or
    const s2_run_t *run;
    uint64_t ord;              // ordinal of the next record in the run
    s2_batch_t *mb;            // the in-memory final batch
    size_t mpos, mstart, mend;
    s2_buf_t buf;
    s2_key_t k;
    s2_rec_t *rec;
    const s2_bound_t *upper;   // range end, exclusive (NULL = none)
    uint64_t rbytes;           // raw bytes read from the run
    int done;
} s2_cursor_t;

typedef struct {
    s2_cursor_t *c;
    enum s2_order order;
} merge_ctx_t;

static int cursor_next(const s2_state_t *st, s2_cursor_t *c) {
    if (c->mb) {
        s2_batch_t *b = c->mb;
        if (c->mpos >= c->mend) {
            c->done = 1;
            return 0;
        }
        c->k = batch_key(st, b, c->mpos);
        c->rec = batch_rec(b, b->sorted[c->mpos++].ref);
    } else {
        uint32_t hdr[2];
        if (c->ord >= c->run->nrec) {   // end of this run or segment
            c->done = 1;
            return 0;
        }
        if (bgzf_read(c->fp, hdr, 8) != 8)
            return -1;
        c->ord++;
        c->rbytes += 8 + (uint64_t) hdr[0];
        size_t body = hdr[0], kpad = ALIGN8((size_t) hdr[1]);
        if (body < kpad + REC_DATA_OFF)
            return -1;
        c->buf.len = 0;
        if (buf_reserve(&c->buf, body) < 0)
            return -1;
        if (bgzf_read(c->fp, c->buf.p, body) != (ssize_t) body)
            return -1;
        c->rec = (s2_rec_t *) (c->buf.p + kpad);
        if (kpad + REC_DATA_OFF + c->rec->l_data != body)
            return -1;
        if (order_numeric(st->order)) {
            if (hdr[1] != 0)
                return -1;
            num_key(st->order, &c->rec->core, &c->k.k_hi, &c->k.k64,
                    &c->k.k_x);
        } else {
            c->k.p = c->buf.p;
            c->k.len = hdr[1];
        }
    }
    if (c->upper && bcmp_src(st->order, &c->k, c->src, c->upper) >= 0)
        c->done = 1;
    return 0;
}

/* Position a cursor at its first record at or after lower (NULL = start),
 * stopping before upper. */
static int cursor_seek(const s2_state_t *st, s2_cursor_t *c,
                       const s2_bound_t *lower, const s2_bound_t *upper) {
    c->done = 0;
    c->upper = upper;
    if (c->mb) {
        size_t lo = c->mstart, hi = c->mend;
        while (lower && lo < hi) {   // first position not below lower
            size_t mid = lo + (hi - lo) / 2;
            s2_key_t k = batch_key(st, c->mb, mid);
            if (bcmp_src(st->order, &k, c->src, lower) < 0)
                lo = mid + 1;
            else
                hi = mid;
        }
        c->mpos = lo;
        return cursor_next(st, c);
    }

    const s2_run_t *run = c->run;
    if (run->nidx == 0) {
        c->done = 1;
        return 0;
    }
    size_t lo = 0, hi = run->nidx;
    while (lower && lo < hi) {       // chunks whose first key is below lower
        size_t mid = lo + (hi - lo) / 2;
        if (bcmp_src(st->order, &run->idx[mid].k, c->src, lower) < 0)
            lo = mid + 1;
        else
            hi = mid;
    }
    size_t j = lo ? lo - 1 : 0;
    if (bgzf_seek(c->fp, (int64_t) (run->idx[j].coff << 16), SEEK_SET) < 0)
        return -1;
    c->ord = run->idx[j].ord;
    if (cursor_next(st, c) < 0)
        return -1;
    while (!c->done && lower
           && bcmp_src(st->order, &c->k, c->src, lower) < 0)
        if (cursor_next(st, c) < 0)
            return -1;
    return 0;
}

static int merge_less(void *vp, int i, int j) {
    merge_ctx_t *m = vp;
    const s2_cursor_t *a = &m->c[i], *b = &m->c[j];
    if (a->done) return 0;
    if (b->done) return 1;
    int r = kcmp(m->order, &a->k, &b->k);
    if (r)
        return r < 0;
    return i < j;              // earlier run = earlier input
}

static inline void rec_view(bam1_t *v, s2_rec_t *r);

/* A record as it is written to the output: in minimiser order, unplaced
 * reads lose the minimiser cached in pos/mpos/isize (as samtools sort). */
static inline void rec_view_out(const s2_state_t *st, bam1_t *v,
                                s2_rec_t *r) {
    rec_view(v, r);
    if (st->order == ORDER_MINHASH && v->core.tid == -1) {
        v->core.pos = -1;
        v->core.mpos = -1;
        v->core.isize = 0;
    }
}

static inline void rec_view(bam1_t *v, s2_rec_t *r) {
    v->core = r->core;
    v->l_data = r->l_data;
    v->m_data = r->l_data;
    v->data = (uint8_t *) r + REC_DATA_OFF;
}

static void cursors_close(s2_state_t *st, s2_cursor_t *c, int k) {
    uint64_t rb = 0;
    for (int i = 0; i < k; i++)
        rb += c[i].rbytes;
    pthread_mutex_lock(&st->lock);
    st->tmp_uread += rb;
    pthread_mutex_unlock(&st->lock);
    for (int i = 0; i < k; i++) {
        if (c[i].fp)
            bgzf_close(c[i].fp);
        free(c[i].buf.p);
    }
    free(c);
}

/* Open cursors over runs[0..nf) plus the in-memory batch if given, each
 * positioned at its first record.  Run readers get the htslib pool if pool
 * is set (attached after any seek, so read-ahead starts in the right
 * place). */
static s2_cursor_t *cursors_open(s2_state_t *st, s2_run_t *runs, int nf,
                                 s2_batch_t *mem, hts_tpool *pool) {
    int k = nf + (mem ? 1 : 0);
    s2_cursor_t *c = calloc(k ? k : 1, sizeof(*c));
    if (!c)
        return NULL;
    int qsize = nf > 0 ? 256 / nf : 0;
    if (qsize < 2)
        qsize = 2;
    for (int i = 0; i < nf; i++) {
        c[i].src = i;
        c[i].run = &runs[i];
        c[i].fp = bgzf_open(runs[i].fn, "r");
        char magic[4];
        int ok = c[i].fp && bgzf_read(c[i].fp, magic, 4) == 4
            && memcmp(magic, S2_RUN_MAGIC, 4) == 0;
        // Position at the first chunk, just after the magic block.
        if (ok && runs[i].nidx > 0
            && bgzf_tell(c[i].fp) != (int64_t) (runs[i].idx[0].coff << 16))
            ok = bgzf_seek(c[i].fp, (int64_t) (runs[i].idx[0].coff << 16),
                           SEEK_SET) >= 0;
        if (!ok) {
            print_error_errno("sort2", "couldn't read temporary file \"%s\"",
                              runs[i].fn);
            cursors_close(st, c, k);
            return NULL;
        }
        c[i].ord = 0;
        if (pool)
            bgzf_thread_pool(c[i].fp, pool, qsize);
    }
    if (mem) {
        c[nf].src = nf;
        c[nf].mb = mem;
        c[nf].mstart = 0;
        c[nf].mend = mem->n;
    }
    return c;
}

/* Serial merge of runs[0..nf) (+ mem).  Output goes to out if set,
 * otherwise to a new run returned in *run_out. */
static int merge_serial(s2_state_t *st, s2_run_t *runs, int nf,
                        s2_batch_t *mem, samFile *out, sam_hdr_t *oh,
                        s2_run_t *run_out) {
    int k = nf + (mem ? 1 : 0), ret = -1, lt_ok = 0;
    s2_cursor_t *c = cursors_open(st, runs, nf, mem, st->pool);
    // Every run is open now, so the names are no longer needed.
    if (c)
        for (int i = 0; i < nf; i++)
            run_unlink(&runs[i]);
    s2_ltree_t lt;
    cw_t w;
    int w_open = 0;
    if (!c)
        return -1;
    // cursors_open() positioned every source at its first record.
    for (int i = 0; i < k; i++)
        if (cursor_next(st, &c[i]) < 0)
            goto read_err;
    if (!out) {
        if (cw_open(st, &w) < 0)
            goto out;
        w_open = 1;
    }
    merge_ctx_t mc = { c, st->order };
    if (s2_ltree_init(&lt, k, merge_less, &mc) < 0)
        goto out;
    lt_ok = 1;

    bam1_t view;
    memset(&view, 0, sizeof(view));
    for (;;) {
        s2_cursor_t *ci = &c[s2_ltree_top(&lt)];
        if (ci->done)
            break;
        if (out) {
            rec_view_out(st, &view, ci->rec);
            if (sam_write1(out, oh, &view) < 0) {
                print_error_errno("sort2", "failed writing output");
                goto out;
            }
        } else if (cw_put(&w, &ci->k, ci->rec) < 0) {
            print_error_errno("sort2", "failed writing temporary file");
            goto out;
        }
        if (cursor_next(st, ci) < 0)
            goto read_err;
        s2_ltree_replay(&lt);
    }
    ret = 0;
    goto out;

 read_err:
    print_error("sort2", "error reading temporary file");
 out:
    if (lt_ok)
        s2_ltree_destroy(&lt);
    cursors_close(st, c, k);
    if (w_open && cw_close(&w, ret == 0, run_out) < 0)
        ret = -1;
    return ret;
}

/* Too many runs: merge them all into one, which keeps its place at the
 * front of the (input-ordered) list. */
static int merge_intermediate(s2_state_t *st) {
    double t0 = now();
    s2_run_t merged;
    int n = st->n_runs;
    if (merge_serial(st, st->runs, n, NULL, NULL, NULL, &merged) < 0)
        return -1;
    for (int i = 0; i < n; i++)
        run_free(&st->runs[i], 1);
    st->runs[0] = merged;
    st->n_runs = 1;
    st->n_intermediate++;
    st->t_inter += now() - t0;
    return 0;
}

/* ---------------------------------------------------------------------- */
/* Parallel BAM output */

typedef struct {
    s2_state_t *st;
    s2_run_t *runs;
    s2_batch_t *mem;
    int nf;
    s2_bound_t *split;         // nsplit splitters -> nsplit+1 ranges
    size_t nsplit;
    hFILE *hf;                 // the real output, positioned after the header
    sam_hdr_t *hdr;
    size_t drain_recs;         // head-of-line worker drains this often
    pthread_mutex_t m;
    pthread_cond_t cv;
    size_t next_claim, next_write;
    size_t n_written;
    uint64_t drained;          // bytes copied from scratch to the output
    hts_idx_t *idx;            // --write-index (BAM): built in output order
    void *ixpend;              // the record waiting for its end offset
    int opened, nworkers;      // workers that have opened the runs, of all
    int fail;
} pm_t;

/* Scratch output for one merge worker: an htsFile in the output's format,
 * written with sam_write1() and htslib's own encoders, over a memory file
 * (memfd on Linux, else an unlinked temporary file).  After each range it
 * is flushed, copied to the real output in order, and truncated, so memory
 * holds only finished ranges waiting their turn. */
static int scratch_fd(const s2_state_t *st) {
#if defined(__linux__) && defined(SYS_memfd_create)
    int mfd = (int) syscall(SYS_memfd_create, "samtools-sort2", 0);
    if (mfd >= 0)
        return mfd;
#endif
    char *fn = malloc(strlen(st->tmp_prefix) + 32);
    if (!fn)
        return -1;
    sprintf(fn, "%s.out.XXXXXX", st->tmp_prefix);
    int fd = mkstemp(fn);
    if (fd >= 0)
        unlink(fn);
    free(fn);
    return fd;
}

static int scratch_flush(htsFile *w, hFILE *hf) {
    if (w->is_bgzf && bgzf_flush(w->fp.bgzf) < 0)
        return -1;
    return hflush(hf) < 0 ? -1 : 0;
}

static int scratch_reset(hFILE *hf, int fd) {
    if (hseek(hf, 0, SEEK_SET) < 0)
        return -1;
    return ftruncate(fd, 0) < 0 ? -1 : 0;
}

/* Copy everything in the scratch file to the real output and empty it.
 * Only called by the worker holding the output turn (next_write == its
 * range), so no lock is needed: nobody else writes until it passes the
 * turn on.  The scratch file is mapped, so the bytes are copied once. */
static int scratch_drain(pm_t *pm, htsFile *w, hFILE *shf, int sfd,
                         int64_t *out_base) {
    if (scratch_flush(w, shf) < 0)
        return -1;
    off_t len = htell(shf);
    if (len < 0)
        return -1;
    // Where this range starts in the output, for the index.
    if (*out_base < 0 && (*out_base = htell(pm->hf)) < 0)
        return -1;
    if (len > 0) {
        void *p = mmap(NULL, (size_t) len, PROT_READ, MAP_SHARED, sfd, 0);
        if (p == MAP_FAILED)
            return -1;
        ssize_t wr = hwrite(pm->hf, p, (size_t) len);
        munmap(p, (size_t) len);
        if (wr != (ssize_t) len)
            return -1;
    }
    pthread_mutex_lock(&pm->m);
    pm->drained += (uint64_t) len;
    pthread_mutex_unlock(&pm->m);
    return scratch_reset(shf, sfd);
}

/* Open a worker's scratch htsFile.  The header is written only to put
 * htslib into its normal state, then discarded.  On failure *wp may still
 * be set, for the caller to close. */
static int scratch_open(const s2_state_t *st, sam_hdr_t *hdr, int *fdp,
                        hFILE **hfp, htsFile **wp) {
    int fd = scratch_fd(st);
    if (fd < 0)
        return -1;
    hFILE *hf = hdopen(fd, "w");
    if (!hf) {
        close(fd);
        return -1;
    }
    htsFile *w = hts_hopen(hf, st->out_fn, st->out_mode);
    if (!w) {
        if (hclose(hf) < 0) { /* already failing */ }
        return -1;
    }
    *fdp = fd;
    *hfp = hf;
    *wp = w;
    if (st->out_opts && hts_opt_apply(w, (hts_opt *) st->out_opts) < 0)
        return -1;
    if (sam_hdr_write(w, hdr) < 0 || scratch_flush(w, hf) < 0
        || scratch_reset(hf, fd) < 0)
        return -1;
    return 0;
}

/* One output record's index entry, with virtual offsets in the worker's
 * scratch file (block addresses counted from the scratch htsFile's start). */
typedef struct {
    int32_t tid;
    int mapped;
    hts_pos_t beg, end;
    uint64_t voff_beg, voff_end;
} s2_ixrec_t;

/* Push a finished range's index entries, translating scratch offsets to
 * the output: block address A in the scratch file is (A - a0) bytes after
 * the range's start, which is at byte out_base in the output.  Called by
 * the worker holding the output turn, so pushes are in output order.
 *
 * htslib's own indexed BAM writing pushes each record with the offset after
 * it, then amends that to the next record's start (hts_idx_amend_last(),
 * which isn't public API).  Here each record is held back until the next
 * one's start is known and pushed with that; the last is pushed by
 * index_flush_pending() with its own end. */
typedef struct {
    int have;
    s2_ixrec_t r;              // voff_end already in output terms
} s2_ixpend_t;

static int push_index(hts_idx_t *idx, void *vpend, const s2_ixrec_t *ix,
                      size_t n, uint64_t a0, int64_t out_base) {
    s2_ixpend_t *pend = vpend;
    for (size_t i = 0; i < n; i++) {
        uint64_t b = (((uint64_t) out_base + (ix[i].voff_beg >> 16) - a0) << 16)
                     | (ix[i].voff_beg & 0xffff);
        uint64_t e = (((uint64_t) out_base + (ix[i].voff_end >> 16) - a0) << 16)
                     | (ix[i].voff_end & 0xffff);
        if (pend->have && hts_idx_push(idx, pend->r.tid, pend->r.beg,
                                       pend->r.end, b, pend->r.mapped) < 0)
            return -1;
        pend->r = ix[i];
        pend->r.voff_end = e;
        pend->have = 1;
    }
    return 0;
}

static int index_flush_pending(hts_idx_t *idx, s2_ixpend_t *pend) {
    if (!pend->have)
        return 0;
    pend->have = 0;
    return hts_idx_push(idx, pend->r.tid, pend->r.beg, pend->r.end,
                        pend->r.voff_end, pend->r.mapped);
}

static void pm_worker(void *vp, int t, int nt) {
    (void) t; (void) nt;
    pm_t *pm = vp;
    s2_state_t *st = pm->st;
    int k = pm->nf + (pm->mem ? 1 : 0), lt_ok = 0, err = 0;
    size_t R = pm->nsplit + 1, written = 0;
    s2_ltree_t lt;
    // Each worker decompresses its own runs: no shared htslib pool.
    s2_cursor_t *cur = cursors_open(st, pm->runs, pm->nf, pm->mem, NULL);
    // Once every worker has opened every run, the names can go.  The last
    // worker to get here removes them, outside the lock.
    pthread_mutex_lock(&pm->m);
    int last_opener = ++pm->opened == pm->nworkers && !pm->fail;
    pthread_mutex_unlock(&pm->m);
    if (last_opener)
        for (int i = 0; i < pm->nf; i++)
            run_unlink(&pm->runs[i]);
    merge_ctx_t mc = { cur, st->order };
    int sfd = -1;
    hFILE *shf = NULL;
    htsFile *w = NULL;
    bam1_t view;
    memset(&view, 0, sizeof(view));
    s2_ixrec_t *ix = NULL;     // this range's index entries (reused)
    size_t mix = 0;
    if (!cur) {
        err = 1;
    } else if (scratch_open(st, pm->hdr, &sfd, &shf, &w) < 0) {
        print_error_errno("sort2", "couldn't open scratch output");
        err = 1;
    }

    while (!err) {
        pthread_mutex_lock(&pm->m);
        size_t r = pm->fail ? R : pm->next_claim++;
        pthread_mutex_unlock(&pm->m);
        if (r >= R)
            break;
        const s2_bound_t *lower = r > 0 ? &pm->split[r - 1] : NULL;
        const s2_bound_t *upper = r < R - 1 ? &pm->split[r] : NULL;

        for (int i = 0; i < k && !err; i++)
            if (cursor_seek(st, &cur[i], lower, upper) < 0) {
                print_error("sort2", "error reading temporary file");
                err = 1;
            }
        if (!err) {
            if (lt_ok)
                s2_ltree_destroy(&lt);
            lt_ok = s2_ltree_init(&lt, k, merge_less, &mc) == 0;
            err = !lt_ok;
        }

        // Index bookkeeping for this range: the scratch file's block
        // address at the start, and where the range lands in the output.
        size_t nix = 0;
        uint64_t a0 = w ? (uint64_t) (bgzf_tell(w->fp.bgzf) >> 16) : 0;
        int64_t out_base = -1;
        int head = 0;          // our range is the next to be written
        size_t since = 0;
        while (!err) {
            s2_cursor_t *ci = &cur[s2_ltree_top(&lt)];
            if (ci->done)
                break;
            rec_view_out(st, &view, ci->rec);
            uint64_t vb = 0;
            if (pm->idx) {
                // As bam_write1(): start a new block if the record won't
                // fit, then the record starts here.
                size_t sz = 4 + 32 + view.l_data - view.core.l_extranul
                    + (view.core.n_cigar > 0xffff ? 16 : 0);
                if (bgzf_flush_try(w->fp.bgzf, sz) < 0) {
                    err = 1;
                    break;
                }
                vb = bgzf_tell(w->fp.bgzf);
            }
            if (sam_write1(w, pm->hdr, &view) < 0) {
                print_error_errno("sort2", "failed writing output");
                err = 1;
                break;
            }
            if (pm->idx) {
                if (nix == mix) {
                    size_t m = mix ? mix * 2 : 4096;
                    s2_ixrec_t *x = realloc(ix, m * sizeof(*x));
                    if (!x) {
                        err = 1;
                        break;
                    }
                    ix = x;
                    mix = m;
                }
                s2_ixrec_t *e = &ix[nix++];
                e->tid = view.core.tid;
                e->beg = view.core.pos;
                e->end = bam_endpos(&view);
                e->mapped = !(view.core.flag & BAM_FUNMAP);
                e->voff_beg = vb;
                e->voff_end = bgzf_tell(w->fp.bgzf);
            }
            written++;
            // The worker whose range is next streams straight to the
            // output as it goes, rather than buffering the whole range.
            if (++since >= pm->drain_recs) {
                since = 0;
                if (!head) {
                    pthread_mutex_lock(&pm->m);
                    head = pm->next_write == r;
                    pthread_mutex_unlock(&pm->m);
                }
                if (head && scratch_drain(pm, w, shf, sfd, &out_base) < 0) {
                    print_error_errno("sort2", "failed writing output");
                    err = 1;
                    break;
                }
            }
            if (cursor_next(st, ci) < 0) {
                print_error("sort2", "error reading temporary file");
                err = 1;
                break;
            }
            s2_ltree_replay(&lt);
        }

        // Wait for the output turn, write without holding the lock, then
        // pass the turn on.
        pthread_mutex_lock(&pm->m);
        while (!err && !pm->fail && pm->next_write != r)
            pthread_cond_wait(&pm->cv, &pm->m);
        int go = !err && !pm->fail;
        pthread_mutex_unlock(&pm->m);
        if (go && scratch_drain(pm, w, shf, sfd, &out_base) < 0) {
            print_error_errno("sort2", "failed writing output");
            err = 1;
        }
        // Still holding the turn: index this range, in output order.
        if (go && !err && pm->idx && nix
            && push_index(pm->idx, pm->ixpend, ix, nix, a0, out_base) < 0) {
            print_error("sort2", "failed building the index");
            err = 1;
        }
        pthread_mutex_lock(&pm->m);
        if (err)
            pm->fail = 1;
        else if (go)
            pm->next_write++;
        pthread_cond_broadcast(&pm->cv);
        pthread_mutex_unlock(&pm->m);
        if (!go)
            break;
    }

    pthread_mutex_lock(&pm->m);
    if (err)
        pm->fail = 1;
    pm->n_written += written;
    pthread_cond_broadcast(&pm->cv);
    pthread_mutex_unlock(&pm->m);
    if (lt_ok)
        s2_ltree_destroy(&lt);
    if (cur)
        cursors_close(st, cur, k);
    // Closing writes an EOF marker into the scratch file, which is dropped.
    if (w)
        hts_close(w);
    else if (shf)
        if (hclose(shf) < 0) { /* nothing more to do */ }
    free(ix);
}

static enum s2_order g_split_order;   // for qsort; set before each use
static int bound_qcmp(const void *a, const void *b) {
    const s2_bound_t *x = a, *y = b;
    int c = kcmp(g_split_order, &x->k, &y->k);
    return c ? c : (x->src > y->src) - (x->src < y->src);
}

/* Choose range splitters at quantiles of all run index keys plus samples
 * of the in-memory batch.  Returns the number of splitters (may be 0). */
static size_t choose_splitters(s2_state_t *st, s2_run_t *runs, int nf,
                               s2_batch_t *mem, int nt, s2_bound_t **out) {
    size_t npool = 0, mstep = 0;
    uint64_t total = 0;
    for (int i = 0; i < nf; i++) {
        npool += runs[i].nidx;
        total += runs[i].ubytes;
    }
    if (mem && mem->n) {
        mstep = chunk_records(st, mem);
        npool += (mem->n + mstep - 1) / mstep;
        total += mem->used;
    }
    *out = NULL;
    if (npool < 2)
        return 0;
    s2_bound_t *pool = malloc(npool * sizeof(*pool));
    if (!pool)
        return 0;
    size_t n = 0;
    for (int i = 0; i < nf; i++)
        for (size_t j = 0; j < runs[i].nidx; j++)
            pool[n++] = (s2_bound_t){ runs[i].idx[j].k, i };
    if (mstep)
        for (size_t j = 0; j < mem->n; j += mstep)
            pool[n++] = (s2_bound_t){ batch_key(st, mem, j), nf };
    g_split_order = st->order;
    qsort(pool, n, sizeof(*pool), bound_qcmp);

    /* Range size.  Starting a range costs each run about one chunk of
     * decompression and skipping, so ranges must be large compared with
     * nf * S2_CHUNK_BYTES; aim for under ~10% overhead.  Within that,
     * prefer smaller ranges (about 1/4 of -m buffered across nt workers)
     * for load balance and less buffered output. */
    uint64_t target = st->mem_total / 4 / nt;
    uint64_t seek_min = (uint64_t) nf * S2_CHUNK_BYTES * 10;
    if (target < seek_min)
        target = seek_min;
    if (target < ((uint64_t) 8 << 20))
        target = (uint64_t) 8 << 20;
    size_t R = total / target + 1;
    /* Whole rounds of ranges.  With a part-filled last round, most workers
     * sit idle while the last few ranges finish (e.g. 65 ranges on 32
     * workers costs three rounds of time for two rounds of work).  At least
     * four rounds, so the tail is small even when ranges vary in cost; the
     * ranges just get smaller, which also reduces buffered output. */
    size_t rounds = (R + nt - 1) / nt;
    if (rounds < 4)
        rounds = 4;
    R = rounds * (size_t) nt;
    if (R > n)
        R = n;

    size_t ns = 0;
    for (size_t j = 1; j < R; j++) {
        const s2_bound_t *b = &pool[j * n / R];
        if (ns == 0 || bound_qcmp(b, &pool[ns - 1]) > 0)
            pool[ns++] = *b;   // compacts in place: j*n/R >= ns always
    }
    if (ns == 0) {
        free(pool);
        return 0;
    }
    *out = pool;
    return ns;
}

/* Merge everything into a BAM file, in parallel.  out must have had its
 * header written and must have no thread pool. */
static int merge_parallel(s2_state_t *st, s2_run_t *runs, int nf,
                          s2_batch_t *mem, samFile *out, sam_hdr_t *hdr,
                          int nt, size_t expect) {
    pm_t pm;
    memset(&pm, 0, sizeof(pm));
    pm.st = st;
    pm.runs = runs;
    pm.mem = mem;
    pm.nf = nf;
    pm.nsplit = choose_splitters(st, runs, nf, mem, nt, &pm.split);
    pm.hdr = hdr;
    pm.drain_recs = st->drain_recs;
    if (out->is_bgzf) {
        // Push the header out so raw blocks can follow it.
        if (bgzf_flush(out->fp.bgzf) < 0) {
            free(pm.split);
            return -1;
        }
        pm.hf = out->fp.bgzf->fp;
    } else {
        pm.hf = out->fp.hfile;  // uncompressed SAM
    }
    pthread_mutex_init(&pm.m, NULL);
    pthread_cond_init(&pm.cv, NULL);
    size_t R = pm.nsplit + 1;
    if (st->idx_fn) {
        // As sam_idx_init() for CSI: levels to cover the longest reference.
        int min_shift = 14, n_lvls = 0;
        int64_t max_len = 0, sz;
        for (int i = 0; i < sam_hdr_nref(hdr); i++)
            if (sam_hdr_tid2len(hdr, i) > max_len)
                max_len = sam_hdr_tid2len(hdr, i);
        max_len += 256;
        for (sz = (int64_t) 1 << min_shift; max_len > sz; n_lvls++, sz <<= 3)
            ;
        off_t hdr_end = htell(pm.hf);   // header already flushed
        if (hdr_end < 0 || !(pm.idx = hts_idx_init(sam_hdr_nref(hdr),
                                                    HTS_FMT_CSI,
                                                    (uint64_t) hdr_end << 16,
                                                    min_shift, n_lvls))) {
            print_error("sort2", "couldn't initialise the index");
            free(pm.split);
            return -1;
        }
    }
    s2_ixpend_t ixpend = { 0 };
    pm.ixpend = &ixpend;
    pm.nworkers = R < (size_t) nt ? (int) R : nt;
    s2_parallel(pm.nworkers, pm_worker, &pm);
    if (pm.idx) {
        off_t end = htell(pm.hf);       // end of the records, before EOF
        if (!pm.fail && index_flush_pending(pm.idx, &ixpend) < 0)
            pm.fail = 1;
        if (!pm.fail && (end < 0
                         || hts_idx_finish(pm.idx, (uint64_t) end << 16) < 0)) {
            print_error("sort2", "failed finishing the index");
            pm.fail = 1;
        }
        if (pm.fail)
            hts_idx_destroy(pm.idx);
        else
            st->idx = pm.idx;         // saved once the output is closed
    }
    pthread_mutex_destroy(&pm.m);
    pthread_cond_destroy(&pm.cv);
    free(pm.split);
    st->out_ranges = R;
    st->scratch_bytes += pm.drained;
    st->out_threads = R < (size_t) nt ? (int) R : nt;
    if (!pm.fail && pm.n_written != expect) {
        print_error("sort2", "internal error: wrote %zu of %zu records",
                    pm.n_written, expect);
        return -1;
    }
    return pm.fail ? -1 : 0;
}

static int output_batch_serial(const s2_state_t *st, s2_batch_t *b,
                               samFile *out, sam_hdr_t *oh) {
    bam1_t view;
    memset(&view, 0, sizeof(view));
    for (size_t i = 0; i < b->n; i++) {
        if (i + 16 < b->n)
            S2_PREFETCH(b->arena + b->off[b->sorted[i + 16].ref]);
        rec_view_out(st, &view, batch_rec(b, b->sorted[i].ref));
        if (sam_write1(out, oh, &view) < 0) {
            print_error_errno("sort2", "failed writing output");
            return -1;
        }
    }
    return 0;
}

/* ---------------------------------------------------------------------- */
/* Pipeline */

static s2_batch_t *get_free_batch(s2_state_t *st) {
    double t0 = now();
    pthread_mutex_lock(&st->lock);
    while (!st->free_mask && !st->error)
        pthread_cond_wait(&st->cond, &st->lock);
    s2_batch_t *b = NULL;
    if (!st->error) {
        int i = (st->free_mask & 1) ? 0 : 1;
        st->free_mask &= ~(1 << i);
        b = &st->batch[i];
    }
    pthread_mutex_unlock(&st->lock);
    st->t_read_stall += now() - t0;
    if (b)
        batch_reset(b);
    return b;
}

static int hand_off(s2_state_t *st, s2_batch_t *b) {
    double t0 = now();
    pthread_mutex_lock(&st->lock);
    while (st->full && !st->error)
        pthread_cond_wait(&st->cond, &st->lock);
    int err = st->error;
    if (!err) {
        st->full = b;
        pthread_cond_broadcast(&st->cond);
    }
    pthread_mutex_unlock(&st->lock);
    st->t_read_stall += now() - t0;
    return err ? -1 : 0;
}

static void *sorter_thread(void *vp) {
    s2_state_t *st = vp;
    for (;;) {
        double t0 = now();
        pthread_mutex_lock(&st->lock);
        while (!st->full && !st->error)
            pthread_cond_wait(&st->cond, &st->lock);
        s2_batch_t *b = st->full;
        st->full = NULL;
        int err = st->error;
        pthread_cond_broadcast(&st->cond);
        pthread_mutex_unlock(&st->lock);
        st->t_sort_idle += now() - t0;
        if (err || !b)
            return NULL;

        st->n_batches++;
        if (batch_sort(st, b) < 0) {
            set_error(st);
            return NULL;
        }
        if (b->eof) {
            st->final = b;
            st->phase1_done = 1;
            return NULL;
        }
        if (write_run(st, b) < 0) {
            set_error(st);
            return NULL;
        }
        if (st->n_runs >= st->max_runs && merge_intermediate(st) < 0) {
            set_error(st);
            return NULL;
        }

        pthread_mutex_lock(&st->lock);
        st->free_mask |= 1 << (int) (b - st->batch);
        pthread_cond_broadcast(&st->cond);
        pthread_mutex_unlock(&st->lock);
    }
}

/* ---------------------------------------------------------------------- */
/* Command line */

static size_t parse_mem(const char *s) {
    char *end;
    errno = 0;
    double v = strtod(s, &end);
    if (errno || end == s || v < 0)
        return 0;
    switch (*end) {
    case 'k': case 'K': v *= 1024; end++; break;
    case 'm': case 'M': v *= 1024 * 1024; end++; break;
    case 'g': case 'G': v *= 1024.0 * 1024 * 1024; end++; break;
    }
    if (*end)
        return 0;
    return (size_t) v;
}

static void sort2_usage(FILE *fp) {
    fprintf(fp,
"Usage: samtools sort2 [options...] [in.bam]\n"
"\n"
"Experimental pipelined, cache-aware sort.  Output order is intended to\n"
"match \"samtools sort\" for the modes it supports.\n"
"\n"
"Options:\n"
"  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)\n"
"  -u         Output uncompressed data (equivalent to -l 0)\n"
"  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]\n"
"  -n         Sort by read name (natural order)\n"
"  -N         Sort by read name (byte order)\n"
"  -t TAG     Sort by value of TAG. Uses position as secondary index (or read\n"
"             name if -n or -N is set)\n"
"  -M         Use minimiser for clustering unaligned/unplaced reads\n"
"  -R         Do not use reverse strand (only compatible with -M)\n"
"  -K INT     Kmer size to use for minimiser [20]\n"
"  -I FILE    Order minimisers by their position in FILE FASTA\n"
"  -w INT     Window size for minimiser indexing via -I ref.fa [100]\n"
"  -H         Squash homopolymers when computing minimiser\n"
"  -o FILE    Write final output to FILE rather than standard output\n"
"  -T PREFIX  Write temporary files to PREFIX.nnnn.s2r\n"
"  --tmp-level INT\n"
"             Compression level for temporary files [1]\n"
"  --max-tmp-files INT\n"
"             With plain runs, merge them once there are this many [auto]\n"
"  --serial-merge\n"
"             Merge on one thread through sam_write1(), with htslib's\n"
"             threaded output (as samtools sort does)\n"

"  --stats    Report per-phase timings on stderr\n"
"  --no-PG    Do not add a PG line\n"
"  --template-coordinate\n"
"             Sort by template-coordinate\n");
    sam_global_opt_help(fp, "-.O..@-.");
    fprintf(fp,
"\n"
"SAM and BAM output are merged in parallel; CRAM uses a serial merge.\n"
"--write-index writes a CSI index (CRAI for CRAM), as samtools sort does.\n"
"All samtools sort orders are supported.  Positions beyond 32 bits are\n"
"supported, but as with \"samtools sort\" only SAM/CRAM output can hold them.\n");
}

static char *default_tmp_prefix(const char *fnout, const char *opt) {
    char *p;
    struct stat sb;
    if (opt && stat(opt, &sb) == 0 && S_ISDIR(sb.st_mode)) {
        p = malloc(strlen(opt) + 64);
        if (p)
            sprintf(p, "%s/samtools.%d.sort2.tmp", opt, (int) getpid());
    } else if (opt) {
        p = strdup(opt);
    } else if (strcmp(fnout, "-") != 0) {
        p = malloc(strlen(fnout) + 16);
        if (p)
            sprintf(p, "%s.s2tmp", fnout);
    } else {
        p = malloc(64);
        if (p)
            sprintf(p, "samtools.%d.sort2.tmp", (int) getpid());
    }
    return p;
}

/* Descriptors available for run files: the soft RLIMIT_NOFILE less a
 * reserve for stdio, input, output, reference and so on.  If fewer than
 * want are available, the soft limit is raised towards the hard limit. */
#define S2_FD_RESERVE 64
static long fd_budget(long want) {
    struct rlimit rl;
    if (getrlimit(RLIMIT_NOFILE, &rl) != 0)
        return 1024 - S2_FD_RESERVE;
    rlim_t need = (rlim_t) want + S2_FD_RESERVE;
    if (rl.rlim_cur != RLIM_INFINITY && need > rl.rlim_cur
        && (rl.rlim_max == RLIM_INFINITY || rl.rlim_cur < rl.rlim_max)) {
        struct rlimit nr = rl;
        nr.rlim_cur = (rl.rlim_max != RLIM_INFINITY && need > rl.rlim_max)
            ? rl.rlim_max : need;
        if (setrlimit(RLIMIT_NOFILE, &nr) == 0)
            rl = nr;
    }
    if (rl.rlim_cur == RLIM_INFINITY)
        return want;
    long avail = (long) rl.rlim_cur - S2_FD_RESERVE;
    return avail > 0 ? avail : 0;
}

/* Merge runs early once there are this many.  Only the serial merge opens
 * every run at once here, so the cap is about memory for per-run read
 * buffers rather than descriptors. */
static int auto_max_runs(void) {
    long lim = fd_budget(512);
    if (lim > 512)
        lim = 512;
    if (lim < 2)
        lim = 2;
    return (int) lim;
}

/* @HD tags as samtools sort sets them: SO, and GO and SS where the order
 * has them (removed otherwise). */
/* @HD tags as samtools develop's set_sort_order(): SO, and GO and SS where
 * the order has them (removed otherwise). */
static int set_header_order(sam_hdr_t *h, const s2_state_t *st) {
    const char *so, *go = NULL, *ss = NULL;
    char sstag[64];
    switch (st->order) {
    case ORDER_COORD:     so = "coordinate"; break;
    case ORDER_NAME_NAT:  so = "queryname"; ss = "queryname:natural"; break;
    case ORDER_NAME_BYTE: so = "queryname"; ss = "queryname:lexicographical"; break;
    case ORDER_MINHASH:
        so = st->placed ? "coordinate" : "unsorted";
        ss = st->placed ? "coordinate:minhash" : "unsorted:minhash";
        break;
    case ORDER_TAG_COORD:
    case ORDER_TAG_NAME:
        so = "unsorted";
        snprintf(sstag, sizeof(sstag), "unsorted:%c%c:%s",
                 st->sort_tag[0], st->sort_tag[1],
                 st->order == ORDER_TAG_COORD ? "coordinate"
                 : st->tag_natural ? "queryname:natural"
                 : "queryname:lexicographical");
        ss = sstag;
        break;
    case ORDER_TEMPLATE:
        so = "unsorted"; go = "query"; ss = "unsorted:template-coordinate";
        break;
    default:              so = "unknown"; break;
    }
    if (sam_hdr_count_lines(h, "HD") <= 0) {
        if (sam_hdr_add_line(h, "HD", "VN", SAM_FORMAT_VERSION, "SO", so,
                             NULL) < 0)
            return -1;
    } else if (sam_hdr_update_hd(h, "SO", so) < 0) {
        return -1;
    }
    if (go ? sam_hdr_update_hd(h, "GO", go) < 0
           : sam_hdr_remove_tag_hd(h, "GO") < 0)
        return -1;
    if (ss ? sam_hdr_update_hd(h, "SS", ss) < 0
           : sam_hdr_remove_tag_hd(h, "SS") < 0)
        return -1;
    return 0;
}

int bam_sort2(int argc, char *argv[]) {
    enum {
        OPT_TMP_LEVEL = 1000, OPT_MAX_TMP, OPT_STATS, OPT_NO_PG, OPT_SERIAL_MERGE,
        OPT_TEMPLATE
    };
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        {"tmp-level", required_argument, NULL, OPT_TMP_LEVEL},
        {"max-tmp-files", required_argument, NULL, OPT_MAX_TMP},
        {"stats", no_argument, NULL, OPT_STATS},
        {"no-PG", no_argument, NULL, OPT_NO_PG},
        {"serial-merge", no_argument, NULL, OPT_SERIAL_MERGE},
        {"template-coordinate", no_argument, NULL, OPT_TEMPLATE},
        {NULL, 0, NULL, 0}
    };
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    s2_state_t st;
    memset(&st, 0, sizeof(st));
    st.order = ORDER_COORD;
    st.tmp_level = 1;
    st.kmer = 20;
    st.try_rev = 1;
    st.no_squash = 1;
    st.window = 100;
    const char *minimiser_ref = NULL;
    const char *sort_tag = NULL;
    int level = -1, no_pg = 0, c, ret = EXIT_FAILURE, sorter_started = 0;
    size_t mem = 0;
    char *fnout = "-", *tmp_opt = NULL, *arg_list = NULL;
    samFile *in = NULL, *out = NULL;
    sam_hdr_t *hin = NULL, *hout = NULL;
    bam1_t *b = NULL;
    pthread_t sorter;
    int lock_init = 0;
    double t_start = now(), t_closed = 0;
    s2_usage_t u0, u1, u2;
    usage_now(&u0);

    while ((c = getopt_long(argc, argv, "l:m:nNo:O:T:@:ut:MK:I:Rw:H", lopts,
                            NULL)) >= 0) {
        switch (c) {
        case 'l': level = atoi(optarg); break;
        case 'u': level = 0; break;
        case 'm':
            if (!(mem = parse_mem(optarg))) {
                print_error("sort2", "invalid -m value \"%s\"", optarg);
                goto end;
            }
            break;
        case 'n': st.order = ORDER_NAME_NAT; break;
        case 'N': st.order = ORDER_NAME_BYTE; break;
        // As samtools sort, the last of -n, -M and --template-coordinate
        // wins, and -t then makes it a tag sort.
        case 'M': st.order = ORDER_MINHASH; break;
        case 'I': st.order = ORDER_MINHASH; minimiser_ref = optarg; break;
        case 'R': st.try_rev = 0; break;
        case 'H': st.no_squash = 0; break;
        case 'w': {
            char *e;
            long w = strtol(optarg, &e, 10);
            if (*e || e == optarg || w < 1 || w > INT_MAX) {
                print_error("sort2", "invalid window \"%s\"", optarg);
                goto end;
            }
            st.window = (int) w;
            break;
        }
        case OPT_TEMPLATE: st.order = ORDER_TEMPLATE; break;
        case 't': sort_tag = optarg; break;
        case 'K':
            st.kmer = atoi(optarg);
            if (st.kmer < 1)
                st.kmer = 1;
            else if (st.kmer > 31)
                st.kmer = 31;
            break;
        case 'o': fnout = optarg; break;
        case 'T': tmp_opt = optarg; break;
        case OPT_TMP_LEVEL: st.tmp_level = atoi(optarg); break;
        case OPT_MAX_TMP: st.max_runs = atoi(optarg); break;
        case OPT_STATS: st.stats = 1; break;
        case OPT_NO_PG: no_pg = 1; break;
        case OPT_SERIAL_MERGE: st.serial_merge = 1; break;
        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0)
                break;
            /* else fall through */
        case '?':
            sort2_usage(stderr);
            goto end;
        }
    }
    if (minimiser_ref) {
        fprintf(stderr, "Building index ... ");
        if (s2_build_minhash_index(&st, minimiser_ref)) {
            fprintf(stderr, "\n");
            print_error("sort2", "couldn't index \"%s\"", minimiser_ref);
            goto end;
        }
        fprintf(stderr, "done\n");
    }
    if (sort_tag) {
        if (!sort_tag[0] || !sort_tag[1]) {
            print_error("sort2", "-t needs a two-character tag");
            goto end;
        }
        st.sort_tag[0] = sort_tag[0];
        st.sort_tag[1] = sort_tag[1];
        st.tag_natural = st.order != ORDER_NAME_BYTE;
        st.order = st.order == ORDER_NAME_NAT || st.order == ORDER_NAME_BYTE
                   ? ORDER_TAG_NAME : ORDER_TAG_COORD;
    }
    if (argc - optind > 1) {
        sort2_usage(stderr);
        goto end;
    }
    if (argc == 1 && isatty(STDIN_FILENO)) {
        sort2_usage(stdout);
        ret = EXIT_SUCCESS;
        goto end;
    }
    const char *fnin = optind < argc ? argv[optind] : "-";

    st.nthreads = ga.nthreads > 0 ? ga.nthreads : 1;
    // -m is per thread, as in "samtools sort".
    mem = (mem ? mem : S2_DEFAULT_MEM_PER_THREAD) * st.nthreads;
    if (mem < ((size_t)2 << 20))
        mem = (size_t)2 << 20;
    st.mem_total = mem;
    if (ga.write_index) {
        if (!order_numeric(st.order)) {
            // As samtools sort: a warning, not an error.
            fprintf(stderr, "[W::bam_sort2] Ignoring --write-index as it "
                    "only works for position sorted files.\n");
            ga.write_index = 0;
        }
        if (ga.write_index && strcmp(fnout, "-") == 0) {
            print_error("sort2", "--write-index needs an output file (-o)");
            goto end;
        }
    }
    st.drain_recs = 16384;
    if (getenv("SORT2_TEST_DRAIN"))   // test hook: exercise streaming
        st.drain_recs = (size_t) atol(getenv("SORT2_TEST_DRAIN")) + 1;
    st.batch_limit = mem / 2;
    if (st.max_runs < 2)
        st.max_runs = auto_max_runs();
    st.tmp_prefix = default_tmp_prefix(fnout, tmp_opt);
    if (!st.tmp_prefix)
        goto end;

    if (!no_pg && !(arg_list = stringify_argv(argc + 1, argv - 1))) {
        print_error("sort2", "failed to create arg_list");
        goto end;
    }

    if (ga.nthreads > 0 && !(st.pool = hts_tpool_init(ga.nthreads))) {
        print_error("sort2", "couldn't create thread pool");
        goto end;
    }
    htsThreadPool tp = { st.pool, 0 };

    in = sam_open_format(fnin, "r", &ga.in);
    if (!in) {
        print_error_errno("sort2", "can't open \"%s\"", fnin);
        goto end;
    }
    if (st.pool)
        hts_set_opt(in, HTS_OPT_THREAD_POOL, &tp);
    if (!(hin = sam_hdr_read(in))) {
        print_error("sort2", "failed to read header from \"%s\"", fnin);
        goto end;
    }
    // Library per read group, needed by the template-coordinate keys.
    if (st.order == ORDER_TEMPLATE && !(st.libs = template_libraries(hin))) {
        print_error("sort2", "couldn't read the read groups' libraries");
        goto end;
    }

    if (pthread_mutex_init(&st.lock, NULL) != 0
        || pthread_cond_init(&st.cond, NULL) != 0)
        goto end;
    lock_init = 1;
    for (int i = 0; i < 2; i++) {
        st.batch[i].arena = big_alloc(st.batch_limit);
        st.batch[i].cap = st.batch[i].arena ? st.batch_limit : 0;
    }
    st.free_mask = 3;
    if (pthread_create(&sorter, NULL, sorter_thread, &st) != 0) {
        print_error("sort2", "couldn't start sorter thread");
        goto end;
    }
    sorter_started = 1;

    // ---- Phase 1: read, with sorting/run writing in the background
    if (!(b = bam_init1()))
        goto end_pipeline;
    s2_batch_t *cur = get_free_batch(&st);
    int r = 0;
    while (cur && (r = sam_read1(in, hin, b)) >= 0) {
        st.placed |= b->core.tid >= 0;
        int a = batch_add(&st, cur, b);
        if (a == 0) {
            if (hand_off(&st, cur) < 0) {
                cur = NULL;
                break;
            }
            cur = get_free_batch(&st);
            a = cur ? batch_add(&st, cur, b) : -1;
        }
        if (a < 0) {
            if (cur)
                print_error("sort2", "couldn't allocate memory for records");
            set_error(&st);
            cur = NULL;
            break;
        }
        st.n_records++;
    }
    if (cur && r < -1) {
        print_error("sort2", "truncated or corrupt input file \"%s\"", fnin);
        set_error(&st);
        cur = NULL;
    }
    if (cur) {
        cur->eof = 1;
        hand_off(&st, cur);
    }

 end_pipeline:
    if (!b)
        set_error(&st);
    pthread_join(sorter, NULL);
    sorter_started = 0;
    usage_now(&u1);
    double t_p1 = now() - t_start;
    if (st.error || !st.phase1_done)
        goto end;

    // ---- Phase 2: write output
    double t0 = now();
    if (!(hout = sam_hdr_dup(hin)) || set_header_order(hout, &st) < 0) {
        print_error("sort2", "failed to update header");
        goto end;
    }
    if (!no_pg && sam_hdr_add_pg(hout, "samtools", "VN", samtools_version(),
                                 arg_list ? "CL" : NULL, arg_list, NULL) < 0) {
        print_error("sort2", "failed to add PG line");
        goto end;
    }
    char modeout[12];
    strcpy(modeout, "wb");
    sam_open_mode(modeout + 1, fnout, NULL);
    if (level >= 0)
        sprintf(strchr(modeout, '\0'), "%d", level > 9 ? 9 : level);
    if (!(out = sam_open_format(fnout, modeout, &ga.out))) {
        print_error_errno("sort2", "can't open \"%s\"", fnout);
        goto end;
    }

    // Parallel merge needs BGZF-compressed BAM, records BAM can hold
    // without special handling, a little-endian host, and enough file
    // descriptors for every worker to open every run.
    int nt_out = st.nthreads;
    if (st.n_runs > 0) {
        // Every merge worker opens every run.
        long by_fd = fd_budget((long) st.nthreads * st.n_runs) / st.n_runs;
        if (by_fd < nt_out)
            nt_out = (int) by_fd;
    }
    // htslib opens "wb" as binary_format and only turns it into bam when
    // the header is written, so accept either (likewise text_format/sam).
    const htsFormat *of = hts_get_format(out);
    int is_bam = (of->format == bam || of->format == binary_format)
        && of->compression == bgzf;
    int is_sam = of->format == sam || of->format == text_format;
    st.merge_kind = MERGE_HTS;
    if (st.serial_merge)
        st.serial_why = "--serial-merge";
    else if (st.nthreads < 2)
        st.serial_why = "fewer than 2 threads";
    else if (nt_out < 2)
        st.serial_why = "too few file descriptors for 2 merge workers";
    else if (!is_bam && !is_sam)
        st.serial_why = "the output format (only SAM and BAM merge in parallel)";
    else if (ga.write_index && !is_bam)
        st.serial_why = "--write-index with non-BAM output";
    if (st.serial_why)
        st.merge_kind = MERGE_SERIAL;
    int parallel_out = st.merge_kind != MERGE_SERIAL;
    st.out_fn = fnout;
    st.out_mode = modeout;
    st.out_opts = ga.out.specific;

    if (st.pool && !parallel_out)
        hts_set_opt(out, HTS_OPT_THREAD_POOL, &tp);
    if (ga.write_index) {
        // As samtools sort --write-index: CSI (min_shift 14), or CRAI.
        const char *ext = hts_get_format(out)->format == cram ? ".crai" : ".csi";
        if (!(st.idx_fn = malloc(strlen(fnout) + 8)))
            goto end;
        sprintf(st.idx_fn, "%s%s", fnout, ext);
    }
    if (sam_hdr_write(out, hout) < 0) {
        print_error_errno("sort2", "failed writing header to \"%s\"", fnout);
        goto end;
    }
    if (st.idx_fn && !parallel_out
        && sam_idx_init(out, hout, 14, st.idx_fn) < 0) {
        print_error_errno("sort2", "couldn't start the index \"%s\"",
                          st.idx_fn);
        goto end;
    }

    st.out_threads = parallel_out ? nt_out : 1;
    if (parallel_out) {
        if (merge_parallel(&st, st.runs, st.n_runs, st.final, out, hout,
                           nt_out, st.n_records) < 0)
            goto end;
    } else if (st.n_runs == 0) {
        // Everything fitted in one batch: no temporary files at all.
        if (output_batch_serial(&st, st.final, out, hout) < 0)
            goto end;
    } else if (merge_serial(&st, st.runs, st.n_runs, st.final, out, hout,
                            NULL) < 0) {
        goto end;
    }
    if (st.idx_fn && !parallel_out && sam_idx_save(out) < 0) {
        print_error_errno("sort2", "couldn't write the index \"%s\"",
                          st.idx_fn);
        goto end;
    }
    double tc0 = now();
    int cret = sam_close(out);
    t_closed = now();
    st.t_close = t_closed - tc0;
    out = NULL;
    if (cret < 0) {
        print_error_errno("sort2", "error closing \"%s\"", fnout);
        goto end;
    }
    if (st.idx) {
        double ti = now();
        int r = hts_idx_save_as(st.idx, fnout, st.idx_fn, HTS_FMT_CSI);
        hts_idx_destroy(st.idx);
        st.idx = NULL;
        st.t_idx_save = now() - ti;
        if (r < 0) {
            print_error_errno("sort2", "couldn't write the index \"%s\"",
                              st.idx_fn);
            goto end;
        }
    }
    st.t_final = now() - t0;
    usage_now(&u2);
    ret = EXIT_SUCCESS;

    if (st.stats) {
        static const char *kind[] = {
            "serial merge", "parallel merge via htsFile buffers"
        };
        fprintf(stderr,
                "[sort2] records %zu, batches %zu\n"
                "[sort2] reader stalled on sorter %.2fs; sorter idle %.2fs\n"
                "[sort2] keys %.2fs, sort %.2fs, temp writes %.2fs\n",
                st.n_records, st.n_batches, st.t_read_stall, st.t_sort_idle,
                st.t_keys, st.t_sort, st.t_runs);
        if (st.run_seq)
            fprintf(stderr, "[sort2] runs: %d written, %zu intermediate "
                    "merges (%.2fs); deleting temp files took %.2fs\n",
                    st.run_seq, st.n_intermediate, st.t_inter, g_t_unlink);
        else
            fprintf(stderr, "[sort2] runs: none, sorted in memory\n");
        const double MB = 1048576.0;
        fprintf(stderr, "[sort2] temp: wrote %.0f MB (%.0f MB raw); merges "
                "read back %.0f MB raw (%.2fx)\n",
                st.tmp_cbytes / MB, st.tmp_ubytes / MB, st.tmp_uread / MB,
                st.tmp_ubytes ? (double) st.tmp_uread / st.tmp_ubytes : 0.0);
        const char *pn[2] = { "phase 1 (read, sort, temp writes)",
                              "phase 2 (merge, output)" };
        s2_usage_t *ua[3] = { &u0, &u1, &u2 };
        double wall[2] = { t_p1, st.t_final };
        for (int p = 0; p < 2; p++) {
            fprintf(stderr, "[sort2] %s: %.2fs wall, %.2fs user, %.2fs sys",
                    pn[p], wall[p], ua[p + 1]->user - ua[p]->user,
                    ua[p + 1]->sys - ua[p]->sys);
            if (ua[p]->rchar >= 0 && ua[p + 1]->rchar >= 0)
                fprintf(stderr, ", syscall read %.0f MB, write %.0f MB",
                        (ua[p + 1]->rchar - ua[p]->rchar) / MB,
                        (ua[p + 1]->wchar - ua[p]->wchar) / MB);
            fputc('\n', stderr);
        }
        if (st.scratch_bytes)
            fprintf(stderr, "[sort2] phase 2 syscall I/O includes %.0f MB "
                    "written to and read back from scratch files\n",
                    st.scratch_bytes / MB);
        fprintf(stderr, "[sort2] final output %.2fs: %s%s%s, %d thread%s",
                st.t_final, kind[st.merge_kind],
                st.serial_why ? " because of " : "",
                st.serial_why ? st.serial_why : "",
                st.out_threads, st.out_threads == 1 ? "" : "s");
        if (!st.serial_why && st.out_ranges)
            fprintf(stderr, ", %zu ranges", st.out_ranges);
        fprintf(stderr, "\n[sort2] output close (incl. flush/sync) %.2fs",
                st.t_close);
        if (st.idx_fn)
            fprintf(stderr, "; index %s built %s", st.idx_fn,
                    st.serial_why ? "by htslib" : "during the merge");
        fprintf(stderr, "\n[sort2] total %.2fs\n", now() - t_start);
    }

 end:
    if (sorter_started) {
        set_error(&st);
        pthread_join(sorter, NULL);
    }
    if (out)
        sam_close(out);
    for (int i = 0; i < st.n_runs; i++)
        run_free(&st.runs[i], 1);
    free(st.runs);
    for (int i = 0; i < 2; i++)
        batch_free(&st.batch[i], st.nthreads);
    if (lock_init) {
        pthread_mutex_destroy(&st.lock);
        pthread_cond_destroy(&st.cond);
    }
    if (b)
        bam_destroy1(b);
    if (hin)
        sam_hdr_destroy(hin);
    if (hout)
        sam_hdr_destroy(hout);
    if (in)
        sam_close(in);
    if (st.pool)
        hts_tpool_destroy(st.pool);
    free(st.tmp_prefix);
    free(st.idx_fn);
    template_libraries_free(st.libs);
    s2_minhash_index_free(st.kmer_h);
    if (st.idx)
        hts_idx_destroy(st.idx);
    free(arg_list);
    sam_global_args_free(&ga);
    if (st.stats && ret == EXIT_SUCCESS && t_closed > 0)
        fprintf(stderr, "[sort2] cleanup after closing the output: %.2fs "
                "(input close, frees, thread pool shutdown)\n",
                now() - t_closed);
    return ret;
}
