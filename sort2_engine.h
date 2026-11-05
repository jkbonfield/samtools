/*  sort2_engine.h -- htslib-independent sorting engine for samtools sort2.

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
 * Design summary
 * --------------
 * Every sort order is reduced to a byte-comparable key, so a single engine
 * serves all modes:
 *
 *  - Coordinate order: one 64-bit integer, (tid' << 33) | ((pos+1) << 1) | rev
 *    where unmapped (tid -1) maps to the largest tid.  Sorted with a parallel
 *    stable LSD radix sort that skips constant digits.
 *
 *  - Name order: the query name is encoded so that memcmp() order equals
 *    strnum_cmp() order (natural) or strcmp() order (byte), followed by a
 *    0x00 terminator and a byte holding (flag & 0xC0) >> 6.  Sorted with an
 *    MSD radix sort taking 8 key bytes per level, starting after the
 *    batch-wide longest common prefix.
 *
 * The sort only ever moves 16-byte {key, ref} pairs; the records themselves
 * are touched once, in a final gather pass.  All sorts are stable with
 * respect to the initial order of the kv array, which the caller sets to
 * input order, so ties keep input order.
 */

#ifndef SORT2_ENGINE_H
#define SORT2_ENGINE_H

#include <stddef.h>
#include <stdint.h>

#if defined(__GNUC__) || defined(__clang__)
#define S2_PREFETCH(p) __builtin_prefetch(p)
#else
#define S2_PREFETCH(p) ((void)0)
#endif

typedef struct {
    uint64_t key;   // sort key (coordinate) or 8-byte key window (name)
    uint64_t ref;   // index of the record in its batch
} s2_kv_t;

/* Run fn(arg, t, nthreads) for t in [0, nthreads) on separate threads and
 * wait for all to finish.  nthreads <= 1 runs inline. */
typedef void (*s2_task_fn)(void *arg, int t, int nthreads);
void s2_parallel(int nthreads, s2_task_fn fn, void *arg);

/* Contiguous chunk [*lo, *hi) of n items for worker t of nt. */
static inline void s2_chunk(size_t n, int t, int nt, size_t *lo, size_t *hi) {
    // 64-bit intermediate: n * t cannot overflow for any realistic batch,
    // including on 32-bit platforms.
    *lo = (size_t)((uint64_t)n * t / nt);
    *hi = (size_t)((uint64_t)n * (t + 1) / nt);
}

/* Stable LSD radix sort on kv.key.  tmp must hold n elements.  Returns a
 * pointer to whichever of a or tmp holds the result. */
s2_kv_t *s2_radix_sort(s2_kv_t *a, s2_kv_t *tmp, size_t n, int nthreads);

/* Upper bound on the encoded key length for a name of name_len bytes. */
static inline size_t s2_name_key_bound(size_t name_len) {
    return 3 * name_len + 8;
}

/* Encode a name.  natural != 0 gives samtools strnum_cmp ordering,
 * otherwise strcmp ordering.  flagbyte is appended after a 0x00 terminator
 * so names that compare equal are then ordered by it.  Returns bytes
 * written (at most s2_name_key_bound(name_len)). */
size_t s2_name_key(const char *name, size_t name_len, int natural,
                   uint8_t flagbyte, uint8_t *out);

/* Compare two encoded keys: memcmp order, shorter-is-smaller on a common
 * prefix. */
int s2_key_cmp(const uint8_t *a, uint32_t la, const uint8_t *b, uint32_t lb);

/* Stable sort of a[0..n) by the encoded keys key[a[i].ref] / len[a[i].ref].
 * The result is left in a; tmp is scratch of n elements.  a[i].key is
 * overwritten. */
void s2_name_sort(s2_kv_t *a, s2_kv_t *tmp, size_t n,
                  const uint8_t *const *key, const uint32_t *len,
                  int nthreads);

/* Loser tree for k-way merging.  less(ctx, i, j) must be a strict total
 * order over sources, treating exhausted sources as +infinity and breaking
 * ties by source index. */
typedef int (*s2_less_fn)(void *ctx, int i, int j);
typedef struct {
    int k, kp;          // sources, leaves (power of 2)
    int *loser;         // internal nodes [1, kp)
    int winner;
    s2_less_fn less;
    void *ctx;
} s2_ltree_t;

int  s2_ltree_init(s2_ltree_t *t, int k, s2_less_fn less, void *ctx);
static inline int s2_ltree_top(const s2_ltree_t *t) { return t->winner; }
/* Call after the current winner's source has advanced to its next item. */
void s2_ltree_replay(s2_ltree_t *t);
void s2_ltree_destroy(s2_ltree_t *t);

#endif
