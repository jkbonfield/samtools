/*  sort2_engine.c -- htslib-independent sorting engine for samtools sort2.

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

#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include "sort2_engine.h"

// Below this many items the radix sort runs single-threaded.
#define RADIX_PAR_MIN 65536
// Name groups at or below this size are finished with insertion sort.
#define INSERTION_MAX 24

/* ---------------------------------------------------------------------- */
/* Fork/join helper */

typedef struct {
    s2_task_fn fn;
    void *arg;
    int t, nt;
} s2_job_t;

static void *s2_job_run(void *vp) {
    s2_job_t *j = vp;
    j->fn(j->arg, j->t, j->nt);
    return NULL;
}

void s2_parallel(int nt, s2_task_fn fn, void *arg) {
    if (nt <= 1) {
        fn(arg, 0, 1);
        return;
    }
    pthread_t *tid = malloc(nt * sizeof(*tid));
    s2_job_t *job = malloc(nt * sizeof(*job));
    char *started = calloc(nt, 1);
    if (!tid || !job || !started) {
        // Degrade to serial; each chunk is still processed exactly once.
        for (int t = 0; t < nt; t++)
            fn(arg, t, nt);
        free(tid); free(job); free(started);
        return;
    }
    for (int t = 0; t < nt; t++) {
        job[t] = (s2_job_t){fn, arg, t, nt};
        if (t > 0 && pthread_create(&tid[t], NULL, s2_job_run, &job[t]) == 0)
            started[t] = 1;
    }
    // Worker 0 runs on the calling thread, as does any worker whose
    // thread failed to start.
    for (int t = 0; t < nt; t++)
        if (!started[t])
            fn(arg, t, nt);
    for (int t = 1; t < nt; t++)
        if (started[t])
            pthread_join(tid[t], NULL);
    free(tid); free(job); free(started);
}

/* ---------------------------------------------------------------------- */
/* Parallel stable LSD radix sort, 8-bit digits, constant digits skipped */

typedef struct {
    s2_kv_t *src, *dst;
    size_t n;
    int shift;
    size_t *hist;                // nt * 256, per-thread counts then offsets
    uint64_t *or_v, *and_v;      // per-thread OR / AND of keys
} radix_arg_t;

static void radix_orand(void *vp, int t, int nt) {
    radix_arg_t *r = vp;
    size_t lo, hi;
    s2_chunk(r->n, t, nt, &lo, &hi);
    uint64_t o = 0, a = ~(uint64_t)0;
    for (size_t i = lo; i < hi; i++) {
        o |= r->src[i].key;
        a &= r->src[i].key;
    }
    r->or_v[t] = o;
    r->and_v[t] = a;
}

static void radix_hist(void *vp, int t, int nt) {
    radix_arg_t *r = vp;
    size_t lo, hi, *h = r->hist + (size_t)t * 256;
    s2_chunk(r->n, t, nt, &lo, &hi);
    memset(h, 0, 256 * sizeof(*h));
    const int sh = r->shift;
    const s2_kv_t *s = r->src;
    for (size_t i = lo; i < hi; i++)
        h[(s[i].key >> sh) & 0xff]++;
}

static void radix_scatter(void *vp, int t, int nt) {
    radix_arg_t *r = vp;
    size_t lo, hi, off[256];
    s2_chunk(r->n, t, nt, &lo, &hi);
    memcpy(off, r->hist + (size_t)t * 256, sizeof(off));
    const int sh = r->shift;
    const s2_kv_t *s = r->src;
    s2_kv_t *d = r->dst;
    for (size_t i = lo; i < hi; i++)
        d[off[(s[i].key >> sh) & 0xff]++] = s[i];
}

s2_kv_t *s2_radix_sort(s2_kv_t *a, s2_kv_t *tmp, size_t n, int nt) {
    if (n < 2)
        return a;
    if (nt < 1 || n < RADIX_PAR_MIN)
        nt = 1;

    size_t hist_stack[256];
    uint64_t or_stack, and_stack;
    radix_arg_t r = { a, tmp, n, 0, hist_stack, &or_stack, &and_stack };
    if (nt > 1) {
        r.hist  = malloc((size_t)nt * 256 * sizeof(size_t));
        r.or_v  = malloc(nt * sizeof(uint64_t));
        r.and_v = malloc(nt * sizeof(uint64_t));
        if (!r.hist || !r.or_v || !r.and_v) {
            free(r.hist); free(r.or_v); free(r.and_v);
            nt = 1;
            r.hist = hist_stack; r.or_v = &or_stack; r.and_v = &and_stack;
        }
    }

    // Digits where every key has the same byte need no pass.
    s2_parallel(nt, radix_orand, &r);
    uint64_t o = 0, an = ~(uint64_t)0;
    for (int t = 0; t < nt; t++) {
        o |= r.or_v[t];
        an &= r.and_v[t];
    }
    uint64_t varying = o ^ an;

    for (int d = 0; d < 8; d++) {
        r.shift = d * 8;
        if (((varying >> r.shift) & 0xff) == 0)
            continue;
        s2_parallel(nt, radix_hist, &r);
        // Bucket-major, thread-minor prefix sum keeps the sort stable.
        size_t sum = 0;
        for (int b = 0; b < 256; b++) {
            for (int t = 0; t < nt; t++) {
                size_t c = r.hist[(size_t)t * 256 + b];
                r.hist[(size_t)t * 256 + b] = sum;
                sum += c;
            }
        }
        s2_parallel(nt, radix_scatter, &r);
        s2_kv_t *x = r.src; r.src = r.dst; r.dst = x;
    }

    if (nt > 1) {
        free(r.hist); free(r.or_v); free(r.and_v);
    }
    return r.src;
}

/* ---------------------------------------------------------------------- */
/* Name keys */

size_t s2_name_key(const char *name, size_t n, int natural,
                   uint8_t flagbyte, uint8_t *out) {
    const unsigned char *s = (const unsigned char *) name;
    uint8_t *p = out;
    if (!natural) {
        memcpy(p, s, n);
        p += n;
    } else {
        size_t i = 0;
        while (i < n) {
            unsigned c = s[i];
            if (c < '0' || c > '9') {
                *p++ = c;
                i++;
                continue;
            }
            // A digit run.  Leading zeros are ignored by strnum_cmp.
            // Marker 0x30 sorts correctly against any non-digit, because
            // strnum_cmp compares a digit to a non-digit as raw bytes and
            // all digits lie in 0x30..0x39.  Then an order-preserving
            // length, then the significant digits.
            while (i < n && s[i] == '0')
                i++;
            size_t st = i;
            while (i < n && s[i] >= '0' && s[i] <= '9')
                i++;
            size_t len = i - st;
            *p++ = 0x30;
            if (len < 255) {
                *p++ = (uint8_t) len;
            } else {
                *p++ = 0xff;
                *p++ = (uint8_t)(len >> 24);
                *p++ = (uint8_t)(len >> 16);
                *p++ = (uint8_t)(len >> 8);
                *p++ = (uint8_t) len;
            }
            memcpy(p, s + st, len);
            p += len;
        }
    }
    // Terminator sorts before any further name byte (names contain no NUL,
    // and the natural encoding never emits 0x00 at a token start).
    *p++ = 0;
    *p++ = flagbyte;
    return p - out;
}

int s2_key_cmp(const uint8_t *a, uint32_t la, const uint8_t *b, uint32_t lb) {
    uint32_t m = la < lb ? la : lb;
    int c = memcmp(a, b, m);
    if (c)
        return c;
    return (la > lb) - (la < lb);
}

/* ---------------------------------------------------------------------- */
/* MSD name sort */

typedef struct {
    const uint8_t *const *key;
    const uint32_t *len;
} nkeys_t;

static inline uint64_t load_be64(const uint8_t *p) {
    uint64_t v;
    memcpy(&v, p, 8);
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__ \
    && (defined(__GNUC__) || defined(__clang__))
    v = __builtin_bswap64(v);
#elif !(defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
    v = 0;
    for (int i = 0; i < 8; i++)
        v = (v << 8) | p[i];
#endif
    return v;
}

// 8 key bytes from off as a big-endian integer, zero padded past the end.
static inline uint64_t key_window(const nkeys_t *k, uint64_t ref,
                                  uint32_t off) {
    uint32_t l = k->len[ref];
    const uint8_t *p = k->key[ref];
    if (l >= off + 8)
        return load_be64(p + off);
    uint64_t v = 0;
    for (uint32_t i = off; i < l; i++)
        v |= (uint64_t) p[i] << (56 - 8 * (i - off));
    return v;
}

// Compare two keys known to be equal before off.
static inline int key_cmp_from(const nkeys_t *k, uint64_t a, uint64_t b,
                               uint32_t off) {
    uint32_t la = k->len[a], lb = k->len[b];
    uint32_t m = la < lb ? la : lb;
    if (off < m) {
        int c = memcmp(k->key[a] + off, k->key[b] + off, m - off);
        if (c)
            return c;
    }
    return (la > lb) - (la < lb);
}

static void insertion_sort(const nkeys_t *k, s2_kv_t *a, size_t n,
                           uint32_t off) {
    for (size_t i = 1; i < n; i++) {
        s2_kv_t x = a[i];
        size_t j = i;
        // Strictly greater: equal keys keep their order (stable).
        while (j > 0 && key_cmp_from(k, a[j-1].ref, x.ref, off) > 0) {
            a[j] = a[j-1];
            j--;
        }
        a[j] = x;
    }
}

// Serial MSD sort of a group whose keys agree on bytes [0, off).
static void sort_group(const nkeys_t *k, s2_kv_t *a, s2_kv_t *t, size_t n,
                       uint32_t off) {
    if (n <= INSERTION_MAX) {
        insertion_sort(k, a, n, off);
        return;
    }
    int more = 0;
    for (size_t i = 0; i < n; i++) {
        a[i].key = key_window(k, a[i].ref, off);
        more |= k->len[a[i].ref] > off + 8;
    }
    s2_kv_t *r = s2_radix_sort(a, t, n, 1);
    if (r != a)
        memcpy(a, r, n * sizeof(*a));

    for (size_t i = 0; i < n; ) {
        size_t j = i + 1;
        while (j < n && a[j].key == a[i].key)
            j++;
        if (j - i > 1) {
            if (more)
                sort_group(k, a + i, t + i, j - i, off + 8);
            else
                // All keys end within this window and agree on it, so only
                // lengths can differ; insertion sort is O(n) when sorted.
                insertion_sort(k, a + i, j - i, off);
        }
        i = j;
    }
}

typedef struct {
    size_t start, n;
    uint32_t off;
} ngroup_t;

typedef struct {
    const nkeys_t *k;
    s2_kv_t *a, *t;
    // the large group currently being processed
    ngroup_t g;
    s2_kv_t *res;            // radix output for g (a+g.start or t+g.start)
    size_t big;              // groups above this are split in parallel
    // pending large groups
    pthread_mutex_t lock;
    ngroup_t *stack;
    size_t nstack, mstack;
    int oom;
    // lcp computation
    uint32_t *lcp;
    // per-worker [bounds[t], bounds[t+1]) for ns_split
    size_t *bounds;
} nsort_t;

static void ns_lcp(void *vp, int t, int nt) {
    nsort_t *s = vp;
    size_t lo, hi;
    s2_chunk(s->g.n, t, nt, &lo, &hi);
    uint64_t r0 = s->a[0].ref;
    const uint8_t *k0 = s->k->key[r0];
    uint32_t l = s->k->len[r0];
    for (size_t i = lo; i < hi && l; i++) {
        uint64_t r = s->a[i].ref;
        uint32_t li = s->k->len[r], m = li < l ? li : l, j = 0;
        const uint8_t *ki = s->k->key[r];
        while (j < m && ki[j] == k0[j])
            j++;
        l = j;
    }
    s->lcp[t] = l;
}

static void ns_fill(void *vp, int t, int nt) {
    nsort_t *s = vp;
    size_t lo, hi;
    s2_chunk(s->g.n, t, nt, &lo, &hi);
    s2_kv_t *a = s->a + s->g.start;
    for (size_t i = lo; i < hi; i++)
        a[i].key = key_window(s->k, a[i].ref, s->g.off);
}

static void ns_copyback(void *vp, int t, int nt) {
    nsort_t *s = vp;
    size_t lo, hi;
    s2_chunk(s->g.n, t, nt, &lo, &hi);
    memcpy(s->a + s->g.start + lo, s->res + lo, (hi - lo) * sizeof(s2_kv_t));
}

static int ns_push(nsort_t *s, ngroup_t g) {
    pthread_mutex_lock(&s->lock);
    if (s->nstack == s->mstack) {
        size_t m = s->mstack ? s->mstack * 2 : 64;
        ngroup_t *x = realloc(s->stack, m * sizeof(*x));
        if (!x) {
            s->oom = 1;
            pthread_mutex_unlock(&s->lock);
            return -1;
        }
        s->stack = x;
        s->mstack = m;
    }
    s->stack[s->nstack++] = g;
    pthread_mutex_unlock(&s->lock);
    return 0;
}

// Each worker scans a chunk of the radix-sorted large group, with chunk
// edges moved to subgroup boundaries.  Small subgroups are finished
// inline; large ones are queued for another parallel level.
static void ns_split(void *vp, int t, int nt) {
    (void) nt;
    nsort_t *s = vp;
    s2_kv_t *a = s->a + s->g.start, *tt = s->t + s->g.start;
    // Bounds were moved to subgroup starts before any worker began, so
    // no worker reads keys that another one is rewriting.
    size_t lo = s->bounds[t], hi = s->bounds[t+1];
    uint32_t off = s->g.off;
    for (size_t i = lo; i < hi; ) {
        size_t j = i + 1;
        int more = s->k->len[a[i].ref] > off + 8;
        while (j < hi && a[j].key == a[i].key) {
            more |= s->k->len[a[j].ref] > off + 8;
            j++;
        }
        if (j - i > 1) {
            if (more && j - i > s->big) {
                ngroup_t g = { s->g.start + i, j - i, off + 8 };
                if (ns_push(s, g) < 0)  // could not queue it; do it here
                    sort_group(s->k, a + i, tt + i, j - i, off + 8);
            } else if (more) {
                sort_group(s->k, a + i, tt + i, j - i, off + 8);
            } else {
                insertion_sort(s->k, a + i, j - i, off);
            }
        }
        i = j;
    }
}

void s2_name_sort(s2_kv_t *a, s2_kv_t *tmp, size_t n,
                  const uint8_t *const *key, const uint32_t *len, int nt) {
    if (n < 2)
        return;
    if (nt < 1)
        nt = 1;
    nkeys_t k = { key, len };

    nsort_t s;
    memset(&s, 0, sizeof(s));
    s.k = &k;
    s.a = a;
    s.t = tmp;
    s.big = n / (8 * (size_t) nt);
    if (s.big < RADIX_PAR_MIN)
        s.big = RADIX_PAR_MIN;
    pthread_mutex_init(&s.lock, NULL);

    // Skip the batch-wide common prefix (e.g. instrument:run:flowcell).
    uint32_t lcp_buf[nt];
    size_t bounds[nt + 1];
    s.lcp = lcp_buf;
    s.bounds = bounds;
    s.g = (ngroup_t){0, n, 0};
    int lnt = n < RADIX_PAR_MIN ? 1 : nt;
    s2_parallel(lnt, ns_lcp, &s);
    uint32_t lcp = lcp_buf[0];
    for (int t = 1; t < lnt; t++)
        if (lcp_buf[t] < lcp)
            lcp = lcp_buf[t];

    if (nt == 1 || n <= s.big) {
        sort_group(&k, a, tmp, n, lcp);
        pthread_mutex_destroy(&s.lock);
        return;
    }

    if (ns_push(&s, (ngroup_t){0, n, lcp}) < 0) {
        sort_group(&k, a, tmp, n, lcp);
        pthread_mutex_destroy(&s.lock);
        return;
    }
    while (s.nstack > 0 && !s.oom) {
        s.g = s.stack[--s.nstack];
        s2_parallel(nt, ns_fill, &s);
        s.res = s2_radix_sort(a + s.g.start, tmp + s.g.start, s.g.n, nt);
        if (s.res != a + s.g.start)
            s2_parallel(nt, ns_copyback, &s);
        // Chunk edges, each moved forward to the start of a subgroup.
        s2_kv_t *ga = a + s.g.start;
        for (int t = 0; t <= nt; t++) {
            size_t lo, hi;
            s2_chunk(s.g.n, t < nt ? t : nt - 1, nt, &lo, &hi);
            size_t b = t < nt ? lo : hi;
            if (t > 0 && b < bounds[t-1])
                b = bounds[t-1];
            while (b > 0 && b < s.g.n && ga[b].key == ga[b-1].key)
                b++;
            bounds[t] = b;
        }
        s2_parallel(nt, ns_split, &s);
    }
    // If the work stack could not grow, finish what remains serially.
    while (s.nstack > 0) {
        ngroup_t g = s.stack[--s.nstack];
        sort_group(&k, a + g.start, tmp + g.start, g.n, g.off);
    }
    free(s.stack);
    pthread_mutex_destroy(&s.lock);
}

/* ---------------------------------------------------------------------- */
/* Loser tree */

static inline int lt_less(s2_ltree_t *t, int i, int j) {
    if (i >= t->k) return 0;     // padding leaves are +infinity
    if (j >= t->k) return 1;
    return t->less(t->ctx, i, j);
}

static int lt_build(s2_ltree_t *t, int node) {
    if (node >= t->kp)
        return node - t->kp;
    int l = lt_build(t, 2 * node), r = lt_build(t, 2 * node + 1);
    if (lt_less(t, r, l)) {
        t->loser[node] = l;
        return r;
    }
    t->loser[node] = r;
    return l;
}

int s2_ltree_init(s2_ltree_t *t, int k, s2_less_fn less, void *ctx) {
    memset(t, 0, sizeof(*t));
    if (k < 1)
        return -1;
    int kp = 1;
    while (kp < k)
        kp *= 2;
    t->k = k;
    t->kp = kp;
    t->less = less;
    t->ctx = ctx;
    t->loser = malloc((kp > 1 ? kp : 2) * sizeof(int));
    if (!t->loser)
        return -1;
    t->winner = kp > 1 ? lt_build(t, 1) : 0;
    return 0;
}

void s2_ltree_replay(s2_ltree_t *t) {
    int w = t->winner;
    for (int node = (w + t->kp) / 2; node >= 1; node /= 2) {
        int l = t->loser[node];
        if (lt_less(t, l, w)) {
            t->loser[node] = w;
            w = l;
        }
    }
    t->winner = w;
}

void s2_ltree_destroy(s2_ltree_t *t) {
    free(t->loser);
    t->loser = NULL;
}
