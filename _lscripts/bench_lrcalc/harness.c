// Exhaustive Schubert-product timing harness for liblrcalc (1.2 API or lrcalc-rs/2.x API).
// Usage: harness n start_u end_u outfile
// Enumerates S_n in lexicographic order; for u index in [start_u, end_u) and every v index >= u,
// times mult_schubert(u, v, 0) and appends a 24-byte record:
//   uint32 ui, uint32 vi, uint32 time_ns, uint32 nterms, uint64 hash
// hash = sum over terms of (fnv1a64(trimmed perm) * coeff) mod 2^64 (order independent).
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef LRCALC_RS
#include <lrcalc/ivector.h>
#include <lrcalc/ivlincomb.h>
#include <lrcalc/schublib.h>
typedef ivector PVEC;
#define PV_NEW(n) iv_new(n)
#define PV_FREE(v) iv_free(v)
#define PV_ELEM(v, i) iv_elem(v, i)
#else
#include <lrcalc/vector.h>
#include <lrcalc/hashtab.h>
#include <lrcalc/schublib.h>
typedef vector PVEC;
#define PV_NEW(n) v_new(n)
#define PV_FREE(v) v_free(v)
#define PV_ELEM(v, i) v_elem(v, i)
#endif

static uint64_t fnv_perm(const int32_t *a, int len) {
    while (len > 0 && a[len - 1] == len) len--;  // strip trailing fixed points
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < len; i++) {
        h ^= (uint64_t)(uint32_t)a[i];
        h *= 1099511628211ULL;
    }
    return h;
}

static int next_perm(int *a, int n) {
    int i = n - 2;
    while (i >= 0 && a[i] > a[i + 1]) i--;
    if (i < 0) return 0;
    int j = n - 1;
    while (a[j] < a[i]) j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) { t = a[l]; a[l] = a[r]; a[r] = t; }
    return 1;
}

typedef struct { uint32_t ui, vi, time_ns, nterms; uint64_t hash; } rec_t;

int main(int argc, char **argv) {
    if (argc != 5) { fprintf(stderr, "usage: %s n start_u end_u outfile\n", argv[0]); return 2; }
    int n = atoi(argv[1]);
    long start_u = atol(argv[2]), end_u = atol(argv[3]);
    FILE *out = fopen(argv[4], "wb");
    if (!out) { perror("fopen"); return 1; }

    // enumerate all perms of S_n once
    long total = 1; for (int i = 2; i <= n; i++) total *= i;
    int *perms = malloc(sizeof(int) * n * total);
    int cur[16]; for (int i = 0; i < n; i++) cur[i] = i + 1;
    long idx = 0;
    do { memcpy(perms + idx * n, cur, sizeof(int) * n); idx++; } while (next_perm(cur, n));

    PVEC *u = PV_NEW(n), *v = PV_NEW(n);
    rec_t buf[4096]; int nb = 0;
    struct timespec t0, t1;
    for (long ui = start_u; ui < end_u && ui < total; ui++) {
        for (int i = 0; i < n; i++) PV_ELEM(u, i) = perms[ui * n + i];
        for (long vi = ui; vi < total; vi++) {
            for (int i = 0; i < n; i++) PV_ELEM(v, i) = perms[vi * n + i];
            clock_gettime(CLOCK_MONOTONIC, &t0);
#ifdef LRCALC_RS
            ivlincomb *res = mult_schubert(u, v, 0);
#else
            hashtab *res = mult_schubert(u, v, 0);
#endif
            clock_gettime(CLOCK_MONOTONIC, &t1);
            uint64_t h = 0; uint32_t nt = 0;
#ifdef LRCALC_RS
            ivlc_iter itr;
            for (ivlc_first(res, &itr); ivlc_good(&itr); ivlc_next(&itr)) {
                ivector *k = ivlc_key(&itr);
                h += fnv_perm(k->array, (int)iv_length(k)) * (uint64_t)(int64_t)ivlc_value(&itr);
                nt++;
            }
            ivlc_free_all(res);
#else
            hash_itr itr;
            for (hash_first(res, itr); hash_good(itr); hash_next(itr)) {
                vector *k = (vector *)hash_key(itr);
                h += fnv_perm(k->array, (int)v_length(k)) * (uint64_t)(int64_t)hash_intvalue(itr);
                nt++;
            }
            for (hash_first(res, itr); hash_good(itr); hash_next(itr)) v_free((vector *)hash_key(itr));
            hash_free(res);
#endif
            long ns = (t1.tv_sec - t0.tv_sec) * 1000000000L + (t1.tv_nsec - t0.tv_nsec);
            rec_t r = { (uint32_t)ui, (uint32_t)vi, (uint32_t)(ns > 0xFFFFFFFFL ? 0xFFFFFFFFL : ns), nt, h };
            buf[nb++] = r;
            if (nb == 4096) { fwrite(buf, sizeof(rec_t), nb, out); nb = 0; }
        }
    }
    if (nb) fwrite(buf, sizeof(rec_t), nb, out);
    fclose(out);
    return 0;
}
