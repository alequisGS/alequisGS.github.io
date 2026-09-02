#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

/*
 * Exact verification of the composite-multiplicity criterion
 *
 *   (m/q - 1) (A - floor(A/q)) >= A (m-p),
 *
 * where p=P(A+1), for every composite p<m<=A+1 and at least one prime
 * divisor q|m. The default/full advertised run is 1 <= A <= 100000000.
 */

typedef struct { uint64_t A, m, p; } Failure;

static uint16_t *build_lpf(uint32_t limit) {
    uint16_t *lpf = (uint16_t *)calloc((size_t)limit + 1u, sizeof(uint16_t));
    if (!lpf) return NULL;
    for (uint32_t i = 2; (uint64_t)i * i <= limit; ++i) {
        if (lpf[i] != 0) continue; /* i has a smaller factor */
        for (uint64_t j = (uint64_t)i * i; j <= limit; j += i) {
            if (lpf[j] == 0) lpf[j] = (uint16_t)i;
        }
    }
    return lpf;
}

static int criterion_holds(uint64_t A, uint64_t m, uint64_t p,
                           const uint16_t *lpf) {
    uint64_t x = m;
    while (x > 1) {
        uint64_t q = lpf[x] ? (uint64_t)lpf[x] : x;
        uint64_t left = (m / q - 1u) * (A - A / q);
        uint64_t right = A * (m - p);
        if (left >= right) return 1;
        do { x /= q; } while (x > 1 && x % q == 0);
    }
    return 0;
}

int main(int argc, char **argv) {
    uint64_t maxA = 100000000ULL;
    if (argc > 2) {
        fprintf(stderr, "usage: %s [max_A]\n", argv[0]);
        return 2;
    }
    if (argc == 2) maxA = strtoull(argv[1], NULL, 10);
    if (maxA < 1 || maxA > 100000000ULL) {
        fprintf(stderr, "max_A must lie in 1..100000000\n");
        return 2;
    }

    uint32_t limit = (uint32_t)(maxA + 1u);
    uint16_t *lpf = build_lpf(limit);
    if (!lpf) {
        fprintf(stderr, "allocation failure for least-prime-factor table\n");
        return 2;
    }

    uint64_t p = 0, previous_prime = 0;
    uint64_t maxgap = 0, gap_after = 0;
    uint64_t pairs = 0, failures = 0;
    Failure first[16];

    for (uint64_t A = 1; A <= maxA; ++A) {
        uint64_t n = A + 1u;
        if (n >= 2 && lpf[n] == 0) {
            previous_prime = p;
            p = n;
            if (previous_prime && p - previous_prime > maxgap) {
                maxgap = p - previous_prime;
                gap_after = previous_prime;
            }
        }
        if (p < 2) continue;

        for (uint64_t m = p + 1u; m <= n; ++m) {
            ++pairs;
            if (!criterion_holds(A, m, p, lpf)) {
                if (failures < 16) {
                    first[failures].A=A; first[failures].m=m; first[failures].p=p;
                }
                ++failures;
                printf("CRITERION_FAILURE A=%" PRIu64 " m=%" PRIu64 " p=%" PRIu64 "\n", A,m,p);
            }
        }
    }

    printf("CRITERION_SUMMARY max_A=%" PRIu64 " pairs_checked=%" PRIu64
           " failures=%" PRIu64 " max_consecutive_prime_gap=%" PRIu64
           " after_prime=%" PRIu64 "\n",
           maxA,pairs,failures,maxgap,gap_after);

    int ok = 1;
    if (maxA == 100000000ULL) {
        const Failure e[4]={{3,4,3},{8,9,7},{9,9,7},{9,10,7}};
        if (pairs != 1449077567ULL || failures != 4 ||
            maxgap != 220 || gap_after != 47326693ULL) ok=0;
        if (failures == 4) {
            for (int i=0;i<4;i++)
                if (first[i].A!=e[i].A || first[i].m!=e[i].m || first[i].p!=e[i].p) ok=0;
        }
    }
    printf("CRITERION_ASSERTIONS status=%s\n", ok ? "PASS" : "FAIL");
    free(lpf);
    return ok ? 0 : 1;
}
