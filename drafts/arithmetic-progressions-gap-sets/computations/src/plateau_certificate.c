#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/*
 * Exhaustive finite Kunz certificate for A=9 and m=8,9,10.
 *
 * The search uses only
 *   (i) 1 <= k_i <= 9,
 *  (ii) all Kunz inequalities, and
 * (iii) the explicit divisor-orbit constraints stated in the manuscript.
 * It does not construct gap sets and does not compute AP(S) directly.
 */

typedef struct {
    int m;
    int k[11];          /* coordinates 1,...,m-1; 0 means unassigned */
    int best;
    uint64_t best_count;
    int best_vec[11];
    uint64_t nodes;
    uint64_t leaves;
    uint64_t kunz_pruned;
    uint64_t orbit_pruned;
    uint64_t bound_pruned;
} Search;

static int kunz_partial_ok(const Search *s) {
    const int m = s->m;
    for (int i = 1; i < m; ++i) {
        if (s->k[i] == 0) continue;
        for (int j = i; j < m; ++j) {
            if (s->k[j] == 0) continue;
            const int sum = i + j;
            if (sum == m) continue;
            if (sum < m) {
                if (s->k[sum] != 0 && s->k[i] + s->k[j] < s->k[sum])
                    return 0;
            } else {
                const int t = sum - m;
                if (s->k[t] != 0 && s->k[i] + s->k[j] + 1 < s->k[t])
                    return 0;
            }
        }
    }
    return 1;
}

/* Check an orbit c[0],...,c[len-1] once all its coordinates are assigned. */
static int run_orbit_ok(const Search *s, const int *c, int len, int A) {
    int mu = 1000;
    int rho = -1;
    for (int j = 0; j < len; ++j) {
        if (s->k[c[j]] == 0) return 1; /* not yet decidable */
        if (s->k[c[j]] < mu) {
            mu = s->k[c[j]];
            rho = j;
        }
    }
    return len * mu + rho <= A;
}

static int orbit_partial_ok(const Search *s) {
    if (s->m == 10) {
        /* d=5: nonzero 2-cycles (r,r+5), r=1,...,4. */
        for (int r = 1; r <= 4; ++r) {
            if (s->k[r] && s->k[r+5] &&
                s->k[r] > 4 && s->k[r+5] > 4)
                return 0;
        }
        /* d=2: odd orbit (1,3,5,7,9), length 5, hence minimum 1. */
        const int odd[5] = {1,3,5,7,9};
        if (!run_orbit_ok(s, odd, 5, 9)) return 0;
    } else if (s->m == 9) {
        const int o1[3] = {1,4,7};
        const int o2[3] = {2,5,8};
        if (!run_orbit_ok(s, o1, 3, 9)) return 0;
        if (!run_orbit_ok(s, o2, 3, 9)) return 0;
    } else if (s->m == 8) {
        /* d=4: nonzero 2-cycles (i,i+4), i=1,2,3. */
        for (int i = 1; i <= 3; ++i) {
            if (s->k[i] && s->k[i+4] &&
                s->k[i] > 4 && s->k[i+4] > 4)
                return 0;
        }
        /* d=2: odd orbit (1,3,5,7), retaining the positional term. */
        const int odd[4] = {1,3,5,7};
        if (!run_orbit_ok(s, odd, 4, 9)) return 0;
    }
    return 1;
}

static void dfs(Search *s, int idx, int sum) {
    s->nodes++;
    if (idx == s->m) {
        s->leaves++;
        if (sum > s->best) {
            s->best = sum;
            s->best_count = 1;
            for (int i = 1; i < s->m; ++i) s->best_vec[i] = s->k[i];
        } else if (sum == s->best) {
            s->best_count++;
        }
        return;
    }

    const int remaining_after = (s->m - 1) - idx;
    for (int v = 9; v >= 1; --v) {
        if (sum + v + 9 * remaining_after < s->best) {
            s->bound_pruned++;
            break;
        }
        s->k[idx] = v;
        if (!kunz_partial_ok(s)) {
            s->kunz_pruned++;
            s->k[idx] = 0;
            continue;
        }
        if (!orbit_partial_ok(s)) {
            s->orbit_pruned++;
            s->k[idx] = 0;
            continue;
        }
        dfs(s, idx + 1, sum + v);
        s->k[idx] = 0;
    }
}

static int vector_equal(const Search *s, const int *expected) {
    for (int i = 1; i < s->m; ++i)
        if (s->best_vec[i] != expected[i-1]) return 0;
    return 1;
}

static void print_vector(const Search *s) {
    putchar('(');
    for (int i = 1; i < s->m; ++i) {
        if (i > 1) putchar(',');
        printf("%d", s->best_vec[i]);
    }
    putchar(')');
}

static int run_one(int m, int expected_best, const int *expected_vec) {
    Search s;
    memset(&s, 0, sizeof(s));
    s.m = m;
    s.best = 0;
    dfs(&s, 1, 0);

    uint64_t theoretical = 1;
    for (int i = 1; i < m; ++i) theoretical *= 9;

    printf("KUNZ_RESULT m=%d max_genus=%d maximizers=%llu vector=",
           m, s.best, (unsigned long long)s.best_count);
    print_vector(&s);
    printf(" theoretical=%llu evaluated_leaves=%llu nodes=%llu kunz_pruned=%llu orbit_pruned=%llu bound_pruned=%llu\n",
           (unsigned long long)theoretical,
           (unsigned long long)s.leaves,
           (unsigned long long)s.nodes,
           (unsigned long long)s.kunz_pruned,
           (unsigned long long)s.orbit_pruned,
           (unsigned long long)s.bound_pruned);

    return s.best == expected_best && s.best_count == 1 && vector_equal(&s, expected_vec);
}

int main(void) {
    const int e8[]  = {7,9,2,9,4,4,9};
    const int e9[]  = {3,3,6,6,9,9,9,9};
    const int e10[] = {9,9,9,8,6,4,4,3,1};

    int ok = 1;
    ok &= run_one(8, 44, e8);
    ok &= run_one(9, 54, e9);
    ok &= run_one(10, 53, e10);

    printf("KUNZ_ASSERTIONS status=%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
