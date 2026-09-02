#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/*
 * Direct numerical-semigroup search.
 *
 * This program does NOT use the Run Lemma or divisor-orbit loss bounds.
 * It enumerates bounded Kunz vectors, checks the Kunz inequalities, constructs
 * their gap sets, and computes arithmetic progressions directly in those gaps.
 */

typedef struct {
    int A, m, n;
    int *k;
    int *order;
    int best;
    uint64_t nodes;
    uint64_t complete_candidates;
    uint64_t direct_checks;
} Search;

static int is_prime_int(int n) {
    if (n < 2) return 0;
    for (int d = 2; d * d <= n; ++d) if (n % d == 0) return 0;
    return 1;
}

static int largest_prime_leq(int n) {
    for (int p = n; p >= 2; --p) if (is_prime_int(p)) return p;
    return 0;
}

static int smallest_prime_factor(int n) {
    for (int p = 2; p * p <= n; ++p) if (n % p == 0) return p;
    return n;
}

static int kunz_partial_ok(const Search *s) {
    int m = s->m;
    for (int i = 1; i < m; ++i) {
        if (!s->k[i]) continue;
        for (int j = i; j < m; ++j) {
            if (!s->k[j]) continue;
            int z = i + j;
            if (z == m) continue;
            if (z < m) {
                if (s->k[z] && s->k[i] + s->k[j] < s->k[z]) return 0;
            } else {
                z -= m;
                if (s->k[z] && s->k[i] + s->k[j] + 1 < s->k[z]) return 0;
            }
        }
    }
    return 1;
}

static int frobenius_from_kunz(int m, const int *k) {
    int F = -1;
    for (int r = 1; r < m; ++r) {
        if (!k[r]) continue;
        int x = r + m * (k[r] - 1);
        if (x > F) F = x;
    }
    return F;
}

static uint8_t *gap_array_from_kunz(int m, const int *k, int *F_out) {
    int F = frobenius_from_kunz(m, k);
    if (F < 0) F = 0;
    uint8_t *gap = (uint8_t *)calloc((size_t)F + 1u, 1u);
    if (!gap) return NULL;
    for (int r = 1; r < m; ++r) {
        for (int h = 0; h < k[r]; ++h) {
            int x = r + h * m;
            if (x <= F) gap[x] = 1;
        }
    }
    *F_out = F;
    return gap;
}

/* Returns the maximum AP length; if stop_after>=0, may return immediately
 * once a run longer than stop_after is found. */
static int direct_ap_from_gap(const uint8_t *gap, int F, int stop_after,
                              int *max_diffs, int *ndiffs, int capdiffs) {
    int best = 0;
    int count = 0;
    if (F <= 0) {
        if (ndiffs) *ndiffs = 0;
        return 0;
    }
    for (int d = 1; d <= F; ++d) {
        int best_d = 0;
        /* A maximal run starts at a gap x for which x-d is not a positive gap. */
        for (int x = 1; x <= F; ++x) {
            if (!gap[x]) continue;
            if (x - d >= 1 && gap[x-d]) continue;
            int len = 0;
            for (int y = x; y <= F && gap[y]; y += d) ++len;
            if (len > best_d) best_d = len;
            if (stop_after >= 0 && len > stop_after) {
                if (ndiffs) *ndiffs = 0;
                return len;
            }
        }
        if (best_d > best) {
            best = best_d;
            count = 0;
            if (max_diffs && count < capdiffs) max_diffs[count] = d;
            count = 1;
        } else if (best_d == best && best_d > 0) {
            if (max_diffs && count < capdiffs) max_diffs[count] = d;
            ++count;
        }
    }
    if (ndiffs) *ndiffs = count;
    return best;
}

static int direct_ap_kunz(int m, const int *k, int stop_after) {
    int F = 0;
    uint8_t *gap = gap_array_from_kunz(m, k, &F);
    if (!gap) {
        fprintf(stderr, "allocation failure while constructing gap set\n");
        exit(2);
    }
    int ans = direct_ap_from_gap(gap, F, stop_after, NULL, NULL, 0);
    free(gap);
    return ans;
}

/* Direct partial pruning: assigned coordinates already determine genuine gaps.
 * If those guaranteed gaps contain an (A+1)-term AP, every completion fails. */
static int partial_direct_ok(Search *s) {
    int assigned = 0;
    for (int i=1;i<s->m;i++) if (s->k[i]) ++assigned;
    if (assigned < 2) return 1;
    s->direct_checks++;
    return direct_ap_kunz(s->m, s->k, s->A) <= s->A;
}

static void make_order(Search *s) {
    int *used = (int *)calloc((size_t)s->m, sizeof(int));
    if (!used) { fprintf(stderr, "allocation failure while building search order\n"); exit(2); }
    int pos = 0;
    int q = smallest_prime_factor(s->m);
    if (q < s->m) {
        int d = s->m / q;
        /* Put complete nonzero d-orbits first. This changes only search order;
           admissibility is still decided by direct gap AP checks. */
        for (int r = 1; r < d; ++r) {
            for (int x = r; x < s->m; x += d) {
                if (!used[x]) { s->order[pos++] = x; used[x] = 1; }
            }
        }
        for (int x = d; x < s->m; x += d) {
            if (!used[x]) { s->order[pos++] = x; used[x] = 1; }
        }
    }
    for (int r=1;r<s->m;r++) if(!used[r]) s->order[pos++]=r;
    s->n = pos;
    free(used);
}

static void dfs(Search *s, int depth, int sum) {
    s->nodes++;
    if (depth == s->n) {
        s->complete_candidates++;
        if (sum <= s->best) return;
        s->direct_checks++;
        if (direct_ap_kunz(s->m, s->k, s->A) <= s->A) s->best = sum;
        return;
    }
    int left = s->n - depth - 1;
    int idx = s->order[depth];
    for (int v = s->A; v >= 1; --v) {
        if (sum + v + s->A * left <= s->best) break;
        s->k[idx] = v;
        if (!kunz_partial_ok(s)) { s->k[idx]=0; continue; }
        if (!partial_direct_ok(s)) { s->k[idx]=0; continue; }
        dfs(s, depth+1, sum+v);
        s->k[idx]=0;
    }
}

static int compute_M(int A, uint64_t *nodes, uint64_t *complete, uint64_t *checks) {
    int p = largest_prime_leq(A+1);
    int best = A * (p-1); /* verified lower bound from the uniform prime vector */
    uint64_t tn=0,tc=0,td=0;

    /* For m<=p, sum k_i <= A(m-1) <= A(p-1), so only m>p can improve. */
    for (int m = p+1; m <= A+1; ++m) {
        Search s; memset(&s,0,sizeof(s));
        s.A=A; s.m=m; s.best=best;
        s.k=(int *)calloc((size_t)m,sizeof(int));
        s.order=(int *)calloc((size_t)(m-1),sizeof(int));
        if(!s.k || !s.order){fprintf(stderr,"allocation failure in direct search\n");free(s.k);free(s.order);exit(2);}
        make_order(&s);
        dfs(&s,0,0);
        if (s.best > best) best=s.best;
        tn+=s.nodes; tc+=s.complete_candidates; td+=s.direct_checks;
        free(s.k); free(s.order);
    }
    if (nodes) *nodes = tn;
    if (complete) *complete = tc;
    if (checks) *checks = td;
    return best;
}

static int test_vector(const char *name,int m,const int *v,int n,int eg,int eap){
    int *k=(int *)calloc((size_t)m,sizeof(int)); int g=0;
    if(!k){fprintf(stderr,"allocation failure in unit test\n");exit(2);}
    for(int i=0;i<n;i++){k[i+1]=v[i];g+=v[i];}
    int ap=direct_ap_kunz(m,k,-1);
    free(k);
    int ok=(g==eg && ap==eap);
    printf("UNIT_TEST name=%s genus=%d AP=%d status=%s\n",name,g,ap,ok?"PASS":"FAIL");
    return ok;
}

static int mode_tests(void){
    const int a[]={1,1,1,1,1};
    const int b[]={3,3,1};
    const int c[]={4,4,4,4};
    const int d[]={3,3,6,6,9,9,9,9};
    int ok=1;
    ok &= test_vector("ordinary_m6",6,a,5,5,5);
    ok &= test_vector("generators_4_7_13",4,b,3,7,3);
    ok &= test_vector("uniform_prime_m5_A4",5,c,4,16,4);
    ok &= test_vector("extremizer_m9_A9",9,d,8,54,9);
    printf("UNIT_TEST_SUMMARY status=%s\n",ok?"PASS":"FAIL");
    return ok?0:1;
}

static int mode_extremizer9(void){
    const int v[]={3,3,6,6,9,9,9,9};
    int *k=(int *)calloc(9,sizeof(int)); int g=0; if(!k)return 2;
    for(int i=0;i<8;i++){k[i+1]=v[i];g+=v[i];}
    int F=0; uint8_t *gap=gap_array_from_kunz(9,k,&F); free(k); if(!gap)return 2;
    int diffs[256],nd=0; int ap=direct_ap_from_gap(gap,F,-1,diffs,&nd,256); free(gap);
    printf("EXTREMIZER_RESULT m=9 genus=%d F=%d AP=%d maximizing_differences=",g,F,ap);
    for(int i=0;i<nd;i++){if(i)putchar(',');printf("%d",diffs[i]);} putchar('\n');
    int ok=(g==54 && F==80 && ap==9 && nd==2 && diffs[0]==3 && diffs[1]==9);
    printf("EXTREMIZER_ASSERTIONS status=%s\n",ok?"PASS":"FAIL");
    return ok?0:1;
}

static int mode_a8(void){
    uint64_t n=0,c=0,d=0; int M=compute_M(8,&n,&c,&d);
    printf("DIRECT_RESULT A=8 M=%d nodes=%llu complete_candidates=%llu direct_AP_checks=%llu\n",
      M,(unsigned long long)n,(unsigned long long)c,(unsigned long long)d);
    int ok=(M==48); printf("DIRECT_A8_ASSERTIONS status=%s\n",ok?"PASS":"FAIL"); return ok?0:1;
}

static int mode_table(int maxA){
    if(maxA<1 || maxA>20){fprintf(stderr,"table max_A must lie in 1..20\n");return 2;}
    int vals[21]={0}; int ok=1;
    const int expected12[]={0,1,4,7,16,20,36,42,48,54,100,110,144};
    for(int A=1;A<=maxA;A++) vals[A]=compute_M(A,NULL,NULL,NULL);
    printf("TABLE_RESULT max_A=%d values=",maxA);
    for(int A=1;A<=maxA;A++){if(A>1)putchar(',');printf("%d:%d",A,vals[A]);} putchar('\n');
    if(maxA==12) for(int A=1;A<=12;A++) if(vals[A]!=expected12[A]) ok=0;
    printf("TABLE_ASSERTIONS status=%s\n",ok?"PASS":"FAIL");
    return ok?0:1;
}

int main(int argc,char **argv){
    if(argc<2){fprintf(stderr,"usage: %s tests|extremizer9|a8|table N\n",argv[0]);return 2;}
    if(strcmp(argv[1],"tests")==0)return mode_tests();
    if(strcmp(argv[1],"extremizer9")==0)return mode_extremizer9();
    if(strcmp(argv[1],"a8")==0)return mode_a8();
    if(strcmp(argv[1],"table")==0 && argc==3)return mode_table(atoi(argv[2]));
    fprintf(stderr,"unknown mode\n"); return 2;
}
