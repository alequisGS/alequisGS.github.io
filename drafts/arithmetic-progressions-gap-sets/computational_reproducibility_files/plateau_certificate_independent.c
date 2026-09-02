#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/*
 * Independent finite cross-check for the A=9 plateau.
 *
 * This implementation uses a different coordinate order and encodes the
 * orbit restrictions directly as forbidden completed patterns, rather than
 * calling the Run-Lemma evaluator used in plateau_certificate.c.
 */

typedef struct {
    int m, n;
    int order[10];
    int k[11];
    int best;
    uint64_t best_count;
    int best_vec[11];
    uint64_t nodes, leaves, kunz_rejects, orbit_rejects, bound_prunes;
} Search;

static int kunz_ok_assigned(const Search *s) {
    int m = s->m;
    for (int i = 1; i < m; ++i) {
        if (!s->k[i]) continue;
        for (int j = i; j < m; ++j) {
            if (!s->k[j]) continue;
            int t = i + j;
            if (t == m) continue;
            if (t < m) {
                if (s->k[t] && s->k[i] + s->k[j] < s->k[t]) return 0;
            } else {
                t -= m;
                if (s->k[t] && s->k[i] + s->k[j] + 1 < s->k[t]) return 0;
            }
        }
    }
    return 1;
}

static int all_set(const Search *s, const int *a, int n) {
    for (int i = 0; i < n; ++i) if (!s->k[a[i]]) return 0;
    return 1;
}

static int orbit_ok(const Search *s) {
    if (s->m == 10) {
        for (int r = 1; r <= 4; ++r)
            if (s->k[r] && s->k[r+5] && s->k[r] >= 5 && s->k[r+5] >= 5)
                return 0;
        const int a[5] = {1,3,5,7,9};
        if (all_set(s,a,5)) {
            int has1 = 0;
            for (int i=0;i<5;i++) if (s->k[a[i]]==1) has1=1;
            if (!has1) return 0;
        }
    } else if (s->m == 9) {
        const int cyc[2][3] = {{1,4,7},{2,5,8}};
        for (int c=0;c<2;c++) if (all_set(s,cyc[c],3)) {
            int mu=10, first=-1;
            for (int j=0;j<3;j++) {
                int v=s->k[cyc[c][j]];
                if (v<mu) { mu=v; first=j; }
            }
            /* Equivalent to 3*mu+rho <= 9. */
            if (mu>3 || (mu==3 && first!=0)) return 0;
        }
    } else if (s->m == 8) {
        for (int r=1;r<=3;r++)
            if (s->k[r] && s->k[r+4] && s->k[r] >=5 && s->k[r+4] >=5)
                return 0;
        const int a[4]={1,3,5,7};
        if (all_set(s,a,4)) {
            int mu=10, first=-1;
            for(int j=0;j<4;j++) {
                int v=s->k[a[j]];
                if(v<mu){mu=v;first=j;}
            }
            /* Equivalent to 4*mu+rho <= 9. */
            if(mu>2 || (mu==2 && first>1)) return 0;
        }
    }
    return 1;
}

static void search(Search *s, int depth, int sum) {
    s->nodes++;
    if (depth == s->n) {
        s->leaves++;
        if (sum > s->best) {
            s->best=sum; s->best_count=1;
            for(int i=1;i<s->m;i++) s->best_vec[i]=s->k[i];
        } else if(sum==s->best) s->best_count++;
        return;
    }
    int left=s->n-depth-1;
    int idx=s->order[depth];
    for(int v=9;v>=1;--v) {
        if(sum+v+9*left < s->best){s->bound_prunes++;break;}
        s->k[idx]=v;
        if(!orbit_ok(s)){s->orbit_rejects++;s->k[idx]=0;continue;}
        if(!kunz_ok_assigned(s)){s->kunz_rejects++;s->k[idx]=0;continue;}
        search(s,depth+1,sum+v);
        s->k[idx]=0;
    }
}

static void set_order(Search *s) {
    /* Deliberately different from natural 1,2,... order. */
    if(s->m==8){int o[]={1,5,3,7,2,6,4}; memcpy(s->order,o,sizeof(o));}
    if(s->m==9){int o[]={1,4,7,2,5,8,3,6}; memcpy(s->order,o,sizeof(o));}
    if(s->m==10){int o[]={1,6,3,8,5,7,9,2,4}; memcpy(s->order,o,sizeof(o));}
}

static int eqv(const Search *s,const int *e){for(int i=1;i<s->m;i++)if(s->best_vec[i]!=e[i-1])return 0;return 1;}
static void pv(const Search *s){putchar('(');for(int i=1;i<s->m;i++){if(i>1)putchar(',');printf("%d",s->best_vec[i]);}putchar(')');}

static int run(int m,int eb,const int *ev){
    Search s; memset(&s,0,sizeof(s)); s.m=m; s.n=m-1; set_order(&s); search(&s,0,0);
    printf("INDEPENDENT_RESULT m=%d max_genus=%d maximizers=%llu vector=",m,s.best,(unsigned long long)s.best_count);pv(&s);
    printf(" nodes=%llu leaves=%llu kunz_rejects=%llu orbit_rejects=%llu bound_prunes=%llu\n",
      (unsigned long long)s.nodes,(unsigned long long)s.leaves,(unsigned long long)s.kunz_rejects,
      (unsigned long long)s.orbit_rejects,(unsigned long long)s.bound_prunes);
    return s.best==eb && s.best_count==1 && eqv(&s,ev);
}

int main(void){
    const int e8[]={7,9,2,9,4,4,9};
    const int e9[]={3,3,6,6,9,9,9,9};
    const int e10[]={9,9,9,8,6,4,4,3,1};
    int ok=run(8,44,e8)&run(9,54,e9)&run(10,53,e10);
    printf("INDEPENDENT_ASSERTIONS status=%s\n",ok?"PASS":"FAIL");
    return ok?0:1;
}
