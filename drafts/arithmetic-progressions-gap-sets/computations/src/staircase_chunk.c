#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

typedef struct { uint64_t A, m, p; } Failure;

static uint16_t *build_lpf(uint32_t limit) {
    uint16_t *lpf=(uint16_t*)calloc((size_t)limit+1u,sizeof(uint16_t));
    if(!lpf)return NULL;
    for(uint32_t i=2;(uint64_t)i*i<=limit;++i){
        if(lpf[i])continue;
        for(uint64_t j=(uint64_t)i*i;j<=limit;j+=i) if(!lpf[j]) lpf[j]=(uint16_t)i;
    }
    return lpf;
}

static int holds(uint64_t A,uint64_t m,uint64_t p,const uint16_t *lpf){
    uint64_t x=m;
    while(x>1){
        uint64_t q=lpf[x]?(uint64_t)lpf[x]:x;
        uint64_t left=(m/q-1u)*(A-A/q);
        uint64_t right=A*(m-p);
        if(left>=right)return 1;
        do{x/=q;}while(x>1 && x%q==0);
    }
    return 0;
}

static uint64_t previous_prime_leq(uint64_t x,const uint16_t *lpf){
    while(x>=2){if(lpf[x]==0)return x;--x;}return 0;
}

int main(int argc,char **argv){
    if(argc!=3){fprintf(stderr,"usage: %s start_A end_A   (inclusive)\n",argv[0]);return 2;}
    uint64_t A0=strtoull(argv[1],NULL,10), A1=strtoull(argv[2],NULL,10);
    if(A0<1 || A1<A0 || A1>100000000ULL){
        fprintf(stderr,"require 1 <= start_A <= end_A <= 100000000\n");return 2;
    }
    uint32_t limit=(uint32_t)(A1+1u);
    uint16_t *lpf=build_lpf(limit);
    if(!lpf){fprintf(stderr,"allocation failure for least-prime-factor table\n");return 2;}

    /* Before processing A0, p must be the largest prime <= A0; the iteration
       itself updates p if A0+1 is prime. */
    uint64_t p=previous_prime_leq(A0,lpf);
    uint64_t pairs=0, failures=0;
    Failure first[16];

    for(uint64_t A=A0;A<=A1;++A){
        uint64_t n=A+1u;
        if(n>=2 && lpf[n]==0) p=n;
        if(p<2)continue;
        for(uint64_t m=p+1u;m<=n;++m){
            ++pairs;
            if(!holds(A,m,p,lpf)){
                if(failures<16){first[failures].A=A;first[failures].m=m;first[failures].p=p;}
                ++failures;
                printf("CHUNK_FAILURE A=%" PRIu64 " m=%" PRIu64 " p=%" PRIu64 "\n",A,m,p);
            }
        }
        if(A==UINT64_MAX)break;
    }

    printf("CHUNK_SUMMARY start=%" PRIu64 " end=%" PRIu64
           " pairs_checked=%" PRIu64 " failures=%" PRIu64 "\n",A0,A1,pairs,failures);

    int ok=1;
    if(A0==1 && (A1==100000ULL || A1==100000000ULL)){
        const Failure e[4]={{3,4,3},{8,9,7},{9,9,7},{9,10,7}};
        uint64_t expected_pairs=(A1==100000ULL)?779997ULL:1449077567ULL;
        if(pairs!=expected_pairs || failures!=4)ok=0;
        if(failures==4)for(int i=0;i<4;i++)
            if(first[i].A!=e[i].A || first[i].m!=e[i].m || first[i].p!=e[i].p)ok=0;
    }
    printf("CHUNK_ASSERTIONS status=%s\n",ok?"PASS":"FAIL");
    free(lpf);return ok?0:1;
}
