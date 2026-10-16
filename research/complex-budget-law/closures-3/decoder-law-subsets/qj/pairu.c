// min |A ∪ B| for A in cands(j=0, span sa), B in cands(j=jb, span sb), sizes <= pmax
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static inline int pc(u128 x){return __builtin_popcountll((uint64_t)x)+__builtin_popcountll((uint64_t)(x>>64));}
int main(int argc,char**argv){
  int sa=atoi(argv[2]), jb=atoi(argv[3]), sb=atoi(argv[4]), pmax=atoi(argv[5]);
  FILE*f=fopen(argv[1],"r"); u128 *A=malloc(sizeof(u128)*4000000),*B=malloc(sizeof(u128)*4000000); long na=0,nb=0;
  int j,s; unsigned long long hi,lo;
  while(fscanf(f,"%d %d %llu %llu",&j,&s,&hi,&lo)==4){u128 m=((u128)hi<<64)|lo; if(pc(m)>pmax)continue; if(j==0&&s==sa)A[na++]=m; if(j==jb&&s==sb)B[nb++]=m;}
  int best=999; long cnt[200]={0};
  for(long i=0;i<na;i++) for(long k=0;k<nb;k++){int u=pc(A[i]|B[k]); if(u<best)best=u; if(u<200)cnt[u]++;}
  printf("sa=%d jb=%d sb=%d na=%ld nb=%ld min_union=%d\n",sa,jb,sb,na,nb,best);
  return 0;
}
