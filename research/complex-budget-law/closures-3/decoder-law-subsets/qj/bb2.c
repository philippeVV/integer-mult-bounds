// max sum of delta (span7 ->1, span6 ->2) over distinct root lines simultaneously low, union of witnesses <= P.
// Exact B&B; root 0 WLOG low. Stops early (reports >LIMIT) if the sum exceeds LIMIT.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static inline int pc(u128 x){return __builtin_popcountll((uint64_t)x)+__builtin_popcountll((uint64_t)(x>>64));}
static u128 *W[120]; static int *D[120]; static int nw[120]; static int P, best, LIMIT; static long nodes;
static void rec(int last,u128 U,int m){ nodes++; if(m>best)best=m; if(best>LIMIT) return;
  for(int b=last+1;b<120;b++) for(int i=0;i<nw[b];i++){u128 V=U|W[b][i]; if(pc(V)<=P) {rec(b,V,m+D[b][i]); if(best>LIMIT)return;}}
}
int main(int argc,char**argv){ P=atoi(argv[2]); LIMIT=atoi(argv[3]);
  FILE*f=fopen(argv[1],"r"); int j,s; unsigned long long hi,lo; for(int a=0;a<120;a++){W[a]=malloc(sizeof(u128)*40000);D[a]=malloc(sizeof(int)*40000);nw[a]=0;}
  while(fscanf(f,"%d %d %llu %llu",&j,&s,&hi,&lo)==4){u128 m=((u128)hi<<64)|lo; if(pc(m)<=P){W[j][nw[j]]=m; D[j][nw[j]++]=8-s;}}
  best=0; for(int i=0;i<nw[0]&&best<=LIMIT;i++) rec(0,W[0][i],D[0][i]);
  printf("P=%d witnesses(root0)=%d max_sum_delta=%d%s (limit %d) nodes=%ld\n",P,nw[0],best,best>LIMIT?"+ EXCEEDS":"",LIMIT,nodes); return 0; }
