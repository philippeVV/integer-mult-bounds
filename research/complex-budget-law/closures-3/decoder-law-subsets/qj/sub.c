// For every subspace X of F2^9 (as kernel of a c-dim functional space, c=1..8, RREF enumeration)
// and every input set S_j of lines: best[j][c] = max |S_j ∩ Lines(X)|.  Also counts subspaces.
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static int NL; static int lab[128]; static u128 H[512]; static int K; static u128 *S; static int *best; static long cnt[10];
static int c_target;
static inline int pc(u128 x){return __builtin_popcountll((uint64_t)x)+__builtin_popcountll((uint64_t)(x>>64));}
static void rec(int nrows, int *piv, int startpiv, u128 mask){
  if(nrows==c_target){ cnt[c_target]++; for(int j=0;j<K;j++){int v=pc(mask&S[j]); if(v>best[j*10+c_target])best[j*10+c_target]=v;} return; }
  // choose next pivot q>=startpiv; but rows must be built with free entries only in non-pivot cols > pivot.
  // We build rows in order of pivots; a row with pivot q has free entries in columns > q that are not pivots
  // of LATER rows. To handle that, enumerate pivot sets first.
}
// simpler: enumerate pivot set, then all fillings
static int pivs[9];
static void fill(int r, int c, u128 mask, int *rowsfree, int nfree_total){}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); if(fscanf(f,"%d",&NL)!=1)return 1; for(int i=0;i<NL;i++) if(fscanf(f,"%d",&lab[i])!=1) return 1;
  if(fscanf(f,"%d",&K)!=1)return 1; S=malloc(sizeof(u128)*K); best=calloc(K*10,sizeof(int));
  for(int j=0;j<K;j++){unsigned long long hi,lo; if(fscanf(f,"%llu %llu",&hi,&lo)!=2)return 1; S[j]=((u128)hi<<64)|lo;}
  for(int g=1;g<512;g++){u128 m=0; for(int i=0;i<NL;i++) if(!(__builtin_popcount(lab[i]&g)&1)) m|=((u128)1)<<i; H[g]=m;}
  for(int c=1;c<=8;c++){
    for(int P=0;P<512;P++){ if(__builtin_popcount(P)!=c) continue;
      int pv[9],np=0; for(int b=0;b<9;b++) if(P>>b&1) pv[np++]=b;
      // row i has pivot pv[i] (its lowest set bit), free bits at non-pivot positions > pv[i]
      int fr[9][9],nf[9],tot=0; for(int i=0;i<c;i++){nf[i]=0; for(int b=pv[i]+1;b<9;b++) if(!(P>>b&1)) fr[i][nf[i]++]=b; tot+=nf[i];}
      for(long code=0; code<(1L<<tot); code++){
        long cc=code; u128 m=~(u128)0;
        for(int i=0;i<c;i++){int row=1<<pv[i]; for(int t=0;t<nf[i];t++){ if(cc&1) row|=1<<fr[i][t]; cc>>=1;} m&=H[row];}
        cnt[c]++;
        for(int j=0;j<K;j++){int v=pc(m&S[j]); if(v>best[j*10+c])best[j*10+c]=v;}
      }
    }
  }
  for(int c=1;c<=8;c++) fprintf(stderr,"codim %d: %ld subspaces\n",c,cnt[c]);
  for(int j=0;j<K;j++){ for(int c=1;c<=8;c++) printf("%d ",best[j*10+c]); printf("\n"); }
  return 0;
}
