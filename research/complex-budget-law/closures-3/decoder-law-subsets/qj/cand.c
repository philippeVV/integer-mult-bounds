// For each input set N_j and each span s (dim X = s, s=1..8) dump candidate deletion sets N_j \ Lines(X)
// with size <= PMAX (deduplicated by caller). Output: j s hi lo
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
typedef unsigned __int128 u128;
static int lab[128]; static u128 H[512];
static inline int pc(u128 x){return __builtin_popcountll((uint64_t)x)+__builtin_popcountll((uint64_t)(x>>64));}
int main(int argc,char**argv){
  int PMAX=atoi(argv[2]); FILE*f=fopen(argv[1],"r"); int NL,K; if(fscanf(f,"%d",&NL)!=1)return 1;
  for(int i=0;i<NL;i++) if(fscanf(f,"%d",&lab[i])!=1)return 1; if(fscanf(f,"%d",&K)!=1)return 1;
  u128 S[300]; for(int j=0;j<K;j++){unsigned long long hi,lo; if(fscanf(f,"%llu %llu",&hi,&lo)!=2)return 1; S[j]=((u128)hi<<64)|lo;}
  for(int g=1;g<512;g++){u128 m=0; for(int i=0;i<NL;i++) if(!(__builtin_popcount(lab[i]&g)&1)) m|=((u128)1)<<i; H[g]=m;}
  for(int c=1;c<=8;c++){
    for(int P=0;P<512;P++){ if(__builtin_popcount(P)!=c) continue;
      int pv[9],np=0; for(int b=0;b<9;b++) if(P>>b&1) pv[np++]=b;
      int fr[9][9],nf[9],tot=0; for(int i=0;i<c;i++){nf[i]=0; for(int b=pv[i]+1;b<9;b++) if(!(P>>b&1)) fr[i][nf[i]++]=b; tot+=nf[i];}
      for(long code=0; code<(1L<<tot); code++){
        long cc=code; u128 m=~(u128)0;
        for(int i=0;i<c;i++){int row=1<<pv[i]; for(int t=0;t<nf[i];t++){ if(cc&1) row|=1<<fr[i][t]; cc>>=1;} m&=H[row];}
        for(int j=0;j<K;j++){u128 A=S[j]&~m; if(pc(A)<=PMAX) printf("%d %d %llu %llu\n",j,9-c,(unsigned long long)(A>>64),(unsigned long long)A);}
      }
    }
  }
  return 0;
}
