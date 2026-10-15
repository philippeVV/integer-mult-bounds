/* Node profiles for every subspace X of F_2^9 (dims 1..8) w.r.t. the 112 ports and 8 centre normals (vecs.txt).
   For X (content side, ports only): Pp = #ports in X, C1 = max #ports in one hyperplane of X, C2 = max #ports in the union of two.
   For X^perp (reach side, V = ports + normals): QV = #V in X^perp, Qp = #ports in X^perp,
   R1/R2 = max #V in one / union of two hyperplanes of X^perp  (an out-edge X->Y has reach in Y^perp, a hyperplane of X^perp).
   Output: aggregated distinct types per dim: "d Pp C1 C2 QV Qp R1 R2 count". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef unsigned __int128 u128;
static int isport[512], isnorm[512];
static int pc128(u128 x){return __builtin_popcountll((unsigned long long)x)+__builtin_popcountll((unsigned long long)(x>>64));}
/* stats of span(b[0..e-1]) for set mask: which=1 ports, which=3 ports|normals */
static void stat(int*b,int e,int which,int*cnt,int*c1,int*c2){
  static int el[256]; static int co[256]; int n=0, N=1<<e;
  el[0]=0; for(int c=1;c<N;c++){int lb=__builtin_ctz(c); el[c]=el[c&(c-1)]^b[lb];}
  for(int c=0;c<N;c++){int v=el[c]; if(isport[v]||((which&2)&&isnorm[v])) co[n++]=c;}
  *cnt=n; *c1=0; if(c2)*c2=0; if(e==0||n==0||!c2) return;
  static u128 m[256]; int W=N;
  for(int w=1;w<W;w++){u128 x=0; for(int i=0;i<n;i++) if(!(__builtin_popcount(co[i]&w)&1)) x|=((u128)1)<<i; m[w]=x; int p=pc128(x); if(p>*c1)*c1=p;}
  if(*c1==n){*c2=n;return;}
  int best=*c1;
  for(int a=1;a<W;a++){ if(pc128(m[a])+*c1<=best) continue; for(int b2=a+1;b2<W;b2++){int p=pc128(m[a]|m[b2]); if(p>best){best=p; if(best==n) goto done;}}}
  done: *c2=best;
}
#define HS (1<<22)
static unsigned long long hk[HS]; static long long hv[HS];
static void add(unsigned long long k){unsigned long long h=(k*0x9E3779B97F4A7C15ULL)>>42; while(hk[h]&&hk[h]!=k) h=(h+1)&(HS-1); hk[h]=k; hv[h]++;}
int main(int argc,char**argv){
  int D=atoi(argv[1]); FILE*f=fopen("vecs.txt","r"); int v,i=0; while(fscanf(f,"%d",&v)==1){ if(i<112) isport[v]=1; else isnorm[v]=1; i++;}
  int piv[9]; long long nsub=0;
  for(int pm=0;pm<512;pm++){ if(__builtin_popcount(pm)!=D) continue;
    int k=0; for(int c=0;c<9;c++) if(pm>>c&1) piv[k++]=c;
    /* free positions: row r may have bits at nonpivot columns c > piv[r] */
    int fr[64],fc[64],F=0; for(int r=0;r<D;r++) for(int c=piv[r]+1;c<9;c++) if(!(pm>>c&1)){fr[F]=r;fc[F]=c;F++;}
    for(long long fm=0;fm<(1LL<<F);fm++){
      int row[9]; for(int r=0;r<D;r++) row[r]=1<<piv[r];
      for(int t=0;t<F;t++) if(fm>>t&1) row[fr[t]]|=1<<fc[t];
      int pb[9],e=0; for(int j=0;j<9;j++) if(!(pm>>j&1)){int x=1<<j; for(int r=0;r<D;r++) if(row[r]>>j&1) x|=1<<piv[r]; pb[e++]=x;}
      int Pp,C1,C2,QV,R1,R2,Qp,d1,d2; stat(row,D,1,&Pp,&C1,&C2); stat(pb,e,3,&QV,&R1,&R2); stat(pb,e,1,&Qp,&d1,NULL);
      unsigned long long key=((unsigned long long)Pp)|((unsigned long long)C1<<7)|((unsigned long long)C2<<14)|((unsigned long long)QV<<21)|((unsigned long long)Qp<<28)|((unsigned long long)R1<<35)|((unsigned long long)R2<<42)|(1ULL<<60);
      add(key); nsub++;
    }
  }
  for(long long h=0;h<HS;h++) if(hk[h]){unsigned long long k=hk[h]; printf("%d %llu %llu %llu %llu %llu %llu %llu %lld\n",D,k&127,(k>>7)&127,(k>>14)&127,(k>>21)&127,(k>>28)&127,(k>>35)&127,(k>>42)&127,hv[h]);}
  fprintf(stderr,"dim %d subspaces %lld\n",D,nsub); return 0;
}
