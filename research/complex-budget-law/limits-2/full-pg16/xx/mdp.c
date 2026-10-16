/* Exact DP over all merges of the cyclic point sequence (rank a = m i mod N) and line sequence (rank b = (m l - c) mod N),
   with the class orientation O (out-classes P->L).  F_class = sum_a f_a(h_a), h nondecreasing = #lines before point a.
   f_a(h) = #{d in O: (a+y_d) mod N < h} + #{d notin O: (a+y_d) mod N >= h},  y_d = (m d - c) mod N.
   usage: mdp mode [m c]   mode 0: scan m (all units), c, O = 9 cyclically consecutive in y-sorted order (17 windows)
                            mode 1: given m c, all C(17,9) out-sets  ; prints best and writes merge order to mdp_best.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 273
int D[17];
int gcd(int a,int b){while(b){int t=a%b;a=b;b=t;}return a;}
static int dp[N+1],nd[N+1],cnt[N+1],choice[N][N+1];
int solve(int *y,int *out,int keep){
  /* dp over a: dp[h] = min cost of points 0..a with h_a = h */
  for(int h=0;h<=N;h++)dp[h]=0;
  for(int a=0;a<N;a++){
    /* f_a(h): start h=0: all in-classes count (pos>=0), out none. as h increases past value v=(a+y_d)%N: out d adds 1 at h=v+1, in d removes 1 */
    memset(cnt,0,sizeof cnt); int base=0;
    for(int k=0;k<17;k++){int v=(a+y[k])%N; if(out[k]) cnt[v+1]++; else {base++; cnt[v+1]--;}}
    int f=base, run=1<<30, arg=0;
    for(int h=0;h<=N;h++){ f+= (h? cnt[h]:0); if(dp[h]<run){run=dp[h];arg=h;} nd[h]=run+f; if(keep)choice[a][h]=arg; }
    memcpy(dp,nd,sizeof dp);
  }
  int best=1<<30; for(int h=0;h<=N;h++)if(dp[h]<best)best=dp[h]; return best;
}
int main(int argc,char**argv){
  FILE*f=fopen("D16.txt","r");for(int k=0;k<17;k++)fscanf(f,"%d",&D[k]);fclose(f);
  int mode=atoi(argv[1]); int y[17],out[17];
  if(mode==0){
    int gbest=1<<30;
    for(int m=1;m<N;m++){ if(gcd(m,N)!=1)continue;
      for(int c=0;c<N;c++){
        int ord[17]; for(int k=0;k<17;k++){y[k]=((m*D[k]-c)%N+N)%N; ord[k]=k;}
        for(int i=0;i<17;i++)for(int j=i+1;j<17;j++)if(y[ord[j]]<y[ord[i]]){int t=ord[i];ord[i]=ord[j];ord[j]=t;}
        for(int w=0;w<17;w++){ memset(out,0,sizeof out); for(int t=0;t<9;t++)out[ord[(w+t)%17]]=1;
          int v=solve(y,out,0); if(v<=gbest){gbest=v; printf("%d m %d c %d w %d\n",v,m,c,w); fflush(stdout);} }
      }
    }
  } else {
    int m=atoi(argv[2]),c=atoi(argv[3]); for(int k=0;k<17;k++)y[k]=((m*D[k]-c)%N+N)%N;
    int best=1<<30,bm=0;
    for(int mask=0;mask<(1<<17);mask++){ if(__builtin_popcount(mask)!=9)continue; for(int k=0;k<17;k++)out[k]=mask>>k&1;
      int v=solve(y,out,0); if(v<best){best=v;bm=mask;printf("%d mask %x\n",v,mask);fflush(stdout);} }
    for(int k=0;k<17;k++)out[k]=bm>>k&1; int v=solve(y,out,1);
    /* backtrack h */
    int h[N]; int hb=0,bv=1<<30; for(int x=0;x<=N;x++)if(dp[x]<bv){bv=dp[x];hb=x;}
    for(int a=N-1;a>=0;a--){h[a]=hb; hb=choice[a][hb];}
    /* build order: point rank a -> i = a * m^{-1}; line rank b -> l with (m l - c) = b */
    int minv=0; for(int t=1;t<N;t++)if((m*t)%N==1)minv=t;
    FILE*g=fopen("mdp_best.txt","w"); int b=0;
    for(int a=0;a<N;a++){ while(b<h[a]){int l=((long)(b+c)*minv)%N; fprintf(g,"%d\n",N+l); b++;} fprintf(g,"%d\n",(a*minv)%N); }
    while(b<N){int l=((long)(b+c)*minv)%N; fprintf(g,"%d\n",N+l); b++;}
    fclose(g); printf("best %d (check %d) written mdp_best.txt\n",best,v);
  }
  return 0;}
