/* two-level lex Singer orders on PG(2,16): point i key (m1 i, m2 i, 0), line l key (m1 l - c1, m2 l - c2, 1) mod N.
   Boot cost with k_P=8,k_L=9. random sampling. usage: lex samples seed mode(0=nonunit m1,1=any)  prints best list */
#include <stdio.h>
#include <stdlib.h>
#define N 273
int D[17],nbv[2*N][17];
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long xr(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
int gcd(int a,int b){while(b){int t=a%b;a=b;b=t;}return a;}
int cost(int m1,int c1,int m2,int c2){
  static int key[2*N],pos[2*N],cnt[2*N*N*2/1+2];
  /* key ordering via counting sort on composite key up to 2*N*N */
  static int bucket[2*N*N+1];
  for(int i=0;i<N;i++){key[i]=(((m1*i)%N)*N+((m2*i)%N))*2; key[N+i]=((((m1*i-c1)%N+N)%N)*N+(((m2*i-c2)%N+N)%N))*2+1;}
  /* keys are distinct? points distinct only if (m1,m2) separate; ties broken by index (stable) */
  static int idx[2*N];
  for(int v=0;v<2*N;v++)idx[v]=v;
  /* simple sort: qsort with comparator */
  int cmp(const void*a,const void*b){int x=*(int*)a,y=*(int*)b; return key[x]!=key[y]?key[x]-key[y]:x-y;}
  qsort(idx,2*N,sizeof(int),cmp);
  for(int t=0;t<2*N;t++)pos[idx[t]]=t;
  int c=0;
  for(int v=0;v<2*N;v++){int j=0;for(int k=0;k<17;k++)j+=pos[nbv[v][k]]<pos[v]; int K=v<N?8:9; if(j<K)c+=K-j;}
  (void)cnt;(void)bucket;
  return c;
}
int main(int argc,char**argv){
  FILE*f=fopen("D16.txt","r");for(int k=0;k<17;k++)fscanf(f,"%d",&D[k]);fclose(f);
  for(int i=0;i<N;i++)for(int k=0;k<17;k++){int l=(i+D[k])%N;nbv[i][k]=N+l;}
  for(int l=0;l<N;l++)for(int k=0;k<17;k++){nbv[N+l][k]=((l-D[k])%N+N)%N;}
  long long S=atoll(argv[1]); rs^=atoll(argv[2])*0x9E3779B97F4A7C15ULL; int mode=atoi(argv[3]);
  int nu=0,units[N],nn=0,non[N]; for(int m=0;m<N;m++){int g=gcd(m,N); if(g==1)units[nu++]=m; else if(m) non[nn++]=m;}
  int best=1<<30;
  for(long long s=0;s<S;s++){
    int m1=mode? (int)(xr()%N): non[xr()%nn]; int m2=units[xr()%nu]; int c1=xr()%N,c2=xr()%N;
    int c=cost(m1,c1,m2,c2);
    if(c<=best+2){ if(c<best)best=c; printf("%d %d %d %d %d\n",c,m1,c1,m2,c2); fflush(stdout);}
  }
  return 0;
}
