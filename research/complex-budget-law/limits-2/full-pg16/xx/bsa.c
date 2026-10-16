/* Boot-cost SA over orders, K=8 on v<N (points), 9 on lines. insertion moves.  usage: bsa graph orderfile|rand iters T0 T1 seed out */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MV 600
int V,d,nb[MV][17],ord_[MV],pos[MV],j[MV],K[MV];
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long xr(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double ur(){return (xr()>>11)*(1.0/9007199254740992.0);}
static inline int pen(int v,int x){return x<K[v]?K[v]-x:0;}
int full(){int c=0;for(int v=0;v<V;v++){int s=0;for(int k=0;k<d;k++)s+=pos[nb[v][k]]<pos[v];j[v]=s;c+=pen(v,s);}return c;}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r");fscanf(f,"%d %d",&V,&d);for(int v=0;v<V;v++)for(int k=0;k<d;k++)fscanf(f,"%d",&nb[v][k]);fclose(f);
  for(int v=0;v<V;v++)K[v]=v<V/2?8:9;
  rs^=atoll(argv[6])*0x9E3779B97F4A7C15ULL;xr();
  if(!strcmp(argv[2],"rand")){for(int i=0;i<V;i++)ord_[i]=i;for(int i=V-1;i>0;i--){int r=xr()%(i+1),t=ord_[i];ord_[i]=ord_[r];ord_[r]=t;}}
  else{f=fopen(argv[2],"r");for(int i=0;i<V;i++)fscanf(f,"%d",&ord_[i]);fclose(f);}
  for(int t=0;t<V;t++)pos[ord_[t]]=t;
  long long it_=atoll(argv[3]);double T0=atof(argv[4]),T1=atof(argv[5]);
  int cur=full(),best=cur,bo[MV];memcpy(bo,ord_,sizeof bo);printf("start %d\n",cur);
  for(long long it=0;it<it_;it++){
    double T=T0*pow(T1/T0,(double)it/it_);
    int v=xr()%V,p=pos[v];int span=(xr()&1)?V:30;int np=p+(int)(xr()%(2*span+1))-span;
    if(np<0||np>=V||np==p)continue;
    int lo=p<np?p:np,hi=p<np?np:p;int dl=0,jv=j[v];
    for(int k=0;k<d;k++){int u=nb[v][k],pu=pos[u];if(pu<lo||pu>hi)continue;
      if(np>p){ /* u was after v, now before: jv+1, ju-1 */ dl+=pen(u,j[u]-1)-pen(u,j[u]); jv++; }
      else { jv--; dl+=pen(u,j[u]+1)-pen(u,j[u]); } }
    dl+=pen(v,jv)-pen(v,j[v]);
    if(dl<=0||ur()<exp(-dl/T)){
      for(int k=0;k<d;k++){int u=nb[v][k],pu=pos[u];if(pu<lo||pu>hi)continue; if(np>p)j[u]--;else j[u]++;}
      j[v]=jv;
      if(np>p){for(int t=p;t<np;t++){ord_[t]=ord_[t+1];pos[ord_[t]]=t;}}else{for(int t=p;t>np;t--){ord_[t]=ord_[t-1];pos[ord_[t]]=t;}}
      ord_[np]=v;pos[v]=np;cur+=dl;
      if(cur<best){best=cur;memcpy(bo,ord_,sizeof bo);}
    }
    if((it&((1<<25)-1))==0){printf("it %lld T %.3f cur %d best %d\n",it,T,cur,best);fflush(stdout);}
  }
  for(int t=0;t<V;t++)pos[bo[t]]=t; printf("final best %d check %d\n",best,full());
  f=fopen(argv[7],"w");for(int t=0;t<V;t++)fprintf(f,"%d\n",bo[t]);fclose(f);return 0;}
