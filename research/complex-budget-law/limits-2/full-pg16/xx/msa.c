/* ILS on exact F(sigma) = min #back arcs over orientations with in-degree K (8 points, 9 lines), PG(2,16).
   round: kick (relocate nk random vertices), quench with joint (order,orientation) SA at low T (fixed in-degrees),
   then exact primal-dual min-cost flow re-optimises the orientation for the new order. accept if F <= cur (or rarely worse).
   usage: ils graph seedfile rounds quenchIters Tq nk seed out */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MV 560
#define MD 17
int V,d,nb[MV][MD],idx_[MV][MD],K[MV];
int pos[MV],ord_[MV],dir[MV][MD];
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long xr(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double ur(){return (xr()>>11)*(1.0/9007199254740992.0);}
/* ---------- min-cost flow (primal-dual, unit capacities on edge arcs) ---------- */
#define MN (MV+2)
#define MA (2*(MV*MD/2)+4*MV+10)
int hd[MN],nx[2*MA],to[2*MA],cp[2*MA],cs[2*MA],na;
void addarc(int a,int b,int c,int w){to[na]=b;cp[na]=c;cs[na]=w;nx[na]=hd[a];hd[a]=na++;to[na]=a;cp[na]=0;cs[na]=-w;nx[na]=hd[b];hd[b]=na++;}
int dist[MN],inq[MN],qu[MN*64],it_[MN],edgearc[MV][MD];
int vis[MN];
int dfs(int v,int t,int f){ if(v==t)return f; vis[v]=1; for(int *p=&it_[v];*p!=-1;*p=nx[*p]){int a=*p,b=to[a]; if(cp[a]>0&&!vis[b]&&dist[b]==dist[v]+cs[a]){int g=dfs(b,t,f<cp[a]?f:cp[a]); if(g){cp[a]-=g;cp[a^1]+=g;return g;}}} return 0;}
int Fopt(int setdir){ /* returns F for current pos; if setdir, writes optimal orientation into dir */
  int s=V,t=V+1,n=V+2; na=0; for(int v=0;v<n;v++)hd[v]=-1; int j[MV];
  for(int v=0;v<V;v++){j[v]=0;for(int i=0;i<d;i++)if(pos[nb[v][i]]<pos[v])j[v]++;}
  for(int v=0;v<V;v++)for(int i=0;i<d;i++){int u=nb[v][i]; if(pos[u]<pos[v]){ /* forward u->v; flip = flow v->u */ edgearc[v][i]=na; addarc(v,u,1,1);} else edgearc[v][i]=-1;}
  int tot=0,base=0;
  for(int v=0;v<V;v++){int sp=j[v]-K[v]; if(sp>0){addarc(s,v,sp,0);tot+=sp;} else if(sp<0) addarc(v,t,-sp,0);}
  int flow=0,cost=0;
  while(flow<tot){
    for(int v=0;v<n;v++){dist[v]=1<<29;inq[v]=0;} dist[s]=0; int qh=0,qt=0; qu[qt++]=s; inq[s]=1;
    while(qh!=qt){int v=qu[qh];qh=(qh+1)%(MN*64);inq[v]=0; for(int a=hd[v];a!=-1;a=nx[a]) if(cp[a]>0&&dist[v]+cs[a]<dist[to[a]]){dist[to[a]]=dist[v]+cs[a]; if(!inq[to[a]]){inq[to[a]]=1;qu[qt]=to[a];qt=(qt+1)%(MN*64);}}}
    if(dist[t]>=1<<29){printf("infeasible\n");exit(1);}
    for(int v=0;v<n;v++)it_[v]=hd[v];
    int f; while(1){memset(vis,0,sizeof(int)*n); f=dfs(s,t,1<<29); if(!f)break; flow+=f;cost+=f*dist[t];}
  }
  (void)base;
  if(setdir){ for(int v=0;v<V;v++)for(int i=0;i<d;i++){int u=nb[v][i]; if(pos[u]<pos[v]){int a=edgearc[v][i]; int flipped=(cp[a]==0); /* forward u->v unless flipped */
        dir[v][i]=flipped?1:0; dir[u][idx_[v][i]]=flipped?0:1;}}}
  return cost;
}
int backcount(){int c=0;for(int v=0;v<V;v++)for(int i=0;i<d;i++)if(dir[v][i]&&pos[v]>pos[nb[v][i]])c++;return c;}
void relocate(int v,int np){int p=pos[v]; if(np>p){for(int t=p;t<np;t++){ord_[t]=ord_[t+1];pos[ord_[t]]=t;}}else{for(int t=p;t>np;t--){ord_[t]=ord_[t-1];pos[ord_[t]]=t;}} ord_[np]=v;pos[v]=np;}
int quench(long long iters,double T,int cur){
  for(long long it=0;it<iters;it++){
    if((xr()&3)!=0){int v=xr()%V,p=pos[v];int span=(xr()&1)?V:40;int np=p+(int)(xr()%(2*span+1))-span; if(np<0||np>=V||np==p)continue;
      int lo=p<np?p:np,hi=p<np?np:p,dl=0; for(int i=0;i<d;i++){int u=nb[v][i],pu=pos[u]; if(pu<lo||pu>hi)continue; int bb=dir[v][i]?(p>pu):(pu>p); dl+=bb?-1:1;}
      if(dl<=0||ur()<exp(-dl/T)){relocate(v,np);cur+=dl;}
    } else { int u=xr()%V,i=xr()%d; if(!dir[u][i])continue; int v=nb[u][i]; if(pos[u]<pos[v])continue;
      int L=6+2*(xr()%3),path[12],pi[12],x=v,ok=1; path[0]=u;pi[0]=i;
      for(int s=1;s<L-1;s++){int cnt=0,ch[MD];for(int k=0;k<d;k++)if(dir[x][k])ch[cnt++]=k; if(!cnt){ok=0;break;} int k=ch[xr()%cnt];path[s]=x;pi[s]=k;x=nb[x][k];}
      if(!ok)continue; int kk=-1; for(int k=0;k<d;k++)if(nb[x][k]==u&&dir[x][k])kk=k; if(kk<0)continue; path[L-1]=x;pi[L-1]=kk;
      for(int a=0;a<L&&ok;a++)for(int b=a+1;b<L;b++)if(path[a]==path[b]){ok=0;break;} if(!ok)continue;
      int bb=0; for(int s=0;s<L;s++){int a=path[s],b=nb[a][pi[s]]; if(pos[a]>pos[b])bb++;} int dl=L-2*bb;
      if(dl<=0||ur()<exp(-dl/T)){for(int s=0;s<L;s++){int a=path[s],k=pi[s],b=nb[a][k];dir[a][k]=0;dir[b][idx_[a][k]]=1;} cur+=dl;}
    }
  }
  return cur;
}

/* merge-path SA: points keep their seed relative order, lines keep theirs; moves swap adjacent point/line.
   objective exact F.  usage: msa graph seedorder iters T0 T1 seed out */
int incident(int a,int b){for(int i=0;i<d;i++)if(nb[a][i]==b)return 1;return 0;}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r");fscanf(f,"%d %d",&V,&d);for(int v=0;v<V;v++)for(int i=0;i<d;i++)fscanf(f,"%d",&nb[v][i]);fclose(f);
  for(int v=0;v<V;v++){K[v]=v<V/2?8:9;for(int i=0;i<d;i++){int u=nb[v][i];for(int k=0;k<d;k++)if(nb[u][k]==v)idx_[v][i]=k;}}
  f=fopen(argv[2],"r");for(int t=0;t<V;t++){fscanf(f,"%d",&ord_[t]);pos[ord_[t]]=t;}fclose(f);
  long long iters=atoll(argv[3]);double T0=atof(argv[4]),T1=atof(argv[5]);rs^=atoll(argv[6])*0x9E3779B97F4A7C15ULL;xr();
  int cur=Fopt(0),best=cur,bo[MV];memcpy(bo,ord_,sizeof bo);printf("start F %d\n",cur);fflush(stdout);
  long long evals=0;
  for(long long it=0;it<iters;it++){
    double T=T0*pow(T1/T0,(double)it/iters);
    int t=xr()%(V-1),a=ord_[t],b=ord_[t+1]; if((a<V/2)==(b<V/2))continue;
    ord_[t]=b;ord_[t+1]=a;pos[a]=t+1;pos[b]=t;
    if(!incident(a,b))continue;
    int F=Fopt(0);evals++;
    if(F<=cur||ur()<exp(-(F-cur)/T)){cur=F; if(cur<best){best=cur;memcpy(bo,ord_,sizeof bo);printf("it %lld best %d\n",it,best);fflush(stdout);
        FILE*g=fopen(argv[7],"w");for(int s=0;s<V;s++)fprintf(g,"%d\n",bo[s]);fclose(g);}}
    else {ord_[t]=a;ord_[t+1]=b;pos[a]=t;pos[b]=t+1;}
    if(evals%20000==0){printf("evals %lld T %.3f cur %d best %d\n",evals,T,cur,best);fflush(stdout);}
  }
  printf("final best %d evals %lld\n",best,evals);return 0;}
