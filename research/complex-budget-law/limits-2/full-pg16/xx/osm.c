/* SA over (order sigma, orientation O) with fixed in-degrees: objective = #back arcs |R|.
   moves: vertex relocation (O fixed), directed-cycle reversal (len 6/8/10; in-degrees preserved).
   usage: osm graph seedfile iters T0 T1 seed out mu0 mu1  (objective #back + mu*sum|indeg-k|, k from seed orientation)   */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MV 600
#define MD 17
int V,d,nb[MV][MD],dir[MV][MD]; /* dir[v][i]=1 if arc v->nb[v][i] */
int idx[MV][MD]; /* index of v in nb[nb[v][i]] */
int pos[MV],ord_[MV];
static unsigned long long rs=88172645463325252ULL;
static inline unsigned long long xr(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double ur(){return (xr()>>11)*(1.0/9007199254740992.0);}
int backcount(){int c=0;for(int v=0;v<V;v++)for(int i=0;i<d;i++)if(dir[v][i]&&pos[v]>pos[nb[v][i]])c++;return c;}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); fscanf(f,"%d %d",&V,&d);
  for(int v=0;v<V;v++)for(int i=0;i<d;i++)fscanf(f,"%d",&nb[v][i]); fclose(f);
  for(int v=0;v<V;v++)for(int i=0;i<d;i++){int u=nb[v][i];for(int k=0;k<d;k++)if(nb[u][k]==v)idx[v][i]=k;}
  f=fopen(argv[2],"r"); for(int t=0;t<V;t++){fscanf(f,"%d",&ord_[t]);pos[ord_[t]]=t;}
  int ne;fscanf(f,"%d",&ne); memset(dir,0,sizeof dir);
  for(int e=0;e<ne;e++){int a,b;fscanf(f,"%d %d",&a,&b);for(int i=0;i<d;i++)if(nb[a][i]==b)dir[a][i]=1;} fclose(f);
  long long iters=atoll(argv[3]); double T0=atof(argv[4]),T1=atof(argv[5]); rs^=atoll(argv[6])*0x9E3779B97F4A7C15ULL; xr();
  int cur=backcount(),best=cur; printf("start %d\n",cur); fflush(stdout);
  int bord[MV]; int bdir[MV][MD]; memcpy(bord,ord_,sizeof ord_); memcpy(bdir,dir,sizeof dir);
  int indeg0[MV]; for(int v=0;v<V;v++){indeg0[v]=0;for(int i=0;i<d;i++)if(!dir[v][i])indeg0[v]++;}
  long long acc=0; double mu0=atof(argv[8]),mu1=atof(argv[9]); int indeg[MV]; memcpy(indeg,indeg0,sizeof indeg0);
  int imb=0; double obj=cur; int bestfeas=cur;
  for(long long it=0;it<iters;it++){
    double T=T0*pow(T1/T0,(double)it/iters); double mu=mu0+(mu1-mu0)*(double)it/iters;
    if((xr()&7)==0){ /* single arc reversal a->b becomes b->a: indeg b -1, indeg a +1 */
      int a=xr()%V,i=xr()%d; if(!dir[a][i])continue; int b=nb[a][i];
      int db=(pos[b]>pos[a])?0:0; int wasback=pos[a]>pos[b]; int dback=wasback?-1:1;
      int di=0; di+=abs(indeg[b]-1-indeg0[b])-abs(indeg[b]-indeg0[b]); di+=abs(indeg[a]+1-indeg0[a])-abs(indeg[a]-indeg0[a]);
      double dl=dback+mu*di; (void)db;
      if(dl<=0||ur()<exp(-dl/T)){dir[a][i]=0;dir[b][idx[a][i]]=1;indeg[b]--;indeg[a]++;cur+=dback;imb+=di;}
      if(imb==0&&cur<bestfeas){bestfeas=cur;memcpy(bord,ord_,sizeof ord_);memcpy(bdir,dir,sizeof dir);if(cur<=906){printf("HIT %d\n",cur);fflush(stdout);}}
      continue;}
    if((xr()&3)!=0){ /* relocation */
      int v=xr()%V,p=pos[v]; int span=(xr()&1)?V:40; int np=p+(int)(xr()%(2*span+1))-span;
      if(np<0||np>=V||np==p)continue;
      int lo=p<np?p:np,hi=p<np?np:p; int dl=0;
      for(int i=0;i<d;i++){int u=nb[v][i],pu=pos[u]; if(pu<lo||pu>hi)continue;
        /* order between v,u flips */ int back_before = dir[v][i]? (p>pu):(pu>p); dl+= back_before?-1:1;}
      if(dl<=0||ur()<exp(-dl/T)){
        if(np>p){for(int t=p;t<np;t++){ord_[t]=ord_[t+1];pos[ord_[t]]=t;}}
        else{for(int t=p;t>np;t--){ord_[t]=ord_[t-1];pos[ord_[t]]=t;}}
        ord_[np]=v;pos[v]=np;cur+=dl;acc++;
      }
    } else { /* cycle reversal through a back arc */
      int u=xr()%V,i=xr()%d; if(!dir[u][i])continue; int v=nb[u][i]; if(pos[u]<pos[v])continue; /* u->v back */
      int L=6+2*(xr()%3); int path[12],pi[12]; path[0]=u;pi[0]=i; int x=v,ok=1;
      for(int s=1;s<L-1;s++){ int cnt=0,ch[MD]; for(int k=0;k<d;k++)if(dir[x][k])ch[cnt++]=k; if(!cnt){ok=0;break;} int k=ch[xr()%cnt]; path[s]=x;pi[s]=k;x=nb[x][k];}
      if(!ok)continue; int kk=-1; for(int k=0;k<d;k++)if(nb[x][k]==u&&dir[x][k])kk=k; if(kk<0)continue;
      path[L-1]=x;pi[L-1]=kk;
      /* simple check */ for(int a=0;a<L&&ok;a++)for(int b=a+1;b<L;b++)if(path[a]==path[b]){ok=0;break;} if(!ok)continue;
      int bb=0; for(int s=0;s<L;s++){int a=path[s],b=nb[a][pi[s]]; if(pos[a]>pos[b])bb++;}
      int dl=L-2*bb;
      if(dl<=0||ur()<exp(-dl/T)){
        for(int s=0;s<L;s++){int a=path[s],k=pi[s],b=nb[a][k]; dir[a][k]=0; dir[b][idx[a][k]]=1;}
        cur+=dl;acc++;
      }
    }
    if(imb==0&&cur<bestfeas){bestfeas=cur;best=cur;memcpy(bord,ord_,sizeof ord_);memcpy(bdir,dir,sizeof dir);
      if(best<=906){printf("HIT %d at it %lld\n",best,it);fflush(stdout);} }
    if((it&((1<<24)-1))==0){printf("it %lld T %.3f mu %.3f cur %d imb %d obj %.1f bestfeas %d\n",it,T,mu,cur,imb,cur+mu*imb,bestfeas);fflush(stdout);}
  }
  memcpy(ord_,bord,sizeof ord_);memcpy(dir,bdir,sizeof dir);for(int t=0;t<V;t++)pos[ord_[t]]=t;
  int chk=backcount(); for(int v=0;v<V;v++){int in=0;for(int i=0;i<d;i++)if(!dir[v][i])in++; if(in!=indeg0[v]){printf("INDEG FAIL\n");return 1;}}
  printf("final best %d (check %d)\n",best,chk);
  f=fopen(argv[7],"w"); for(int t=0;t<V;t++)fprintf(f,"%d ",ord_[t]); fprintf(f,"\n%d\n",ne);
  for(int v=0;v<V;v++)for(int i=0;i<d;i++)if(dir[v][i])fprintf(f,"%d %d\n",v,nb[v][i]); fclose(f);
  return 0;
}
