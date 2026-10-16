// SA over linear orders, cost = sum_v (8 - j_v)^+, insertion moves (windowed). osa graph init(file of order|S-file via 'S:'path) seed secs T0 T1 W out
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MAXV 8400
#define MAXD 64
int DEG,K;
int V, nb[MAXV][MAXD], ord[MAXV], pos[MAXV], j[MAXV];
unsigned long long rs;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double urand(){return (rnd()>>11)*(1.0/9007199254740992.0);}
static inline int pen(int x){return x<K?K-x:0;}
int full(){int c=0; for(int v=0;v<V;v++){int s=0; for(int k=0;k<DEG;k++) s+=pos[nb[v][k]]<pos[v]; j[v]=s; c+=pen(s);} return c;}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); fscanf(f,"%d %d",&V,&DEG); K=DEG/2; for(int v=0;v<V;v++) for(int k=0;k<DEG;k++) fscanf(f,"%d",&nb[v][k]); fclose(f);
  rs=88172645463325252ULL ^ (unsigned long long)atoll(argv[3])*0x9E3779B97F4A7C15ULL; for(int i=0;i<10;i++) rnd();
  if(strcmp(argv[2],"rand")==0){ for(int i=0;i<V;i++) ord[i]=i; for(int i=V-1;i>0;i--){int r=rnd()%(i+1); int t=ord[i];ord[i]=ord[r];ord[r]=t;} }
  else { FILE*g=fopen(argv[2],"r"); for(int i=0;i<V;i++) fscanf(g,"%d",&ord[i]); fclose(g);} 
  for(int i=0;i<V;i++) pos[ord[i]]=i;
  double secs=atof(argv[4]), T0=atof(argv[5]), T1=atof(argv[6]); int W=atoi(argv[7]);
  int cur=full(), best=cur; fprintf(stderr,"start %d\n",cur);
  clock_t t0=clock(); long it=0; double T=T0;
  while(1){
    it++; if((it&4095)==0){ double el=(double)(clock()-t0)/CLOCKS_PER_SEC; if(el>secs) break; T=T0*pow(T1/T0, el/secs);}
    int v=rnd()%V; int a=pos[v]; int b;
    if(rnd()%8==0){ int u=nb[v][rnd()%DEG]; b=pos[u]+(int)(rnd()%3)-1; } else b=a+(int)(rnd()%(2*W+1))-W;
    if(b<0) b=0; if(b>=V) b=V-1; if(b==a) continue;
    int d=0, jv=j[v];
    // neighbours u with pos in between change
    int dj=0;
    for(int k=0;k<DEG;k++){int u=nb[v][k]; int p=pos[u];
      if(a<b){ if(p>a && p<=b){ d+=pen(j[u]-1)-pen(j[u]); dj++; } }
      else { if(p>=b && p<a){ d+=pen(j[u]+1)-pen(j[u]); dj--; } } }
    d+=pen(jv+dj)-pen(jv);
    if(d<=0 || urand()<exp(-d/T)){
      for(int k=0;k<DEG;k++){int u=nb[v][k]; int p=pos[u]; if(a<b){ if(p>a&&p<=b) j[u]--; } else { if(p>=b&&p<a) j[u]++; } }
      j[v]=jv+dj;
      if(a<b){ for(int i=a;i<b;i++){ord[i]=ord[i+1]; pos[ord[i]]=i;} } else { for(int i=a;i>b;i--){ord[i]=ord[i-1]; pos[ord[i]]=i;} }
      ord[b]=v; pos[v]=b; cur+=d;
      if(cur<best){ best=cur; fprintf(stderr,"it %ld best %d T %.3f\n",it,best,T); FILE*o=fopen(argv[8],"w"); for(int i=0;i<V;i++) fprintf(o,"%d\n",ord[i]); fclose(o);} }
  }
  int chk=full(); fprintf(stderr,"final best %d cur %d chk %d iters %ld\n",best,cur,chk,it);
  return 0;
}
