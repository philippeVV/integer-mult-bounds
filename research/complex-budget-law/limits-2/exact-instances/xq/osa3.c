// SA over orders, objective = cost + mu * sum_t (c_t - TH)^+, c_t = cut after position t. args: graph init seed secs T0 T1 W mu TH out
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define MAXV 8400
#define MAXD 64
int V,DEG,K, nb[MAXV][MAXD], ord[MAXV], pos[MAXV], j[MAXV], c[MAXV], tmpj[MAXV];
unsigned long long rs;
static inline unsigned long long rnd(){rs^=rs<<13;rs^=rs>>7;rs^=rs<<17;return rs;}
static inline double urand(){return (rnd()>>11)*(1.0/9007199254740992.0);}
static inline int pen(int x){return x<K?K-x:0;}
double mu; int TH;
static inline double over(int x){return x>TH? (double)(x-TH):0.0;}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); fscanf(f,"%d %d",&V,&DEG); K=DEG/2; for(int v=0;v<V;v++) for(int k=0;k<DEG;k++) fscanf(f,"%d",&nb[v][k]); fclose(f);
  rs=88172645463325252ULL ^ (unsigned long long)atoll(argv[3])*0x9E3779B97F4A7C15ULL; for(int i=0;i<10;i++) rnd();
  FILE*g=fopen(argv[2],"r"); for(int i=0;i<V;i++) fscanf(g,"%d",&ord[i]); fclose(g);
  for(int i=0;i<V;i++) pos[ord[i]]=i;
  double secs=atof(argv[4]), T0=atof(argv[5]), T1=atof(argv[6]); int W=atoi(argv[7]); mu=atof(argv[8]); TH=atoi(argv[9]);
  int cost=0; for(int v=0;v<V;v++){int s=0; for(int k=0;k<DEG;k++) s+=pos[nb[v][k]]<pos[v]; j[v]=s; cost+=pen(s);} 
  int acc=0; double ov=0; for(int t=0;t<V;t++){acc+=DEG-2*j[ord[t]]; c[t]=acc; ov+=over(acc);} 
  double cur=cost+mu*ov; double best=1e18; int bestcost=0, bestcw=0;
  fprintf(stderr,"start cost %d over %.0f\n",cost,ov);
  clock_t t0=clock(); long it=0; double T=T0; int seg[MAXV];
  while(1){
    it++; if((it&4095)==0){ double el=(double)(clock()-t0)/CLOCKS_PER_SEC; if(el>secs) break; T=T0*pow(T1/T0, el/secs);}
    int v=rnd()%V; int a=pos[v]; int b;
    if(rnd()%8==0){ int u=nb[v][rnd()%DEG]; b=pos[u]+(int)(rnd()%3)-1; } else b=a+(int)(rnd()%(2*W+1))-W;
    if(b<0) b=0; if(b>=V) b=V-1; if(b==a) continue;
    int lo=a<b?a:b, hi=a<b?b:a;
    // build new segment order and j's
    int dcost=0, dj=0;
    for(int k=0;k<DEG;k++){int u=nb[v][k]; int p=pos[u];
      if(a<b){ if(p>a && p<=b){ dcost+=pen(j[u]-1)-pen(j[u]); dj++; } }
      else { if(p>=b && p<a){ dcost+=pen(j[u]+1)-pen(j[u]); dj--; } } }
    dcost+=pen(j[v]+dj)-pen(j[v]);
    double dov=0;
    if(mu>0){
      int L=0; if(a<b){ for(int i=a+1;i<=b;i++) seg[L++]=ord[i]; seg[L++]=v; } else { seg[L++]=v; for(int i=b;i<a;i++) seg[L++]=ord[i]; }
      int base= lo>0? c[lo-1]:0; int x=base;
      for(int i=0;i<L;i++){ int w=seg[i]; int jw=j[w];
        if(w==v) jw=j[v]+dj; else { int adj=0; for(int k=0;k<DEG;k++) if(nb[w][k]==v){adj=1;break;} if(adj) jw+= (a<b)? -1: +1; }
        x+=DEG-2*jw; dov+=over(x)-over(c[lo+i]); }
    }
    double d=dcost+mu*dov;
    if(d<=0 || urand()<exp(-d/T)){
      for(int k=0;k<DEG;k++){int u=nb[v][k]; int p=pos[u]; if(a<b){ if(p>a&&p<=b) j[u]--; } else { if(p>=b&&p<a) j[u]++; } }
      j[v]+=dj;
      if(a<b){ for(int i=a;i<b;i++){ord[i]=ord[i+1]; pos[ord[i]]=i;} } else { for(int i=a;i>b;i--){ord[i]=ord[i-1]; pos[ord[i]]=i;} }
      ord[b]=v; pos[v]=b; cost+=dcost; ov+=dov; cur+=d;
      int x= lo>0? c[lo-1]:0; for(int i=lo;i<=hi;i++){ x+=DEG-2*j[ord[i]]; c[i]=x; }
      if(cur<best-1e-9){ best=cur; int cw=0; for(int t=0;t<V;t++) if(c[t]>cw) cw=c[t];
        fprintf(stderr,"it %ld obj %.1f cost %d cw %d over %.0f T %.3f\n",it,best,cost,cw,ov,T); FILE*o=fopen(argv[10],"w"); for(int i=0;i<V;i++) fprintf(o,"%d\n",ord[i]); fclose(o);} }
  }
  int chk=0; for(int v=0;v<V;v++){int s=0; for(int k=0;k<DEG;k++) s+=pos[nb[v][k]]<pos[v]; chk+=pen(s);} 
  fprintf(stderr,"final best %.1f cur cost %d chk %d iters %ld\n",best,cost,chk,it);
  return 0;
}
