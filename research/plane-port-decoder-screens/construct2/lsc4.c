// Tabu local search, exact move evaluation from stored independent-set lists.
// Lists L[0]=T-2, L[1]=T-1, L[2]=T (independent sets of exactly that size, whole graph).
// f = W*N_T + W1*N_{T-1}. Move (a->c): new N_k = N_k - freq_k[a] + #{I in L_{k-1}: a notin I, I cap N_c(a) = empty}.
// usage: lsc2 ports.txt start.txt T W W1 seconds seed out.txt nport tlo thi stall kickfrac
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
typedef unsigned __int128 u128;
#define NP 112
#define NC 128
static int Pq[NP], RC[NP][NC], ch[NP];
static u128 O[NP][NC][NP]; // O[a][c][b] bit d
static u128 ALL, nadj[NP], adj[NP];
static inline int ctz(u128 x){uint64_t lo=(uint64_t)x; return lo?__builtin_ctzll(lo):64+__builtin_ctzll((uint64_t)(x>>64));}
static inline u128 B(int i){return (u128)1<<i;}
static int dotp(int a,int b){return __builtin_popcount(a&b)&1;}
static uint64_t rs_=88172645463325252ULL; static inline uint64_t rnd(){rs_^=rs_<<13;rs_^=rs_>>7;rs_^=rs_<<17;return rs_;}
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static u128 nbr(int a,int c){u128 m=0; for(int b=0;b<NP;b++) if(b!=a && ((O[a][c][b]>>ch[b])&1)) m|=B(b); return m;}
static void rebuild(){ for(int a=0;a<NP;a++){adj[a]=nbr(a,ch[a]); nadj[a]=ALL&~adj[a]&~B(a);} }
static long long cnt, cap; static int K; static u128 *store; static long long scap; static u128 base;
static void ex(u128 R,int rs,u128 P){
  if(cnt>=cap) return;
  if(rs==K){ if(store && cnt<scap) store[cnt]=R|base; cnt++; return; }
  int order[NP], col[NP], m=0; u128 U=P; int k=0;
  while(U){ k++; u128 Q=U; while(Q){int v=ctz(Q); Q&=~B(v); U&=~B(v); Q&=~nadj[v]; order[m]=v; col[m]=k; m++;} }
  for(int i=m-1;i>=0;i--){
    if(rs+col[i]<K) return;
    int v=order[i]; ex(R|B(v),rs+1,P&nadj[v]);
    if(cnt>=cap) return;
    P&=~B(v);
  }
}
static long long countk(u128 P,int k,long long cp,u128 *st,long long sc,u128 bs){
  cnt=0; cap=cp; K=k; store=st; scap=sc; base=bs; ex(0,0,P); return cnt;
}
static int best_;
static void mx(u128 R,int rs,u128 P){
  int order[NP], col[NP], m=0; u128 U=P; int k=0;
  while(U){ k++; u128 Q=U; while(Q){int v=ctz(Q); Q&=~B(v); U&=~B(v); Q&=~nadj[v]; order[m]=v; col[m]=k; m++;} }
  for(int i=m-1;i>=0;i--){
    if(rs+col[i]<=best_) return;
    int v=order[i]; u128 P2=P&nadj[v];
    if(!P2){ if(rs+1>best_) best_=rs+1; } else mx(R|B(v),rs+1,P2);
    P&=~B(v);
  }
}
static int alpha(){best_=0; mx(0,0,ALL); return best_;}
static long long MAXS=6000000; static int KCUR=0; static int KTOP=4; static long long CEILM=10, CEILA=200;
typedef struct { int k; u128 *S; long long ns; long long freq[NP]; } Lev;
static Lev L[3]; static long long Wt[3]; static int LO=0;
static void refreq(Lev*V){ memset(V->freq,0,sizeof(V->freq)); for(long long i=0;i<V->ns;i++){u128 x=V->S[i]; while(x){int v=ctz(x); x&=x-1; V->freq[v]++;}} }
static int full_enum(){
  for(int l=LO;l<3;l++){ long long c=countk(ALL,L[l].k,MAXS,L[l].S,MAXS,0); if(c>=MAXS) return 0; L[l].ns=c; refreq(&L[l]); }
  return 1;
}
static long long fval(){return Wt[2]*L[2].ns+Wt[1]*L[1].ns;}
// evaluate all candidates c of port a: delta[c]
static u128 Oa[NP]; // Oa[b] bit c: (a,c) orth (b,ch[b])
static void eval_port(int a,long long *delta){
  for(int b=0;b<NP;b++){ u128 m=0; if(b!=a) for(int c=0;c<NC;c++) if((O[a][c][b]>>ch[b])&1) m|=B(c); Oa[b]=m; }
  static long long cntc[2][NC];
  for(int l=0;l<2;l++){ memset(cntc[l],0,sizeof(cntc[l]));
    if(Wt[l+1]==0) continue;
    Lev*V=&L[l]; long long *cc=cntc[l];
    for(long long i=0;i<V->ns;i++){ u128 x=V->S[i]; if((x>>a)&1) continue;
      u128 kill=0; while(x){int v=ctz(x); x&=x-1; kill|=Oa[v];}
      u128 s=~kill; while(s){int c=ctz(s); s&=s-1; cc[c]++;} }
  }
  for(int c=0;c<NC;c++){ delta[c]=0; for(int l=0;l<2;l++) if(Wt[l+1]) delta[c]+=Wt[l+1]*(cntc[l][c]-L[l+1].freq[a]); }
}
static int apply_move(int a,int c){
  u128 Nn=nbr(a,c);
  for(int b=0;b<NP;b++){ u128 m=0; if(b!=a) for(int cc=0;cc<NC;cc++) if((O[a][cc][b]>>ch[b])&1) m|=B(cc); Oa[b]=m; }
  // top-down: levels 2,1 get additions from lower list (sets without a)
  for(int l=2;l>=LO;l--){
    Lev*V=&L[l]; long long j=0;
    for(long long i=0;i<V->ns;i++){ u128 x=V->S[i]; if(!((x>>a)&1)) V->S[j++]=x; }
    V->ns=j;
  }
  ch[a]=c; adj[a]=Nn; nadj[a]=ALL&~Nn&~B(a);
  for(int b=0;b<NP;b++) if(b!=a){ if((Nn>>b)&1){adj[b]|=B(a); nadj[b]&=~B(a);} else {adj[b]&=~B(a); nadj[b]|=B(a);} }
  for(int l=2;l>=LO+1;l--){
    Lev*V=&L[l], *U=&L[l-1]; long long n0=V->ns;
    for(long long i=0;i<U->ns;i++){ if(!(U->S[i]&Nn)){ if(V->ns>=MAXS) return 0; V->S[V->ns++]=U->S[i]|B(a);} }
    (void)n0;
  }
  { Lev*V=&L[LO]; long long room=MAXS-V->ns; long long c2=countk(ALL&~Nn&~B(a),V->k-1,room,V->S+V->ns,room,B(a)); if(c2>=room) return 0; V->ns+=c2; }
  for(int l=LO;l<3;l++) refreq(&L[l]);
  return 1;
}
int main(int argc,char**argv){
  if(argc<14){fprintf(stderr,"usage: ports start T W W1 secs seed out nport tlo thi stall kickfrac\n");return 1;}
  FILE*f=fopen(argv[1],"r"); for(int i=0;i<NP;i++) if(fscanf(f,"%d",&Pq[i])!=1) return 2; fclose(f);
  int r0[NP]; f=fopen(argv[2],"r"); for(int i=0;i<NP;i++) if(fscanf(f,"%d",&r0[i])!=1) return 3; fclose(f);
  int T=atoi(argv[3]); Wt[2]=atoll(argv[4]); Wt[1]=atoll(argv[5]); Wt[0]=0; double TL=atof(argv[6]);
  rs_^=(uint64_t)atoll(argv[7])*0x9E3779B97F4A7C15ULL; for(int i=0;i<10;i++) rnd();
  const char*outp=argv[8]; int nport=atoi(argv[9]), tlo=atoi(argv[10]), thi=atoi(argv[11]), stall=atoi(argv[12]); double kf=atof(argv[13]);
  if(getenv("MAXS")) MAXS=atoll(getenv("MAXS")); if(Wt[1]==0) LO=1; if(getenv("KTOP")) KTOP=atoi(getenv("KTOP")); if(getenv("KCUR")) KCUR=atoi(getenv("KCUR")); if(getenv("CEILM")) CEILM=atoll(getenv("CEILM")); if(getenv("CEILA")) CEILA=atoll(getenv("CEILA"));
  ALL=0; for(int i=0;i<NP;i++) ALL|=B(i);
  for(int a=0;a<NP;a++){int m=0; for(int r=0;r<512;r++) if(!dotp(r,Pq[a]) && (__builtin_popcount(r)&1)) RC[a][m++]=r; if(m!=NC) return 4;}
  for(int a=0;a<NP;a++) for(int c=0;c<NC;c++){int q=Pq[a],r=RC[a][c]; for(int b=0;b<NP;b++){u128 m=0; if(b!=a && !dotp(q,Pq[b])) for(int d=0;d<NC;d++){int s=RC[b][d]; if(!dotp(q,s)&&!dotp(r,Pq[b])&&!dotp(r,s)) m|=B(d);} O[a][c][b]=m;}}
  for(int a=0;a<NP;a++){int c; for(c=0;c<NC;c++) if(RC[a][c]==r0[a]) break; if(c==NC) return 5; ch[a]=c;}
  rebuild();
  for(int l=LO;l<3;l++){L[l].k=T-2+l; L[l].S=malloc(sizeof(u128)*MAXS);} L[0].ns=0;
  if(!full_enum()){fprintf(stderr,"too many sets at start\n");return 6;}
  long long fcur=fval(), fbest=fcur; int bestch[NP]; memcpy(bestch,ch,sizeof ch);
  long long *tabu=calloc(NP*NC,sizeof(long long)); long long it=0, lastimp=0; double t0=now(), tlast=t0; int kicks=0;
  printf("start alpha %d f %lld N%d %lld N%d %lld N%d %lld\n",alpha(),fcur,T,L[2].ns,T-1,L[1].ns,T-2,L[0].ns); fflush(stdout);
  long long delta[NC];
  while(now()-t0<TL){
    it++;
    long long bd=(long long)4e18; int ba=-1,bc=-1, nties=0;
    long long tot=0; long long wt[NP]; for(int a=0;a<NP;a++){wt[a]=L[2].freq[a]*16+L[1].freq[a]/8+1; tot+=wt[a];}
    for(int s=0;s<nport;s++){
      int a;
      if(nport>=NP) a=s; else { long long x=rnd()%tot; a=0; while(x>=wt[a]){x-=wt[a]; a++;} }
      eval_port(a,delta);
      for(int c=0;c<NC;c++){ if(c==ch[a]) continue;
        long long d=delta[c]; int tb=tabu[a*NC+c]>it; if(tb && !(fcur+d<fbest)) continue;
        if(d<bd){bd=d; ba=a; bc=c; nties=1;} else if(d==bd){ nties++; if(rnd()%nties==0){ba=a; bc=c;} } }
    }
    if(ba<0) continue;
    int oc=ch[ba];
    if(!apply_move(ba,bc)){ memcpy(ch,bestch,sizeof ch); rebuild(); full_enum(); fcur=fval(); continue; }
    long long fn=fval(); if(fn!=fcur+bd){fprintf(stderr,"delta mismatch %lld %lld\n",fn,fcur+bd);} fcur=fn;
    tabu[ba*NC+oc]=it+tlo+rnd()%(thi-tlo+1);
    if(fcur<fbest){fbest=fcur; memcpy(bestch,ch,sizeof ch); lastimp=it;
      FILE*o=fopen(outp,"w"); for(int a=0;a<NP;a++) fprintf(o,"%d ",RC[a][ch[a]]); fprintf(o,"\n"); fclose(o);
      printf("%.0fs it %lld NEW BEST f %lld N%d %lld N%d %lld\n",now()-t0,it,fbest,T,L[2].ns,T-1,L[1].ns); fflush(stdout);
      if(L[2].ns==0){ printf("TARGET REACHED: alpha <= %d (alpha=%d)\n",T-1,alpha()); fflush(stdout); }
    }
    if(fcur>CEILM*fbest+CEILA) lastimp=it-stall-1;
    if(it-lastimp>stall){
      if(rnd()%100>=KCUR || fcur>2*fbest+20){memcpy(ch,bestch,sizeof ch); rebuild(); full_enum();}
      u128 sset = L[2].ns>0 ? L[2].S[rnd()%L[2].ns] : L[1].S[rnd()%L[1].ns];
      memset(tabu,0,NP*NC*sizeof(long long));
      { int vs[NP], nv=0; u128 x=sset; while(x){int v=ctz(x); x&=x-1; vs[nv++]=v;}
        int nk=0;
        for(int i=0;i<nv;i++){ if((rnd()%1000)>=kf*1000 && !(i==nv-1 && nk==0)) continue;
          int v=vs[i]; eval_port(v,delta); int top[KTOP]; long long tv[KTOP]; int nt=0;
          for(int c=0;c<NC;c++){ if(c==ch[v]) continue; long long d=delta[c]; int full=(nt==KTOP); if(full && d>=tv[KTOP-1]) continue; int j=full?KTOP-1:nt++; while(j>0 && tv[j-1]>d){tv[j]=tv[j-1]; top[j]=top[j-1]; j--;} tv[j]=d; top[j]=c; }
          int c=top[rnd()%nt]; int oc=ch[v];
          if(!apply_move(v,c)){ memcpy(ch,bestch,sizeof ch); rebuild(); full_enum(); break; }
          tabu[v*NC+oc]=it+thi; nk++; }
      }
      fcur=fval(); lastimp=it; kicks++;
      printf("%.0fs it %lld kick %d f %lld N%d %lld (best %lld)\n",now()-t0,it,kicks,fcur,T,L[2].ns,fbest); fflush(stdout);
    }
    if(now()-tlast>30){tlast=now(); printf("%.0fs it %lld f %lld N%d %lld N%d %lld N%d %lld best %lld\n",now()-t0,it,fcur,T,L[2].ns,T-1,L[1].ns,T-2,L[0].ns,fbest); fflush(stdout);}
  }
  printf("final best f %lld\n",fbest);
  return 0;
}
