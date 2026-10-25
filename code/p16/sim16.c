// Ogata thinning for a stationary-ish multi-exponential Hawkes truth; also tracks model states S_k (rates beta_k) from the observed path.
// usage: sim16 seed nev warm mu M u1 w1 ... uM wM K b1..bK out.bin
// out record per event i (after warm): w, lamstar(t_i-), s_1..s_K  where s_k = S_k(t_{i-1}+) (state just after previous event)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static uint64_t st[2];
static uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static uint64_t nx(void){uint64_t s0=st[0],s1=st[1],r=s0+s1;s1^=s0;st[0]=rotl(s0,55)^s1^(s1<<14);st[1]=rotl(s1,36);return r;}
static double U(void){return ((nx()>>11)+0.5)*(1.0/9007199254740992.0);}
int main(int c,char**v){
  int a=1; uint64_t seed=strtoull(v[a++],0,10); long nev=atol(v[a++]); long warm=atol(v[a++]); double mu=atof(v[a++]);
  int M=atoi(v[a++]); double u[256],w[256],x[256]; for(int j=0;j<M;j++){u[j]=atof(v[a++]);w[j]=atof(v[a++]);x[j]=0;}
  int K=atoi(v[a++]); double b[16],s[16]; for(int k=0;k<K;k++){b[k]=atof(v[a++]);s[k]=0;}
  FILE*f=fopen(v[a],"wb");
  st[0]=seed*0x9E3779B97F4A7C15ULL+1; st[1]=seed^0xD1B54A32D192ED03ULL; for(int i=0;i<20;i++)nx();
  double t=0; long n=0; double tprev=0;
  while(n<nev+warm){
    double lam=mu; for(int j=0;j<M;j++) lam+=w[j]*x[j];
    double dt=-log(U())/lam; // candidate, upper bound = current lam (nonincreasing between events)
    double lam2=mu; for(int j=0;j<M;j++) lam2+=w[j]*x[j]*exp(-u[j]*dt);
    t+=dt; for(int j=0;j<M;j++) x[j]*=exp(-u[j]*dt); for(int k=0;k<K;k++) s[k]*=exp(-b[k]*dt);
    if(U()*lam<=lam2){ // event at t
      if(n>=warm){ double rec[20]; // s currently = S_k(t-) ; need S_k(tprev+) -> recompute
        // store w and lam2 and S_k(t-) (we convert later: s(tprev+) = S(t-) e^{+b w})
        rec[0]=t-tprev; rec[1]=lam2; for(int k=0;k<K;k++) rec[2+k]=s[k]*exp(b[k]*(t-tprev));
        fwrite(rec,sizeof(double),2+K,f);}
      for(int j=0;j<M;j++) x[j]+=1; for(int k=0;k<K;k++) s[k]+=1; n++; tprev=t;
    }
  }
  fclose(f); return 0;}
