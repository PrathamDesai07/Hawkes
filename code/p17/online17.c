// Online projected Newton (EF or OI preconditioner) on a simulated stationary Hawkes truth (multi-exponential).
// usage: online17 cfg.txt out.bin
// cfg: nrep seed0 nev warm mu M [u w]*M K [beta]*K eps mu_min precond rho0 [theta0]*(K+1) ncp [cp]*ncp
// out per replicate, per checkpoint: theta(d), nactive_so_far, dist_to_face
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#define DM 6
static uint64_t st[2];
static uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static uint64_t nx(void){uint64_t s0=st[0],s1=st[1],r=s0+s1;s1^=s0;st[0]=rotl(s0,55)^s1^(s1<<14);st[1]=rotl(s1,36);return r;}
static double U(void){return ((nx()>>11)+0.5)*(1.0/9007199254740992.0);}
static int d,K; static double betas[8], cvec[DM], eps_, mumin;
static int solve_lin(int n,double A[][2*DM],double*bv,double*x){ // Gaussian elimination, A n x n, b n
  double M[2*DM][2*DM+1]; for(int i=0;i<n;i++){for(int j=0;j<n;j++)M[i][j]=A[i][j];M[i][n]=bv[i];}
  for(int c=0;c<n;c++){int p=c;for(int r=c+1;r<n;r++)if(fabs(M[r][c])>fabs(M[p][c]))p=r; if(fabs(M[p][c])<1e-14)return 0;
    if(p!=c)for(int j=0;j<=n;j++){double t=M[p][j];M[p][j]=M[c][j];M[c][j]=t;}
    for(int r=0;r<n;r++)if(r!=c){double f=M[r][c]/M[c][c];for(int j=c;j<=n;j++)M[r][j]-=f*M[c][j];}}
  for(int i=0;i<n;i++)x[i]=M[i][n]/M[i][i]; return 1;}
// constraints G x >= h : m rows
static int mcon; static double G[DM+2][DM], hh[DM+2];
static void build_cons(void){mcon=0; memset(G,0,sizeof G);
  G[mcon][0]=1;hh[mcon]=mumin;mcon++;
  for(int k=0;k<K;k++){G[mcon][1+k]=1;hh[mcon]=0;mcon++;}
  if(eps_>0){for(int k=0;k<K;k++)G[mcon][1+k]=-1.0/betas[k];hh[mcon]=-(1-eps_);mcon++;}}
static int feasible(const double*x,double tol){for(int j=0;j<mcon;j++){double s=0;for(int i=0;i<d;i++)s+=G[j][i]*x[i];if(s<hh[j]-tol)return 0;}return 1;}
// project y onto polytope in F-norm; returns 1 if changed
static int project(double F[][DM],const double*y,double*x){
  if(feasible(y,0)){memcpy(x,y,d*sizeof(double));return 0;}
  for(int sz=1;sz<=d;sz++)for(int mask=1;mask<(1<<mcon);mask++){
    if(__builtin_popcount(mask)!=sz)continue; int idx[DM+2],na=0; for(int j=0;j<mcon;j++)if(mask>>j&1)idx[na++]=j;
    int n=d+na; double A[2*DM][2*DM]; double bv[2*DM],sol[2*DM]; memset(A,0,sizeof A);
    for(int i=0;i<d;i++){for(int j=0;j<d;j++)A[i][j]=F[i][j]; double s=0;for(int j=0;j<d;j++)s+=F[i][j]*y[j]; bv[i]=s;
      for(int a=0;a<na;a++)A[i][d+a]=-G[idx[a]][i];}
    for(int a=0;a<na;a++){for(int i=0;i<d;i++)A[d+a][i]=G[idx[a]][i]; bv[d+a]=hh[idx[a]];}
    if(!solve_lin(n,A,bv,sol))continue; int ok=1; for(int a=0;a<na;a++)if(sol[d+a]<-1e-10)ok=0;
    if(ok&&feasible(sol,1e-10)){memcpy(x,sol,d*sizeof(double));return 1;}}
  memcpy(x,y,d*sizeof(double)); return 1;}
static void chol_solve(double F[][DM],const double*g,double*x){double L[DM][DM]={{0}};
  for(int i=0;i<d;i++)for(int j=0;j<=i;j++){double s=F[i][j];for(int k=0;k<j;k++)s-=L[i][k]*L[j][k];if(i==j)L[i][j]=sqrt(s>1e-300?s:1e-300);else L[i][j]=s/L[j][j];}
  double z[DM];for(int i=0;i<d;i++){double s=g[i];for(int k=0;k<i;k++)s-=L[i][k]*z[k];z[i]=s/L[i][i];}
  for(int i=d-1;i>=0;i--){double s=z[i];for(int k=i+1;k<d;k++)s-=L[k][i]*x[k];x[i]=s/L[i][i];}}
int main(int c,char**v){
  FILE*f=fopen(v[1],"r"); int nrep; long seed0,nev,warm; double mu; int M; double u[256],w[256],x[256];
  fscanf(f,"%d %ld %ld %ld %lf %d",&nrep,&seed0,&nev,&warm,&mu,&M); for(int j=0;j<M;j++)fscanf(f,"%lf %lf",&u[j],&w[j]);
  fscanf(f,"%d",&K); for(int k=0;k<K;k++)fscanf(f,"%lf",&betas[k]); d=K+1;
  int precond; double rho0,kappa,th0[DM]; fscanf(f,"%lf %lf %d %lf %lf",&eps_,&mumin,&precond,&rho0,&kappa); for(int i=0;i<d;i++)fscanf(f,"%lf",&th0[i]);
  int ncp; long cp[64]; fscanf(f,"%d",&ncp); for(int i=0;i<ncp;i++)fscanf(f,"%ld",&cp[i]); fclose(f);
  build_cons(); FILE*fo=fopen(v[2],"wb");
  double gx[10],gw[10]; { // Gauss-Legendre 10 on [-1,1] via Newton
    for(int i=0;i<10;i++){double z=cos(M_PI*(i+0.75)/10.5),pp; for(int it=0;it<100;it++){double p1=1,p2=0;for(int j=1;j<=10;j++){double p3=p2;p2=p1;p1=((2*j-1)*z*p2-(j-1)*p3)/j;}pp=10*(z*p1-p2)/(z*z-1);double z1=z;z=z1-p1/pp;if(fabs(z-z1)<1e-15)break;}
      double p1=1,p2=0;for(int j=1;j<=10;j++){double p3=p2;p2=p1;p1=((2*j-1)*z*p2-(j-1)*p3)/j;}pp=10*(z*p1-p2)/(z*z-1);gx[i]=z;gw[i]=2/((1-z*z)*pp*pp);}}
  for(int rep=0;rep<nrep;rep++){
    uint64_t seed=seed0+rep; st[0]=seed*0x9E3779B97F4A7C15ULL+1; st[1]=seed^0xD1B54A32D192ED03ULL; for(int i=0;i<20;i++)nx();
    double s[8],sp[8]; for(int j=0;j<M;j++)x[j]=0; for(int k=0;k<K;k++){s[k]=0;sp[k]=0;}
    double th[DM],F[DM][DM]; for(int i=0;i<d;i++){th[i]=th0[i];for(int j=0;j<d;j++)F[i][j]=(i==j)?rho0*kappa:0;}
    long n=0,na=0; double t=0,tprev=0; int icp=0; long ev=0;
    while(icp<ncp){
      double lam=mu; for(int j=0;j<M;j++)lam+=w[j]*x[j];
      double dt=-log(U())/lam; double lam2=mu; for(int j=0;j<M;j++)lam2+=w[j]*x[j]*exp(-u[j]*dt);
      t+=dt; for(int j=0;j<M;j++)x[j]*=exp(-u[j]*dt); for(int k=0;k<K;k++)s[k]*=exp(-betas[k]*dt);
      if(U()*lam<=lam2){
        double wg=t-tprev;
        if(ev>=warm){ // algorithm step using gap
          double a[DM],b[DM],g[DM],dF[DM][DM]; a[0]=1;b[0]=wg; for(int k=0;k<K;k++){a[1+k]=s[k];b[1+k]=sp[k]*(1-exp(-betas[k]*wg))/betas[k];}
          double lm=0;for(int i=0;i<d;i++)lm+=a[i]*th[i]; for(int i=0;i<d;i++)g[i]=b[i]-a[i]/lm;
          memset(dF,0,sizeof dF);
          if(precond==1){for(int i=0;i<d;i++)for(int j=0;j<d;j++)dF[i][j]=a[i]*a[j]/(lm*lm);}
          else{for(int q=0;q<10;q++){double uu=0.5*wg*(gx[q]+1),ww=0.5*wg*gw[q];double ph[DM];ph[0]=1;for(int k=0;k<K;k++)ph[1+k]=sp[k]*exp(-betas[k]*uu);
              double l=0;for(int i=0;i<d;i++)l+=ph[i]*th[i];for(int i=0;i<d;i++)for(int j=0;j<d;j++)dF[i][j]+=ww*ph[i]*ph[j]/l;}}
          for(int i=0;i<d;i++)for(int j=0;j<d;j++)F[i][j]+=kappa*dF[i][j];
          double step[DM],y[DM],xn[DM]; chol_solve(F,g,step); for(int i=0;i<d;i++)y[i]=th[i]-step[i];
          int ch=project(F,y,xn); na+=ch; for(int i=0;i<d;i++)th[i]=xn[i]; n++;
          if(n==cp[icp]){double rec[DM+2]; for(int i=0;i<d;i++)rec[i]=th[i]; rec[d]=(double)na; double dist=1-eps_; for(int k=0;k<K;k++)dist-=th[1+k]/betas[k]; rec[d+1]=dist; fwrite(rec,sizeof(double),d+2,fo); icp++;}
        }
        for(int j=0;j<M;j++)x[j]+=1; for(int k=0;k<K;k++){s[k]+=1;sp[k]=s[k];} tprev=t; ev++;
      }
    }
  }
  fclose(fo); return 0;}
