#define EPS 1e-15

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

namespace nummat{
   typedef vector<double> VECTOR;
   typedef vector<vector<double>> MATRIX;
   void read(VECTOR& v){
      int n=v.size();
      for(int i=0;i<n;i++){
         scanf("%lf",&v[i]);
      }
   }
   void read(MATRIX& A){
      int n=A.size();
      for(int i=0;i<n;i++){
         read(A[i]);
      }
   }
   void print(VECTOR& v){
      int n=v.size();
      for(int i=0;i<n;i++){
         printf("%lf ",v[i]);
      }
      printf("\n");
   }
   void print(MATRIX& A){
      int n=A.size();
      for(int i=0;i<n;i++){
         print(A[i]);
      }
      printf("\n");
   }

   void mul(VECTOR& x,double C){
      int n=x.size();
      for(int i=0;i<n;i++){
         x[i]*=C;
      }
   }

   void mul(VECTOR& x,VECTOR& y, double& xy){
      int n=x.size();
      xy=0;
      for(int i=0;i<n;i++){
         xy+=x[i]*y[i];
      }
   }

   void mul(MATRIX& A,VECTOR& x, VECTOR& Ax){
      int n=A.size();
      for(int i=0;i<n;i++){
         mul(A[i],x,Ax[i]);
      }
   }
   double dist(VECTOR& v, VECTOR& w){
      int n=v.size();
      double ans=0.0;
      for(int i=0;i<n;i++){
         ans+=(v[i]-w[i])*(v[i]-w[i]);
      }
      return sqrt(ans);
   }
   double normalize(VECTOR& v){//returns the original length
      double nv;
      mul(v,v,nv);
      nv=sqrt(nv);
      mul(v,1.0/nv);
      return nv;
   }
}


int main(){
   using namespace nummat;
   int n,maxit;
   double tol;
   scanf("%d%lf%d",&n,&tol,&maxit);
   MATRIX A(n);
   for(int i=0;i<n;i++){
      A[i].resize(n);
   }
   VECTOR v(n);

   read(A);
   read(v);

   
   if(normalize(v)<EPS){
      printf("zero starting vector\n");
   }else{
      VECTOR Av(n);
      mul(A,v,Av);
      double lam1,lam0;
      mul(Av,v,lam1);
   
      int i=2;
      while(i<=maxit){
         lam0=lam1;
         v=Av;
         normalize(v);
         mul(A,v,Av);
         mul(Av,v,lam1);
         if(abs(lam0-lam1)<tol){
            break;
         }
         i+=1;
      }
      
      if(i>maxit){
         printf("maxit\n");
      }else{
         mul(v,lam1);
         if(dist(v,Av)<tol*abs(lam1)){
            printf("success\n");
         }else{
            printf("fail\n");
         }
      }
      
   }

   return 0;
}
