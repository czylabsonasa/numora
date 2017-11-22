#define EPS 1e-15

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
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


	typedef function<double(double)> valos1;

    double gyok2(double x){// [0,2]-> f(sqrt(2))=0
        return x*x-2.0;
    }
    double Dgyok2(double x){// [0,2]-> f(sqrt(2))=0
        return 2.0*x;
    }



}// end of ns

using namespace nummat;


double felez(valos1 f, double a, double b, int it){
    double fa=f(a);
    double fb=f(b);
    while(it--){
        double mid=0.5*(a+b);
        double fm=f(mid);
        if(fm*fa<=0.0){
            b=mid;
            fb=fm;
        }else{
            a=mid;
            fa=fm;
        }
    }
    return 0.5*(a+b);
}

double nyuton(valos1 f, valos1 df, double x0, int it){
    while(it--){
		x0-=f(x0)/df(x0);
    }
    return x0;
}


int main(){

	int it=1;
	double a=0,b=2;

	int maxi=5;
	for(int i=1;i<=maxi;i++){
		it*=2;
		double Fans=felez(gyok2,a,b,it);
		double Nans=nyuton(gyok2,Dgyok2,0.5*(a+b),it);
		printf("%05d iteráció után: Fhiba:%20e vs. Nhiba:%20e\n",it,Fans-sqrt(2.0),Nans-sqrt(2.0));   
	}


   return 0;
}
