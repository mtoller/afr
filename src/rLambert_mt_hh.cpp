#include "stdlib.h"
#include "math.h"
#include "cmath"
#include <vector>
#include <stdio.h>
#include <float.h>
// #include <gmpxx.h>
// #include <Rcpp.h>
std::vector<double> result(3);

// [[Rcpp::export]]
double LambertW0(const double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p,e,t,w;
  if (z<-em1) { 
    fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z);
    return NAN;
  }
  if (0.0==z) return 0.0;
  if (z<-em1+1e-4) { // series near -em1 in sqrt(q)
    double q=z+em1,r=sqrt(q),q2=q*q,q3=q2*q;
    return 
      -1.0
    +2.331643981597124203363536062168*r
    -1.812187885639363490240191647568*q
    +1.936631114492359755363277457668*r*q
    -2.353551201881614516821543561516*q2
    +3.066858901050631912893148922704*r*q2
    -4.175335600258177138854984177460*q3
    +5.858023729874774148815053846119*r*q3
    -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }
  /* initial approx for iteration... */
  if (z<1.0) { /* series near 0 */
  p=sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else 
    w=log(z)-log(log(z)); /* asymptotic */
  for (i=0; i<10; i++) { /* Halley iteration */
  e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  //fprintf(stderr,"LambertW: No convergence at z=%g, exiting.\n",z); 
  return 0;
}
// [[Rcpp::export]]
double LambertW1(const double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p=1.0,e,t,w,l1,l2;
  if (z<-em1 || z>=0.0) { 
    fprintf(stderr,"LambertW1: bad argument %g, exiting.\n",z);
    return NAN;
  }
  /* initial approx for iteration... */
  if (z<-1e-6) { /* series about -1/e */
  p=-sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));
    w=-1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); 
  } else { /* asymptotic near zero */
  l1=log(-z);
    l2=log(-l1);
    w=l1-l2+l2/l1;
  }
  if (fabs(p)<1e-4) return w;
  for (i=0; i<10; i++) { /* Halley iteration */
  e=exp(w); 
    t=w*e-z;
    p=w+1.0;
    t/=e*p-0.5*(p+1.0)*t/p; 
    w-=t;
    if (fabs(t)<eps*(1.0+fabs(w))) return w; /* rel-abs error */
  }
  /* should never get here */
  //fprintf(stderr,"LambertW1: No convergence at z=%g, exiting.\n",z);
  return 0;
}

const double EPS = 1e-12;

double f(const double y, const double r){
  return y * exp(y) + r * y;
}

double bs(double a, double b, const double x, const double r, const bool increasing){
  double mid, x_;
  int i = 0;
  while(abs(a-b) > EPS){
    i++;
    if (i > 100000){
      break;
    }
    mid = 0.5 * (a+b);
    x_ = f(mid, r);
    if(increasing)
      if(x_ < x) a = mid; else b = mid;
    else
      if(x_ > x) a = mid; else b = mid;
  }

  return a;
}
// [[Rcpp::export]]
std::vector<double> rLambert(const double x, const double r){
  result[0] = NAN;
  result[1] = NAN;
  result[2] = NAN;
  const double e = exp( 1. );
  const double em2 = exp ( -2. );
  if (x == 0. | x == -0.){
    result[0] = 0;
    result[1] = 0;
    result[2] = 0;
    return result;
  }
  
  if (r == 0.){
    // printf("Case0\n");
    result[1] = LambertW1(x);
    result[2] = LambertW0(x);
    return result;
  }
  
  if (r >= em2){ //Case 1 & Case2 are the same
    // printf("Case1\n");
    double a = DBL_MIN / 2.;
    double b = DBL_MAX / 2.;
    if (x > 0.) a = 0; else b = 0;
    result[2] = bs(a,b,x,r,true);
  }
  else if (r > 0 && r < em2){ //Case 3
    // printf("Case3\n");
    double alpha = LambertW1(-r*e) - 1;
    double beta = LambertW0(-r*e) - 1;
    
    double f_alpha = f(alpha,r);
    double f_beta = f(beta,r);
    
    if (x < f_alpha){
      result[0] = bs(-DBL_MAX/ 2.,alpha,x,r,true);
    }
    else if (x >= f_alpha && x <= f_beta){
      result[0] = bs(-DBL_MAX/ 2.,f_alpha,x,r,true);
      result[1] = bs(alpha,beta,x,r,false);
      result[2] = bs(alpha,DBL_MAX /2.,x,r,true);
    }
    else if (x > f_beta){
      result[2] = bs(alpha,DBL_MAX /2.,x,r,true);
    }
    else{
      abort();
    }
  }
  else if (r < 0){ // Case 4
    // printf("Case4\n");
    double gamma = LambertW0(-r*e) - 1;
    double f_gamma = f(gamma,r);
    
    // printf("gamma=%.*e\n",DECIMAL_DIG,gamma);
    // printf("f_gamma=%.*e\n",DECIMAL_DIG,f_gamma);
    if (x >= f_gamma){
      result[1] = bs(-DBL_MAX / 2., gamma,x,r,false);
      result[2] = bs(gamma,DBL_MAX / 2.,x,r,true);
    }
  }
  
  return result;
}

// [[Rcpp::export]]
std::vector<double> theLambertSolutionsC(const double x2_bar, const double x_bar, const double mu,
                                         const double a, const double b)
{
  double alpha;// = x2_bar - mu*x_bar + (mu-x_bar)*a;
  alpha = fma(mu,a,fma(-x_bar,a,fma(-x_bar,mu,x2_bar)));
  double beta;// = x2_bar - mu*x_bar + (mu-x_bar)*b;
  beta = fma(mu,b,fma(-x_bar,b,fma(-x_bar,mu,x2_bar)));
  std::vector<double> result2(3);
  
  if (abs(beta) < 1e-10){
    beta = 0.0;
    for (int i=0;i<=2;i++)result2[i]=0;
    return result2;
  }
  
  double m;// = 0.5 * (pow(a-mu,2) - pow(b-mu,2));
  m = fma(a-b,(a+b)/2.,(b-a)*mu);
  double r = -alpha / beta;
  double e_power = exp(-m/beta);
  // double inp = m * ((alpha - beta)/pow(beta,2)) * e_power;
  double inp = m * ((alpha - beta)/pow(beta,2)) * e_power;

  // printf("mu=%f\n",mu);
  // printf("x_bar=%f\n",x_bar);
  // printf("x2_bar=%f\n",x2_bar);
  // printf("a=%f\n",a);
  // printf("b=%f\n",b);
  // printf("alpha=%f\n",alpha);
  // printf("beta=%.*e\n",DECIMAL_DIG,beta);
  // printf("-alpha/beta=%.*e\n",DECIMAL_DIG,-alpha/beta);
  // printf("(alpha-beta)/beta^2=%.*e\n",DECIMAL_DIG,(alpha-beta)/pow(beta,2));
  // printf("-m/beta=%.*e\n",DECIMAL_DIG,-m/beta);
  // printf("m=%.*e\n",DECIMAL_DIG,m);
  // printf("inp=%.*e\n",DECIMAL_DIG,inp);
  // printf("e_power=%.*e\n",DECIMAL_DIG,e_power);
  // printf("r*e_power=%.*e\n",DECIMAL_DIG,r*e_power);

  std::vector<double> lamb = rLambert(inp,r*e_power);
  
  double d;
  int i;
  
  for (i=0;i<=2;i++){
    // printf("lamb[i]=%f\n",lamb[i]);
    // printf("lamb[i]=%.*e\n",DECIMAL_DIG,lamb[i]);
    if (isnan(lamb[i])){
      result2[i] = NAN;
    }
    else{
      d = fma(lamb[i],beta,m);
      // printf("beta=%.*e\n",DECIMAL_DIG,beta);
      // printf("m=%.*e\n",DECIMAL_DIG,m);
      // printf("d=%.*e\n",DECIMAL_DIG,d);
      result2[i] = beta * m / d;
      // printf("sigma2[i]=%f\n",result2[i]);
    }
  }
  
  
  return result2;
}


int main(){
  
  std::vector<double> res(3);
  res = theLambertSolutionsC(1264.25,33.9375,27.,20.,50.);  
  printf("res0=%f\nres1=%f\nres2=%f\n",res[0],res[1],res[2]);
  return 0;
}