// Original code by Istvan Mezo
// Modification by (anonymous)
// rLambert.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"
#include "stdlib.h"
#include "math.h"
#include "cmath"
#include <vector>
#include <stdio.h>
// #include <Rcpp.h>
std::vector<double> result(3);
// [[Rcpp::export]]
double LambertW(const double z) {
  int i; 
  const double eps=4.0e-16, em1=0.3678794411714423215955237701614608; 
  double p,e,t,w;
  if (z<-em1) { 
    fprintf(stderr,"LambertW: bad argument %g, exiting.\n",z); return 0; 
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
    fprintf(stderr,"LambertW1: bad argument %g, exiting.\n",z); return 0; 
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
double F ( const double x, const double r )
{
	return ( x * exp( x ) + r * x );
}

double DF ( const double x, const double r )
{
	return ( (1. + x) * exp( x ) + r );
}

double DDF ( const double x )
{
	return ( (2. + x) * exp( x ) );
}

// [[Rcpp::export]]
std::vector<double> rLambert( const double x, const double r, const int prec = 126 )
{
  result[0] = NAN;
  result[1] = NAN;
  result[2] = NAN;
	if ( x == 0. ){
	  result[0] = 0;
	  result[1] = 0;
	  result[2] = 0;
	  return result; // W(x,r)=0 always
	} 
	if ( r == 0. )
	{
		if ( x < 0 ) {result[1] = LambertW1(x);}
		result[0] = LambertW(x);
		result[1] = result[0];
		result[2] = result[0];
		return result;
	}

	const double e = exp( 1. );
	const double em2 = exp ( -2. );

	double w, wprev; // intermediate values in the Halley method
	const unsigned int maxiter = 20; // max number of iterations. This eliminates problems with convergence
	unsigned int iter; // iteration counter
	unsigned int iter_outer=0;
	const unsigned int maxiter_outer=2000;
	

	// If r >= exp(-2), there is just one branch of W(x,r):
	if ( r >= em2 )
	{
		if ( (r == em2) && (x == -4. * em2) ){
		  result[2] = -2;
		  return result;
		}  // At this x W(x,em2) is not differentiable, so
		// Halley's method does not work. But it can be calculated that W(x,em2) = -2.
		// If x is close but not equal to -4*em2, there is no problem.
		// For example, if x= -4. * exp( -2. ) +/- .00000000001,
		// the program still gives the correct result

		// begin the iteration up to prec precision
		// initial value for the Halley method
		if ( x > 1. )  w = log( x ) - log( log( x ) );
		if ( x < -1.) w = 1. / r * x;
		if ( abs(x) <= 1 ) w = 0.;
		wprev = w;

		iter = 0;
		do
		{
			wprev = w;
			w -= 2.*((F(w,r)-x) * DF(w,r)) /
				(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
			iter++;
		} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
		result[2] = w;
		return result;
	} // if ( r >= em2 )
	
	//here comes the calculation when W(x,r) has three branches (i.e. 0 < r < e^(-2))
	else if ( 0 < r && r < em2 )
	{
		double alpha  = LambertW1( -r * e ) - 1.; // left  branch point
		double beta   = LambertW ( -r * e ) - 1.; // right branch point
		double falpha = F( alpha, r );
		double fbeta  = F( beta,  r );

		if ( x < fbeta ) // leftmost branch
		{
			if ( x < -40.54374295204823 ){
			  result[2] = x/ r;
			  return result; // because x*e^x < 1E-16 as x < -40.5437...
			} 
			wprev = w = x / r;

			iter = 0;
			do
			{
				wprev = w;
				w -= 2.*((F(w,r)-x) * DF(w,r)) /
					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
				iter++;
			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
			result[2] = w;
			return result;
		} // if ( x < fbeta )

		if ( x >= fbeta && x <= falpha ) // leftmost and inner branches
		{
			double wm2 = x / r; // initial value for the leftmost branch W-2(x,r)
			double wm1 = -3.; // initial value for the inner branch W-1(x,r)
			double wm0 = -1.;  // initial value for the rightmost branch W0(x,r)
			
			bool condition1, condition2, condition3, condition4;
			double rounded1, rounded2, rounded3;
			iter_outer=0;
			do
			{
			  iter = 0;
			  do
			  {
			    wprev = wm2;
			    wm2 -= 2.*((F(wm2,r)-x) * DF(wm2,r)) /
			      (2.*DF(wm2,r)*DF(wm2,r) - (F(wm2,r)-x)*DDF(wm2));
			    iter++;
			  } while ( (abs( wm2 - wprev ) > pow( 10., -prec )) && iter < maxiter );
			  result[0] = wm2;
			  
			  // Halley iteration for the inner branch
			  iter = 0;
			  do
			  {
			    wprev = wm1;
			    wm1 -= 2.*((F(wm1,r)-x) * DF(wm1,r)) /
			      (2.*DF(wm1,r)*DF(wm1,r) - (F(wm1,r)-x)*DDF(wm1));
			    iter++;
			  } while ( (abs( wm1 - wprev ) > pow( 10., -prec )) && iter < maxiter );
			  result[1] = wm1;
			  // Halley iteration for the rightmost branch
			  iter = 0;
			  do
			  {
			    wprev = wm0;
			    wm0 -= 2.*((F(wm0,r)-x) * DF(wm0,r)) /
			      (2.*DF(wm0,r)*DF(wm0,r) - (F(wm0,r)-x)*DDF(wm0));
			    iter++;
			  } while ( (abs( wm0 - wprev ) > pow( 10., -prec )) && iter < maxiter );
			  result[2] = wm0;
			  
			  condition1 = result[0] <= alpha;
			  condition2 = result[1] <= beta && result[1] >= alpha;
			  condition3 = result[2] >= beta;
			  rounded1 = round(result[0]*1000.0);
			  rounded2 = round(result[1]*1000.0);
			  rounded3 = round(result[2]*1000.0);
			  condition4 = rounded1 != rounded2 && rounded2 != rounded3 && rounded1 != rounded3;
			  if (condition1 && condition2 && condition3 && condition4){
			    break;
			  }
			  iter_outer++;
			  if (iter_outer >= maxiter_outer){
			    break;
			  }
			  
			  
			  wm2 = x / r * iter_outer;
			  wm1 = (alpha + beta)/2 + (iter_outer*1.0) / (maxiter*1.0) *(alpha+beta)/2;
			  wm0 = beta + log(iter_outer);
			} while (true);
			// Halley iteration for the leftmost branch
			
			return result;
		} // if ( x >= fbeta && x <= falpha ) 

		if ( x > falpha ) // rightmost branch
		{
			wprev = w = ( x > 1. ) ? log ( x ) - log( log ( x ) ) : 0.; // initial value

			iter = 0;
			do
			{
				wprev = w;
				w -= 2.*((F(w,r)-x) * DF(w,r)) /
					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
				iter++;
			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
			result[2] = w;
			return result;
		} // if ( x > falpha )
			
	} // else if
	else if ( r < 0 ) // two branches separated by W(-r*e)-1
	{
	  // printf("In case 4\n");
	  printf("LambertW( -r * e ) - 1. is %-20.15lf\n",LambertW( -r * e ) - 1.);
	  printf("F( LambertW( -r * e ) - 1., r ) is %-20.15lf\n",F( LambertW( -r * e ) - 1., r ));
	  // printf("F( LambertW( -r * e ) - 1., r ) - x is %-20.15lf\n",F( LambertW( -r * e ) - 1., r )-x);
		// minimum of F: W(-r*e)-1, zeros: 0 and log(-r)
		double val = F( LambertW( -r * e ) - 1., r )-x;
		val = round(val* pow(10,10)) / pow(10,10);
		double gam = LambertW( -r * e ) - 1.;
		if ( x < F( gam, r ) && val != 0.0)
		{
      // printf( "%s\n", "First argument of LambertW(x,r) is out of domain" );
		  // printf("F( LambertW( -r * e ) - 1., r ) - x is %-20.15lf\n",F( LambertW( -r * e ) - 1., r )-x);
			return result;
		}
		if ( x == log( -r ) ){
		  result[1] = 0;
		  result[2] = 0;
		  return result;
		} 
	
	
	  double wl,wr;
		bool condition1, condition2,condition3,condition4;
		double rounded1, rounded2;
		iter_outer = 0;
		if (x < 0)
		{
		  wl = LambertW( -r * e ) - 2;
		  wr = LambertW( -r * e );
		}
		else
		{
		  wl = (r > -1 ? log(-r) : 0)-1;
		  wr = (r > -1 ? 0 : log(-r))+1;
		}
		do
		{
		  w = wl;
		  iter = 0;
		  do
		  {
		    wprev = w;
		    w -= 2.*((F(w,r)-x) * DF(w,r)) /
		      (2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
		    iter++;
		  } while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
		  result[1] = w;
		  
		  w = wr;
		  iter = 0;
		  do
		  {
		    wprev = w;
		    w -= 2.*((F(w,r)-x) * DF(w,r)) /
		      (2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
		    iter++;
		  } while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
		  result[2] = w;
		  // return result;
		  
		  rounded1 = round(result[1]*100.0);
		  rounded2 = round(result[2]*100.0);
		  condition1 = result[1] <= gam;
		  condition2 = result[2] >= gam;
		  condition3 = rounded1 != rounded2;
		  condition4 = !isnan(result[1]) && !isnan(result[2]);
		  // printf("Gam is %f\n",gam);
		  // printf("Condition1 is %d\nCondition2 is %d\nCondition3 is %d\nCondition4 is %d\n",condition1,condition2,condition3,condition4);
		  // printf("iter_outer=%d\n",iter_outer);
		  // printf("wl=%f\nwr=%f\n",wl,wr);
		  // printf("result1=%f\nresult2=%f\n",result[1],result[2]);
		  // printf("condition1=%d\ncondition2=%d\ncondition3=%d\n",condition1,condition2,condition3);
		  if (condition1 && condition2 && condition3 && condition4)
		  {
		    break;  
		  }
		  iter_outer++;
		  if (iter_outer >= maxiter_outer)
		  {
		    printf("Warning: No convergence!\n");
		    result[1] = NAN;
		    result[2] = NAN;
		    break;
		  }
		  float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		  float r2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		  if (!condition3){
		    if (r2 > 0.5){
		      wl = wl*r1;
		    }
		    // else{
		    //   wl = wl/r2;
		    // }
		  }
		  wl = wl*1.5;
		  wr = wr*1.5;
		    
		  if (isnan(result[1])){
		    wl = wl/4;
		  }
		  if (isnan(result[2])){
		    wr = wr/4;
		  }
		  
		  if (iter_outer % 50 == 0){
		    wl = -1;
		    wr = 1;
		  }
		} while (true);
		return result;
// 		if ( x < 0 ) // Two initial values, one less than the minimum of f, the other is one greater
// 		{
// 			// left branch
// 			wprev = w = LambertW( -r * e ) - 2;
// 	
// 			iter = 0;
// 			do
// 			{
// 				wprev = w;
// 				w -= 2.*((F(w,r)-x) * DF(w,r)) /
// 					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
// 				iter++;
// 			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
// 			result[1] = w;
// 
// 			// right branch
// 			wprev = w = LambertW( -r * e );
// 			iter = 0;
// 			do
// 			{
// 				wprev = w;
// 				w -= 2.*((F(w,r)-x) * DF(w,r)) /
// 					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
// 				iter++;
// 			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
// 			result[2] = w;
// 			return result;
// 			}
// 
// 		if ( x > 0 )
// 		{
// 			double lzero = r > -1 ? log(-r) : 0;
// 			double rzero = r > -1 ? 0 : log(-r);
// 
// 			// left branch
// 			wprev = w = lzero - 1;
// 			iter = 0;
// 			do
// 			{
// 				wprev = w;
// 				w -= 2.*((F(w,r)-x) * DF(w,r)) /
// 					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
// 				iter++;
// 			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
//       result[1] = w;
//       
// 			// right branch
// 			wprev = w = rzero + 1;
// 			iter = 0;
// 			do
// 			{
// 				wprev = w;
// 				w -= 2.*((F(w,r)-x) * DF(w,r)) /
// 					(2.*DF(w,r)*DF(w,r) - (F(w,r)-x)*DDF(w));
// 				iter++;
// 			} while ( (abs( w - wprev ) > pow( 10., -prec )) && iter < maxiter );
// 			result[2] = w;
// 			return result;
// 		} // if ( x < 0 )
		} // else if ( r < 0 )
	/* should never get here */
	//fprintf(stderr, "LambertW1: No convergence at x=%g, exiting.\n", x);
	return result;
} // rLambertW

// [[Rcpp::export]]
std::vector<double> theLambertSolutionsC(const double x2_bar, const double x_bar, const double mu,
                                        const double a, const double b)
{
  double alpha = x2_bar - mu*x_bar + (mu-x_bar)*a;
  double beta;// = x2_bar - mu*x_bar + (mu-x_bar)*b;
  beta = fma(mu,b,fma(-x_bar,b,fma(-x_bar,mu,x2_bar)));
  double m;// = 0.5 * (pow(a-mu,2) - pow(b-mu,2));
  m = fma(a-b,(a+b)/2.,(b-a)*mu);
  double r = -alpha / beta;
  double e_power = exp(-m/beta);
  double inp = m * (alpha / pow(beta,2) - 1/beta) * e_power;
  std::vector<double> lamb = rLambert(inp,r*e_power,16);
  std::vector<double> result2(3);
  double d;
  int i;
  
  printf("alpha=%f\n",alpha);
  printf("beta=%f\n",beta);
  printf("m=%f\n",m);
  printf("inp=%f\n",inp);
  printf("r=%f\n",r);
  printf("e_power=%f\n",e_power);
  printf("r*e_power=%f\n",r*e_power);
  for (i=0;i<=2;i++){
    printf("lamb[i]=%f\n",lamb[i]);
    if (isnan(lamb[i])){
      result2[i] = NAN;
    }
    else{
      d = fma(lamb[i],beta,m);
      result2[i] = beta * m / d;
    }
  }
  return result2;
}

int main(){
  std::vector<double> result1(3);
  result1 = theLambertSolutionsC(784.,98.,97,49.5,146.5);
  printf("result1%f\nresult2%f\nresult3%f\n",result1[0],result1[1],result1[2]);
  return 0;
}