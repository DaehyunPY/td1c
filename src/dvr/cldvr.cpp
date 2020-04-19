////////////////////////////////////////////////////////////////////////
// DVR: basic
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
#include "radau_alglib.hpp"
////////////////////////////////////////////////////////////////////////
cldvr::cldvr()
{
}
////////////////////////////////////////////////////////////////////////
cldvr::cldvr(long n, long m0, long m1, double x0, double x1)
{
  gen(n, m0, m1, x0, x1);
}
////////////////////////////////////////////////////////////////////////
cldvr::~cldvr()
{
  xpt.resize(0);
  wpt.resize(0);
  dshape.resize(0);
}
////////////////////////////////////////////////////////////////////////
void cldvr::init()
{
}
////////////////////////////////////////////////////////////////////////
void cldvr::final()
{
}
////////////////////////////////////////////////////////////////////////
void cldvr::clear()
{
//  xpt.clear();
//  wpt.clear();
//  dshape.clear();
}
////////////////////////////////////////////////////////////////////////
void cldvr::print() const
{
  printf("# Grid nodes and weights:\n");

  for(long m = 0; m <= ndvr; m ++) {
    printf("%5ld%20.10f%20.10f\n", m, xpt[m], wpt[m]);
  }

  printf("# Lobatto shape functions:\n");

  for(long m = 0; m <= ndvr; m ++) {
    for(long k = 0; k <= ndvr; k ++) {
      long km = (ndvr + 1) * m + k;
      printf( "%5ld%5ld%20.10f\n", m, k, dshape[km]);
    }
  }
}
////////////////////////////////////////////////////////////////////////
void cldvr::gen(long n, long m0, long m1, double x0, double x1)
{
  ndvr = n;
  mmin = m0;
  mmax = m1;
  xmin = x0;
  xmax = x1;

  xpt.resize(ndvr + 1);
  wpt.resize(ndvr + 1);
  dshape.resize((ndvr + 1)*(ndvr + 1));

  gen_xpt();
  gen_wpt();
  double lenx = HALF * (xmax - xmin);
  double midx = HALF * (xmax + xmin);
  for( long k = 0; k <= ndvr; k ++ ) {
    xpt[k] = xpt[k] * lenx + midx;
    wpt[k] = wpt[k] * lenx;
  }

  gen_dshape();
  //  print();
}
////////////////////////////////////////////////////////////////////////
void cldvr::gen_xpt()
{
  const long nled = ndvr - 1;
  //BUG  std::vector<double> amat(nled*nled);
  //BUG  std::vector<double> evec(nled*nled);
  std::vector<double> amat(nled*nled, ZERO);
  std::vector<double> evec(nled*nled, ZERO);

  xpt[0]    = - ONE;
  xpt[ndvr] = + ONE;

  if (ndvr > 1) {

    //BUG    amat.clear();
    //BUG    evec.clear();
    for (long m = 2; m <= nled; m ++) {
      double a, b, xm;
      xm = static_cast<double>( m );
      a = (xm + ONE) / (TWO * xm + ONE);
      b = (xm - ONE) / (TWO * xm - ONE);
      //old      amat(m - 2, m - 1) = sqrt(a * b);
      amat[nled * (m - 1) + m - 2] = sqrt(a * b);
    }
    for (long m = 0; m < nled - 1; m ++) {
      //old      amat(m + 1, m) = amat(m, m + 1);
      amat[nled * m + m + 1] = amat[nled * (m + 1) + m];
    }
    //debug
    //for (long m1 = 0; m1 < nled; m1 ++) {
    //  for (long m2 = 0; m2 < nled; m2 ++) {
    //    printf("%10.5f", amat(m1, m2));
    //  }
    //  printf("\n");
    //}
    //debug

    //  printf("before lapack_dsyev:\n"); 
    //  abort();
    //old    lapack_dsyev_(&nled, &amat(0,0), &evec(0,0));
    lapack_dsyev_(&nled, &amat[0], &evec[0]);
    //  printf("after lapack_dsyev:\n"); 
    //  abort();

    for(long m = 1; m <= ndvr - 1; m ++ ) {
      //old      xpt[m] = amat(m - 1, m - 1);
      xpt[m] = amat[nled * (m - 1) + m - 1];
    }
  }

  //debug
  //  for(long m = 0; m <= ndvr; m ++ ) {
  //    printf("%10ld%20.10f\n", m, xpt[m]);
  //  }
  //  printf("end of cldvr_gen_xpt.\n");
  //  abort();
  //debug
}
////////////////////////////////////////////////////////////////////////
void cldvr::gen_wpt()
{
  const double xn = static_cast<double>( ndvr );
  const double scale = TWO / (xn * (xn + ONE));

  wpt[0]    = scale;
  wpt[ndvr] = scale;

  if (ndvr > 1) {
    double p, absg;
    for (long k = 1; k < ndvr; k ++) {
      absg = get_abs( xpt[ k ] );
      get_legendre( absg, p );
      wpt[ k ] = scale / (p * p);
    }
  }
// DEBUG
//  printf("# %ld's order Gauss-Lobatto weights:\n", ndvr + 1);
//  for(long k = 0; k <= ndvr; k ++ ) {
//    printf( "%5ld %20.10f\n", k, wpt[ k ]);
//  }
// DEBUG
}
////////////////////////////////////////////////////////////////////////
void cldvr::get_legendre(double& x, double& p)
{
  if ( ndvr == 0 ) {
    p = ONE;
  } else if ( ndvr == 1 ) {
    p = x;
  } else {
//    double p0 = ONE;
//    double p1 = x;
//    double p2 = ZERO;
//    double c0, c1, _Rk;
    long double p0 = ONE;
    long double p1 = x;
    long double p2 = ZERO;
    long double c0, c1, rk;

    for ( long k = 2; k <= ndvr; k ++ ) {
      rk = static_cast<double>( k );
      c0 = - (rk - ONE) / rk;
      c1 = (TWO*rk - ONE) / rk;
      p2 = c0 * p0 + c1 * x * p1;
      p0 = p1;
      p1 = p2;
    }
    p = p2;
  }
}
////////////////////////////////////////////////////////////////////////
void cldvr::gen_dshape()
//
// dshape(k, m) = df_m/dx (x_k)
//
{
  long km;
  double dshx, diffm, diffk;

  for ( long m = 0; m <= ndvr; m ++ ) {
    for ( long k = 0; k <= ndvr; k ++ ) {
      km = (ndvr + 1) * m + k;
      if ( m != k ) {
	dshx = ONE / ( xpt[ m ] - xpt[ k ] );
	for ( long l = 0; l <= ndvr; l ++ ) {
	  if ( l != m && l != k ) {
	    diffm = xpt[ m ] - xpt[ l ];
	    diffk = xpt[ k ] - xpt[ l ];
	    dshx *= diffk / diffm;
	  }
	}
	dshape[km] = dshx;
      } else {
	dshape[km] = ZERO;
      }
    }
  }

  dshape[0]                   = -HALF / wpt[ 0 ];
  dshape[(ndvr+1)*(ndvr+1)-1] =  HALF / wpt[ ndvr ];
}
////////////////////////////////////////////////////////////////////////
double cldvr::get_val(long m, double rval) const
//
// val_m(r) = \prod_{k \ne m} \frac{r - r_m}{r_k - r_m}
//
{
  double val = ONE;
  for ( long k = 0; k <= ndvr; k ++ ) {
    if ( m != k ) {
      val *= (rval - xpt[ k ] ) / ( xpt[ m ] - xpt[ k ] );
    }
  }
  return val;
}
////////////////////////////////////////////////////////////////////////
void cldvr::gen_radau(long n, long m0, long m1,
		      double x0, double x1, double d_exp_factor)
{
  ndvr = n;
  mmin = m0;
  mmax = m1;
  xmin = x0;
  xmax = x1;
  exp_alpha = d_exp_factor;

  xpt.resize(ndvr + 1);
  wpt.resize(ndvr + 1);
  dshape.resize((ndvr + 1)*(ndvr + 1));
  
  gen_xw_radau();
  
//   double lenx = HALF * (xmax - xmin);
//   double midx = HALF * (xmax + xmin);
//   for( long k = 0; k <= ndvr; k ++ ) {
//     xpt[k] = xpt[k] * lenx + midx;
//     wpt[k] = wpt[k] * lenx;
//   }

  gen_dshape_radau();
  
  
}
////////////////////////////////////////////////////////////////////////
void cldvr::gen_xw_radau()
{
  //recurrence relation                                                                                                                               
  alglib::real_1d_array alpha;
  alglib::real_1d_array beta;
  double mu0;
  double a;
  alglib::ae_int_t nnum;
  //output                                                                                                                                            
  alglib::ae_int_t info;
  alglib::real_1d_array x;
  alglib::real_1d_array w;

  nnum = ndvr;
  a = xmin;
  
  alpha.setlength(nnum-1);
  beta.setlength(nnum);
  x.setlength(nnum);
  w.setlength(nnum);
  
  //alpha(0) = 0, beta(0) = arbitrary                                                                                                                 
  double di;
  for(int i = 0; i < nnum-1; i++){
    di = static_cast<double>( i );
    alpha[i] = (2*di+1) / exp_alpha + a;
    beta[i] = di*di/(exp_alpha*exp_alpha);
  }
  beta[nnum-1] = static_cast<double>( (nnum-1)*(nnum-1) )/(exp_alpha*exp_alpha);

  //0th moment of weighting function mu0 = integral[a,infinite]exp(-x)dx 
  mu0 = 1.0/exp_alpha;
  
  
  alglib::gqgenerategaussradaurec(alpha,beta,mu0,a,nnum,info,x,w);
  if(info != 1){
    printf("alglib::gqgenerategaussradaurec error : info = %d", info);
    abort();
  }
  for(int i = 0; i < nnum; i++){
    printf("x[%3d] = %15.10E : w[%2d] = %15.10E \n", i, x[i], i, w[i]);
    xpt[i] = static_cast<double>(x[i]);
    //BUG wpt[i] = w[i] * exp(x[i]);
    wpt[i] = static_cast<double>(w[i]);
  }

  xpt[0] = a;  
  //actually, wpt[nnum] and xpt[nnum] have nothing to do with the result
  //so anything must cause no-problem.
  wpt[nnum] = 1;
  xpt[nnum] = xpt[nnum-1] + 1;
  
}
////////////////////////////////////////////////////////////////////////
void cldvr::gen_dshape_radau()
//
// dshape(k, m) = df_m/dx (x_k)
//in last element, radau means x_k|k=ndvr  = infinity 
{
  //bug dshape.clear();

  long km;
  double dshx, diffm, diffk;
  //  int shape = 1; //to compemsate ndvr = ndvr -1
  //not to go to last point :ndvr -> ndvr-1
  for ( long m = 0; m <= ndvr-1; m ++ ) {
    for ( long k = 0; k <= ndvr-1; k ++ ) {
      km = (ndvr + 1) * m + k;
      if ( m != k ) {
	dshx = ONE / ( xpt[ m ] - xpt[ k ] );
	for ( long l = 0; l <= ndvr-1; l ++ ) {
	  if ( l != m && l != k ) {
	    diffm = xpt[ m ] - xpt[ l ];
	    diffk = xpt[ k ] - xpt[ l ];
	    dshx *= diffk / diffm;
	    //BUG
	    //BUG if(l != ndvr) dshx *= diffk / diffm;
	    //BUG else if(l == ndvr) dshx *= ONE;
	    //BUG
	  }
	}
	dshape[km] = dshx;
      } else {
	// the case of  m == k
	for(long l = 0; l <= ndvr-1; l++){
	  if(l == m) continue;
	  dshape[km] += 1/(xpt[m] - xpt[l]);
	}
	dshape[km] += - HALF * exp_alpha;
	
      }
      //BUG if(k == ndvr) dshape[km] = ZERO;
    }
  }

  //printf("%10.8e",dshape[0]);
  //printf("%10.8e", -1.0/(2*wpt[0]));
  //abort();

  //dshape[0] += -HALF / wpt[ 0 ];
  //dshape[(ndvr+1) * (ndvr+1) -1] +=  HALF / wpt[ ndvr ];
  //for(long l = 0; l <= ndvr-2; l++){
  //  dshape[0] += 1/(xpt[0] - xpt[l+1]);
  //  dshape[(ndvr+1) * (ndvr-1) + ndvr-1] += 1/(xpt[ndvr-1] - xpt[l]);  
  //}


//   std::cout<< "+++++++++++++++++++++" <<std::endl;
//   std::cout << wpt[ndvr]  << std::endl;
//   std::cout << dshape[0]  << std::endl;
//   std::cout << dshape[(ndvr+1)*(ndvr+1)-1] << std::endl;
}
