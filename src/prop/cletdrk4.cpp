////////////////////////////////////////////////////////////////////////
// ETD-RK4
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
int cletdrk4::max_dim_numer = 8;
int cletdrk4::max_dim_denom = 3;
////////////////////////////////////////////////////////////////////////
cletdrk4::cletdrk4()
{
  std::cout << "cletdrk4" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cletdrk4::~cletdrk4()
{
  std::cout << "~cletdrk4" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cletdrk4::cletdrk4(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		   const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		   const clfield& Field, const clhprod& HPW)
{
  etd_dt1 = Field.dtime;
  etd_dt2 = Field.dtime * HALF;

  IO.read_info("etd_dim_numer", 2, dim_numer);
  IO.read_info("etd_dim_denom", 2, dim_denom);
  IO.read_info("etd_pade_old", false, pade_old);

  gen_exp(MPIP, IO, Bas);
  gen_gfun(MPIP, IO, Bas);
  gen_f0fun(MPIP, IO, Bas);
  gen_f1fun(MPIP, IO, Bas);
  gen_f2fun(MPIP, IO, Bas);

  //DEBUG
  printf("cletdrk4::exp1:\n"); exp1.print();
  printf("cletdrk4::exp2:\n"); exp2.print();
  printf("cletdrk4::gfun:\n"); gfun.print();
  printf("cletdrk4::f0fun:\n"); f0fun.print();
  printf("cletdrk4::f1fun:\n"); f1fun.print();
  printf("cletdrk4::f2fun:\n"); f2fun.print();
  //  std::cout << "abort for debug in cletdrk4::gen." << std::endl;
  //  abort();
  //DEBUG

  Wfn0.gen(MPIP, IO, Bas);
  Wfn1.gen(MPIP, IO, Bas);
  Wfn2.gen(MPIP, IO, Bas);
  Wfn3.gen(MPIP, IO, Bas);
  hWfn0.gen(MPIP, IO, Bas);
  hWfn1.gen(MPIP, IO, Bas);
  hWfn2.gen(MPIP, IO, Bas);
  hWfn3.gen(MPIP, IO, Bas);
  tWfn.gen(MPIP, IO, Bas);
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		    double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  std::cout << "cletdrk4::prop nyi." << std::endl;
  abort();
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::prop(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double ttime;
  double lfield[9];
  double thalf = Field.time + etd_dt2;

  // step (0)
  HPW.copy(Proc, Bas, Wfn, Wfn0);

  // step (1)
  ttime = Field.time;
  Field.get_value(ttime, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  prod1(Proc, Bas, HPW, thalf, Field, gfun, gfunrk, hWfn0, Wfn1);
  prod1(Proc, Bas, HPW, thalf, Field, exp2, exp2rk, Wfn0, tWfn);
  HPW.xpyz(Proc, Bas, tWfn, Wfn1, Wfn1);

  // step (2)
  ttime = Field.time + etd_dt2;
  Field.get_value(ttime, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
  prod1(Proc, Bas, HPW, thalf, Field, gfun, gfunrk, hWfn1, Wfn2);
  HPW.xpyz(Proc, Bas, tWfn, Wfn2, Wfn2);

  // step (3)
  ttime = Field.time + etd_dt2;
  Field.get_value(ttime, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn2, hWfn2);
  HPW.axpbyz(Proc, Bas, CTWO, hWfn2, -RUNIT, hWfn0, tWfn);
  prod1(Proc, Bas, HPW, thalf, Field, gfun, gfunrk, tWfn, Wfn3);
  prod1(Proc, Bas, HPW, thalf, Field, exp2, exp2rk, Wfn1, tWfn);
  HPW.xpyz(Proc, Bas, tWfn, Wfn3, Wfn3);

  // step (4)
  ttime = Field.time + etd_dt1;
  Field.get_value(ttime, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn3, hWfn3);

  // summation
  HPW.xpyz(Proc, Bas, hWfn1, hWfn2, hWfn2);
  prod1(Proc, Bas, HPW, thalf, Field, exp1, exp1rk, Wfn0, Wfn);
  prod2(Proc, Bas, HPW, thalf, Field, f0fun, f0funrk, hWfn0, Wfn);
  prod2(Proc, Bas, HPW, thalf, Field, f1fun, f1funrk, hWfn2, Wfn);
  prod2(Proc, Bas, HPW, thalf, Field, f2fun, f2funrk, hWfn3, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::prod1(const clmpi& Proc, const clbas& Bas, clhprod& HPW, 
		     double time, const clfield& Field, const clh1rat& FacOrb, 
		     dcomplex FacCI, const clwfn& WIn, clwfn& WOut)
{
  HPW.copy(Proc, Bas, WIn, WOut);
  HPW.scalc(Proc, Bas, FacCI, WOut);
  FacOrb.prod(Proc, Bas, time, Field, HPW, WOut);
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::prod2(const clmpi& Proc, const clbas& Bas, clhprod& HPW, 
		     double time, const clfield& Field, const clh1rat& FacOrb, 
		     dcomplex FacCI, clwfn& WIn, clwfn& WOut)
{
  HPW.scalc(Proc, Bas, FacCI, WIn);
  FacOrb.prod(Proc, Bas, time, Field, HPW, WIn);
  HPW.xpyz(Proc, Bas, WIn, WOut, WOut);
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::gen_exp(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  exp1rk = RUNIT;
  exp2rk = RUNIT;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> exp_n(max_dim_numer + 1);
  std::vector<dcomplex> exp_d0(max_dim_denom);
  std::vector<dcomplex> exp_d1(max_dim_denom);

  if (dim_denom == 0) {
    // 8-th order Taylor expansion of exp(x) generated by Maxima.
    //             1     8     1     7     1    6     1    5    1    4    1   3    1   2        
    // exp(x) = (-----) x  + (----) x  + (---) x  + (---) x  + (--) x  + (-) x  + (-) x  + x + 1
    //           40320        5040        720        120        24        6        2            
    exp_n[8] = 1.0 / 40320.0;
    exp_n[7] = 1.0 / 5040.0;
    exp_n[6] = 1.0 / 720.0;
    exp_n[5] = 1.0 / 120.0;
    exp_n[4] = 1.0 / 24.0;
    exp_n[3] = 1.0 / 6.0;
    exp_n[2] = 1.0 / 2.0;
    exp_n[1] = 1.0;
    exp_n[0] = 1.0;
  } else {
    if (dim_numer == 1 && dim_denom == 1) {
      //             x + 2
      // exp(x) = [- -----]
      //             x - 2
      // x0 = 2
      dfac = -1.0;
      exp_n[0] = 2.0;
      exp_n[1] = 1.0;
      roots[0] = 2.0;
    } else if (dim_numer == 2 && dim_denom == 2) {
      //            2            
      //           x  + 6 x + 12  
      // exp(x) = [-------------] 
      //            2             
      //           x  - 6 x + 12
      // x0 =  3 - sqrt(3) %i
      // x1 =  3 + sqrt(3) %i
      dfac = 1.0;
      exp_n[0] = 12.0;
      exp_n[1] = 6.0;
      exp_n[2] = 1.0;
      roots[0] =  3.0 - sqrt(3.0) * IUNIT;
      roots[1] =  3.0 + sqrt(3.0) * IUNIT;
    } else if (dim_numer == 3 && dim_denom == 3) {
      //              3       2              
      //             x  + 12 x  + 60 x + 120 
      // exp(x) = [- -----------------------]
      //              3       2              
      //             x  - 12 x  + 60 x - 120 
      // x0 = +4.644370709252173768 +0.000000000000000000 %i
      // x1 = +3.677814645373918445 +3.508761919567442433 %i
      // x2 = +3.677814645373918445 -3.508761919567442433 %i
      dfac = -1.0;
      exp_n[0] = 120.0;
      exp_n[1] = 60.0;
      exp_n[2] = 12.0;
      exp_n[3] = 1.0;
      roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
      roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
      roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
    } else {
      std::cout << "cletdrk4:gen_exp. bad dim_numer or dim_denom." << std::endl;
      abort();
    }

    exp_d1[0] = +dfac;
    exp_d0[0] = -dfac * roots[0];
    for (int idim = 1; idim < dim_denom; idim ++) {
      exp_d1[idim] = +RUNIT;
      exp_d0[idim] = -roots[idim];
    }
  }

  dcomplex tunit, tfac, ttmp;
  std::vector<dcomplex> cn;
  std::vector<dcomplex> cd0;
  std::vector<dcomplex> cd1;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // ##### exp(-i*h*dt) #####
  cn = exp_n;
  cd0 = exp_d0;
  cd1 = exp_d1;
  tfac = tunit * etd_dt1;
  ttmp = RUNIT; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  exp1.gen(MPIP, IO, Bas, etd_dt1, -1, dim_numer, dim_denom, cn, cd0, cd1);

  // ##### exp(-i*h*dt/2) #####
  cn = exp_n;
  cd0 = exp_d0;
  cd1 = exp_d1;
  tfac = tunit * etd_dt2;
  ttmp = RUNIT; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  exp2.gen(MPIP, IO, Bas, etd_dt2, -1, dim_numer, dim_denom, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("# exp: real\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = exp(x);
  //   vala = test_pade(exp_n, exp_d0, exp_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# exp: imag\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = exp(x);
  //   vala = test_pade(exp_n, exp_d0, exp_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::gen_gfun(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  gfunrk = etd_dt2;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> gfun_n(max_dim_numer + 1);
  std::vector<dcomplex> gfun_d0(max_dim_denom);
  std::vector<dcomplex> gfun_d1(max_dim_denom);

  if (dim_denom == 0) {
    //  
    // 8-th order Taylor expansion of gfun(x) = (exp(x)-1)/x generated by Maxima.
    //               1      8      1     7     1     6     1    5     1    4    1    3    1   2    1       
    // gfun(x) =  (------) x  + (-----) x  + (----) x  + (---) x  + (---) x  + (--) x  + (-) x  + (-) x + 1
    //             362880        40320        5040        720        120        24        6        2       
    // 
    gfun_n[8] = 1.0 / 362880.0;
    gfun_n[7] = 1.0 / 40320.0; 
    gfun_n[6] = 1.0 / 5040.0;  
    gfun_n[5] = 1.0 / 720.0;   
    gfun_n[4] = 1.0 / 120.0;   
    gfun_n[3] = 1.0 / 24.0;    
    gfun_n[2] = 1.0 / 6.0;     
    gfun_n[1] = 1.0 / 2.0;     
    gfun_n[0] = 1.0;
  } else {
    if (dim_numer == 1 && dim_denom == 1) {
      //               x + 6  
      // gfun(x) = [- -------]
      //              2 x - 6 
      // x0 = 3
      dfac = -2.0;
      gfun_n[0] = 6.0;
      gfun_n[1] = 1.0;
      roots[0] = 3.0;
    } else if (dim_numer == 2 && dim_denom == 2 && pade_old) {
      //                 12
      // gfun(x) = [-------------]
      //             2
      //            x  - 6 x + 12
      // x0 =  3 - sqrt(3) %i
      // x1 =  3 + sqrt(3) %i
      dfac = 1.0;
      gfun_n[0] = 12.0;
      gfun_n[1] = 0.0;
      gfun_n[2] = 0.0;
      roots[0] = 3.0 - sqrt(3.0) * IUNIT;
      roots[1] = 3.0 + sqrt(3.0) * IUNIT;
    } else if (dim_numer == 2 && dim_denom == 2) {
      //               2              
      //              x  + 6 x + 60    
      // gfun(x) =  [----------------] 
      //                2              
      //             3 x  - 24 x + 60 
      // x0 = 4 - 2 %i,
      // x1 = 4 + 2 %i
      dfac = 3.0;
      gfun_n[0] = 60.0;
      gfun_n[1] = 6.0;
      gfun_n[2] = 1.0;
      roots[0] =  4.0 - 2.0 * IUNIT;
      roots[1] =  4.0 + 2.0 * IUNIT;
    } else if (dim_numer == 3 && dim_denom == 3 && pade_old) {
      //(%i101) pade(taylor(gfun33(x),x,0,6),3,3);
      //                      2
      //                   2 x  + 120
      // gfun(x) = [- -----------------------]
      //               3       2
      //              x  - 12 x  + 60 x - 120
      // x0 = +4.644370709252173768 +0.000000000000000000 %i
      // x1 = +3.677814645373918445 +3.508761919567442433 %i
      // x2 = +3.677814645373918445 -3.508761919567442433 %i
      dfac = -1.0;
      gfun_n[0] = 120.0;
      gfun_n[1] = 0.0;
      gfun_n[2] = 2.0;
      gfun_n[3] = 0.0;
      roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
      roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
      roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
    } else if (dim_numer == 3 && dim_denom == 3) {
      //                3       2                
      //               x  + 20 x  + 60 x + 840   
      // gfun(x) = [- --------------------------]
      //                 3       2               
      //              4 x  - 60 x  + 360 x - 840 
      // x0 = +5.648485971016884655 +0.000000000000000000
      // x1 = +4.675757014491553676 +3.913489560603714335
      // x2 = +4.675757014491553676 -3.913489560603714335
      dfac = -4.0;
      gfun_n[0] = 840.0;
      gfun_n[1] = 60.0;
      gfun_n[2] = 20.0;
      gfun_n[3] = 1.0;
      roots[0] = +5.648485971016884655 +0.000000000000000000 * IUNIT;
      roots[1] = +4.675757014491553676 +3.913489560603714335 * IUNIT;
      roots[2] = +4.675757014491553676 -3.913489560603714335 * IUNIT;
    } else {
      std::cout << "cletdrk4:gen_exp. bad dim_numer or dim_denom." << std::endl;
      abort();
    }

    gfun_d1[0] = +dfac;
    gfun_d0[0] = -dfac * roots[0];
    for (int idim = 1; idim < dim_denom; idim ++) {
      gfun_d1[idim] = +RUNIT;
      gfun_d0[idim] = -roots[idim];
    }
  }

  dcomplex tunit, tfac, ttmp;
  std::vector<dcomplex> cn;
  std::vector<dcomplex> cd0;
  std::vector<dcomplex> cd1;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // ##### dt/2 * gfun(-i*h*dt/2) #####
  cn = gfun_n;
  cd0 = gfun_d0;
  cd1 = gfun_d1;
  tfac = tunit * etd_dt2;
  ttmp = RUNIT * etd_dt2; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  gfun.gen(MPIP, IO, Bas, etd_dt2, -1, dim_numer, dim_denom, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# gfun: real\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(gfun_n, gfun_d0, gfun_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# gfun: imag\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(gfun_n, gfun_d0, gfun_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::gen_f0fun(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  f0funrk = etd_dt1 / SIX;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> f0fun_n(max_dim_numer + 1);
  std::vector<dcomplex> f0fun_d0(max_dim_denom);
  std::vector<dcomplex> f0fun_d1(max_dim_denom);

  if (dim_denom == 0) {
    // 8-th order Taylor expansion of f0fun(x) = (-4-x+exp(x)*(4-3x+x^2))/x^3 generated by Maxima.
    //               1      8      1     7      7     6     1     5     5     4    1    3    3    2    1      1
    // f0fun(x) = (------) x  + (-----) x  + (-----) x  + (----) x  + (----) x  + (--) x  + (--) x  + (-) x + -
    //             492800        56700        51840        1120        1008        45        40        6      6
    f0fun_n[8] = 1.0 / 492800.0;
    f0fun_n[7] = 1.0 / 56700.0; 
    f0fun_n[6] = 7.0 / 51840.0; 
    f0fun_n[5] = 1.0 / 1120.0;  
    f0fun_n[4] = 5.0 / 1008.0;  
    f0fun_n[3] = 1.0 / 45.0;    
    f0fun_n[2] = 3.0 / 40.0;    
    f0fun_n[1] = 1.0 / 6.0;     
    f0fun_n[0] = 1.0 / 6.0;     
  } else {
    if (dim_numer == 1 && dim_denom == 1) {
      //               11 x + 20  
      // f0fun(x) = [- ----------]
      //               54 x - 120 
      // x0 = 120 / 54
      dfac = -54.0;
      f0fun_n[0] = 20.0;
      f0fun_n[1] = 11.0;
      roots[0] = 120.0 / 54.0;
    } else if (dim_numer == 2 && dim_denom == 2 && pade_old) {
      //                 x + 2
      // f0fun(x) = [-------------]
      //              2
      //             x  - 6 x + 12
      // x0 =  3 - sqrt(3) %i
      // x1 =  3 + sqrt(3) %i
      dfac = 1.0;
      f0fun_n[0] = 2.0;
      f0fun_n[1] = 1.0;
      f0fun_n[2] = 0.0;
      roots[0] = 3.0 - sqrt(3.0) * IUNIT;
      roots[1] = 3.0 + sqrt(3.0) * IUNIT;
    } else if (dim_numer == 2 && dim_denom == 2) {
      //                    2                    
      //              2657 x  + 19620 x + 34860  
      // f0fun(x) = [---------------------------] 
      //                    2                    
      //             13260 x  - 91440 x + 209160 
      //      762 - sqrt(189762) %i 
      // x0 = ---------------------,
      //               221          
      //      762 + sqrt(189762) %i 
      // x1 = ---------------------
      //               221          
      dfac = 13260.0;
      f0fun_n[0] = 34860.0;
      f0fun_n[1] = 19620.0;
      f0fun_n[2] = 2657.0;
      roots[0] = (762.0 - sqrt(189762.0) * IUNIT) / 221.0;
      roots[1] = (762.0 + sqrt(189762.0) * IUNIT) / 221.0;
    } else if (dim_numer == 3 && dim_denom == 3 && pade_old) {
      //(%i102) pade(taylor(f0fun33(x),x,0,6),3,3);
      //                  2
      //                 x  + 10 x + 20
      // f0fun(x) = [- -----------------------]
      //                3       2
      //               x  - 12 x  + 60 x - 120
      // x0 = +4.644370709252173768 +0.000000000000000000 %i
      // x1 = +3.677814645373918445 +3.508761919567442433 %i
      // x2 = +3.677814645373918445 -3.508761919567442433 %i
      dfac = -1.0;
      f0fun_n[0] = 20.0;
      f0fun_n[1] = 10.0;
      f0fun_n[2] = 1.0;
      f0fun_n[3] = 0.0;
      roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
      roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
      roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
    } else if (dim_numer == 3 && dim_denom == 3) {
      //                      3           2                           
      //               68401 x  + 986832 x  + 6549900 x + 11788560    
      // f0fun(x) = [- ----------------------------------------------]
      //                      3            2                          
      //              399900 x  - 5523840 x  + 31431960 x - 70731360  
      // x0 = +5.264996797311573218 +0.000000000000000000 %i
      // x1 = +4.274028233002134414 +3.914938394690479662 %i
      // x2 = +4.274028233002134414 -3.914938394690479662 %i
      dfac = -399900.0;
      f0fun_n[0] = 11788560.0;
      f0fun_n[1] = 6549900.0;
      f0fun_n[2] = 986832.0;
      f0fun_n[3] = 68401.0;
      roots[0] = +5.264996797311573218 +0.000000000000000000 * IUNIT;
      roots[1] = +4.274028233002134414 +3.914938394690479662 * IUNIT;
      roots[2] = +4.274028233002134414 -3.914938394690479662 * IUNIT;
    } else {
      std::cout << "cletdrk4:gen_exp. bad dim_numer or dim_denom." << std::endl;
      abort();
    }

    f0fun_d1[0] = +dfac;
    f0fun_d0[0] = -dfac * roots[0];
    for (int idim = 1; idim < dim_denom; idim ++) {
      f0fun_d1[idim] = +RUNIT;
      f0fun_d0[idim] = -roots[idim];
    }
  }

  dcomplex tunit, tfac, ttmp;
  std::vector<dcomplex> cn;
  std::vector<dcomplex> cd0;
  std::vector<dcomplex> cd1;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // ##### dt * f0fun(-i*h*dt) #####
  cn = f0fun_n;
  cd0 = f0fun_d0;
  cd1 = f0fun_d1;
  tfac = tunit * etd_dt1;
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  f0fun.gen(MPIP, IO, Bas, etd_dt1, -1, dim_numer, dim_denom, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# f0fun: real\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (-4.0-x+exp(x)*(4.0-3.0*x+x*x))/(x*x*x);
  //   vala = test_pade(f0fun_n, f0fun_d0, f0fun_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# f0fun: imag\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (-4.0-x+exp(x)*(4.0-3.0*x+x*x))/(x*x*x);
  //   vala = test_pade(f0fun_n, f0fun_d0, f0fun_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::gen_f1fun(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  f1funrk = etd_dt1 / THREE;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> f1fun_n(max_dim_numer + 1);
  std::vector<dcomplex> f1fun_d0(max_dim_denom);
  std::vector<dcomplex> f1fun_d1(max_dim_denom);

  if (dim_denom == 0) {
    // 8-th order Taylor expansion of f1fun(x) = 2*(2+x+exp(x)*(x-2))/x^3 generated by Maxima.
    //                1      8      1      7      1     6     1     5     1    4    1    3    1    2    1      1
    // f1fun(x) = (-------) x  + (------) x  + (-----) x  + (----) x  + (---) x  + (--) x  + (--) x  + (-) x + -
    //             2217600        226800        25920        3360        504        90        20        6      3                                   
    f1fun_n[8] = 1.0 / 2217600.0;
    f1fun_n[7] = 1.0 / 226800.0; 
    f1fun_n[6] = 1.0 / 25920.0;  
    f1fun_n[5] = 1.0 / 3360.0;   
    f1fun_n[4] = 1.0 / 504.0;    
    f1fun_n[3] = 1.0 / 90.0;     
    f1fun_n[2] = 1.0 / 20.0;     
    f1fun_n[1] = 1.0 / 6.0;      
    f1fun_n[0] = 1.0 / 3.0;      
  } else {
    if (dim_numer == 1 && dim_denom == 1) {
      //               2 x + 10 
      // f1fun(x) = [- --------]
      //               9 x - 30 
      // x0 = 30 / 9
      dfac = -9.0;
      f1fun_n[0] = 10.0;
      f1fun_n[1] = 2.0;
      roots[0] = 30.0 / 9.0;
    } else if (dim_numer == 2 && dim_denom == 2 && pade_old) {
      //                   4
      // f1fun(x) = [-------------]
      //              2
      //             x  - 6 x + 12
      // x0 =  3 - sqrt(3) %i
      // x1 =  3 + sqrt(3) %i
      dfac = 1.0;
      f1fun_n[0] = 4.0;
      f1fun_n[1] = 0.0;
      f1fun_n[2] = 0.0;
      roots[0] = 3.0 - sqrt(3.0) * IUNIT;
      roots[1] = 3.0 + sqrt(3.0) * IUNIT;
    } else if (dim_numer == 2 && dim_denom == 2) {
      //                  2                 
      //              41 x  + 450 x + 2940  
      // f1fun(x) = [----------------------] 
      //                  2                 
      //             330 x  - 3060 x + 8820 
      //      51 - sqrt(633) %i 
      // x0 = -----------------,
      //             11         
      //      51 + sqrt(633) %i 
      // x1 = ----------------- 
      //             11         
      dfac = 330.0;
      f1fun_n[0] = 2940.0;
      f1fun_n[1] = 450.0;
      f1fun_n[2] = 41.0;
      roots[0] = (51.0 - sqrt(633.0) * IUNIT) / 11.0;
      roots[1] = (51.0 + sqrt(633.0) * IUNIT) / 11.0;
    } else if (dim_numer == 3 && dim_denom == 3 && pade_old) {
      //(%i103) pade(taylor(f1fun33(x),x,0,6),3,3);
      //                       40
      // f1fun(x) = [- -----------------------]
      //                3       2
      //               x  - 12 x  + 60 x - 120
      // x0 = +4.644370709252173768 +0.000000000000000000 %i
      // x1 = +3.677814645373918445 +3.508761919567442433 %i
      // x2 = +3.677814645373918445 -3.508761919567442433 %i
      dfac = -1.0;
      f1fun_n[0] = 40.0;
      f1fun_n[1] = 0.0;
      f1fun_n[2] = 0.0;
      f1fun_n[3] = 0.0;
      roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
      roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
      roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
    } else if (dim_numer == 3 && dim_denom == 3) {
      //                     3         2                     
      //                 73 x  + 1386 x  + 10920 x + 93240   
      // f1fun(x) = [- -------------------------------------]
      //                    3          2                     
      //               915 x  - 15750 x  + 107100 x - 279720 
      // x0 = +6.374187319550069120 +0.000000000000000000 %i
      // x1 = +5.419463717274148884 +4.311524038457500652 %i
      // x2 = +5.419463717274148884 -4.311524038457500652 %i
      dfac = -915.0;
      f1fun_n[0] = 93240.0;
      f1fun_n[1] = 10920.0;
      f1fun_n[2] = 1386.0;
      f1fun_n[3] = 73.0;
      roots[0] = +6.374187319550069120 +0.000000000000000000 * IUNIT;
      roots[1] = +5.419463717274148884 +4.311524038457500652 * IUNIT;
      roots[2] = +5.419463717274148884 -4.311524038457500652 * IUNIT;
    } else {
      std::cout << "cletdrk4:gen_exp. bad dim_numer or dim_denom." << std::endl;
      abort();
    }

    f1fun_d1[0] = +dfac;
    f1fun_d0[0] = -dfac * roots[0];
    for (int idim = 1; idim < dim_denom; idim ++) {
      f1fun_d1[idim] = +RUNIT;
      f1fun_d0[idim] = -roots[idim];
    }
  }

  dcomplex tunit, tfac, ttmp;
  std::vector<dcomplex> cn;
  std::vector<dcomplex> cd0;
  std::vector<dcomplex> cd1;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // ##### dt * f1fun(-i*h*dt) #####
  cn = f1fun_n;
  cd0 = f1fun_d0;
  cd1 = f1fun_d1;
  tfac = tunit * etd_dt1;
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  f1fun.gen(MPIP, IO, Bas, etd_dt1, -1, dim_numer, dim_denom, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# f1fun: real\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = 2.0*(2.0+x+exp(x)*(x-2.0))/(x*x*x);
  //   vala = test_pade(f1fun_n, f1fun_d0, f1fun_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# f1fun: imag\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = 2.0*(2.0+x+exp(x)*(x-2.0))/(x*x*x);
  //   vala = test_pade(f1fun_n, f1fun_d0, f1fun_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrk4::gen_f2fun(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  f2funrk = etd_dt1 / SIX;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> f2fun_n(max_dim_numer + 1);
  std::vector<dcomplex> f2fun_d0(max_dim_denom);
  std::vector<dcomplex> f2fun_d1(max_dim_denom);

  if (dim_denom == 0) {
    // 8-th order Taylor expansion of f2fun(x) = (-4-3*x-x^2+exp(x)*(4-x))/x^3 generated by Maxima.
    //               - 1     8     - 1     7     - 1    6     - 1    5    - 1    4    - 1   3    - 1   2   1
    // f2fun(x) = (-------) x  + (------) x  + (-----) x  + (-----) x  + (----) x  + (---) x  + (---) x  + -
    //             5702400        604800        72576        10080        1680        360        120       6
    f2fun_n[8] = -1.0 / 5702400.0;
    f2fun_n[7] = -1.0 / 604800.0; 
    f2fun_n[6] = -1.0 / 72576.0;  
    f2fun_n[5] = -1.0 / 10080.0;  
    f2fun_n[4] = -1.0 / 1680.0;   
    f2fun_n[3] = -1.0 / 360.0;    
    f2fun_n[2] = -1.0 / 120.0;    
    f2fun_n[1] =  0.0;              
    f2fun_n[0] =  1.0 / 6.0;      
  } else {
    if (dim_numer == 1 && dim_denom == 1) {
      //                            0 x + 1
      // f2fun(x) : impossible ==> [-------]
      //                            0 x + 6 
      std::cout << "Pade(1,1) is impossible for f2fun. Use zero'th order approximant 1/6" << std::endl;
      dfac = 0.0;
      f2fun_n[0] = 1.0 / 6.0;
      f2fun_n[1] = 0.0;
      roots[0] = 0.0;
    } else if (dim_numer == 2 && dim_denom == 2 && pade_old) {
      //                   x - 2
      // f2fun(x) = [- -------------]
      //                2
      //               x  - 6 x + 12
      // x0 =  3 - sqrt(3) %i
      // x1 =  3 + sqrt(3) %i
      dfac = 1.0;
      f2fun_n[0] = 2.0;
      f2fun_n[1] = -1.0;
      f2fun_n[2] = 0.0;
      roots[0] = 3.0 - sqrt(3.0) * IUNIT;
      roots[1] = 3.0 + sqrt(3.0) * IUNIT;
    } else if (dim_numer == 2 && dim_denom == 2) {
      //                    2                 
      //                13 x  + 420 x - 1260  
      // f2fun(x) = [- ----------------------] 
      //                    2                 
      //               300 x  - 2520 x + 7560 
      //     21 - 3 sqrt(21) %i 
      // x = ------------------,
      //             5          
      //     21 + 3 sqrt(21) %i 
      // x = ------------------]
      //             5
      dfac = 300.0;
      f2fun_n[0] = 1260.0;
      f2fun_n[1] = -420.0;
      f2fun_n[2] = -13.0;
      roots[0] = (21.0 - 3.0 * sqrt(21.0) * IUNIT) / 5.0;
      roots[1] = (21.0 + 3.0 * sqrt(21.0) * IUNIT) / 5.0;
    } else if (dim_numer == 3 && dim_denom == 3 && pade_old) {
      //(%i104) pade(taylor(f2fun33(x),x,0,6),3,3);
      //                  2
      //                 x  - 10 x + 20
      // f2fun(x) = [- -----------------------]
      //                3       2
      //               x  - 12 x  + 60 x - 120
      // x0 = +4.644370709252173768 +0.000000000000000000 %i
      // x1 = +3.677814645373918445 +3.508761919567442433 %i
      // x2 = +3.677814645373918445 -3.508761919567442433 %i
      dfac = -1.0;
      f2fun_n[0] = 20.0;
      f2fun_n[1] = -10.0;
      f2fun_n[2] = 1.0;
      f2fun_n[3] = 0.0;
      roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
      roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
      roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
    } else if (dim_numer == 3 && dim_denom == 3) {
      //                   3        2                     
      //               25 x  - 336 x  + 13860 x - 35280   
      // f2fun(x) = [------------------------------------]
      //                  3          2                    
      //             780 x  - 12600 x  + 83160 x - 211680 
      // x0 = +5.865168732584368882 +0.000000000000000000 %i
      // x1 = +5.144338710630889544 +4.450430848227815872 %i
      // x2 = +5.144338710630889544 -4.450430848227815872 %i
      dfac = 780.0;
      f2fun_n[0] = -35280.0;
      f2fun_n[1] = 13860.0;
      f2fun_n[2] = -336.0;
      f2fun_n[3] = 25.0;
      roots[0] = +5.865168732584368882 +0.000000000000000000 * IUNIT;
      roots[1] = +5.144338710630889544 +4.450430848227815872 * IUNIT;
      roots[2] = +5.144338710630889544 -4.450430848227815872 * IUNIT;
    } else {
      std::cout << "cletdrk4:gen_exp. bad dim_numer or dim_denom." << std::endl;
      abort();
    }

    f2fun_d1[0] = +dfac;
    f2fun_d0[0] = -dfac * roots[0];
    for (int idim = 1; idim < dim_denom; idim ++) {
      f2fun_d1[idim] = +RUNIT;
      f2fun_d0[idim] = -roots[idim];
    }
  }

  dcomplex tunit, tfac, ttmp;
  std::vector<dcomplex> cn;
  std::vector<dcomplex> cd0;
  std::vector<dcomplex> cd1;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // ##### dt * f2fun(-i*h*dt) #####
  cn = f2fun_n;
  cd0 = f2fun_d0;
  cd1 = f2fun_d1;
  tfac = tunit * etd_dt1;
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  f2fun.gen(MPIP, IO, Bas, etd_dt1, -1, dim_numer, dim_denom, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# f2fun: real\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (-4.0-3.0*x-x*x+exp(x)*(4.0-x))/(x*x*x);
  //   vala = test_pade(f2fun_n, f2fun_d0, f2fun_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# f2fun: imag\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (-4.0-3.0*x-x*x+exp(x)*(4.0-x))/(x*x*x);
  //   vala = test_pade(f2fun_n, f2fun_d0, f2fun_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
dcomplex cletdrk4::test_pade(const std::vector<dcomplex>& pnum, 
			     const std::vector<dcomplex>& pd0, 
			     const std::vector<dcomplex>& pd1, 
			     dcomplex x) const
{
  dcomplex valp, tmp, tmp_tay;
  tmp = pnum[0];
  tmp_tay = RUNIT;
  for (int inum = 1; inum <= dim_numer; inum ++) {
    tmp_tay *= x;
    tmp += pnum[inum] * tmp_tay;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    tmp /= (pd0[idenom] + pd1[idenom] * x);
  }
  valp = tmp;
  return valp;
}
////////////////////////////////////////////////////////////////////////
