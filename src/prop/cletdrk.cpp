////////////////////////////////////////////////////////////////////////
// Exponential Runge-Kuta with Pade approximation for phi functions
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
int cletdrk::max_dim_numer = 3;
int cletdrk::max_dim_denom = 3;
////////////////////////////////////////////////////////////////////////
cletdrk::cletdrk()
{
  std::cout << "cletdrk" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cletdrk::~cletdrk()
{
  std::cout << "~cletdrk" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cletdrk::cletdrk(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		   const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void cletdrk::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		   const clfield& Field, const clhprod& HPW)
{
  etd_dt1 = Field.dtime;
  etd_dt2 = Field.dtime * HALF;

  IO.read_info("etd_rk_order", 4, rk_order);
  IO.read_info("etd_dim_numer", 2, dim_numer);
  IO.read_info("etd_dim_denom", 3, dim_denom);
  IO.read_info("etd_pade_equiv", false, pade_equiv);
  IO.read_info("etd_rosenbrock", false, rosenbrock);

  if (rk_order < 1 || rk_order > 4) {
    std::cout << "cletdrk::gen bad rk_order." << std::endl;
    abort();
  }

  gen_phi1(MPIP, IO, Bas);
  gen_phi2(MPIP, IO, Bas);
  gen_phi3(MPIP, IO, Bas);

  //DEBUG
  printf("cletdrk::phi1dt1:\n"); phi1dt1.print();
  printf("cletdrk::phi1dt2:\n"); phi1dt2.print();
  printf("cletdrk::phi2dt1:\n"); phi2dt1.print();
  printf("cletdrk::phi2dt2:\n"); phi2dt2.print();
  printf("cletdrk::phi3dt1:\n"); phi3dt1.print();
  printf("cletdrk::phi3dt2:\n"); phi3dt2.print();
  //  std::cout << "abort for debug in cletdrk::gen." << std::endl;
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
void cletdrk::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  std::cout << "cletdrk::prop nyi." << std::endl;
  abort();
}
////////////////////////////////////////////////////////////////////////
void cletdrk::prop(const clmpi& Proc, const clbas& Bas,
		   const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  if (rk_order == 1) {
    prop1(Proc, Bas, Field, HPW, Wfn);
  } else if (rk_order == 2) {
    prop2(Proc, Bas, Field, HPW, Wfn);
  } else if (rk_order == 3) {
    prop3(Proc, Bas, Field, HPW, Wfn);
  } else if (rk_order == 4) {
    prop4(Proc, Bas, Field, HPW, Wfn);
  } else {
    std::cout << "cletdrk::prop bad rk_order." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void cletdrk::prop1(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  double lfield[9];
  double thalf = Field.time + etd_dt2;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // stage 0: -i*h*u0 + N[u0]
  HPW.copy(Proc, Bas, Wfn, Wfn0);
  Field.get_value(thalf, lfield);
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn0, tWfn);
  Field.get_value(Field.time, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfn0, tWfn);
  prod2(Proc, Bas, HPW, thalf, Field, phi1dt1, phi1rk, tWfn, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cletdrk::prop2(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  double lfield[9];
  double thalf = Field.time + etd_dt2;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // stage 0: -i*h*u0 + N[u0]
  HPW.copy(Proc, Bas, Wfn, Wfn0);
  Field.get_value(thalf, lfield);
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn0, tWfn);
  Field.get_value(Field.time, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfn0, tWfn);
  // u0 + phi1dt1(-i*h*u0 + N[u0]) ==> Wfn1, Wfn
  prod1(Proc, Bas, HPW, thalf, Field, phi1dt1, phi1rk, tWfn, Wfn1);
  HPW.axpyz(Proc, Bas, RUNIT, Wfn1, Wfn0, Wfn1);
  HPW.copy(Proc, Bas, Wfn1, Wfn);

  // stage 1:
  Field.get_value(Field.time + etd_dt1, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
  HPW.xmyz(Proc, Bas, hWfn1, hWfn0, tWfn);
  // phi2dt1(N[u1] - N[u0]) +=> Wfn
  prod2(Proc, Bas, HPW, thalf, Field, phi2dt1, phi2rk, tWfn, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cletdrk::prop3(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  double lfield[9];
  double thalf = Field.time + etd_dt2;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // tWfn = -i*h*u0 + N[u0]
  HPW.copy(Proc, Bas, Wfn, Wfn0);
  Field.get_value(thalf, lfield);
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn0, tWfn);
  Field.get_value(Field.time, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfn0, tWfn);

  // Wfn1 = u0 + phi1dt2(-i*h*u0 + N[u0])/2
  prod1(Proc, Bas, HPW, thalf, Field, phi1dt2, phi1rk, tWfn, Wfn1);
  HPW.axpyz(Proc, Bas, CHALF, Wfn1, Wfn0, Wfn1);
  // Wfn2 = Wfn = u0 + phi1dt1(-i*h*u0 + N[u0])
  prod1(Proc, Bas, HPW, thalf, Field, phi1dt1, phi1rk, tWfn, Wfn2);
  HPW.axpyz(Proc, Bas, RUNIT, Wfn2, Wfn0, Wfn2);
  HPW.copy(Proc, Bas, Wfn2, Wfn);

  // tWfn = Wfn3 = (N[u1] - N[u0]) * 2
  Field.get_value(Field.time + etd_dt2, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
  HPW.xmyz(Proc, Bas, hWfn1, hWfn0, tWfn);
  HPW.scal(Proc, Bas, CTWO, tWfn);
  HPW.copy(Proc, Bas, tWfn, Wfn3); // Wfn3=tWfn is a scratch

  // Wfn2 += phi2dt2(N[u1] - N[u0])*2
  prod2(Proc, Bas, HPW, thalf, Field, phi2dt2, phi2rk, tWfn, Wfn2);
  // Wfn2 += phi2dt1(N[u1] - N[u0])*2
  prod2(Proc, Bas, HPW, thalf, Field, phi2dt1, phi2rk, Wfn3, Wfn2);

  // hWfn2 = (-5N[u0] + 4N[u1] + N[u2]) / 3
  Field.get_value(Field.time + etd_dt1, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn2, hWfn2);
  HPW.axpy(Proc, Bas, -CFIVE, hWfn0, hWfn2);
  HPW.axpy(Proc, Bas, +CFOUR, hWfn1, hWfn2);
  HPW.scal(Proc, Bas, CTHIRD, hWfn2);

  // Wfn += phi2dt1(-5N[u0] + 4N[u1] + N[u2])/3
  prod2(Proc, Bas, HPW, thalf, Field, phi2dt1, phi2rk, hWfn2, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cletdrk::prop4(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  double lfield[9];
  double thalf = Field.time + etd_dt2;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // stage 0: -i*h*u0 + N[u0]
  HPW.copy(Proc, Bas, Wfn, Wfn0);
  Field.get_value(thalf, lfield);
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn0, tWfn);
  Field.get_value(Field.time, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  if (rosenbrock) {
    Field.get_value(Field.time, thalf, ONE, -ONE, lfield);
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn0, hWfn0);
  }
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfn0, tWfn);

  // ==> Wfn1, Wfn2
  prod1(Proc, Bas, HPW, thalf, Field, phi1dt2, phi1rk, tWfn, Wfn1);
  HPW.axpyz(Proc, Bas, CHALF, Wfn1, Wfn0, Wfn1);
  HPW.copy(Proc, Bas, Wfn1, Wfn2);
  // ==> Wfn3, Wfn(t+dt)
  prod1(Proc, Bas, HPW, thalf, Field, phi1dt1, phi1rk, tWfn, Wfn3);
  HPW.xpyz(Proc, Bas, Wfn0, Wfn3, Wfn3);
  HPW.copy(Proc, Bas, Wfn3, Wfn);

  // stage 1
  Field.get_value(Field.time + etd_dt2, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
  HPW.xmyz(Proc, Bas, hWfn1, hWfn0, tWfn);
  prod2(Proc, Bas, HPW, thalf, Field, phi2dt2, phi2rk, tWfn, Wfn2);

  // stage 2
  Field.get_value(Field.time + etd_dt2, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn2, hWfn2);
  HPW.xmyz(Proc, Bas, hWfn2, hWfn0, tWfn);
  HPW.scal(Proc, Bas, CTWO, tWfn);
  prod2(Proc, Bas, HPW, thalf, Field, phi2dt1, phi2rk, tWfn, Wfn3);

  // stage 3
  Field.get_value(Field.time + etd_dt1, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn3, hWfn3);
  if (rosenbrock) {
    Field.get_value(Field.time + etd_dt1, thalf, ONE, -ONE, lfield);
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn3, hWfn3);
  }

  // summation
  HPW.axpbyz(Proc, Bas, -CTHREE, hWfn0, CTWO, hWfn1, tWfn);
  HPW.axpy(Proc, Bas, CTWO, hWfn2, tWfn);
  HPW.axpy(Proc, Bas, -RUNIT, hWfn3, tWfn);
  prod2(Proc, Bas, HPW, thalf, Field, phi2dt1, phi2rk, tWfn, Wfn);

  HPW.axpbyz(Proc, Bas, CFOUR, hWfn0, -CFOUR, hWfn1, tWfn);
  HPW.axpy(Proc, Bas, -CFOUR, hWfn2, tWfn);
  HPW.axpy(Proc, Bas, +CFOUR, hWfn3, tWfn);
  prod2(Proc, Bas, HPW, thalf, Field, phi3dt1, phi3rk, tWfn, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cletdrk::prod1(const clmpi& Proc, const clbas& Bas, clhprod& HPW, 
		    double time, const clfield& Field, const clh1rat& FacOrb, 
		    dcomplex FacCI, const clwfn& WIn, clwfn& WOut)
{
  HPW.copy(Proc, Bas, WIn, WOut);
  HPW.scalc(Proc, Bas, FacCI, WOut);
  FacOrb.prod(Proc, Bas, time, Field, HPW, WOut);
}
////////////////////////////////////////////////////////////////////////
void cletdrk::prod2(const clmpi& Proc, const clbas& Bas, clhprod& HPW, 
		     double time, const clfield& Field, const clh1rat& FacOrb, 
		     dcomplex FacCI, clwfn& WIn, clwfn& WOut)
{
  HPW.scalc(Proc, Bas, FacCI, WIn);
  FacOrb.prod(Proc, Bas, time, Field, HPW, WIn);
  HPW.xpyz(Proc, Bas, WIn, WOut, WOut);
}
////////////////////////////////////////////////////////////////////////
void cletdrk::gen_phi1(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  phi1rk = etd_dt1;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi1_n(max_dim_numer + 1);
  std::vector<dcomplex> phi1_d0(max_dim_denom);
  std::vector<dcomplex> phi1_d1(max_dim_denom);

  if (dim_numer == 1 && dim_denom == 1) {
    //               x + 6  
    // phi1(x) = [- -------]
    //              2 x - 6 
    // x0 = 3
    dfac = -2.0;
    phi1_n[0] = 6.0;
    phi1_n[1] = 1.0;
    roots[0] = 3.0;
  } else if (dim_numer == 2 && dim_denom == 2) {
    //               2              
    //              x  + 6 x + 60    
    // phi1(x) =  [----------------] 
    //                2              
    //             3 x  - 24 x + 60 
    dfac = 3.0;
    phi1_n[0] = 60.0;
    phi1_n[1] = 6.0;
    phi1_n[2] = 1.0;
    roots[0] =  4.0 - 2.0 * IUNIT;
    roots[1] =  4.0 + 2.0 * IUNIT;
  } else if (dim_numer == 2 && dim_denom == 3) {
    //                      2		    
    //                   2 x  + 120	    
    // phi1(x) = - -----------------------
    //              3       2		    
    //             x  - 12 x  + 60 x - 120
    dfac = -1.0;
    phi1_n[0] = 120.0;
    phi1_n[1] = 0.0;
    phi1_n[2] = 2.0;
    roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
    roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
    roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
  } else if (dim_numer == 3 && dim_denom == 3 && pade_equiv) {
    //                      2		    
    //                   2 x  + 120	    
    // phi1(x) = - -----------------------
    //              3       2		    
    //             x  - 12 x  + 60 x - 120
    dfac = -1.0;
    phi1_n[0] = 120.0;
    phi1_n[1] = 0.0;
    phi1_n[2] = 2.0;
    phi1_n[3] = 0.0;
    roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
    roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
    roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
  } else if (dim_numer == 3 && dim_denom == 3) {
    //                3       2                
    //               x  + 20 x  + 60 x + 840   
    // phi1(x) = [- --------------------------]
    //                 3       2               
    //              4 x  - 60 x  + 360 x - 840 
    dfac = -4.0;
    phi1_n[0] = 840.0;
    phi1_n[1] = 60.0;
    phi1_n[2] = 20.0;
    phi1_n[3] = 1.0;
    roots[0] = +5.648485971016884655 +0.000000000000000000 * IUNIT;
    roots[1] = +4.675757014491553676 +3.913489560603714335 * IUNIT;
    roots[2] = +4.675757014491553676 -3.913489560603714335 * IUNIT;
  } else {
    std::cout << "cletdrk:gen_exp. bad dim_numer or dim_denom." << std::endl;
    abort();
  }

  phi1_d1[0] = +dfac;
  phi1_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < dim_denom; idim ++) {
    phi1_d1[idim] = +RUNIT;
    phi1_d0[idim] = -roots[idim];
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

  // ##### dt * phi1(-i*h*dt) #####
  cn = phi1_n;
  cd0 = phi1_d0;
  cd1 = phi1_d1;
  tfac = tunit * etd_dt1; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi1dt1.gen(MPIP, IO, Bas, etd_dt1, -1, dim_numer, dim_denom, cn, cd0, cd1);

  // ##### dt * phi1(-i*h*dt/2) #####
  cn = phi1_n;
  cd0 = phi1_d0;
  cd1 = phi1_d1;
  tfac = tunit * etd_dt2; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi1dt2.gen(MPIP, IO, Bas, etd_dt2, -1, dim_numer, dim_denom, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# phi1: real\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi1_n, phi1_d0, phi1_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# phi1: imag\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi1_n, phi1_d0, phi1_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrk::gen_phi2(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  phi2rk = etd_dt1 / TWO;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi2_n(max_dim_numer + 1);
  std::vector<dcomplex> phi2_d0(max_dim_denom);
  std::vector<dcomplex> phi2_d1(max_dim_denom);

  if (dim_numer == 1 && dim_denom == 1) {
    //               x + 12  
    // phi2(x) = [- --------]
    //              6 x - 24 
    dfac = -6.0;
    phi2_n[0] = 12.0;
    phi2_n[1] = 1.0;
    roots[0] = 4.0;
  } else if (dim_numer == 2 && dim_denom == 2) {
    //                  2	       
    //                 x  + 180           
    // phi2(x) = [-------------------] 
    //                2	           
    //            12 x  - 120 x + 360 
    dfac = 12.0;
    phi2_n[0] = 180.0;
    phi2_n[1] = 0.0;
    phi2_n[2] = 1.0;
    roots[0] = +5.000000000000000000 +2.236067977499789805 * IUNIT;
    roots[1] = +5.000000000000000000 -2.236067977499789805 * IUNIT;
  } else if (dim_numer == 2 && dim_denom == 3) {
    //                    2		     
    //                 5 x  - 40 x + 420     	    
    // phi2(x) = - --------------------------
    //                3       2		     	    
    //             4 x  - 60 x  + 360 x - 840
    dfac = -4.0;
    phi2_n[0] = 420.0;
    phi2_n[1] = -40.0;
    phi2_n[2] = 5.0;
    roots[0] = +5.648485971016884655 +0.000000000000000000 * IUNIT;
    roots[1] = +4.675757014491553676 +3.913489560603714335 * IUNIT;
    roots[2] = +4.675757014491553676 -3.913489560603714335 * IUNIT;
  } else if (dim_numer == 3 && dim_denom == 3 && pade_equiv) {
    //                  2		    
    //                 x  - 10 x + 60	    
    // phi2(x) = - -----------------------
    //              3       2		    
    //             x  - 12 x  + 60 x - 120
    dfac = -1.0;
    phi2_n[0] = 60.0;
    phi2_n[1] = -10.0;
    phi2_n[2] = 1.0;
    phi2_n[3] = 0.0;
    roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
    roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
    roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
  } else if (dim_numer == 3 && dim_denom == 3) {
    //                 3       2		   
    //                x  + 40 x  - 140 x + 3360	   
    // phi2(x) = [- ------------------------------]
    //                  3        2		   
    //              20 x  - 360 x  + 2520 x - 6720 
    dfac = -20.0;
    phi2_n[0] = 3360.0;
    phi2_n[1] = -140.0;
    phi2_n[2] = 40.0;
    phi2_n[3] = 1.0;
    roots[0] = +6.651316808950019421 +0.000000000000000000 * IUNIT;
    roots[1] = +5.674341595524996507 +4.279971984629763249 * IUNIT;
    roots[2] = +5.674341595524996507 -4.279971984629763249 * IUNIT;
  } else {
    std::cout << "cletdrk:gen_exp. bad dim_numer or dim_denom." << std::endl;
    abort();
  }

  phi2_d1[0] = +dfac;
  phi2_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < dim_denom; idim ++) {
    phi2_d1[idim] = +RUNIT;
    phi2_d0[idim] = -roots[idim];
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

  // ##### dt * phi2(-i*h*dt) #####
  cn = phi2_n;
  cd0 = phi2_d0;
  cd1 = phi2_d1;
  tfac = tunit * etd_dt1; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi2dt1.gen(MPIP, IO, Bas, etd_dt1, -1, dim_numer, dim_denom, cn, cd0, cd1);

  // ##### dt * phi2(-i*h*dt/2) #####
  cn = phi2_n;
  cd0 = phi2_d0;
  cd1 = phi2_d1;
  tfac = tunit * etd_dt2; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi2dt2.gen(MPIP, IO, Bas, etd_dt2, -1, dim_numer, dim_denom, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# phi2: real\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi2_n, phi2_d0, phi2_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# phi2: imag\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi2_n, phi2_d0, phi2_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrk::gen_phi3(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  phi3rk = etd_dt1 / SIX;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi3_n(max_dim_numer + 1);
  std::vector<dcomplex> phi3_d0(max_dim_denom);
  std::vector<dcomplex> phi3_d1(max_dim_denom);

  if (dim_numer == 1 && dim_denom == 1) {
    //                x + 20   
    // phi3(x) = [- ----------]
    //              24 x - 120 
    dfac = -24.0;
    phi3_n[0] = 20.0;
    phi3_n[1] = 1.0;
    roots[0] = 5.0;
  } else if (dim_numer == 2 && dim_denom == 2) {
    //               2		       
    //              x  - 15 x + 420           
    // phi3(x) = [--------------------] 
    //                2		       
    //            60 x  - 720 x + 2520 
    dfac = 60.0;
    phi3_n[0] = 420.0;
    phi3_n[1] = -15.0;
    phi3_n[2] = 1.0;
    roots[0] = +6.000000000000000000 +2.449489742783177881 * IUNIT;
    roots[1] = +6.000000000000000000 -2.449489742783177881 * IUNIT;
  } else if (dim_numer == 2 && dim_denom == 3) {
    //                      2		         
    //                  11 x  - 140 x + 1120     
    // phi3(x) = - ------------------------------
    //                 3        2		 	     
    //             20 x  - 360 x  + 2520 x - 6720
    dfac = -20.0;
    phi3_n[0] = 1120.0;
    phi3_n[1] = -140.0;
    phi3_n[2] = 11.0;
    roots[0] = +6.651316808950019421 +0.000000000000000000 * IUNIT;
    roots[1] = +5.674341595524996507 +4.279971984629763249 * IUNIT;
    roots[2] = +5.674341595524996507 -4.279971984629763249 * IUNIT;
  } else if (dim_numer == 3 && dim_denom == 3 && pade_equiv) {
    //                    2		     
    //                   x  - 10 x + 40	     
    // phi3(x) = - --------------------------
    //                3       2		     
    //             2 x  - 24 x  + 120 x - 240
    dfac = -2.0;
    phi3_n[0] = 40.0;
    phi3_n[1] = -10.0;
    phi3_n[2] = 1.0;
    phi3_n[3] = 0.0;
    roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
    roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
    roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
  } else if (dim_numer == 3 && dim_denom == 3) {
    //                   3       2		       
    //                  x  + 84 x  - 840 x + 10080     	   
    // phi3(x) = [- ----------------------------------]
    //                   3         2		       
    //              120 x  - 2520 x  + 20160 x - 60480 
    dfac = -120.0;
    phi3_n[0] = 10080.0;
    phi3_n[1] = -840.0;
    phi3_n[2] = 84.0;
    phi3_n[3] = 1.0;
    roots[0] = +7.653383973560514697 +0.000000000000000000 * IUNIT;
    roots[1] = +6.673308013219741319 +4.617378364686957504 * IUNIT;
    roots[2] = +6.673308013219741319 -4.617378364686957504 * IUNIT;
  } else {
    std::cout << "cletdrk:gen_exp. bad dim_numer or dim_denom." << std::endl;
    abort();
  }

  phi3_d1[0] = +dfac;
  phi3_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < dim_denom; idim ++) {
    phi3_d1[idim] = +RUNIT;
    phi3_d0[idim] = -roots[idim];
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

  // ##### dt * phi3(-i*h*dt) #####
  cn = phi3_n;
  cd0 = phi3_d0;
  cd1 = phi3_d1;
  tfac = tunit * etd_dt1; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi3dt1.gen(MPIP, IO, Bas, etd_dt1, -1, dim_numer, dim_denom, cn, cd0, cd1);

  // ##### dt * phi3(-i*h*dt/2) #####
  cn = phi3_n;
  cd0 = phi3_d0;
  cd1 = phi3_d1;
  tfac = tunit * etd_dt2; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi3dt2.gen(MPIP, IO, Bas, etd_dt2, -1, dim_numer, dim_denom, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# phi3: real\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi3_n, phi3_d0, phi3_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# phi3: imag\n");
  // for (int ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi3_n, phi3_d0, phi3_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
