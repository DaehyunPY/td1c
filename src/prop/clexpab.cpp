////////////////////////////////////////////////////////////////////////
// ETD Adams-Bashforth with Pade approximation for phi functions
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
int clexpab::max_dim_numer = 4;
int clexpab::max_dim_denom = 4;
////////////////////////////////////////////////////////////////////////
clexpab::clexpab()
{
  std::cout << "clexpab" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clexpab::~clexpab()
{
  std::cout << "~clexpab" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clexpab::clexpab(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		 const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void clexpab::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		  const clfield& Field, const clhprod& HPW)
{
  std::cout << "clexpab::gen" << std::endl;
  clh2prop::gen_basic(MPIP, IO, Bas, Field, HPW);

  num_step = 0;
  expab_dt = Field.dtime;
  expab_cdt = expab_dt;

  IO.read_info("expab_pc", false, expab_pc);
  IO.read_info("expab_order", 4, expab_order);
  IO.read_info("expab_ndpade", 12, expab_ndpade);
  IO.read_info("expab_h1type", -1, expab_h1type);

  if (expab_order < 1 || expab_order > 5) {
    std::cout << "clexpab::gen bad expab_order." << std::endl;
    abort();
  }

  if (expab_ndpade == 11) {
    expab_npade = 1;
    expab_dpade = 1;
  } else if (expab_ndpade == 12) {
    expab_npade = 1;
    expab_dpade = 2;
  } else if (expab_ndpade == 13) {
    expab_npade = 2;
    expab_dpade = 3;
  } else if (expab_ndpade == 14) {
    expab_npade = 3;
    expab_dpade = 4;
  } else if (expab_ndpade == 22) {
    expab_npade = 2;
    expab_dpade = 2;
  } else if (expab_ndpade == 23) {
    expab_npade = 2;
    expab_dpade = 3;
  } else if (expab_ndpade == 24) {
    expab_npade = 3;
    expab_dpade = 4;
  } else if (expab_ndpade == 33) {
    expab_npade = 3;
    expab_dpade = 3;
  } else if (expab_ndpade == 34) {
    expab_npade = 3;
    expab_dpade = 4;
  } else if (expab_ndpade == 44) {
    expab_npade = 4;
    expab_dpade = 4;
  } else {
    std::cout << "clexpab::gen bad expab_ndpade." << std::endl;
    abort();
  }

  if (expab_order == 1) {
    gen_phi1(MPIP, IO, Bas);
  } else if (expab_order == 2) {
    //debug
    //    gen_phi2_phi0(MPIP, IO, Bas);
    //    gen_phi2_phi1(MPIP, IO, Bas);
    //debug
    gen_phi2(MPIP, IO, Bas);
  } else if (expab_order == 3) {
    gen_phi3(MPIP, IO, Bas);
  } else if (expab_order == 4) {
    gen_phi4(MPIP, IO, Bas);
  } else if (expab_order == 5) {
    gen_phi5(MPIP, IO, Bas);
  } else if (expab_order == 6) {
  }
  //DEBUG
  printf("clexpab::phi:\n"); expab_phi.print();
  //DEBUG

  timeP.resize(expab_order);
  tWfn.gen(MPIP, IO, Bas);
  WfnB.resize(expab_order);
  WfnP.resize(expab_order);
  hWfnP.resize(expab_order);
  for (int is = 0; is < expab_order; is ++) {
    WfnB[is].gen(MPIP, IO, Bas);
    WfnP[is].gen(MPIP, IO, Bas);
    hWfnP[is].gen(MPIP, IO, Bas);
  }
}
////////////////////////////////////////////////////////////////////////
void clexpab::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  std::cout << "clexpab::prop (1) nyi." << std::endl;
  abort();
}
////////////////////////////////////////////////////////////////////////
void clexpab::prop(const clmpi& Proc, const clbas& Bas,
		   const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  if (expab_order == 1) {
    prop1(Proc, Bas, Field, HPW, Wfn);
    //prop_general(Proc, Bas, Field, HPW, Wfn);
  } else if (expab_order == 2) {
    //    std::cout << "calling prop2_0..." << std::endl;
    //    prop2_0(Proc, Bas, Field, HPW, Wfn);
    //    std::cout << "calling prop2_1..." << std::endl;
    //    prop2_1(Proc, Bas, Field, HPW, Wfn);
    //    std::cout << "calling prop2..." << std::endl;
    //    prop2(Proc, Bas, Field, HPW, Wfn);
    prop_general(Proc, Bas, Field, HPW, Wfn);
  } else {
    prop_general(Proc, Bas, Field, HPW, Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void clexpab::prop1(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  double lfield[9];
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // stage 0: -i*h*u0 + N[u0]
  Field.get_value(lfield);
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn, tWfn);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[0]);
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfnP[0], tWfn);
  expab_phi.prod(Proc, Bas, Field.time, Field, HPW, tWfn);
  HPW.xpyz(Proc, Bas, tWfn, Wfn, Wfn);
}
////////////////////////////////////////////////////////////////////////
void clexpab::prop2(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // adjust previous steps
  double lfield[9];
  if (num_step == 0) {
    timeP[0] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[0]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[0]);
    timeP[1] = timeP[0];
    HPW.copy(Proc, Bas, WfnP[0], WfnP[1]);
    HPW.copy(Proc, Bas, hWfnP[0], hWfnP[1]);
  } else {
    timeP[1] = timeP[0];
    HPW.copy(Proc, Bas, WfnP[0], WfnP[1]);
    HPW.copy(Proc, Bas, hWfnP[0], hWfnP[1]);
    timeP[0] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[0]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[0]);
  }
  num_step ++;

  // WfnB(t) = -i{vext(t)-vext(t_n)}*WfnP(t) + hWfnP(t)
  double lfieldP[9];
  HPW.copy(Proc, Bas, hWfnP[0], WfnB[0]);
  HPW.copy(Proc, Bas, hWfnP[1], WfnB[1]);

  Field.get_value(timeP[1], lfieldP);
  for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
  HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[1], WfnB[1]);

  // tWfn = -i*h*u0 + N[u0]
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn, tWfn);
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfnP[0], tWfn);

  dcomplex rjfact;
  rjfact = expab_dt;
  HPW.axpy(Proc, Bas, rjfact, tWfn, Wfn);

  rjfact = tunit*expab_dt;
  const int K1 = expab_order-1;
  HPW.copy(Proc, Bas, tWfn, WfnP[K1]);
  HPW.h1tot(Proc, Bas, ONE, lfield, WfnP[K1], tWfn);
  HPW.scalo(Proc, Bas, rjfact, tWfn);

  HPW.axpy(Proc, Bas, +RUNIT, WfnB[0], tWfn);
  HPW.axpy(Proc, Bas, -RUNIT, WfnB[1], tWfn);

  // Wfn += phi_k(dt*\alpha)*w_k
  rjfact = expab_dt/2.0;
  HPW.scalc(Proc, Bas, rjfact, tWfn);
  expab_phi.prod(Proc, Bas, Field.time, Field, HPW, tWfn);
  HPW.xpy(Proc, Bas, tWfn, Wfn);

  if (expab_pc) {
    HPW.copy(Proc, Bas, WfnB[0], WfnB[1]);

    double ttime = timeP[0]+expab_dt;
    Field.get_value(ttime, lfieldP);
    HPW.htot(Proc, Bas, ONE, lfieldP, Wfn, WfnB[0]);
    for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
    HPW.v1ext(Proc, Bas, tunit, lfieldP, Wfn, WfnB[0]);

    HPW.copy(Proc, Bas, WfnP[0], Wfn);
    HPW.copy(Proc, Bas, WfnP[K1], tWfn);
    rjfact = expab_dt;
    HPW.axpy(Proc, Bas, rjfact, tWfn, Wfn);

    rjfact = tunit*expab_dt;
    HPW.copy(Proc, Bas, tWfn, WfnP[K1]);
    HPW.h1tot(Proc, Bas, ONE, lfield, WfnP[K1], tWfn);
    HPW.scalo(Proc, Bas, rjfact, tWfn);

    HPW.axpy(Proc, Bas, +RUNIT, WfnB[0], tWfn);
    HPW.axpy(Proc, Bas, -RUNIT, WfnB[1], tWfn);
    
    rjfact = expab_dt/2.0;
    HPW.scalc(Proc, Bas, rjfact, tWfn);
    expab_phi.prod(Proc, Bas, Field.time, Field, HPW, tWfn);
    HPW.xpy(Proc, Bas, tWfn, Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void clexpab::prop2_0(const clmpi& Proc, const clbas& Bas,
		      const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // adjust previous steps
  double lfield[9];
  if (num_step == 0) {
    timeP[0] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[0]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[0]);
    timeP[1] = timeP[0];
    HPW.copy(Proc, Bas, WfnP[0], WfnP[1]);
    HPW.copy(Proc, Bas, hWfnP[0], hWfnP[1]);
  } else {
    timeP[1] = timeP[0];
    HPW.copy(Proc, Bas, WfnP[0], WfnP[1]);
    HPW.copy(Proc, Bas, hWfnP[0], hWfnP[1]);
    timeP[0] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[0]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[0]);
  }
  num_step ++;

  // WfnB(t) = -i{vext(t)-vext(t_n)}*WfnP(t) + hWfnP(t)
  double lfieldP[9];
  HPW.copy(Proc, Bas, hWfnP[0], WfnB[0]);
  HPW.copy(Proc, Bas, hWfnP[1], WfnB[1]);
  Field.get_value(timeP[1], lfieldP);
  for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
  HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[1], WfnB[1]);

  // Wfn = phi_0(dt*\alpha)*Wfn
  expab_phi0.prod(Proc, Bas, Field.time, Field, HPW, Wfn);

  // Wfn += dt*phi_1(dt*\alpha)*B0
  dcomplex rjfact;
  rjfact = expab_dt;
  HPW.copy(Proc, Bas, WfnB[0], tWfn);
  HPW.scalc(Proc, Bas, rjfact, tWfn);
  expab_phi1.prod(Proc, Bas, Field.time, Field, HPW, tWfn);
  HPW.xpy(Proc, Bas, tWfn, Wfn);

  // Wfn += dt*phi_2(dt*\alpha)*B1
  rjfact = expab_dt/2.0;
  HPW.copy(Proc, Bas, WfnB[0], tWfn);
  HPW.axpy(Proc, Bas, -RUNIT, WfnB[1], tWfn);
  expab_phi.prod(Proc, Bas, Field.time, Field, HPW, tWfn);
  HPW.xpy(Proc, Bas, tWfn, Wfn);

  //std::cout << "prop2_0 (5)..." << std::endl;
}
////////////////////////////////////////////////////////////////////////
void clexpab::prop2_1(const clmpi& Proc, const clbas& Bas,
		      const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // adjust previous steps
  double lfield[9];
  if (num_step == 0) {
    timeP[0] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[0]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[0]);
    timeP[1] = timeP[0];
    HPW.copy(Proc, Bas, WfnP[0], WfnP[1]);
    HPW.copy(Proc, Bas, hWfnP[0], hWfnP[1]);
  } else {
    timeP[1] = timeP[0];
    HPW.copy(Proc, Bas, WfnP[0], WfnP[1]);
    HPW.copy(Proc, Bas, hWfnP[0], hWfnP[1]);
    timeP[0] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[0]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[0]);
  }
  num_step ++;

  // WfnB(t) = -i{vext(t)-vext(t_n)}*WfnP(t) + hWfnP(t)
  double lfieldP[9];
  HPW.copy(Proc, Bas, hWfnP[0], WfnB[0]);
  HPW.copy(Proc, Bas, hWfnP[1], WfnB[1]);
  Field.get_value(timeP[1], lfieldP);
  for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
  HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[1], WfnB[1]);

  // tWfn = -i*h*u0 + N[u0]
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn, tWfn);
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfnP[0], tWfn);

  // Wfn += dt*phi_1(dt*\alpha)*B0
  dcomplex rjfact;
  rjfact = expab_dt;
  HPW.scalc(Proc, Bas, rjfact, tWfn);
  expab_phi1.prod(Proc, Bas, Field.time, Field, HPW, tWfn);
  HPW.xpy(Proc, Bas, tWfn, Wfn);

  // Wfn += dt*phi_2(dt*\alpha)*B1
  rjfact = expab_dt/2.0;
  HPW.copy(Proc, Bas, WfnB[0], tWfn);
  HPW.axpy(Proc, Bas, -RUNIT, WfnB[1], tWfn);
  expab_phi.prod(Proc, Bas, Field.time, Field, HPW, tWfn);
  HPW.xpy(Proc, Bas, tWfn, Wfn);

  //std::cout << "prop2_0 (5)..." << std::endl;
}
////////////////////////////////////////////////////////////////////////
void clexpab::prop_general(const clmpi& Proc, const clbas& Bas,
			   const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // adjust previous steps
  double lfield[9];
  if (num_step == 0) {
    timeP[0] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[0]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[0]);
    for (int is = 1; is < expab_order; is ++) {
      //timeP[is] = timeP[is-1] - Field.dtime;
      timeP[is] = timeP[0];
      HPW.copy(Proc, Bas, WfnP[0], WfnP[is]);
      HPW.copy(Proc, Bas, hWfnP[0], hWfnP[is]);
    }
  } else {
    for (int is = expab_order-1; is > 0; is --) {
      timeP[is] = timeP[is-1];
      HPW.copy(Proc, Bas, WfnP[is-1], WfnP[is]);
      HPW.copy(Proc, Bas, hWfnP[is-1], hWfnP[is]);
    }
    timeP[0] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[0]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[0]);
  }
  num_step ++;

  // WfnB(t) = -i{vext(t)-vext(t_n)}*WfnP(t) + hWfnP(t)
  double lfieldP[9];
  for (int is = 0; is < expab_order; is ++) {
    HPW.copy(Proc, Bas, hWfnP[is], WfnB[is]);  
  }
  for (int is = 1; is < expab_order; is ++) {
    Field.get_value(timeP[is], lfieldP);
    for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
    HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[is], WfnB[is]);
  }

  // determine multistep coefficients
  int jfact;
  std::vector<double> amat(expab_order*expab_order);
  // A(i,0) = (-i)^0/0! = 1
  for (int is = 0; is < expab_order; is ++) {
    amat[is] = 1.0;
  }
  // A(0,j) = 0^j/j! = 0 (j > 0)
  for (int js = 1; js < expab_order; js ++) {
    amat[expab_order*js] = 0.0;
  }
  // A(i,j) = (-i)^j/j!
  jfact = 1;
  for (int js = 1; js < expab_order; js ++) {
    jfact *= js;
    for (int is = 1; is < expab_order; is ++) {
      amat[expab_order*js+is] = pow(-is,js)/jfact;
    }
  }

  //DEBUG
  //std::cout << "clexpab: amat" << std::endl;
  //for (int is = 0; is < expab_order; is ++) {
  //  for (int js = 0; js < expab_order; js ++) {
  //    printf("%15.5f", amat[expab_order*js+is]);
  //  }
  //  printf("\n");
  //}
  //DEBUG
  util_gmatinv_real_(&expab_order, &ZERO, &amat[0], &amat[0]);
  ////DEBUG
  //std::cout << "clexpab: inv_amat" << std::endl;
  //for (int is = 0; is < expab_order; is ++) {
  //  for (int js = 0; js < expab_order; js ++) {
  //    printf("%15.5f", amat[expab_order*js+is]);
  //  }
  //  printf("\n");
  //}
  //abort();
  //DEBUG

  dcomplex rjfact;
  const int K1 = expab_order-1;

  // tWfn = -i*h*u0 + N[u0]
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn, tWfn);
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfnP[0], tWfn);
  HPW.copy(Proc, Bas, tWfn, hWfnP[K1]);

  jfact = 1;
  for (int js = 1; js < expab_order; js ++) {
    // Wfn += dt/j! * w_j
    jfact *= js;
    rjfact = expab_dt/jfact;
    HPW.axpy(Proc, Bas, rjfact, tWfn, Wfn);

    // w_{j+1} = dt*\alpha*w_{j} + sum_i A^{-1}_{ji} beta_{i}
    // WfnP[k-1] is used as a scratch
    HPW.copy(Proc, Bas, tWfn, WfnP[K1]);
    HPW.h1tot(Proc, Bas, ONE, lfield, WfnP[K1], tWfn);
    HPW.scalo(Proc, Bas, tunit*expab_dt, tWfn);
    for (int is = 0; is < expab_order; is ++) {
      dcomplex afac = amat[expab_order*is+js];
      HPW.axpy(Proc, Bas, afac, WfnB[is], tWfn);
    }
  }

  // Wfn += phi_k(dt*\alpha)*w_k
  jfact *= expab_order;
  rjfact = expab_dt/jfact;
  HPW.scalc(Proc, Bas, rjfact, tWfn);
  expab_phi.prod(Proc, Bas, Field.time, Field, HPW, tWfn);
  HPW.xpy(Proc, Bas, tWfn, Wfn);

  // At this moment
  // WfnP[K1] holds f(t_n,y_n)
  if (expab_pc) {
    prop_general_corrector(Proc, Bas, Field, HPW, Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void clexpab::prop_general_corrector(const clmpi& Proc, const clbas& Bas,
				     const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  ////////////
  int jfact;
  std::vector<double> amat(expab_order*expab_order);
  // A(i,0) = (1-i)^0/0! = 1
  for (int is = 0; is < expab_order; is ++) {
    amat[is] = 1.0;
  }
  // A(1,j) = 0^j/j! = 0 (j > 0)
  for (int js = 1; js < expab_order; js ++) {
    amat[expab_order*js+1] = 0.0;
  }
  // A(i,j) = (-i)^j/j!
  jfact = 1;
  for (int js = 1; js < expab_order; js ++) {
    jfact *= js;
    for (int is = 0; is < expab_order; is ++) {
      if (is != 1) {
	amat[expab_order*js+is] = pow(1-is,js)/jfact;
      }
    }
  }
  //DEBUG
  //std::cout << "clexpab: amat" << std::endl;
  //for (int is = 0; is < expab_order; is ++) {
  //  for (int js = 0; js < expab_order; js ++) {
  //    printf("%15.5f", amat[expab_order*js+is]);
  //  }
  //  printf("\n");
  //}
  //DEBUG
  util_gmatinv_real_(&expab_order, &ZERO, &amat[0], &amat[0]);
  ////DEBUG
  //std::cout << "clexpab: inv_amat" << std::endl;
  //for (int is = 0; is < expab_order; is ++) {
  //  for (int js = 0; js < expab_order; js ++) {
  //    printf("%15.5f", amat[expab_order*js+is]);
  //  }
  //  printf("\n");
  //}
  //abort();
  //DEBUG
  ////////////

  for (int is = expab_order-1; is > 0; is --) {
    HPW.copy(Proc, Bas, WfnB[is-1], WfnB[is]);
  }
  double lfield[9], lfield1[9];
  double ttime = timeP[0]+expab_dt;
  Field.get_value(timeP[0], lfield);
  Field.get_value(ttime, lfield1);
  HPW.htot(Proc, Bas, ONE, lfield1, Wfn, WfnB[0]);
  for (int ii = 0; ii < 9; ii++) lfield1[ii] -= lfield[ii];
  HPW.v1ext(Proc, Bas, tunit, lfield1, Wfn, WfnB[0]);

  dcomplex rjfact;
  const int K1 = expab_order-1;
  HPW.copy(Proc, Bas, WfnP[0], Wfn);
  HPW.copy(Proc, Bas, hWfnP[K1], tWfn);

  jfact = 1;
  for (int js = 1; js < expab_order; js ++) {
    // Wfn += dt/j! * w_j
    jfact *= js;
    rjfact = expab_dt/jfact;
    HPW.axpy(Proc, Bas, rjfact, tWfn, Wfn);

    // w_{j+1} = dt*\alpha*w_{j} + sum_i A^{-1}_{ji} beta_{i}
    // WfnP[k-1] is used as a scratch
    HPW.copy(Proc, Bas, tWfn, WfnP[K1]);
    HPW.h1tot(Proc, Bas, ONE, lfield, WfnP[K1], tWfn);
    HPW.scalo(Proc, Bas, tunit*expab_dt, tWfn);
    for (int is = 0; is < expab_order; is ++) {
      //Old      dcomplex gfac = gcoeff[js][is];
      dcomplex gfac = amat[expab_order*is+js];
      HPW.axpy(Proc, Bas, gfac, WfnB[is], tWfn);
    }
  }

  // Wfn += phi_k(dt*\alpha)*w_k
  jfact *= expab_order;
  rjfact = expab_dt/jfact;
  HPW.scalc(Proc, Bas, rjfact, tWfn);
  expab_phi.prod(Proc, Bas, Field.time, Field, HPW, tWfn);
  HPW.xpy(Proc, Bas, tWfn, Wfn);
}
//////////////////////////////////////////////////////////////////////////
void clexpab::gen_phi1(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi1_n(max_dim_numer + 1);
  std::vector<dcomplex> phi1_d0(max_dim_denom);
  std::vector<dcomplex> phi1_d1(max_dim_denom);

  if (expab_ndpade == 11) {
    //                2  
    // phi1(x) = [- ------]
    //              x - 2 
    // x0 = 2
    dfac = -1.0;
    phi1_n[0] = 2.0;
    phi1_n[1] = 0.0;
    roots[0] = 2.0;
  } else if (expab_ndpade == 12) {
    //
    //              - x + 6
    // phi1(x) = [------------] 
    //             2              
    //            x  - 4 x + 6
    // x0 = 2.0000000000000000000 - 1.4142135623730950488i
    // x0 = 2.0000000000000000000 + 1.4142135623730950488i
    dfac = 1.0;
    phi1_n[0] = 6.0;
    phi1_n[1] = -1.0;
    roots[0] = 2.0000000000000000000 - 1.4142135623730950488*IUNIT;
    roots[1] = 2.0000000000000000000 + 1.4142135623730950488*IUNIT;
  } else if (expab_ndpade == 22) {
    //
    //                  12
    // phi1(x) = [-------------] 
    //             2              
    //            x  - 6 x + 12
    // x0 = 3.0000000000000000000 - 1.7320508075688772935i
    // x0 = 3.0000000000000000000 + 1.7320508075688772935i
    dfac = 1.0;
    phi1_n[0] = 12.0;
    phi1_n[1] = 0.0;
    phi1_n[2] = 0.0;
    roots[0] = 3.0000000000000000000-1.7320508075688772935*IUNIT;
    roots[1] = 3.0000000000000000000+1.7320508075688772935*IUNIT;
  } else if (expab_ndpade == 13) {
    // Pade(1,3) approximation of exp(x)
    //                  2
    //               - x  + 6 x - 24
    // phi1(x) = [---------------------]
    //               3      2               
    //            x  - 6 x  + 18 x - 24
    // x0 = 1.687091590520766641994056 - 2.508731754924880510838744i
    // x0 = 1.687091590520766641994056 + 2.508731754924880510838744i
    // x0 = 2.625816818958466716011889
    dfac = 1.0;
    phi1_n[0] = -24.0;
    phi1_n[1] = 6.0;
    phi1_n[2] = -1.0;
    roots[0] = 1.687091590520766641994056 - 2.508731754924880510838744*IUNIT;
    roots[1] = 1.687091590520766641994056 + 2.508731754924880510838744*IUNIT;
    roots[2] = 2.625816818958466716011889;
  } else if (expab_ndpade == 23) {
    // Pade(2,3) of exp(x)
    //                    2
    //                 - x  + 6 x - 60
    // phi1(x) = [-----------------------]
    //             3      2               
    //            x  - 9 x  + 36 x - 60
    // x0 = 2.6810828736277521339 + 3.0504301992474105694i
    // x0 = 2.6810828736277521339 - 3.0504301992474105694i
    // x0 = 3.6378342527444957322
    dfac = 1.0;
    phi1_n[0] = -60.0;
    phi1_n[1] = 6.0;
    phi1_n[2] = -1.0;
    roots[0] = 2.6810828736277521339 + 3.0504301992474105694*IUNIT;
    roots[1] = 2.6810828736277521339 - 3.0504301992474105694*IUNIT;
    roots[2] = 3.6378342527444957322;
  } else if (expab_ndpade == 33) {
    //                       2
    //                    2 x  + 120   
    // phi1(x) = [- -----------------------]
    //               3       2               
    //              x  - 12 x  + 60 x - 120 
    // x0 = 3.6778146453739144071 + 3.5087619195674433219i
    // x0 = 3.6778146453739144071 - 3.5087619195674433219i
    // x0 = 4.6443707092521711858
    dfac = -1.0;
    phi1_n[0] = 120.0;
    phi1_n[1] = 0.0;
    phi1_n[2] = 2.0;
    phi1_n[3] = 0.0;
    roots[0] = 3.6778146453739144071 + 3.5087619195674433219*IUNIT;
    roots[1] = 3.6778146453739144071 - 3.5087619195674433219*IUNIT;
    roots[2] = 4.6443707092521711858;
  } else if (expab_ndpade == 34) {
    //                    3       2
    //                 - x  + 20 x  - 60 x + 840
    // phi1(x) = [ ----------------------------------]
    //              4       3        2       
    //             x  - 16 x  + 120 x  - 480 x + 840
    // x0 = 3.2128068968715339829 - 4.7730874332766424998i
    // x0 = 3.2128068968715339829 + 4.7730874332766424998i
    // x0 = 4.7871931031284660171 - 1.5674764168952081241i
    // x0 = 4.7871931031284660171 + 1.5674764168952081241i
    dfac = 1.0;
    phi1_n[0] = 840.0;
    phi1_n[1] = -60.0;
    phi1_n[2] = 20.0;
    phi1_n[3] = -1.0;
    roots[0] = 3.2128068968715339829 - 4.7730874332766424998*IUNIT;
    roots[1] = 3.2128068968715339829 + 4.7730874332766424998*IUNIT;
    roots[2] = 4.7871931031284660171 - 1.5674764168952081241*IUNIT;
    roots[3] = 4.7871931031284660171 + 1.5674764168952081241*IUNIT;
  } else if (expab_ndpade == 44) {
    //                          2
    //                      40 x  + 1680
    // phi1(x) = [ ----------------------------------]
    //              4       3        2       
    //             x  - 20 x  + 180 x  - 840 x + 1680
    // x0 = 
    // x0 = 
    // x0 = 
    // x0 = 
    dfac = 1.0;
    phi1_n[0] = 1680.0;
    phi1_n[1] = 0.0;
    phi1_n[2] = 40.0;
    phi1_n[3] = 0.0;
    phi1_n[4] = 0.0;
    roots[0] = 4.2075787943592556632 - 5.3148360837135054337*IUNIT;
    roots[1] = 4.2075787943592556632 + 5.3148360837135054337*IUNIT;
    roots[2] = 5.7924212056407443368 - 1.7344682578690075036*IUNIT;
    roots[3] = 5.7924212056407443368 + 1.7344682578690075036*IUNIT;
  } else {
    std::cout << "clexpab:gen_exp. bad expab_npade or expab_dpade." << std::endl;
    abort();
  }

  phi1_d1[0] = +dfac;
  phi1_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < expab_dpade; idim ++) {
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
  tfac = tunit * expab_dt; // BE CAREFUL!
  ttmp = RUNIT * expab_dt; // BE CAREFUL!
  for (int inumer = 0; inumer <= expab_npade; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < expab_dpade; idenom ++) {
    cd1[idenom] *= tfac;
  }
  expab_phi.gen(MPIP, IO, Bas, expab_dt, expab_h1type, expab_npade, expab_dpade, cn, cd0, cd1);

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
//////////////////////////////////////////////////////////////////////////
void clexpab::gen_phi2(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi2_n(max_dim_numer + 1);
  std::vector<dcomplex> phi2_d0(max_dim_denom);
  std::vector<dcomplex> phi2_d1(max_dim_denom);

  if (expab_ndpade == 11) {
    //                1
    // phi2(x) = [- ------]
    //              x - 2 
    // x0 = 2
    dfac = -1.0;
    phi2_n[0] = 1.0;
    phi2_n[1] = 0.0;
    roots[0] = 2.0;
  } else if (expab_ndpade == 12) {
    //
    //              - x + 3
    // phi2(x) = [------------] 
    //             2              
    //            x  - 4 x + 6
    // x0 = 2.0000000000000000000 - 1.4142135623730950488i
    // x0 = 2.0000000000000000000 + 1.4142135623730950488i
    dfac = 1.0;
    phi2_n[0] = 3.0;
    phi2_n[1] = -1.0;
    roots[0] = 2.0000000000000000000 - 1.4142135623730950488*IUNIT;
    roots[1] = 2.0000000000000000000 + 1.4142135623730950488*IUNIT;
  } else if (expab_ndpade == 22) {
    //
    //                6 - x
    // phi2(x) = [-------------] 
    //             2              
    //            x  - 6 x + 12
    // x0 = 3.0000000000000000000 - 1.7320508075688772935i
    // x0 = 3.0000000000000000000 + 1.7320508075688772935i
    dfac = 1.0;
    phi2_n[0] = 6.0;
    phi2_n[1] = -1.0;
    phi2_n[2] = 0.0;
    roots[0] = 3.0-sqrt(3.0)*IUNIT;
    roots[1] = 3.0+sqrt(3.0)*IUNIT;
    //    roots[0] = 3.0000000000000000000-1.7320508075688772935*IUNIT;
    //    roots[1] = 3.0000000000000000000+1.7320508075688772935*IUNIT;
  } else if (expab_ndpade == 13) {
    // Pade(1,3) approximation of exp(x)
    //                  2
    //               - x  + 5 x - 12
    // phi2(x) = [---------------------]
    //               3      2               
    //            x  - 6 x  + 18 x - 24
    // x0 = 1.687091590520766641994056 - 2.508731754924880510838744i
    // x0 = 1.687091590520766641994056 + 2.508731754924880510838744i
    // x0 = 2.625816818958466716011889
    dfac = 1.0;
    phi2_n[0] = -12.0;
    phi2_n[1] = 5.0;
    phi2_n[2] = -1.0;
    roots[0] = 2.625816818958466716011889;
    roots[1] = 1.687091590520766641994056 - 2.508731754924880510838744*IUNIT;
    roots[2] = 1.687091590520766641994056 + 2.508731754924880510838744*IUNIT;
  } else if (expab_ndpade == 23) {
    // Pade(2,3) of exp(x)
    //                    2
    //                 - x  + 8 x - 30
    // phi2(x) = [-----------------------]
    //             3      2               
    //            x  - 9 x  + 36 x - 60
    // x0 = 2.6810828736277521339 - 3.0504301992474105694i
    // x0 = 2.6810828736277521339 + 3.0504301992474105694i
    // x0 = 3.6378342527444957322
    dfac = 1.0;
    phi2_n[0] = -30.0;
    phi2_n[1] = 8.0;
    phi2_n[2] = -1.0;
//    roots[0] = 3.6378342527444957322;
//    roots[1] = 2.6810828736277521339 - 3.0504301992474105694*IUNIT;
//    roots[2] = 2.6810828736277521339 + 3.0504301992474105694*IUNIT;
    roots[0] = 3.0 - pow(3.0,1.0/3.0) + pow(3.0,2.0/3.0);
    roots[1] = 3.0
      + 0.5*pow(3.0,1.0/3.0)*(1.0-pow(3.0,1.0/2.0)*IUNIT)
      - 0.5*pow(3.0,2.0/3.0)*(1.0+pow(3.0,1.0/2.0)*IUNIT);
    roots[2] = 3.0
      - 0.5*pow(3.0,2.0/3.0)*(1.0-pow(3.0,1.0/2.0)*IUNIT)
      + 0.5*pow(3.0,1.0/3.0)*(1.0+pow(3.0,1.0/2.0)*IUNIT);
  } else if (expab_ndpade == 33) {
    // Pade(3,3) of exp(x)
    //                    2
    //                 - x  + 10 x - 60   
    // phi2(x) = [-----------------------]
    //             3       2               
    //            x  - 12 x  + 60 x - 120 
    // x0 = 3.6778146453739144071 + 3.5087619195674433219i
    // x0 = 3.6778146453739144071 - 3.5087619195674433219i
    // x0 = 4.6443707092521711858
    dfac = 1.0;
    phi2_n[0] = -60.0;
    phi2_n[1] = 10.0;
    phi2_n[2] = -1.0;
    phi2_n[3] = 0.0;
    roots[0] = 3.6778146453739144071 + 3.5087619195674433219*IUNIT;
    roots[1] = 3.6778146453739144071 - 3.5087619195674433219*IUNIT;
    roots[2] = 4.6443707092521711858;
  } else if (expab_ndpade == 34) {
    //                    3       2
    //                 - x  + 15 x  - 100 x + 420
    // phi2(x) = [ ----------------------------------]
    //              4       3        2       
    //             x  - 16 x  + 120 x  - 480 x + 840
    // x0 = 3.2128068968715339829 - 4.7730874332766424998i
    // x0 = 3.2128068968715339829 + 4.7730874332766424998i
    // x0 = 4.7871931031284660171 - 1.5674764168952081241i
    // x0 = 4.7871931031284660171 + 1.5674764168952081241i
    dfac = 1.0;
    phi2_n[0] = 420.0;
    phi2_n[1] = -100.0;
    phi2_n[2] = 15.0;
    phi2_n[3] = -1.0;
    roots[0] = 3.2128068968715339829 - 4.7730874332766424998*IUNIT;
    roots[1] = 3.2128068968715339829 + 4.7730874332766424998*IUNIT;
    roots[2] = 4.7871931031284660171 - 1.5674764168952081241*IUNIT;
    roots[3] = 4.7871931031284660171 + 1.5674764168952081241*IUNIT;
  } else if (expab_ndpade == 44) {
    //                    3       2
    //                 - x  + 20 x  - 140 x + 840
    // phi2(x) = [ ----------------------------------]
    //              4       3        2       
    //             x  - 20 x  + 180 x  - 840 x + 1680
    // x0 =  4.2075787943592556632 - 5.3148360837135054337i
    // x0 =  4.2075787943592556632 + 5.3148360837135054337i
    // x0 =  5.7924212056407443368 - 1.7344682578690075036i
    // x0 =  5.7924212056407443368 + 1.7344682578690075036i
    dfac = 1.0;
    phi2_n[0] = 840.0;
    phi2_n[1] = -140.0;
    phi2_n[2] = 20.0;
    phi2_n[3] = -1.0;
    phi2_n[4] = 0.0;
    roots[0] = 4.2075787943592556632 - 5.3148360837135054337*IUNIT;
    roots[1] = 4.2075787943592556632 + 5.3148360837135054337*IUNIT;
    roots[2] = 5.7924212056407443368 - 1.7344682578690075036*IUNIT;
    roots[3] = 5.7924212056407443368 + 1.7344682578690075036*IUNIT;
  } else {
    std::cout << "clexpab:gen_exp. bad expab_npade or expab_dpade." << std::endl;
    abort();
  }

  phi2_d1[0] = +dfac;
  phi2_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < expab_dpade; idim ++) {
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
  tfac = tunit * expab_dt; // BE CAREFUL!
  ttmp = RUNIT * expab_dt; // BE CAREFUL!
  for (int inumer = 0; inumer <= expab_npade; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < expab_dpade; idenom ++) {
    cd1[idenom] *= tfac;
  }
  expab_phi.gen(MPIP, IO, Bas, expab_dt, expab_h1type, expab_npade, expab_dpade, cn, cd0, cd1);
}
//////////////////////////////////////////////////////////////////////////
void clexpab::gen_phi3(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi3_n(max_dim_numer + 1);
  std::vector<dcomplex> phi3_d0(max_dim_denom);
  std::vector<dcomplex> phi3_d1(max_dim_denom);

  if (expab_ndpade == 11) {
    //                  1
    // phi3(x) = [- --------]
    //              2(x - 2)
    // x0 = 2
    dfac = -2.0;
    phi3_n[0] = 1.0;
    phi3_n[1] = 0.0;
    roots[0] = 2.0;
  } else if (expab_ndpade == 12) {
    //
    //                - x + 2
    // phi3(x) = [---------------] 
    //               2              
    //            2(x  - 4 x + 6)
    // x0 = 2.0000000000000000000 - 1.4142135623730950488i
    // x0 = 2.0000000000000000000 + 1.4142135623730950488i
    dfac = 2.0;
    phi3_n[0] = 2.0;
    phi3_n[1] = -1.0;
    roots[0] = 2.0000000000000000000 - 1.4142135623730950488*IUNIT;
    roots[1] = 2.0000000000000000000 + 1.4142135623730950488*IUNIT;
  } else if (expab_ndpade == 22) {
    //
    //                  4 - x
    // phi3(x) = [----------------] 
    //               2              
    //            2(x  - 6 x + 12)
    // x0 = 3.0000000000000000000 - 1.7320508075688772935i
    // x0 = 3.0000000000000000000 + 1.7320508075688772935i
    dfac = 2.0;
    phi3_n[0] = 4.0;
    phi3_n[1] = -1.0;
    phi3_n[2] = 0.0;
    roots[0] = 3.0000000000000000000-1.7320508075688772935*IUNIT;
    roots[1] = 3.0000000000000000000+1.7320508075688772935*IUNIT;
  } else if (expab_ndpade == 13) {
    // Pade(1,3) approximation of exp(x)
    //                    2
    //                 - x  + 4 x - 8
    // phi3(x) = [------------------------]
    //               3      2               
    //            2(x  - 6 x  + 18 x - 24)
    // x0 = 1.687091590520766641994056 - 2.508731754924880510838744i
    // x0 = 1.687091590520766641994056 + 2.508731754924880510838744i
    // x0 = 2.625816818958466716011889
    dfac = 2.0;
    phi3_n[0] = -8.0;
    phi3_n[1] = 4.0;
    phi3_n[2] = -1.0;
    roots[0] = 1.687091590520766641994056 - 2.508731754924880510838744*IUNIT;
    roots[1] = 1.687091590520766641994056 + 2.508731754924880510838744*IUNIT;
    roots[2] = 2.625816818958466716011889;
  } else if (expab_ndpade == 23) {
    // Bad: Pade(2,3) of exp(x)
    //                    2
    //                 - x  + 7 x - 20
    // phi3(x) = [------------------------]
    //               3      2               
    //            2(x  - 9 x  + 36 x - 60)
    // x0 = 2.6810828736277521339 + 3.0504301992474105694i
    // x0 = 2.6810828736277521339 - 3.0504301992474105694i
    // x0 = 3.6378342527444957322
    dfac = 2.0;
    phi3_n[0] = -20.0;
    phi3_n[1] = 7.0;
    phi3_n[2] = -1.0;
    roots[0] = 2.6810828736277521339 + 3.0504301992474105694*IUNIT;
    roots[1] = 2.6810828736277521339 - 3.0504301992474105694*IUNIT;
    roots[2] = 3.6378342527444957322;
  } else if (expab_ndpade == 33) {
    //                    2
    //                 - x  + 10 x - 40
    // phi3(x) = [--------------------------]
    //               3       2               
    //            2(x  - 12 x  + 60 x - 120)
    // x0 = 3.6778146453739144071 + 3.5087619195674433219i
    // x0 = 3.6778146453739144071 - 3.5087619195674433219i
    // x0 = 4.6443707092521711858
    dfac = 2.0;
    phi3_n[0] = -40.0;
    phi3_n[1] = 10.0;
    phi3_n[2] = -1.0;
    phi3_n[3] = 0.0;
    roots[0] = 3.6778146453739144071 + 3.5087619195674433219*IUNIT;
    roots[1] = 3.6778146453739144071 - 3.5087619195674433219*IUNIT;
    roots[2] = 4.6443707092521711858;
  } else if (expab_ndpade == 34) {
    //                    3       2
    //                 - x  + 14 x  - 90 x + 280
    // phi3(x) = [ ------------------------------------]
    //                4       3        2       
    //             2(x  - 16 x  + 120 x  - 480 x + 840)
    // x0 = 3.2128068968715339829 - 4.7730874332766424998i
    // x0 = 3.2128068968715339829 + 4.7730874332766424998i
    // x0 = 4.7871931031284660171 - 1.5674764168952081241i
    // x0 = 4.7871931031284660171 + 1.5674764168952081241i
    dfac = 2.0;
    phi3_n[0] = 280.0;
    phi3_n[1] = -90.0;
    phi3_n[2] = 14.0;
    phi3_n[3] = -1.0;
    roots[0] = 3.2128068968715339829 - 4.7730874332766424998*IUNIT;
    roots[1] = 3.2128068968715339829 + 4.7730874332766424998*IUNIT;
    roots[2] = 4.7871931031284660171 - 1.5674764168952081241*IUNIT;
    roots[3] = 4.7871931031284660171 + 1.5674764168952081241*IUNIT;
  } else if (expab_ndpade == 44) {
    //                    3       2
    //                 - x  + 18 x  - 140 x + 560
    // phi3(x) = [ -------------------------------------]
    //                4       3        2       
    //             2(x  - 20 x  + 180 x  - 840 x + 1680)
    // x0 =  4.2075787943592556632 - 5.3148360837135054337i
    // x0 =  4.2075787943592556632 + 5.3148360837135054337i
    // x0 =  5.7924212056407443368 - 1.7344682578690075036i
    // x0 =  5.7924212056407443368 + 1.7344682578690075036i
    dfac = 2.0;
    phi3_n[0] = 560.0;
    phi3_n[1] = -140.0;
    phi3_n[2] = 18.0;
    phi3_n[3] = -1.0;
    phi3_n[4] = 0.0;
    roots[0] = 4.2075787943592556632 - 5.3148360837135054337*IUNIT;
    roots[1] = 4.2075787943592556632 + 5.3148360837135054337*IUNIT;
    roots[2] = 5.7924212056407443368 - 1.7344682578690075036*IUNIT;
    roots[3] = 5.7924212056407443368 + 1.7344682578690075036*IUNIT;
  } else {
    std::cout << "clexpab:gen_exp. bad expab_npade or expab_dpade." << std::endl;
    abort();
  }

  phi3_d1[0] = +dfac;
  phi3_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < expab_dpade; idim ++) {
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
  tfac = tunit * expab_dt; // BE CAREFUL!
  ttmp = RUNIT * expab_dt; // BE CAREFUL!
  for (int inumer = 0; inumer <= expab_npade; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < expab_dpade; idenom ++) {
    cd1[idenom] *= tfac;
  }
  expab_phi.gen(MPIP, IO, Bas, expab_dt, expab_h1type, expab_npade, expab_dpade, cn, cd0, cd1);
}
//////////////////////////////////////////////////////////////////////////
void clexpab::gen_phi4(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi4_n(max_dim_numer + 1);
  std::vector<dcomplex> phi4_d0(max_dim_denom);
  std::vector<dcomplex> phi4_d1(max_dim_denom);

  if (expab_ndpade == 12) {
    //
    //                - x + 1
    // phi4(x) = [---------------] 
    //               2              
    //            6(x  - 4 x + 6)
    // x0 = 2.0000000000000000000 - 1.4142135623730950488i
    // x0 = 2.0000000000000000000 + 1.4142135623730950488i
    dfac = 6.0;
    phi4_n[0] = 1.0;
    phi4_n[1] = -1.0;
    roots[0] = 2.0000000000000000000 - 1.4142135623730950488*IUNIT;
    roots[1] = 2.0000000000000000000 + 1.4142135623730950488*IUNIT;
  } else if (expab_ndpade == 22) {
    // Bad: Phi4 consistent with Pade(2,2) approximation of exp(x)
    //
    //                  3 - x
    // phi4(x) = [----------------] 
    //               2              
    //            6(x  - 6 x + 12)
    // x0 = 3.0000000000000000000 - 1.7320508075688772935i
    // x0 = 3.0000000000000000000 + 1.7320508075688772935i
    dfac = 6.0;
    phi4_n[0] = 3.0;
    phi4_n[1] = -1.0;
    phi4_n[2] = 0.0;
    roots[0] = 3.0000000000000000000-1.7320508075688772935*IUNIT;
    roots[1] = 3.0000000000000000000+1.7320508075688772935*IUNIT;
  } else if (expab_ndpade == 13) {
    // Good: Phi4 consistent with Pade(1,3) approximation of exp(x)
    //                    2
    //                 - x  + 3 x - 6
    // phi4(x) = [------------------------]
    //               3      2               
    //            6(x  - 6 x  + 18 x - 24)
    // x0 = 1.687091590520766641994056 - 2.508731754924880510838744i
    // x0 = 1.687091590520766641994056 + 2.508731754924880510838744i
    // x0 = 2.625816818958466716011889
    dfac = 6.0;
    phi4_n[0] = -6.0;
    phi4_n[1] = 3.0;
    phi4_n[2] = -1.0;
    roots[0] = 1.687091590520766641994056 - 2.508731754924880510838744*IUNIT;
    roots[1] = 1.687091590520766641994056 + 2.508731754924880510838744*IUNIT;
    roots[2] = 2.625816818958466716011889;
  } else if (expab_ndpade == 23) {
    // Bad: Phi4 consistent with Pade(2,3) approximation of exp(x)
    //                    2
    //                 - x  + 6 x - 15
    // phi4(x) = [------------------------]
    //               3      2               
    //            6(x  - 9 x  + 36 x - 60)
    // x0 = 2.6810828736277521339 + 3.0504301992474105694i
    // x0 = 2.6810828736277521339 - 3.0504301992474105694i
    // x0 = 3.6378342527444957322
    dfac = 6.0;
    phi4_n[0] = -15.0;
    phi4_n[1] = 6.0;
    phi4_n[2] = -1.0;
    roots[0] = 2.6810828736277521339 + 3.0504301992474105694*IUNIT;
    roots[1] = 2.6810828736277521339 - 3.0504301992474105694*IUNIT;
    roots[2] = 3.6378342527444957322;
  } else if (expab_ndpade == 33) {
    // Bad: Phi4 consistent with Pade(3,3) approximation of exp(x)
    //                    2
    //                 - x  + 9 x - 30
    // phi4(x) = [--------------------------]
    //               3       2               
    //            6(x  - 12 x  + 60 x - 120)
    // x0 = 3.6778146453739144071 + 3.5087619195674433219i
    // x0 = 3.6778146453739144071 - 3.5087619195674433219i
    // x0 = 4.6443707092521711858
    dfac = 6.0;
    phi4_n[0] = -30.0;
    phi4_n[1] = 9.0;
    phi4_n[2] = -1.0;
    phi4_n[3] = 0.0;
    roots[0] = 3.6778146453739144071 + 3.5087619195674433219*IUNIT;
    roots[1] = 3.6778146453739144071 - 3.5087619195674433219*IUNIT;
    roots[2] = 4.6443707092521711858;
  } else if (expab_ndpade == 14) {
    // Bad: Phi4 consistent with Pade(1,4) of exp(x)
    //                    3      2
    //                 - x  + 5 x  - 18 x + 30
    // phi4(x) = [ ------------------------------------]
    //                4      3       2       
    //             6(x  - 8 x  + 36 x  - 96 x + 120)
    // x0 = 1.2356535842848997066486550 - 3.4376524936710512909193657i
    // x0 = 1.2356535842848997066486550 + 3.4376524936710512909193657i
    // x0 = 2.7643464157151002933513450 - 1.1623236292832751740121108i
    // x0 = 2.7643464157151002933513450 + 1.1623236292832751740121108i
    dfac = 6.0;
    phi4_n[0] = 30.0;
    phi4_n[1] = -18.0;
    phi4_n[2] = 5.0;
    phi4_n[3] = -1.0;
    roots[0] = 1.2356535842848997066486550 - 3.4376524936710512909193657*IUNIT;
    roots[1] = 1.2356535842848997066486550 + 3.4376524936710512909193657*IUNIT;
    roots[2] = 2.7643464157151002933513450 - 1.1623236292832751740121108*IUNIT;
    roots[3] = 2.7643464157151002933513450 + 1.1623236292832751740121108*IUNIT;
  } else if (expab_ndpade == 24) {
    // Bad: Phi4 consistent with Pade(2,4) of exp(x)
    //                    3      2
    //                 - x  + 9 x  - 42 x + 90
    // phi4(x) = [ ------------------------------------]
    //                4       3       2       
    //             6(x  - 12 x  + 72 x  - 240 x + 360)
    // x0 = 2.2209800329898068974239251 - 4.1603914455069319822284852i
    // x0 = 2.2209800329898068974239251 + 4.1603914455069319822284852i
    // x0 = 3.7790199670101931025760749 - 1.3801765242728430462268849i
    // x0 = 3.7790199670101931025760749 + 1.3801765242728430462268849i
    dfac = 6.0;
    phi4_n[0] = 90.0;
    phi4_n[1] = -42.0;
    phi4_n[2] = 9.0;
    phi4_n[3] = -1.0;
    roots[0] = 2.2209800329898068974239251 - 4.1603914455069319822284852*IUNIT;
    roots[1] = 2.2209800329898068974239251 + 4.1603914455069319822284852*IUNIT;
    roots[2] = 3.7790199670101931025760749 - 1.3801765242728430462268849*IUNIT;
    roots[3] = 3.7790199670101931025760749 + 1.3801765242728430462268849*IUNIT;
  } else if (expab_ndpade == 34) {
    // Bad: Phi4 consistent with Pade(3,4) of exp(x)
    //                    3       2
    //                 - x  + 13 x  - 78 x + 210
    // phi4(x) = [ ------------------------------------]
    //                4       3        2       
    //             6(x  - 16 x  + 120 x  - 480 x + 840)
    // x0 = 3.2128068968715339829 - 4.7730874332766424998i
    // x0 = 3.2128068968715339829 + 4.7730874332766424998i
    // x0 = 4.7871931031284660171 - 1.5674764168952081241i
    // x0 = 4.7871931031284660171 + 1.5674764168952081241i
    dfac = 6.0;
    phi4_n[0] = 210.0;
    phi4_n[1] = -78.0;
    phi4_n[2] = 13.0;
    phi4_n[3] = -1.0;
    roots[0] = 3.2128068968715339829 - 4.7730874332766424998*IUNIT;
    roots[1] = 3.2128068968715339829 + 4.7730874332766424998*IUNIT;
    roots[2] = 4.7871931031284660171 - 1.5674764168952081241*IUNIT;
    roots[3] = 4.7871931031284660171 + 1.5674764168952081241*IUNIT;
  } else if (expab_ndpade == 44) {
    // Bad: Phi4 consistent with Pade(4,4) of exp(x)
    //                    3       2
    //                 - x  + 17 x  - 126 x + 420
    // phi4(x) = [ -------------------------------------]
    //                4       3        2       
    //             6(x  - 20 x  + 180 x  - 840 x + 1680)
    // x0 =  4.2075787943592556632 - 5.3148360837135054337i
    // x0 =  4.2075787943592556632 + 5.3148360837135054337i
    // x0 =  5.7924212056407443368 - 1.7344682578690075036i
    // x0 =  5.7924212056407443368 + 1.7344682578690075036i
    dfac = 6.0;
    phi4_n[0] = 420.0;
    phi4_n[1] = -126.0;
    phi4_n[2] = 17.0;
    phi4_n[3] = -1.0;
    phi4_n[4] = 0.0;
    roots[0] = 4.2075787943592556632 - 5.3148360837135054337*IUNIT;
    roots[1] = 4.2075787943592556632 + 5.3148360837135054337*IUNIT;
    roots[2] = 5.7924212056407443368 - 1.7344682578690075036*IUNIT;
    roots[3] = 5.7924212056407443368 + 1.7344682578690075036*IUNIT;
  } else {
    std::cout << "clexpab:gen_exp. bad expab_npade or expab_dpade." << std::endl;
    abort();
  }

  phi4_d1[0] = +dfac;
  phi4_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < expab_dpade; idim ++) {
    phi4_d1[idim] = +RUNIT;
    phi4_d0[idim] = -roots[idim];
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

  // ##### dt * phi4(-i*h*dt) #####
  cn = phi4_n;
  cd0 = phi4_d0;
  cd1 = phi4_d1;
  tfac = tunit * expab_dt; // BE CAREFUL!
  ttmp = RUNIT * expab_dt; // BE CAREFUL!
  for (int inumer = 0; inumer <= expab_npade; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < expab_dpade; idenom ++) {
    cd1[idenom] *= tfac;
  }
  expab_phi.gen(MPIP, IO, Bas, expab_dt, expab_h1type, expab_npade, expab_dpade, cn, cd0, cd1);
}
//////////////////////////////////////////////////////////////////////////
void clexpab::gen_phi5(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi5_n(max_dim_numer + 1);
  std::vector<dcomplex> phi5_d0(max_dim_denom);
  std::vector<dcomplex> phi5_d1(max_dim_denom);

  if (expab_ndpade == 22) {
    //
    //                  2 - x
    // phi5(x) = [----------------] 
    //               2              
    //           24(x  - 6 x + 12)
    // x0 = 3.0000000000000000000 - 1.7320508075688772935i
    // x0 = 3.0000000000000000000 + 1.7320508075688772935i
    dfac = 24.0;
    phi5_n[0] = 2.0;
    phi5_n[1] = -1.0;
    phi5_n[2] = 0.0;
    roots[0] = 3.0000000000000000000-1.7320508075688772935*IUNIT;
    roots[1] = 3.0000000000000000000+1.7320508075688772935*IUNIT;
  } else if (expab_ndpade == 13) {
    // Pade(1,3) of exp(x)
    //                    2
    //                 - x  + 2 x - 6
    // phi5(x) = [------------------------]
    //               3      2               
    //           24(x  - 6 x  + 18 x - 24)
    // x0 = 1.687091590520766641994056 - 2.508731754924880510838744i
    // x0 = 1.687091590520766641994056 + 2.508731754924880510838744i
    // x0 = 2.625816818958466716011889
    dfac = 24.0;
    phi5_n[0] = -6.0;
    phi5_n[1] = 2.0;
    phi5_n[2] = -1.0;
    roots[0] = 1.687091590520766641994056 - 2.508731754924880510838744*IUNIT;
    roots[1] = 1.687091590520766641994056 + 2.508731754924880510838744*IUNIT;
    roots[2] = 2.625816818958466716011889;
  } else if (expab_ndpade == 23) {
    // Bad: Pade(2,3) of exp(x)
    //                    2
    //                 - x  + 5 x - 12
    // phi5(x) = [------------------------]
    //               3      2               
    //           24(x  - 9 x  + 36 x - 60)
    // x0 = 2.6810828736277521339 + 3.0504301992474105694i
    // x0 = 2.6810828736277521339 - 3.0504301992474105694i
    // x0 = 3.6378342527444957322
    dfac = 24.0;
    phi5_n[0] = -12.0;
    phi5_n[1] = 5.0;
    phi5_n[2] = -1.0;
    roots[0] = 2.6810828736277521339 + 3.0504301992474105694*IUNIT;
    roots[1] = 2.6810828736277521339 - 3.0504301992474105694*IUNIT;
    roots[2] = 3.6378342527444957322;
  } else if (expab_ndpade == 33) {
    //                    2
    //                 - x  + 8 x - 24
    // phi5(x) = [--------------------------]
    //               3       2               
    //           24(x  - 12 x  + 60 x - 120)
    // x0 = 3.6778146453739144071 + 3.5087619195674433219i
    // x0 = 3.6778146453739144071 - 3.5087619195674433219i
    // x0 = 4.6443707092521711858
    dfac = 24.0;
    phi5_n[0] = -24.0;
    phi5_n[1] = 8.0;
    phi5_n[2] = -1.0;
    phi5_n[3] = 0.0;
    roots[0] = 3.6778146453739144071 + 3.5087619195674433219*IUNIT;
    roots[1] = 3.6778146453739144071 - 3.5087619195674433219*IUNIT;
    roots[2] = 4.6443707092521711858;
  } else if (expab_ndpade == 14) {
    // Pade(1,4) of exp(x)
    //                    3      2
    //                 - x  + 4 x  - 16 x + 24
    // phi5(x) = [ ----------------------------------]
    //                4      3       2       
    //            24(x  - 8 x  + 36 x  - 96 x + 120)
    // x0 = 1.2356535842848997066486550 - 3.4376524936710512909193657i
    // x0 = 1.2356535842848997066486550 + 3.4376524936710512909193657i
    // x0 = 2.7643464157151002933513450 - 1.1623236292832751740121108i
    // x0 = 2.7643464157151002933513450 + 1.1623236292832751740121108i
    dfac = 24.0;
    phi5_n[0] = 24.0;
    phi5_n[1] = -16.0;
    phi5_n[2] = 4.0;
    phi5_n[3] = -1.0;
    roots[0] = 1.2356535842848997066486550 - 3.4376524936710512909193657*IUNIT;
    roots[1] = 1.2356535842848997066486550 + 3.4376524936710512909193657*IUNIT;
    roots[2] = 2.7643464157151002933513450 - 1.1623236292832751740121108*IUNIT;
    roots[3] = 2.7643464157151002933513450 + 1.1623236292832751740121108*IUNIT;
  } else if (expab_ndpade == 24) {
    // Pade(2,4) of exp(x)
    //                    3      2
    //                 - x  + 8 x  - 36 x + 72
    // phi5(x) = [ ------------------------------------]
    //                4       3       2       
    //            24(x  - 12 x  + 72 x  - 240 x + 360)
    // x0 = 2.2209800329898068974239251 - 4.1603914455069319822284852i
    // x0 = 2.2209800329898068974239251 + 4.1603914455069319822284852i
    // x0 = 3.7790199670101931025760749 - 1.3801765242728430462268849i
    // x0 = 3.7790199670101931025760749 + 1.3801765242728430462268849i
    dfac = 24.0;
    phi5_n[0] = 72.0;
    phi5_n[1] = -36.0;
    phi5_n[2] = 8.0;
    phi5_n[3] = -1.0;
    roots[0] = 2.2209800329898068974239251 - 4.1603914455069319822284852*IUNIT;
    roots[1] = 2.2209800329898068974239251 + 4.1603914455069319822284852*IUNIT;
    roots[2] = 3.7790199670101931025760749 - 1.3801765242728430462268849*IUNIT;
    roots[3] = 3.7790199670101931025760749 + 1.3801765242728430462268849*IUNIT;
  } else if (expab_ndpade == 34) {
    // Bad: Pade(3,4) of exp(x)
    //                    3       2
    //                 - x  + 12 x  - 68 x + 168
    // phi5(x) = [ ------------------------------------]
    //                4       3        2       
    //            24(x  - 16 x  + 120 x  - 480 x + 840)
    // x0 = 3.2128068968715339829 - 4.7730874332766424998i
    // x0 = 3.2128068968715339829 + 4.7730874332766424998i
    // x0 = 4.7871931031284660171 - 1.5674764168952081241i
    // x0 = 4.7871931031284660171 + 1.5674764168952081241i
    dfac = 24.0;
    phi5_n[0] = 168.0;
    phi5_n[1] = -68.0;
    phi5_n[2] = 12.0;
    phi5_n[3] = -1.0;
    roots[0] = 3.2128068968715339829 - 4.7730874332766424998*IUNIT;
    roots[1] = 3.2128068968715339829 + 4.7730874332766424998*IUNIT;
    roots[2] = 4.7871931031284660171 - 1.5674764168952081241*IUNIT;
    roots[3] = 4.7871931031284660171 + 1.5674764168952081241*IUNIT;
  } else if (expab_ndpade == 44) {
    //                    3       2
    //                 - x  + 16 x  - 112 x + 336
    // phi5(x) = [ -------------------------------------]
    //                4       3        2       
    //            24(x  - 20 x  + 180 x  - 840 x + 1680)
    // x0 =  4.2075787943592556632 - 5.3148360837135054337i
    // x0 =  4.2075787943592556632 + 5.3148360837135054337i
    // x0 =  5.7924212056407443368 - 1.7344682578690075036i
    // x0 =  5.7924212056407443368 + 1.7344682578690075036i
    dfac = 24.0;
    phi5_n[0] = 336.0;
    phi5_n[1] = -112.0;
    phi5_n[2] = 16.0;
    phi5_n[3] = -1.0;
    phi5_n[4] = 0.0;
    roots[0] = 4.2075787943592556632 - 5.3148360837135054337*IUNIT;
    roots[1] = 4.2075787943592556632 + 5.3148360837135054337*IUNIT;
    roots[2] = 5.7924212056407443368 - 1.7344682578690075036*IUNIT;
    roots[3] = 5.7924212056407443368 + 1.7344682578690075036*IUNIT;
  } else {
    std::cout << "clexpab:gen_exp. bad expab_npade or expab_dpade." << std::endl;
    abort();
  }

  phi5_d1[0] = +dfac;
  phi5_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < expab_dpade; idim ++) {
    phi5_d1[idim] = +RUNIT;
    phi5_d0[idim] = -roots[idim];
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

  // ##### dt * phi5(-i*h*dt) #####
  cn = phi5_n;
  cd0 = phi5_d0;
  cd1 = phi5_d1;
  tfac = tunit * expab_dt; // BE CAREFUL!
  ttmp = RUNIT * expab_dt; // BE CAREFUL!
  for (int inumer = 0; inumer <= expab_npade; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < expab_dpade; idenom ++) {
    cd1[idenom] *= tfac;
  }
  expab_phi.gen(MPIP, IO, Bas, expab_dt, expab_h1type, expab_npade, expab_dpade, cn, cd0, cd1);
}
//////////////////////////////////////////////////////////////////////////
void clexpab::gen_phi2_phi0(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi0_n(max_dim_numer + 1);
  std::vector<dcomplex> phi0_d0(max_dim_denom);
  std::vector<dcomplex> phi0_d1(max_dim_denom);

  if (expab_ndpade == 11) {
    //            2 + x
    // phi0(x) = [-----]
    //            2 - x
    // x0 = 2
    dfac = -1.0;
    phi0_n[0] = 2.0;
    phi0_n[1] = 1.0;
    roots[0] = 2.0;
  } else if (expab_ndpade == 12) {
    //
    //               2 x + 6
    // phi0(x) = [------------] 
    //             2              
    //            x  - 4 x + 6
    // x0 = 2.0000000000000000000 - 1.4142135623730950488i
    // x0 = 2.0000000000000000000 + 1.4142135623730950488i
    dfac = 1.0;
    phi0_n[0] = 6.0;
    phi0_n[1] = 2.0;
    roots[0] = 2.0000000000000000000 - 1.4142135623730950488*IUNIT;
    roots[1] = 2.0000000000000000000 + 1.4142135623730950488*IUNIT;
  } else if (expab_ndpade == 22) {
    //             2           
    //            x  + 6 x + 12
    // phi0(x) = [-------------] 
    //             2              
    //            x  - 6 x + 12
    // x0 = 3.0000000000000000000 - 1.7320508075688772935i
    // x0 = 3.0000000000000000000 + 1.7320508075688772935i
    dfac = 1.0;
    phi0_n[0] = 12.0;
    phi0_n[1] = 6.0;
    phi0_n[2] = 1.0;
    roots[0] = 3.0000000000000000000-1.7320508075688772935*IUNIT;
    roots[1] = 3.0000000000000000000+1.7320508075688772935*IUNIT;
  } else if (expab_ndpade == 13) {
    // Pade(1,3) approximation of exp(x)
    //                 
    //                 -6 x - 24
    // phi0(x) = [---------------------]
    //               3      2               
    //            x  - 6 x  + 18 x - 24
    // x0 = 1.687091590520766641994056 - 2.508731754924880510838744i
    // x0 = 1.687091590520766641994056 + 2.508731754924880510838744i
    // x0 = 2.625816818958466716011889
    dfac = 1.0;
    phi0_n[0] = -24.0;
    phi0_n[1] = -6.0;
    phi0_n[2] = 0.0;
    roots[0] = 2.625816818958466716011889;
    roots[1] = 1.687091590520766641994056 - 2.508731754924880510838744*IUNIT;
    roots[2] = 1.687091590520766641994056 + 2.508731754924880510838744*IUNIT;
  } else if (expab_ndpade == 23) {
    // Pade(2,3) of exp(x)
    //                   2
    //              - 3 x  - 24 x - 60
    // phi0(x) = [-----------------------]
    //             3      2               
    //            x  - 9 x  + 36 x - 60
    // x0 = 2.6810828736277521339 - 3.0504301992474105694i
    // x0 = 2.6810828736277521339 + 3.0504301992474105694i
    // x0 = 3.6378342527444957322
    dfac = 1.0;
    phi0_n[0] = -60.0;
    phi0_n[1] = -24.0;
    phi0_n[2] = -3.0;
    roots[0] = 3.6378342527444957322;
    roots[1] = 2.6810828736277521339 - 3.0504301992474105694*IUNIT;
    roots[2] = 2.6810828736277521339 + 3.0504301992474105694*IUNIT;
  } else {
    std::cout << "clexpab:gen_exp. bad expab_npade or expab_dpade." << std::endl;
    abort();
  }

  phi0_d1[0] = +dfac;
  phi0_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < expab_dpade; idim ++) {
    phi0_d1[idim] = +RUNIT;
    phi0_d0[idim] = -roots[idim];
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

  // ##### phi0(-i*h*dt) #####
  cn = phi0_n;
  cd0 = phi0_d0;
  cd1 = phi0_d1;
  tfac = tunit * expab_dt; // BE CAREFUL!
  ttmp = RUNIT;            // BE CAREFUL!
  for (int inumer = 0; inumer <= expab_npade; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < expab_dpade; idenom ++) {
    cd1[idenom] *= tfac;
  }
  expab_phi0.gen(MPIP, IO, Bas, expab_dt, expab_h1type, expab_npade, expab_dpade, cn, cd0, cd1);
}
//////////////////////////////////////////////////////////////////////////
void clexpab::gen_phi2_phi1(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi1_n(max_dim_numer + 1);
  std::vector<dcomplex> phi1_d0(max_dim_denom);
  std::vector<dcomplex> phi1_d1(max_dim_denom);

  if (expab_ndpade == 11) {
    //              2
    // phi1(x) = [-----]
    //            2 - x
    // x0 = 2
    dfac = -1.0;
    phi1_n[0] = 2.0;
    phi1_n[1] = 0.0;
    roots[0] = 2.0;
  } else if (expab_ndpade == 12) {
    //
    //               - x + 6
    // phi1(x) = [------------] 
    //             2              
    //            x  - 4 x + 6
    // x0 = 2.0000000000000000000 - 1.4142135623730950488i
    // x0 = 2.0000000000000000000 + 1.4142135623730950488i
    dfac = 1.0;
    phi1_n[0] = 6.0;
    phi1_n[1] = -1.0;
    roots[0] = 2.0000000000000000000 - 1.4142135623730950488*IUNIT;
    roots[1] = 2.0000000000000000000 + 1.4142135623730950488*IUNIT;
  } else if (expab_ndpade == 22) {
    //            
    //                  12
    // phi1(x) = [-------------] 
    //             2              
    //            x  - 6 x + 12
    // x0 = 3.0000000000000000000 - 1.7320508075688772935i
    // x0 = 3.0000000000000000000 + 1.7320508075688772935i
    dfac = 1.0;
    phi1_n[0] = 12.0;
    phi1_n[1] = 0.0;
    phi1_n[2] = 0.0;
    roots[0] = 3.0000000000000000000-1.7320508075688772935*IUNIT;
    roots[1] = 3.0000000000000000000+1.7320508075688772935*IUNIT;
  } else if (expab_ndpade == 13) {
    // Pade(1,3) approximation of exp(x)
    //                   2
    //                - x  + 6 x - 24
    // phi1(x) = [---------------------]
    //               3      2               
    //            x  - 6 x  + 18 x - 24
    // x0 = 1.687091590520766641994056 - 2.508731754924880510838744i
    // x0 = 1.687091590520766641994056 + 2.508731754924880510838744i
    // x0 = 2.625816818958466716011889
    dfac = 1.0;
    phi1_n[0] = -24.0;
    phi1_n[1] = +6.0;
    phi1_n[2] = -1.0;
    roots[0] = 2.625816818958466716011889;
    roots[1] = 1.687091590520766641994056 - 2.508731754924880510838744*IUNIT;
    roots[2] = 1.687091590520766641994056 + 2.508731754924880510838744*IUNIT;
  } else if (expab_ndpade == 23) {
    // Pade(2,3) of exp(x)
    //                   2
    //                - x  + 6 x - 60
    // phi1(x) = [-----------------------]
    //             3      2               
    //            x  - 9 x  + 36 x - 60
    // x0 = 2.6810828736277521339 - 3.0504301992474105694i
    // x0 = 2.6810828736277521339 + 3.0504301992474105694i
    // x0 = 3.6378342527444957322
    dfac = 1.0;
    phi1_n[0] = -60.0;
    phi1_n[1] = 6.0;
    phi1_n[2] = -1.0;
    roots[0] = 3.6378342527444957322;
    roots[1] = 2.6810828736277521339 - 3.0504301992474105694*IUNIT;
    roots[2] = 2.6810828736277521339 + 3.0504301992474105694*IUNIT;
  } else {
    std::cout << "clexpab:gen_exp. bad expab_npade or expab_dpade." << std::endl;
    abort();
  }

  phi1_d1[0] = +dfac;
  phi1_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < expab_dpade; idim ++) {
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
  tfac = tunit * expab_dt; // BE CAREFUL!
  ttmp = RUNIT * expab_dt; // BE CAREFUL!
  for (int inumer = 0; inumer <= expab_npade; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (int idenom = 0; idenom < expab_dpade; idenom ++) {
    cd1[idenom] *= tfac;
  }
  expab_phi1.gen(MPIP, IO, Bas, expab_dt, expab_h1type, expab_npade, expab_dpade, cn, cd0, cd1);
}
////////////////////////////////////////////////////////////////////////
