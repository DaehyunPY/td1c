////////////////////////////////////////////////////////////////////////
// Lawson Adams-Bashforth with Pade approximation for exponential
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
int cllawab::max_dim_numer = 4;
int cllawab::max_dim_denom = 4;
int cllawab::max_order = 4;
dcomplex cllawab::coeffAB[10][10];
dcomplex cllawab::coeffAM[10][10];
////////////////////////////////////////////////////////////////////////
cllawab::cllawab()
{
  std::cout << "cllawab" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cllawab::~cllawab()
{
  std::cout << "~cllawab" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cllawab::cllawab(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		 const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void cllawab::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		  const clfield& Field, const clhprod& HPW)
{
  std::cout << "cllawab::gen" << std::endl;
  clh2prop::gen_basic(MPIP, IO, Bas, Field, HPW);

  num_step = 0;
  lawab_dt = Field.dtime;
  lawab_cdt = lawab_dt;

  IO.read_info("lawab_order", 4, lawab_order);
  IO.read_info("lawab_ndpade", 13, lawab_ndpade);
  IO.read_info("lawab_h1type", -1, lawab_h1type);
  IO.read_info("lawab_ncorr", 1, lawab_ncorr);

  if (lawab_order < 1 || lawab_order > max_order) {
    std::cout << "cllawab::gen bad lawab_order." << std::endl;
    abort();
  }
  for (int order_ = 1; order_ <= lawab_order; order_ ++) {
    get_coeff(order_, coeffAB[order_], coeffAM[order_]); 
  }

  if (lawab_ndpade == 11) {
    lawab_npade = 1;
    lawab_dpade = 1;
  } else if (lawab_ndpade == 12) {
    lawab_npade = 1;
    lawab_dpade = 2;
  } else if (lawab_ndpade == 22) {
    lawab_npade = 2;
    lawab_dpade = 2;
  } else if (lawab_ndpade == 23) {
    lawab_npade = 2;
    lawab_dpade = 3;
  } else if (lawab_ndpade == 33) {
    lawab_npade = 3;
    lawab_dpade = 3;
  } else if (lawab_ndpade == 34) {
    lawab_npade = 3;
    lawab_dpade = 4;
  } else if (lawab_ndpade == 44) {
    lawab_npade = 4;
    lawab_dpade = 4;
  } else {
    std::cout << "cllawab::gen bad lawab_ndpade." << std::endl;
    abort();
  }

  lawab_exp.resize(lawab_order);
  gen_exp(MPIP, IO, Bas);
  //  printf("cllawab::phi:\n"); lawab_exp.print();

  timeP.resize(lawab_order+1);
  tWfn.gen(MPIP, IO, Bas);
  WfnP.resize(lawab_order+1);
  hWfnP.resize(lawab_order+1);
  for (int is = 0; is <= lawab_order; is ++) {
    WfnP[is].gen(MPIP, IO, Bas);
    hWfnP[is].gen(MPIP, IO, Bas);
  }
}
////////////////////////////////////////////////////////////////////////
void cllawab::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  std::cout << "cllawab::prop (1) nyi." << std::endl;
  abort();
}
////////////////////////////////////////////////////////////////////////
void cllawab::prop(const clmpi& Proc, const clbas& Bas,
		   const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  if (lawab_order == 1) {
    //    prop1(Proc, Bas, Field, HPW, Wfn);
    prop_general(Proc, Bas, 1, Field, HPW, Wfn);

  } else if (lawab_order == 2) {
    //    prop2(Proc, Bas, Field, HPW, Wfn);
    //    prop_general(Proc, Bas, 2, Field, HPW, Wfn);
    if (num_step == 0) {
      prop_general(Proc, Bas, 1, Field, HPW, Wfn);
    } else {
      prop_general(Proc, Bas, 2, Field, HPW, Wfn);
    }

  } else if (lawab_order == 3) {
    if (num_step == 0) {
      prop_general(Proc, Bas, 1, Field, HPW, Wfn);
    } else if (num_step == 1) {
      prop_general(Proc, Bas, 2, Field, HPW, Wfn);
    } else {
      prop_general(Proc, Bas, 3, Field, HPW, Wfn);
    }

  } else if (lawab_order == 4) {
    //    prop_general(Proc, Bas, Field, HPW, Wfn);
    //    prop_general(Proc, Bas, lawab_order, Field, HPW, Wfn);
    if (num_step == 0) {
      prop_general(Proc, Bas, 1, Field, HPW, Wfn);
    } else if (num_step == 1) {
      prop_general(Proc, Bas, 2, Field, HPW, Wfn);
    } else if (num_step == 2) {
      prop_general(Proc, Bas, 3, Field, HPW, Wfn);
    } else {
      prop_general(Proc, Bas, 4, Field, HPW, Wfn);
    }
  } else {
    std::cout << "cllawab:prop. bad lawab_order." << std::endl;
    abort();    
  }
}
////////////////////////////////////////////////////////////////////////
void cllawab::prop1(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
}
////////////////////////////////////////////////////////////////////////
void cllawab::prop2(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
}
////////////////////////////////////////////////////////////////////////
void cllawab::prop_general(const clmpi& Proc, const clbas& Bas,
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
    timeP[1] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[1]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[1]);
    for (int is = 2; is <= lawab_order; is ++) {
      timeP[is] = timeP[1];
      HPW.copy(Proc, Bas, WfnP[1], WfnP[is]);
      HPW.copy(Proc, Bas, hWfnP[1], hWfnP[is]);
    }
  } else {
    for (int is = lawab_order; is >= 2; is --) {
      timeP[is] = timeP[is-1];
      HPW.copy(Proc, Bas, WfnP[is-1], WfnP[is]);
      HPW.copy(Proc, Bas, hWfnP[is-1], hWfnP[is]);
    }
    timeP[1] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[1]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[1]);
  }
  num_step ++;

  dcomplex cAB[10];
  dcomplex cAM[10];
  get_coeff(cAB, cAM);

  HPW.copy(Proc, Bas, WfnP[1], tWfn);
  //  lawab_exp[0].prod(Proc, Bas, timeP[1], Field, HPW, tWfn);
  //  lawab_exp[0].prod(Proc, Bas, timeP[1], lawab_cdt, Field, HPW, tWfn);
  lawab_exp[0].prod(Proc, Bas, timeP[1], RUNIT, Field, HPW, tWfn);
  HPW.copy(Proc, Bas, tWfn, Wfn);
  HPW.copy(Proc, Bas, tWfn, WfnP[0]);

  HPW.copy(Proc, Bas, hWfnP[1], tWfn);
  //  lawab_exp[0].prod(Proc, Bas, timeP[1], Field, HPW, tWfn);
  lawab_exp[0].prod(Proc, Bas, timeP[1], lawab_cdt, Field, HPW, tWfn);
  HPW.axpy(Proc, Bas, cAB[1], tWfn, WfnP[0]);
  HPW.axpy(Proc, Bas, cAM[1], tWfn, Wfn);

  for (int is = 2; is <= lawab_order; is ++) {
    // td-part of h1
    double lfieldP[9];
    HPW.copy(Proc, Bas, hWfnP[is], tWfn);
    Field.get_value(timeP[is], lfieldP);
    for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
    HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[is], tWfn);

    //    lawab_exp[is-1].prod(Proc, Bas, timeP[is], Field, HPW, tWfn);
    lawab_exp[is-1].prod(Proc, Bas, timeP[is], lawab_cdt, Field, HPW, tWfn);
    HPW.axpy(Proc, Bas, cAB[is], tWfn, WfnP[0]);
    HPW.axpy(Proc, Bas, cAM[is], tWfn, Wfn);
  }

  if (lawab_ncorr == 0) {
    HPW.copy(Proc, Bas, WfnP[0], Wfn);
  } else if (lawab_ncorr >= 1) {
    double lfieldP[9];
    timeP[0] = Field.time + Field.dtime;
    HPW.copy(Proc, Bas, Wfn, tWfn);

    Field.get_value(timeP[0], lfieldP);
    HPW.htot(Proc, Bas, ONE, lfieldP, WfnP[0], hWfnP[0]);
    for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
    HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[0], hWfnP[0]);
    HPW.axpy(Proc, Bas, lawab_cdt*cAM[0], hWfnP[0], Wfn);

    for (int icorr = 2; icorr <= lawab_ncorr; icorr++) {
      HPW.copy(Proc, Bas, Wfn, WfnP[0]);
      HPW.copy(Proc, Bas, tWfn, Wfn);
      Field.get_value(timeP[0], lfieldP);
      HPW.htot(Proc, Bas, ONE, lfieldP, WfnP[0], hWfnP[0]);
      for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
      HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[0], hWfnP[0]);
      HPW.axpy(Proc, Bas, lawab_cdt*cAM[0], hWfnP[0], Wfn);
    }
  }
}
////////////////////////////////////////////////////////////////////////
void cllawab::prop_general(const clmpi& Proc, const clbas& Bas, int order,
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
//  if (num_step == 0) {
//    timeP[1] = Field.time;
//    Field.get_value(lfield);
//    HPW.copy(Proc, Bas, Wfn, WfnP[1]);
//    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[1]);
//    for (int is = 2; is <= order; is ++) {
//      timeP[is] = timeP[1];
//      HPW.copy(Proc, Bas, WfnP[1], WfnP[is]);
//      HPW.copy(Proc, Bas, hWfnP[1], hWfnP[is]);
//    }
//  } else {
    for (int is = order; is >= 2; is --) {
      timeP[is] = timeP[is-1];
      HPW.copy(Proc, Bas, WfnP[is-1], WfnP[is]);
      HPW.copy(Proc, Bas, hWfnP[is-1], hWfnP[is]);
    }
    timeP[1] = Field.time;
    Field.get_value(lfield);
    HPW.copy(Proc, Bas, Wfn, WfnP[1]);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfnP[1]);
//  }
  num_step ++;

  HPW.copy(Proc, Bas, WfnP[1], tWfn);
  //  lawab_exp[0].prod(Proc, Bas, timeP[1], Field, HPW, tWfn);
  //  lawab_exp[0].prod(Proc, Bas, timeP[1], lawab_cdt, Field, HPW, tWfn);
  lawab_exp[0].prod(Proc, Bas, timeP[1], RUNIT, Field, HPW, tWfn);
  HPW.copy(Proc, Bas, tWfn, Wfn);
  HPW.copy(Proc, Bas, tWfn, WfnP[0]);

  HPW.copy(Proc, Bas, hWfnP[1], tWfn);
  //  lawab_exp[0].prod(Proc, Bas, timeP[1], Field, HPW, tWfn);
  lawab_exp[0].prod(Proc, Bas, timeP[1], lawab_cdt, Field, HPW, tWfn);
  HPW.axpy(Proc, Bas, coeffAB[order][1], tWfn, WfnP[0]);
  HPW.axpy(Proc, Bas, coeffAM[order][1], tWfn, Wfn);

  for (int is = 2; is <= order; is ++) {
    // td-part of h1
    double lfieldP[9];
    HPW.copy(Proc, Bas, hWfnP[is], tWfn);
    Field.get_value(timeP[is], lfieldP);
    for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
    HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[is], tWfn);

    //    lawab_exp[is-1].prod(Proc, Bas, timeP[is], Field, HPW, tWfn);
    lawab_exp[is-1].prod(Proc, Bas, timeP[is], lawab_cdt, Field, HPW, tWfn);
    HPW.axpy(Proc, Bas, coeffAB[order][is], tWfn, WfnP[0]);
    HPW.axpy(Proc, Bas, coeffAM[order][is], tWfn, Wfn);
  }

  if (lawab_ncorr == 0) {
    HPW.copy(Proc, Bas, WfnP[0], Wfn);
  } else if (lawab_ncorr >= 1) {
    double lfieldP[9];
    timeP[0] = Field.time + Field.dtime;
    HPW.copy(Proc, Bas, Wfn, tWfn);

    Field.get_value(timeP[0], lfieldP);
    HPW.htot(Proc, Bas, ONE, lfieldP, WfnP[0], hWfnP[0]);
    for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
    HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[0], hWfnP[0]);
    HPW.axpy(Proc, Bas, lawab_cdt*coeffAM[order][0], hWfnP[0], Wfn);

    for (int icorr = 2; icorr <= lawab_ncorr; icorr++) {
      HPW.copy(Proc, Bas, Wfn, WfnP[0]);
      HPW.copy(Proc, Bas, tWfn, Wfn);
      Field.get_value(timeP[0], lfieldP);
      HPW.htot(Proc, Bas, ONE, lfieldP, WfnP[0], hWfnP[0]);
      for (int ii = 0; ii < 9; ii++) lfieldP[ii] -= lfield[ii];
      HPW.v1ext(Proc, Bas, tunit, lfieldP, WfnP[0], hWfnP[0]);
      HPW.axpy(Proc, Bas, lawab_cdt*coeffAM[order][0], hWfnP[0], Wfn);
    }
  }
}
//////////////////////////////////////////////////////////////////////////
void cllawab::gen_exp(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi0_n(max_dim_numer + 1);
  std::vector<dcomplex> phi0_d0(max_dim_denom);
  std::vector<dcomplex> phi0_d1(max_dim_denom);

  if (lawab_ndpade == 11) {
    //            2 + x
    // phi0(x) = [-----]
    //            2 - x
    // x0 = 2
    dfac = -1.0;
    phi0_n[0] = 2.0;
    phi0_n[1] = 1.0;
    roots[0] = 2.0;
  } else if (lawab_ndpade == 12) {
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
  } else if (lawab_ndpade == 22) {
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
  } else if (lawab_ndpade == 23) {
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
  } else if (lawab_ndpade == 33) {
    //             3       2             
    //            x  - 12 x  - 60 x - 120 
    // phi0(x) = [-----------------------]
    //             3       2               
    //            x  - 12 x  + 60 x - 120 
    // x0 = 3.6778146453739144071 + 3.5087619195674433219i
    // x0 = 3.6778146453739144071 - 3.5087619195674433219i
    // x0 = 4.6443707092521711858
    dfac = 1.0;
    phi0_n[0] = -120.0;
    phi0_n[1] = -60.0;
    phi0_n[2] = -12.0;
    phi0_n[3] = 1.0;
    roots[0] = 3.6778146453739144071 + 3.5087619195674433219*IUNIT;
    roots[1] = 3.6778146453739144071 - 3.5087619195674433219*IUNIT;
    roots[2] = 4.6443707092521711858;
  } else if (lawab_ndpade == 34) {
    //                    3       2
    //                 4 x  + 60 x  + 360 x + 840
    // phi0(x) = [ ----------------------------------]
    //              4       3        2       
    //             x  - 16 x  + 120 x  - 480 x + 840
    // x0 = 3.2128068968715339829 - 4.7730874332766424998i
    // x0 = 3.2128068968715339829 + 4.7730874332766424998i
    // x0 = 4.7871931031284660171 - 1.5674764168952081241i
    // x0 = 4.7871931031284660171 + 1.5674764168952081241i
    dfac = 1.0;
    phi0_n[0] = 840.0;
    phi0_n[1] = 360.0;
    phi0_n[2] = 60.0;
    phi0_n[3] = 4.0;
    roots[0] = 3.2128068968715339829 - 4.7730874332766424998*IUNIT;
    roots[1] = 3.2128068968715339829 + 4.7730874332766424998*IUNIT;
    roots[2] = 4.7871931031284660171 - 1.5674764168952081241*IUNIT;
    roots[3] = 4.7871931031284660171 + 1.5674764168952081241*IUNIT;
  } else if (lawab_ndpade == 44) {
    //              4       3        2       	     
    //             x  + 20 x  + 180 x  + 840 x + 1680
    // phi0(x) = [ ----------------------------------]
    //              4       3        2       
    //             x  - 20 x  + 180 x  - 840 x + 1680
    // x0 = 
    // x0 = 
    // x0 = 
    // x0 = 
    dfac = 1.0;
    phi0_n[0] = 1680.0;
    phi0_n[1] = 840.0;
    phi0_n[2] = 180.0;
    phi0_n[3] = 20.0;
    phi0_n[4] = 1.0;
    roots[0] = 4.2075787943592556632 - 5.3148360837135054337*IUNIT;
    roots[1] = 4.2075787943592556632 + 5.3148360837135054337*IUNIT;
    roots[2] = 5.7924212056407443368 - 1.7344682578690075036*IUNIT;
    roots[3] = 5.7924212056407443368 + 1.7344682578690075036*IUNIT;
  } else {
    std::cout << "cllawab:gen_exp. bad lawab_npade or lawab_dpade." << std::endl;
    abort();
  }

  phi0_d1[0] = +dfac;
  phi0_d0[0] = -dfac * roots[0];
  for (int idim = 1; idim < lawab_dpade; idim ++) {
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

  // ##### exp(-i*h*dt*js) #####
  for (int js = 1; js <= lawab_order; js ++) {
    double dtj = lawab_dt * js;
    cn = phi0_n;
    cd0 = phi0_d0;
    cd1 = phi0_d1;
    tfac = tunit * dtj; // BE CAREFUL!
    ttmp = RUNIT;       // BE CAREFUL!
    for (int inumer = 0; inumer <= lawab_npade; inumer ++) {
      cn[inumer] *= ttmp;
      ttmp *= tunit;
    }
    for (int idenom = 0; idenom < lawab_dpade; idenom ++) {
      cd1[idenom] *= tfac;
    }
    lawab_exp[js-1].gen(MPIP, IO, Bas, dtj, lawab_h1type, 
			lawab_npade, lawab_dpade, cn, cd0, cd1);
  }
}
////////////////////////////////////////////////////////////////////////
void cllawab::get_coeff(dcomplex* cAB, dcomplex* cAM) const
{
  if (lawab_order == 1) {
    cAB[1] = 1.0;

    cAM[0] = 0.5;
    cAM[1] = 0.5;
  } else if (lawab_order == 2) {
    cAB[1] =  1.5;
    cAB[2] = -0.5;

    cAM[0] =  5.0/12.0;
    cAM[1] =  8.0/12.0;
    cAM[2] = -1.0/12.0;
  } else if (lawab_order == 3) {
    cAB[1] =  23.0/12.0;
    cAB[2] = -16.0/12.0;
    cAB[3] =   5.0/12.0;

    cAM[0] =  9.0/24.0;
    cAM[1] = 19.0/24.0;
    cAM[2] = -5.0/24.0;
    cAM[3] =  1.0/24.0;
  } else if (lawab_order == 4) {
    cAB[1] =  55.0/24.0;
    cAB[2] = -59.0/24.0;
    cAB[3] =  37.0/24.0;
    cAB[4] =  -9.0/24.0;

    cAM[0] =  251.0/720.0;
    cAM[1] =  646.0/720.0;
    cAM[2] = -264.0/720.0;
    cAM[3] =  106.0/720.0;
    cAM[4] =  -19.0/720.0;
  }
}
////////////////////////////////////////////////////////////////////////
void cllawab::get_coeff(int order_, dcomplex* cAB, dcomplex* cAM)
{
  if (order_ == 1) {
    cAB[1] = 1.0;

    cAM[0] = 0.5;
    cAM[1] = 0.5;
  } else if (order_ == 2) {
    cAB[1] =  1.5;
    cAB[2] = -0.5;

    cAM[0] =  5.0/12.0;
    cAM[1] =  8.0/12.0;
    cAM[2] = -1.0/12.0;
  } else if (order_ == 3) {
    cAB[1] =  23.0/12.0;
    cAB[2] = -16.0/12.0;
    cAB[3] =   5.0/12.0;

    cAM[0] =  9.0/24.0;
    cAM[1] = 19.0/24.0;
    cAM[2] = -5.0/24.0;
    cAM[3] =  1.0/24.0;
  } else if (order_ == 4) {
    cAB[1] =  55.0/24.0;
    cAB[2] = -59.0/24.0;
    cAB[3] =  37.0/24.0;
    cAB[4] =  -9.0/24.0;

    cAM[0] =  251.0/720.0;
    cAM[1] =  646.0/720.0;
    cAM[2] = -264.0/720.0;
    cAM[3] =  106.0/720.0;
    cAM[4] =  -19.0/720.0;
  }
}
//////////////////////////////////////////////////////////////////////////
