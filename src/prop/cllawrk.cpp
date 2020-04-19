////////////////////////////////////////////////////////////////////////
// Lawson's exponential RK4
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
cllawrk::cllawrk()
{
}
////////////////////////////////////////////////////////////////////////
cllawrk::~cllawrk()
{
}
////////////////////////////////////////////////////////////////////////
cllawrk::cllawrk(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		 const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void cllawrk::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		  const clfield& Field, const clhprod& HPW)
{
  IO.read_info("lawrk_order", LFOUR, lawrk_order);
  IO.read_info("lawrk_npade", (long) 11, lawrk_npade);

  if (lawrk_npade != 11 &&
      lawrk_npade != 12 &&
      lawrk_npade != 22) {
    std::cout << "cllawrk::gen: bad lawrk_npade." << std::endl;
    abort();
  }
  if (lawrk_order != 1 &&
      lawrk_order != 2 &&
      lawrk_order != 3 &&
      lawrk_order != 4) {
    std::cout << "cllawrk::gen: bad lawrk_order." << std::endl;
    abort();
  }

  long icomp = Field.td_type;
  long isplit = clcontrol::split_type;
  lawrk_dt1 = Field.dtime;
  lawrk_dt2 = Field.dtime * HALF;

  if (lawrk_npade == 11) {
    dPade2p.gen(MPIP, IO, Bas, lawrk_dt2, CTWO, icomp, isplit);
  } else if (lawrk_npade == 12) {
    dcomplex alphp = TWO + IUNIT * sqrt(TWO);
    dcomplex alphm = TWO - IUNIT * sqrt(TWO);
    dPade2p.gen(MPIP, IO, Bas, lawrk_dt2, alphp, icomp, isplit);
    dPade2m.gen(MPIP, IO, Bas, lawrk_dt2, alphm, icomp, isplit);
  } else if (lawrk_npade == 22) {
    dcomplex alphp = THREE + IUNIT * sqrt(THREE);
    dcomplex alphm = THREE - IUNIT * sqrt(THREE);
    dPade2p.gen(MPIP, IO, Bas, lawrk_dt2, alphp, icomp, isplit);
    dPade2m.gen(MPIP, IO, Bas, lawrk_dt2, alphm, icomp, isplit);
  }

  Wfn0.gen(MPIP, IO, Bas);
  Wfn1.gen(MPIP, IO, Bas);
//OLD  Wfn2.gen(MPIP, IO, Bas);
//OLD  Wfn3.gen(MPIP, IO, Bas);
  hWfn0.gen(MPIP, IO, Bas);
  hWfn1.gen(MPIP, IO, Bas);
//OLD  hWfn2.gen(MPIP, IO, Bas);
//OLD  hWfn3.gen(MPIP, IO, Bas);
  v0Wfn.gen(MPIP, IO, Bas);
  vWfn.gen(MPIP, IO, Bas);
  xWfn.gen(MPIP, IO, Bas);

//DEBUG
//  std::cout << "abort for debug @ cllawrk::gen." << std::endl;
//  abort();
//DEBUG
}
////////////////////////////////////////////////////////////////////////
void cllawrk::prop(const clmpi& Proc, const clbas& Bas, 
		   const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  if (lawrk_order == 1) {
    prop1(Proc, Bas, Field, HPW, Wfn);
  } else if (lawrk_order == 2) {
    prop2(Proc, Bas, Field, HPW, Wfn);
  } else if (lawrk_order == 3) {
    prop3(Proc, Bas, Field, HPW, Wfn);
  } else if (lawrk_order == 4) {
    prop4(Proc, Bas, Field, HPW, Wfn);
//old    double lfield[9];
//old    Field.get_value(lfield);
//old    HPW.copy(Proc, Bas, Wfn, Wfn0);
//old    HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
//old    prop4_old(Proc, Bas, Field, HPW, Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void cllawrk::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  if (lawrk_order == 1) {
    prop1(Proc, Bas, Field, HPW, Wfn);
  } else if (lawrk_order == 2) {
    prop2(Proc, Bas, Field, HPW, Wfn);
  } else if (lawrk_order == 3) {
    prop3(Proc, Bas, Field, HPW, Wfn);
  } else if (lawrk_order == 4) {
    prop4(Proc, Bas, Field, HPW, Wfn);
//old    double lfield[9];
//old    Field.get_value(lfield);
//old    HPW.copy(Proc, Bas, Wfn, Wfn0);
//old    HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
//old    prop4_old(Proc, Bas, Field, HPW, Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void cllawrk::prop1(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double lfield[9];

  Field.get_value(lfield);
  HPW.copy(Proc, Bas, Wfn, Wfn0);
  HPW.htot(Proc, Bas, Field.dtime, lfield, Wfn0, hWfn0);
  HPW.xpy(Proc, Bas, hWfn0, Wfn0);

  expwfn(Proc, Bas, Field, HPW, Wfn0, Wfn1);
  expwfn(Proc, Bas, Field, HPW, Wfn1, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cllawrk::prop2(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double lfield[9];
  dcomplex cdt1 = Field.dtime;
  dcomplex cdt2 = Field.dtime * HALF;
  double thalf = Field.time + Field.dtime * HALF;

  HPW.copy(Proc, Bas, Wfn, Wfn0);

  Field.get_value(lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  expwfn(Proc, Bas, Field, HPW, Wfn0, Wfn1);
  expwfn(Proc, Bas, Field, HPW, hWfn0, hWfn1);

  HPW.copy(Proc, Bas, Wfn1, Wfn);
  HPW.axpyz(Proc, Bas, cdt2, hWfn1, Wfn1, Wfn0);

  Field.get_value(thalf, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  HPW.axpyz(Proc, Bas, cdt1, hWfn0, Wfn, Wfn1);
  expwfn(Proc, Bas, Field, HPW, Wfn1, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cllawrk::prop3(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double lfield[9];
  dcomplex cdt1  = Field.dtime;
  dcomplex cdt2  = Field.dtime / TWO;
  dcomplex cdt6  = Field.dtime / SIX;
  dcomplex cdt23 = Field.dtime * TWO / THREE;
  double time0 = Field.time;
  double thalf = Field.time + Field.dtime * HALF;
  double time1 = Field.time + Field.dtime;

  HPW.copy(Proc, Bas, Wfn, Wfn0);

  Field.get_value(time0, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  expwfn(Proc, Bas, Field, HPW, Wfn0, Wfn1);
  expwfn(Proc, Bas, Field, HPW, hWfn0, hWfn1);

  HPW.axpyz(Proc, Bas, cdt6, hWfn1, Wfn1, Wfn);
  HPW.axpyz(Proc, Bas, cdt2, hWfn1, Wfn1, Wfn0);

  Field.get_value(thalf, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  HPW.axpy(Proc, Bas, cdt23, hWfn0, Wfn);
  HPW.axpbyz(Proc, Bas, CTWO, hWfn0, -RUNIT, hWfn1, hWfn1);
  HPW.axpy(Proc, Bas, cdt1, hWfn1, Wfn1);
  expwfn(Proc, Bas, Field, HPW, Wfn1, Wfn0);

  Field.get_value(time1, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  expwfn(Proc, Bas, Field, HPW, Wfn, Wfn1);
  HPW.axpyz(Proc, Bas, cdt6, hWfn0, Wfn1, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cllawrk::prop4(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double lfield[9];
  dcomplex cdt1 = Field.dtime;
  dcomplex cdt2 = Field.dtime / TWO;
  dcomplex cdt3 = Field.dtime / THREE;
  dcomplex cdt6 = Field.dtime / SIX;
  double time0 = Field.time;
  double thalf = Field.time + Field.dtime * HALF;
  double time1 = Field.time + Field.dtime;

  HPW.copy(Proc, Bas, Wfn, Wfn0);

  Field.get_value(time0, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  expwfn(Proc, Bas, Field, HPW, Wfn0, Wfn1);
  expwfn(Proc, Bas, Field, HPW, hWfn0, hWfn1);
  HPW.axpyz(Proc, Bas, cdt6, hWfn1, Wfn1, Wfn);
  HPW.axpyz(Proc, Bas, cdt2, hWfn1, Wfn1, Wfn0);

  Field.get_value(thalf, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  HPW.axpy(Proc, Bas, cdt3, hWfn0, Wfn);
  HPW.axpyz(Proc, Bas, cdt2, hWfn0, Wfn1, Wfn0);

  Field.get_value(thalf, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  HPW.axpy(Proc, Bas, cdt3, hWfn0, Wfn);
  HPW.axpy(Proc, Bas, cdt1, hWfn0, Wfn1);
  expwfn(Proc, Bas, Field, HPW, Wfn1, Wfn0);

  Field.get_value(time1, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  expwfn(Proc, Bas, Field, HPW, Wfn, Wfn1);
  HPW.axpyz(Proc, Bas, cdt6, hWfn0, Wfn1, Wfn);
}
////////////////////////////////////////////////////////////////////////
//oldvoid cllawrk::prop4_old(const clmpi& Proc, const clbas& Bas,
//old			const clfield& Field, clhprod& HPW, clwfn& Wfn)
//old{
//old  double lfield[9];
//old  double ttime;
//old
//old  const dcomplex cdt1 = lawrk_dt1 * RUNIT;
//old  const dcomplex cdt2 = lawrk_dt1 * CHALF;
//old  const dcomplex cdt3 = lawrk_dt1 / THREE;
//old  const dcomplex cdt6 = lawrk_dt1 / SIX;
//old
//old  // ***** step (1) *****
//old  // CI-RK4
//old  HPW.axpyc(Proc, Bas, cdt6, hWfn0, Wfn);
//old  HPW.axpbyzc(Proc, Bas, RUNIT, Wfn0, cdt2, hWfn0, Wfn1);
//old  // MO-OIFS4
//old  expwfn(Proc, Bas, Field, HPW, Wfn0, v0Wfn);
//old  expwfn(Proc, Bas, Field, HPW, hWfn0, vWfn);
//old  HPW.copyo(Proc, Bas, v0Wfn, xWfn);
//old  HPW.axpyo(Proc, Bas, cdt6, vWfn, xWfn);
//old  HPW.axpyzo(Proc, Bas, cdt2, vWfn, v0Wfn, Wfn1);
//old
//old  // ***** step (2) *****
//old  ttime = Field.time + lawrk_dt2;
//old  Field.get_value(ttime, lfield);
//old  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
//old  // CI-RK4
//old  HPW.axpyc(Proc, Bas, cdt3, hWfn1, Wfn);
//old  HPW.axpbyzc(Proc, Bas, RUNIT, Wfn0, cdt2, hWfn1, Wfn1);
//old  // MO-OIFS4
//old  HPW.axpyo(Proc, Bas, cdt3, hWfn1, xWfn);
//old  HPW.axpyzo(Proc, Bas, cdt2, hWfn1, v0Wfn, Wfn1);
//old
//old  // ***** step (3) *****
//old  ttime = Field.time + lawrk_dt2;
//old  Field.get_value(ttime, lfield);
//old  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
//old  // CI-RK4
//old  HPW.axpyc(Proc, Bas, cdt3, hWfn1, Wfn);
//old  HPW.axpbyzc(Proc, Bas, RUNIT, Wfn0, cdt1, hWfn1, Wfn1);
//old  // MO-OIFS4
//old  HPW.axpyo(Proc, Bas, cdt3, hWfn1, xWfn);
//old  HPW.axpyo(Proc, Bas, cdt1, hWfn1, v0Wfn); // v0Wfn is modified here!
//old  expwfn(Proc, Bas, Field, HPW, v0Wfn, Wfn1);
//old
//old  // ***** step (4) *****
//old  ttime = Field.time + lawrk_dt1;
//old  Field.get_value(ttime, lfield);
//old  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
//old  // CI-RK4
//old  HPW.axpyc(Proc, Bas, cdt6, hWfn1, Wfn);
//old  // MO-OIFS4
//old  expwfn(Proc, Bas, Field, HPW, xWfn, Wfn);
//old  HPW.axpyo(Proc, Bas, cdt6, hWfn1, Wfn);
//old}
////////////////////////////////////////////////////////////////////////
void cllawrk::expwfn(const clmpi& Proc, const clbas& Bas, const clfield& Field, 
		     clhprod& HPW, const clwfn& Wfn, clwfn& ExpWfn)
{
  if (clcontrol::split_type == 0) {
    HPW.copy(Proc, Bas, Wfn, ExpWfn);
  } else if (lawrk_npade == 11) {
    expwfn11(Proc, Bas, Field, HPW, Wfn, ExpWfn);
  } else if (lawrk_npade == 12) {
    expwfn12(Proc, Bas, Field, HPW, Wfn, ExpWfn);
  } else if (lawrk_npade == 22) {
    expwfn22(Proc, Bas, Field, HPW, Wfn, ExpWfn);
  }
}
////////////////////////////////////////////////////////////////////////
void cllawrk::expwfn11(const clmpi& Proc, const clbas& Bas, const clfield& Field, 
		       clhprod& HPW, const clwfn& Wfn, clwfn& ExpWfn)
{
  if (dPade2p.dpade_midpt) {
    double thalf = Field.time + lawrk_dt2;
//new
    HPW.copy(Proc, Bas, Wfn, ExpWfn);
//old    ExpWfn.clearo(Bas);
//old    HPW.copyo(Proc, Bas, Wfn, ExpWfn);
    dPade2p.prod(Proc, Bas, thalf, Field, HPW, ExpWfn);
    HPW.axpbyzo(Proc, Bas, -RUNIT, Wfn, CTWO, ExpWfn, ExpWfn);
  } else {
    std::cout << "cllawrk::expwfn11: only midpt" << std::endl;
    abort();

    double lfield[9];
    double time0 = Field.time;
    double time1 = Field.time + lawrk_dt1;
    ExpWfn.clearo(Bas);
    Field.get_value(time0, lfield);
    HPW.copyo(Proc, Bas, Wfn, ExpWfn);
    //*****
    std::cout << "cllawrk::expwfn11. h1tot DOES NOT multiply dt factor!" << std::endl;
    abort();
    HPW.h1tot(Proc, Bas, lawrk_dt2, lfield, Wfn, ExpWfn);
    //*****
    dPade2p.prod(Proc, Bas, time1, Field, HPW, ExpWfn);
  }
}
////////////////////////////////////////////////////////////////////////
void cllawrk::expwfn12(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		       clhprod& HPW, const clwfn& Wfn, clwfn& ExpWfn)
{
  std::cout << "cllawrk::expwfn12 nyi." << std::endl;
  abort();

  dcomplex F13 = RUNIT / THREE;
//new
  HPW.copy(Proc, Bas, Wfn, ExpWfn);
//old  ExpWfn.clearo(Bas);
  //OLD  dPadeT.tprod(Proc, Bas, dt, Wfn, ExpWfn);
  HPW.axpbyzo(Proc, Bas, RUNIT, Wfn, F13, ExpWfn, ExpWfn);
//nyi about time  dPade2m.prod(Proc, Bas, time, Field, HPW, ExpWfn);
//nyi about time  dPade2p.prod(Proc, Bas, time, Field, HPW, ExpWfn);
}
////////////////////////////////////////////////////////////////////////
void cllawrk::expwfn22(const clmpi& Proc, const clbas& Bas, const clfield& Field, 
		       clhprod& HPW, const clwfn& Wfn, clwfn& ExpWfn)
{
  std::cout << "cllawrk::expwfn22 nyi." << std::endl;
  abort();

//new
  HPW.copy(Proc, Bas, Wfn, ExpWfn);
//old  ExpWfn.clearo(Bas);
  //OLD  dPadeT.tprod(Proc, Bas, dt, Wfn, ExpWfn);
//nyi about time  dPade2m.prod(Proc, Bas, time, Field, HPW, ExpWfn);
//nyi about time  dPade2p.prod(Proc, Bas, time, Field, HPW, ExpWfn);
  HPW.xpyzo(Proc, Bas, Wfn, ExpWfn, ExpWfn);
}
////////////////////////////////////////////////////////////////////////
