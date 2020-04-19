////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
clcrnic::clcrnic()
{
}
////////////////////////////////////////////////////////////////////////
clcrnic::~clcrnic()
{
}
////////////////////////////////////////////////////////////////////////
clcrnic::clcrnic(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		 const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void clcrnic::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		  const clfield& Field, const clhprod& HPW)
{
  IO.read_info("crnic_maxcyc", LFOUR, crnic_maxcyc);

  int icomp = Field.td_type;
  int isplit = clcontrol::split_type;
  double dt2 = Field.dtime * HALF;
  dPade.gen(MPIP, IO, Bas, dt2, CTWO, icomp, isplit);

  Wfn0.gen(MPIP, IO, Bas);
  Wfn1.gen(MPIP, IO, Bas);
  hWfn0.gen(MPIP, IO, Bas);
  hWfn1.gen(MPIP, IO, Bas);
}
////////////////////////////////////////////////////////////////////////
void clcrnic::prop(const clmpi& Proc, const clbas& Bas, 
		   const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double res_wfn, lfield[9];
  double dt2 = Field.dtime * HALF;
  dcomplex cdt1 = Field.dtime * RUNIT;
  dcomplex cdt2 = Field.dtime * CHALF;
  double time0 = Field.time;
  double time1 = Field.time + Field.dtime;

  Field.get_value(time0, lfield);
  HPW.copy(Proc, Bas, Wfn, Wfn0);
  HPW.copy(Proc, Bas, Wfn, Wfn1);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn0);
  HPW.h1tot(Proc, Bas, dt2, lfield, Wfn1, Wfn0);

  HPW.axpyz(Proc, Bas, cdt1, hWfn0, Wfn0, Wfn1);
  dPade.prod(Proc, Bas, time1, Field, HPW, Wfn1);

  for (int icyc = 1; icyc < crnic_maxcyc; icyc ++) {
    Field.get_value(time1, lfield);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
    HPW.xpyz(Proc, Bas, hWfn1, hWfn0, Wfn1);
    HPW.axpyz(Proc, Bas, cdt2, Wfn1, Wfn0, Wfn1);
    dPade.prod(Proc, Bas, time1, Field, HPW, Wfn1);
  }

//old  HPW.copy(Proc, Bas, Wfn, Wfn0);
//old  HPW.copy(Proc, Bas, Wfn, Wfn1);
//old
//old  Field.get_value(time0, lfield);
//old  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn0, Wfn1);
//old  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
//old  HPW.xpy(Proc, Bas, hWfn0, Wfn1);
//old  HPW.scal(Proc, Bas, cdt1, Wfn1);
//old
//old  dPade.prod(Proc, Bas, time1, Field, HPW, Wfn1); // dC0 = C1 - C0
//old  res_wfn = get_res(Proc, Bas, Wfn1);
//old  printf(" clcrnic::prop: %10d%20.10e\n", LZERO, res_wfn);
//old  HPW.xpy(Proc, Bas, Wfn0, Wfn1);                 // C1 = C0 + dC0
//old  HPW.copy(Proc, Bas, Wfn1, Wfn0);
//old
//old  for (int icyc = 1; icyc < crnic_maxcyc; icyc ++) {
//old    Field.get_value(time1, lfield);
//old    HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
//old    HPW.xmyz(Proc, Bas, hWfn1, hWfn0, Wfn1);
//old    HPW.scal(Proc, Bas, cdt2, Wfn1);
//old    dPade.prod(Proc, Bas, time1, Field, HPW, Wfn1); // dC(i-1) = C(i) - C(i-1)
//old    res_wfn = get_res(Proc, Bas, Wfn1);
//old    printf(" clcrnic::prop: %10d%20.10e\n", icyc, res_wfn);
//old    HPW.xpy(Proc, Bas, Wfn0, Wfn1);                 // C(i) = C(i-1) + dC(i-1)
//old    HPW.copy(Proc, Bas, Wfn1, Wfn0);
//old    HPW.copy(Proc, Bas, hWfn1, hWfn0);
//old  }

  HPW.copy(Proc, Bas, Wfn1, Wfn);
}
////////////////////////////////////////////////////////////////////////
void clcrnic::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  double res_wfn, lfield[9];
  double dt2 = Field.dtime * HALF;
  dcomplex cdt1 = Field.dtime * RUNIT;
  dcomplex cdt2 = Field.dtime * CHALF;
  double time0 = Field.time;
  double time1 = Field.time + Field.dtime;

  Field.get_value(time0, lfield);
  HPW.copy(Proc, Bas, Wfn, Wfn0);
  HPW.copy(Proc, Bas, Wfn, Wfn1);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn0);
  HPW.h1tot(Proc, Bas, dt2, lfield, Wfn1, Wfn0);

  HPW.axpyz(Proc, Bas, cdt1, hWfn0, Wfn0, Wfn1);
  dPade.prod(Proc, Bas, time1, Field, HPW, Wfn1);

  for (int icyc = 1; icyc < crnic_maxcyc; icyc ++) {
    Field.get_value(time1, lfield);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
    HPW.xpyz(Proc, Bas, hWfn1, hWfn0, Wfn1);
    HPW.axpyz(Proc, Bas, cdt2, Wfn1, Wfn0, Wfn1);
    dPade.prod(Proc, Bas, time1, Field, HPW, Wfn1);
  }

//old  HPW.copy(Proc, Bas, Wfn, Wfn0);
//old  HPW.copy(Proc, Bas, Wfn, Wfn1);
//old
//old  Field.get_value(time0, lfield);
//old  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn0, Wfn1);
//old  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
//old  HPW.xpy(Proc, Bas, hWfn0, Wfn1);
//old  HPW.scal(Proc, Bas, cdt1, Wfn1);
//old
//old  dPade.prod(Proc, Bas, time1, Field, HPW, Wfn1); // dC0 = C1 - C0
//old  res_wfn = get_res(Proc, Bas, Wfn1);
//old  printf(" clcrnic::prop: %10d%20.10e\n", LZERO, res_wfn);
//old  HPW.xpy(Proc, Bas, Wfn0, Wfn1);                 // C1 = C0 + dC0
//old  HPW.copy(Proc, Bas, Wfn1, Wfn0);
//old
//old  for (int icyc = 1; icyc < crnic_maxcyc; icyc ++) {
//old    Field.get_value(time1, lfield);
//old    HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
//old    HPW.xmyz(Proc, Bas, hWfn1, hWfn0, Wfn1);
//old    HPW.scal(Proc, Bas, cdt2, Wfn1);
//old    dPade.prod(Proc, Bas, time1, Field, HPW, Wfn1); // dC(i-1) = C(i) - C(i-1)
//old    res_wfn = get_res(Proc, Bas, Wfn1);
//old    printf(" clcrnic::prop: %10d%20.10e\n", icyc, res_wfn);
//old    HPW.xpy(Proc, Bas, Wfn0, Wfn1);                 // C(i) = C(i-1) + dC(i-1)
//old    HPW.copy(Proc, Bas, Wfn1, Wfn0);
//old    HPW.copy(Proc, Bas, hWfn1, hWfn0);
//old  }

  HPW.copy(Proc, Bas, Wfn1, Wfn);
}
////////////////////////////////////////////////////////////////////////
double clcrnic::get_res(const clmpi& Proc, const clbas& Bas, const clwfn& Wfn)
{
  dcomplex dC;
  zdotc_omp_(&Wfn.size, &Wfn.wfn[0], &Wfn.wfn[0], &dC);
  return sqrt(real(dC));
}
////////////////////////////////////////////////////////////////////////
