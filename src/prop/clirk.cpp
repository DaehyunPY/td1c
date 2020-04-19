////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
clirk::clirk()
{
  std::cout << "clirk" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clirk::~clirk()
{
  std::cout << "~clirk" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clirk::clirk(const clmpi& MPIP, const clio& IO, const clbas& Bas,
	     const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void clirk::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		const clfield& Field, const clhprod& HPW)
{
  std::cout << "clirk::gen" << std::endl;
  clh2prop::gen_basic(MPIP, IO, Bas, Field, HPW);

  IO.read_info("rk_order", LTWO, irk_order);
  IO.read_info("rk_maxcyc", 10, irk_maxcyc);
  IO.read_info("rk_thresh", (double) 1.0E-10, irk_thresh);
  IO.read_info("rk_formula", "trapezoidal", irk_formula);

  Wfn0.gen(MPIP, IO, Bas);
  Wfn1.gen(MPIP, IO, Bas);
  hWfn0.gen(MPIP, IO, Bas);
  hWfn1.gen(MPIP, IO, Bas);
}
////////////////////////////////////////////////////////////////////////
void clirk::prop(const clmpi& Proc, const clbas& Bas, 
		 const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  if (irk_order == 1) {
    prop1(Proc, Bas, Field.time, Field.dtime, Field, HPW, Wfn);
  } else if (irk_order == 2 && irk_formula.compare("midpoint") == 0) {
    prop2m(Proc, Bas, Field.time, Field.dtime, Field, HPW, Wfn);
  } else if (irk_order == 2 && irk_formula.compare("trapezoidal") == 0) {
    prop2t(Proc, Bas, Field.time, Field.dtime, Field, HPW, Wfn);
  } else {
    std::cout << "clirk::gen. bad irk_order or irk_formula." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clirk::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field, 
		 double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  if (irk_order == 1) {
    prop1(Proc, Bas, time, dtime, Field, HPW, Wfn);
  } else if (irk_order == 2 && irk_formula.compare("midpoint") == 0) {
    prop2m(Proc, Bas, time, dtime, Field, HPW, Wfn);
  } else if (irk_order == 2 && irk_formula.compare("trapezoidal") == 0) {
    prop2t(Proc, Bas, time, dtime, Field, HPW, Wfn);
  } else {
    std::cout << "clirk::gen. bad irk_order or irk_formula." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clirk::prop1(const clmpi& Proc, const clbas& Bas, double time, double dtime, 
		  const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double lfield[9];
  double time1 = time + dtime;
  dcomplex dC;

  Field.get_value(time1, lfield);
  HPW.copy(Proc, Bas, Wfn, Wfn0);
  HPW.htot(Proc, Bas, dtime, lfield, Wfn0, hWfn1);
  HPW.xpyz(Proc, Bas, Wfn0, hWfn1, Wfn);

  for (int icyc = 1; icyc <= irk_maxcyc; icyc ++) {
    // update
    Field.get_value(time1, lfield);
    HPW.copy(Proc, Bas, Wfn, Wfn1);
    HPW.htot(Proc, Bas, dtime, lfield, Wfn1, hWfn1);
    HPW.xpyz(Proc, Bas, Wfn0, hWfn1, Wfn);

    // difference
    HPW.xmyz(Proc, Bas, Wfn, Wfn1, Wfn1);
    zdotc_omp_(&Wfn.size, &Wfn1.wfn[0], &Wfn1.wfn[0], &dC);
    printf("clirk::prop1: %10d%20.10f\n", icyc, sqrt(real(dC)));
    if (sqrt(real(dC)) < irk_thresh) break;
  }
}
////////////////////////////////////////////////////////////////////////
void clirk::prop2m(const clmpi& Proc, const clbas& Bas, double time, double dtime, 
		   const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double lfield[9];
  double timeh = time + dtime * HALF;
  dcomplex dC;

  // initial wavefunction
  HPW.copy(Proc, Bas, Wfn, Wfn0);

  // implicit euler to estimate wfn(1/2)
  prop1(Proc, Bas, time, HALF * dtime, Field, HPW, Wfn);

  Field.get_value(timeh, lfield);
  HPW.htot(Proc, Bas, dtime, lfield, Wfn, hWfn1);
  HPW.xpyz(Proc, Bas, Wfn0, hWfn1, Wfn);

  for (int icyc = 1; icyc <= irk_maxcyc; icyc ++) {
    // update
    Field.get_value(timeh, lfield);
    HPW.axpbyz(Proc, Bas, CHALF, Wfn0, CHALF, Wfn, Wfn1);
    HPW.htot(Proc, Bas, dtime, lfield, Wfn1, hWfn1);
    HPW.copy(Proc, Bas, Wfn, Wfn1);
    HPW.xpyz(Proc, Bas, Wfn0, hWfn1, Wfn);

    // difference
    HPW.xmyz(Proc, Bas, Wfn, Wfn1, Wfn1);
    zdotc_omp_(&Wfn.size, &Wfn1.wfn[0], &Wfn1.wfn[0], &dC);
    printf("clirk::prop2m: %10d%20.10f\n", icyc, sqrt(real(dC)));
    if (sqrt(real(dC)) < irk_thresh) break;
  }
}
////////////////////////////////////////////////////////////////////////
void clirk::prop2t(const clmpi& Proc, const clbas& Bas, double time, double dtime, 
		   const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double lfield[9];
  double time0 = time;
  double time1 = time + dtime;

  dcomplex dC;
  dcomplex cdt2 = dtime * CHALF;

  Field.get_value(time0, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn, hWfn1);
  HPW.axpyz(Proc, Bas, cdt2, hWfn1, Wfn, Wfn0);
  HPW.axpyz(Proc, Bas, cdt2, hWfn1, Wfn0, Wfn);

  for (int icyc = 1; icyc <= irk_maxcyc; icyc ++) {
    // update
    Field.get_value(time1, lfield);
    HPW.copy(Proc, Bas, Wfn, Wfn1);
    HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
    HPW.axpyz(Proc, Bas, cdt2, hWfn1, Wfn0, Wfn);

    // difference
    HPW.xmyz(Proc, Bas, Wfn, Wfn1, Wfn1);
    zdotc_omp_(&Wfn.size, &Wfn1.wfn[0], &Wfn1.wfn[0], &dC);
    printf("clirk::prop2t: %10d%20.10f\n", icyc, sqrt(real(dC)));
    if (sqrt(real(dC)) < irk_thresh) break;
  }
}
////////////////////////////////////////////////////////////////////////
