////////////////////////////////////////////////////////////////////////
// ETD-RK2
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
cletdrk2::cletdrk2()
{
  std::cout << "cletdrk2" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cletdrk2::~cletdrk2()
{
  std::cout << "~cletdrk2" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cletdrk2::cletdrk2(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		   const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void cletdrk2::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		   const clfield& Field, const clhprod& HPW)
{
  long icomp = Field.td_type;
  long isplit = clcontrol::split_type;

  Wfn0.gen(MPIP, IO, Bas);
  Wfn1.gen(MPIP, IO, Bas);
  hWfn0.gen(MPIP, IO, Bas);
  hWfn1.gen(MPIP, IO, Bas);
  tWfn.gen(MPIP, IO, Bas);
  dPade1.gen(MPIP, IO, Bas, Field.dtime, CTWO, icomp, isplit);
}
////////////////////////////////////////////////////////////////////////
void cletdrk2::prop(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double ttime;
  double lfield[9];
  double dt1 = Field.dtime;

  HPW.copy(Proc, Bas, Wfn, Wfn0);

  ttime = Field.time;
  Field.get_value(ttime, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  HPW.axpbyz(Proc, Bas, CTWO, Wfn0, dt1 * RUNIT, hWfn0, tWfn);
  if (clcontrol::split_type != 0) dPade1.prod(Proc, Bas, ttime, Field, HPW, tWfn);

  HPW.xmyz(Proc, Bas, tWfn, Wfn0, Wfn1);

  ttime = Field.time + dt1;
  Field.get_value(ttime, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
  HPW.xmyz(Proc, Bas, hWfn1, hWfn0, tWfn);
  HPW.scal(Proc, Bas, dt1 * CHALF, tWfn);
  if (clcontrol::split_type != 0) dPade1.prod(Proc, Bas, ttime, Field, HPW, tWfn);
  
  HPW.xpyz(Proc, Bas, Wfn1, tWfn, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cletdrk2::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		    double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  double ttime;
  double lfield[9];
  double dt1 = dtime;

  HPW.copy(Proc, Bas, Wfn, Wfn0);

  ttime = time;
  Field.get_value(ttime, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);
  HPW.axpbyz(Proc, Bas, CTWO, Wfn0, dt1 * RUNIT, hWfn0, tWfn);
  if (clcontrol::split_type != 0) dPade1.prod(Proc, Bas, ttime, Field, HPW, tWfn);

  HPW.xmyz(Proc, Bas, tWfn, Wfn0, Wfn1);

  ttime = time + dt1;
  Field.get_value(ttime, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1);
  HPW.xmyz(Proc, Bas, hWfn1, hWfn0, tWfn);
  HPW.scal(Proc, Bas, dt1 * CHALF, tWfn);
  if (clcontrol::split_type != 0) dPade1.prod(Proc, Bas, ttime, Field, HPW, tWfn);

  HPW.xpyz(Proc, Bas, Wfn1, tWfn, Wfn);
}
////////////////////////////////////////////////////////////////////////
