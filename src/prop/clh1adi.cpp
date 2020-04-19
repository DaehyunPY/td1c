////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
clh1adi::clh1adi()
{
  std::cout << "clh1adi" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clh1adi::~clh1adi()
{
  std::cout << "clh1adi" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clh1adi::clh1adi(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		 const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void clh1adi::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		  const clfield& Field, const clhprod& HPW)
{
  std::cout << "clh1adi::gen" << std::endl;
  clh1prop::gen_basic(MPIP, IO, Bas, Field, HPW);

  if (clcontrol::split_type == 1) {
    doext = 1;
  } else {
    doext = 0;
  }

  if (projfc) {
    std::cout << "clh1adi::gen: projfc is nyi" << std::endl;
    abort();
  }

  long size;
  size = (2 * Bas.GRad.nmax + 1) * (Bas.GRad.nrad - 1) * (Bas.GAng.lmax1 + 1);
  tadi1.resize(size);
  size = (3 * Bas.GRad.nmax + 1) * (Bas.GRad.nrad - 1) * (Bas.GAng.lmax1 + 1);
  tadi2.resize(size);
  size = (Bas.GRad.nrad - 1) * (Bas.GAng.lmax1 + 1);
  tpiv2.resize(size);

  double dt2 = HALF * Field.dtime;
  adi_gen_tadi_(&Field.td_type, &LONE, &LONE, &dt2, &tadi1[0], &tadi2[0], &tpiv2[0]);
}
////////////////////////////////////////////////////////////////////////
void clh1adi::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		  const clfield& Field, double dtime, const clhprod& HPW)
{
  std::cout << "clh1adi::gen_2 nyi." << std::endl;
  abort();
}
////////////////////////////////////////////////////////////////////////
void clh1adi::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   long STEP, clhprod& HPW, clwfn& Wfn)
{
  if (Field.lgauge) {
    propl(Proc, Bas, Field, STEP, HPW, Wfn);
  } else {
    propv(Proc, Bas, Field, STEP, HPW, Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void clh1adi::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  std::cout << "clh1adi::prop_2 nyi." << std::endl;
  abort();
}
////////////////////////////////////////////////////////////////////////
void clh1adi::propl(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		    long STEP, clhprod& HPW, clwfn& Wfn)
{
  double time, tnow;
  double dt2 = HALF * Field.dtime;
  double dt4 = dt2 * HALF;

  if (STEP == 1) {
    time = Field.time;
  } else {
    time = Field.time + Field.dtime * HALF;
  }

  if (doext == 1) {
    tnow = time;
    laser_lgauge(Proc, Bas, LZERO, tnow, dt4, Field, Wfn);
  }
  t_explicit(Proc, Bas, HPW, Wfn);
  t_implicit(Proc, Bas, HPW, Wfn);
  if (doext == 1) {
    tnow = time + dt2;
    laser_lgauge(Proc, Bas, LONE, tnow, dt4, Field, Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void clh1adi::propv(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		    long STEP, clhprod& HPW, clwfn& Wfn)
{
  double time, tnow;
  double dt2 = HALF * Field.dtime;
  double dt4 = dt2 * HALF;

  if (STEP == 1) {
    time = Field.time;
  } else {
    time = Field.time + Field.dtime * HALF;
  }

  if (doext == 1) {
    tnow = time;
    laser_vgauge(Proc, Bas, LZERO, tnow, dt4, Field, Wfn);
  }
  t_explicit(Proc, Bas, HPW, Wfn);
  t_implicit(Proc, Bas, HPW, Wfn);
  if (doext == 1) {
    tnow = time + dt2;
    laser_vgauge(Proc, Bas, LONE, tnow, dt4, Field, Wfn);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1adi::t_explicit(const clmpi& MPIP, const clbas& Bas, clhprod& HPW, clwfn& Wfn)
{
  adi_t_explicit_(&tadi1[0], &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clh1adi::t_implicit(const clmpi& MPIP, const clbas& Bas, clhprod& HPW, clwfn& Wfn)
{
  adi_t_implicit_(&tadi2[0], &tpiv2[0], &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clh1adi::laser_lgauge(const clmpi& MPIP, const clbas& Bas, long istag, 
			   double time, double dt, const clfield& Field, clwfn& Wfn)
{
  double lfield[9];
  Field.get_value(time, lfield);
  adi_laser_lgauge_(&Field.td_type, &istag, &dt, lfield, &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clh1adi::laser_vgauge(const clmpi& MPIP, const clbas& Bas, long istag, 
			   double time, double dt, const clfield& Field, clwfn& Wfn)
{
  double lfield[9];
  Field.get_value(time, lfield);
  adi_laser_vgauge_(&Field.td_type, &istag, &dt, lfield, &Wfn.wfn[0]);
}
////////////////////////////////////////////////////////////////////////
