////////////////////////////////////////////////////////////////////////
// Exponential Runge-Kuta with Pade approximation for phi functions
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
long cletdrb::max_dim_numer = 3;
long cletdrb::max_dim_denom = 3;
////////////////////////////////////////////////////////////////////////
cletdrb::cletdrb()
{
  std::cout << "cletdrb" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cletdrb::~cletdrb()
{
  std::cout << "~cletdrb" << std::endl;
}
////////////////////////////////////////////////////////////////////////
cletdrb::cletdrb(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		 const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void cletdrb::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		   const clfield& Field, const clhprod& HPW)
{
  etd_dt1 = Field.dtime;
  etd_dt2 = Field.dtime * HALF;

  IO.read_info("etd_rk_order", (long) 4, rk_order);
  IO.read_info("etd_dim_numer1", (long) 2, dim_numer1);
  IO.read_info("etd_dim_numer2", (long) 2, dim_numer2);
  IO.read_info("etd_dim_numer3", (long) 2, dim_numer3);
  IO.read_info("etd_dim_denom1", (long) 3, dim_denom1);
  IO.read_info("etd_dim_denom2", (long) 3, dim_denom2);
  IO.read_info("etd_dim_denom3", (long) 3, dim_denom3);

  IO.read_info("etd_pade_equiv", false, pade_equiv);
  IO.read_info("etd_rosenbrock", true, rosenbrock);
  IO.read_info("etd_time_of_h", (long) 0, time_of_h);
  IO.read_info("etd_cisplit", (long) -1, cisplit);
  if (Bas.ORMAS.nact == 0) cisplit = -1;
  if (cisplit != -1) {
    std::cout << "cletdrb::gen: etd_cisplit = -1 only." << std::endl;
    abort();
  }

  //  if (rk_order < 1 || rk_order > 4) {
  if (rk_order != 4) {
    std::cout << "cletdrb::gen bad rk_order." << std::endl;
    abort();
  }

  gen_phi1(MPIP, IO, Bas);
  gen_phid1(MPIP, IO, Bas);
  gen_phi2(MPIP, IO, Bas);
  gen_phi3(MPIP, IO, Bas);

  //DEBUG
  printf("cletdrb::phi1dt1:\n"); phi1dt1.print();
  printf("cletdrb::phi1dt2:\n"); phi1dt2.print();
  printf("cletdrb::phi2dt1:\n"); phi2dt1.print();
  printf("cletdrb::phi2dt2:\n"); phi2dt2.print();
  printf("cletdrb::phi3dt1:\n"); phi3dt1.print();
  printf("cletdrb::phi3dt2:\n"); phi3dt2.print();
  //  std::cout << "abort for debug in cletdrb::gen." << std::endl;
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
  dWfn.gen(MPIP, IO, Bas);

  //  long NA2 = pow(Bas.ORMAS.nact, 2);
  //  long NA4 = pow(Bas.ORMAS.nact, 4);
  long NA2 = std::max(LONE, Bas.ORMAS.nact * Bas.ORMAS.nact);
  long NA4 = NA2 * NA2;
  Den1.resize(NA2);
  Int1e.resize(NA2);
  Int2e.resize(NA4);
  hDiag.resize(Bas.ORMAS.lcic);
}
////////////////////////////////////////////////////////////////////////
void cletdrb::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  std::cout << "cletdrb::prop nyi." << std::endl;
  abort();
}
////////////////////////////////////////////////////////////////////////
void cletdrb::prop(const clmpi& Proc, const clbas& Bas,
		   const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  if (rk_order == 1) {
    prop1(Proc, Bas, Field, HPW, Wfn);
  } else if (rk_order == 2) {
    prop2(Proc, Bas, Field, HPW, Wfn);
  } else if (rk_order == 3) {
    prop3(Proc, Bas, Field, HPW, Wfn);
  } else if (rk_order == 4) {
    if (cisplit == -1) {
      prop4_orb(Proc, Bas, Field, HPW, Wfn);
    } else if (cisplit == 2) {
      prop4_orbci(Proc, Bas, Field, HPW, Wfn);
    } else {
      std::cout << "cletdrb::prop: bad cisplit." << std::endl;
      abort();
    }
  } else {
    std::cout << "cletdrb::prop bad rk_order." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void cletdrb::prop1(const clmpi& Proc, const clbas& Bas,
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
void cletdrb::prop2(const clmpi& Proc, const clbas& Bas,
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
void cletdrb::prop3(const clmpi& Proc, const clbas& Bas,
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
void cletdrb::prop4(const clmpi& Proc, const clbas& Bas,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  double lfield[9];
  double timeh;
  double time0 = Field.time;
  double time1 = Field.time + etd_dt2;
  double time2 = Field.time + etd_dt2;
  double time3 = Field.time + etd_dt1;
  if (time_of_h == 0) {
    timeh = Field.time;
  } else {
    timeh = Field.time + etd_dt2;
  }

  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  // stage 0
  HPW.copy(Proc, Bas, Wfn, Wfn0); // u0
  Field.get_value(timeh, lfield); // Vext0
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn0, tWfn);
  Field.get_value(time0, lfield); // Vext0
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);  // hWfn0 = W0
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfn0, tWfn); // tWfn  = -i*h0*u0 + W0

  prod1(Proc, Bas, HPW, timeh, Field, phi1dt2, phi1rk, tWfn, Wfn1);
  HPW.axpyz(Proc, Bas, CHALF, Wfn1, Wfn0, Wfn1); // phi1(c1)*(-i*h0*u0 + W0)/2 ==> u1 (c1 = 1/2)
  prod1(Proc, Bas, HPW, timeh, Field, phi1dt1, phi1rk, tWfn, Wfn3);
  HPW.xpyz(Proc, Bas, Wfn0, Wfn3, Wfn3); // phi1(c3)*(-i*h0*u0 + W0) ==> u3 (c3 = 1)
  if (rosenbrock) {
    HPW.clear(Proc, Bas, dWfn);
    Field.get_der(timeh, lfield); // dot(Vext0)
    HPW.v1ext(Proc, Bas, tunit*etd_dt1, lfield, Wfn0, dWfn); // -i*dt*dot(Vext0)*u0
    prod1(Proc, Bas, HPW, timeh, Field, phi2dt2, phi2rk, dWfn, tWfn);
    HPW.axpyz(Proc, Bas, RUNIT/CFOUR, tWfn, Wfn1, Wfn1); // phi2(c1)*du0*c1^2 +=> u1 (c1=0.5)
    prod1(Proc, Bas, HPW, timeh, Field, phi2dt1, phi2rk, dWfn, tWfn);
    HPW.axpyz(Proc, Bas, RUNIT, tWfn, Wfn3, Wfn3); // phi1(c3)*du0*c3^2 +=> u3 (c3=1.0)
  }
  HPW.copy(Proc, Bas, Wfn1, Wfn2); // phi1(c2)*(-i*h0*u0 + W0)*c2 + phi2(c1)*du0*c1^2 ==> u2 (c2=0.5)
  HPW.copy(Proc, Bas, Wfn3, Wfn);  // phi1(1.)*(-i*h0*u0 + W0)*1. + phi2(1.)*du0*c1^2 ==> u(t+dt)

  // stage 1
  Field.get_value(time1, lfield); // Vext1
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1); // W1
  if (rosenbrock) {
    Field.get_value(time1, timeh, ONE, -ONE, lfield); // dVext1 = Vext1 - Vext0
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn1, hWfn1); // W1 - i*dVext1*u1
    HPW.axpyz(Proc, Bas, -CHALF, dWfn, hWfn1, hWfn1); // W1 - i*dVext1*u1 + i*c1*dt*dot(Vext0)*u0
  }
  HPW.xmyz(Proc, Bas, hWfn1, hWfn0, hWfn1); // W1 - W0 - i*dVext1*u1 + i*c1*dt*dot(Vext0)*u0
  HPW.copy(Proc, Bas, hWfn1, tWfn);
  prod2(Proc, Bas, HPW, timeh, Field, phi2dt2, phi2rk, tWfn, Wfn2); // a21*(W1 - W0 - i*dVext1*u1 + i*c1*dt*dot(Vext0)*u0) +=> u2

  // stage 2
  Field.get_value(time2, lfield); // Vext2
  HPW.htot(Proc, Bas, ONE, lfield, Wfn2, hWfn2); // W2
  if (rosenbrock) {
    Field.get_value(time2, timeh, ONE, -ONE, lfield); // dVext2 = Vext2 - Vext0
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn2, hWfn2); // W2 - i*dVext2*u2
    HPW.axpyz(Proc, Bas, -CHALF, dWfn, hWfn2, hWfn2); // W2 - i*dVext2*u2 + i*c2*dt*dot(Vext0)*u0
  }
  HPW.xmyz(Proc, Bas, hWfn2, hWfn0, tWfn); 
  HPW.copy(Proc, Bas, tWfn, hWfn2);        // W2 - W0 - i*dVext2*u2 + i*c2*dt*dot(Vext0)*u0
  HPW.scal(Proc, Bas, CTWO, tWfn);         //(W2 - W0 - i*dVext2*u2 + i*c2*dt*dot(Vext0)*u0)*2
  prod2(Proc, Bas, HPW, timeh, Field, phi2dt1, phi2rk, tWfn, Wfn3); // a32*(W2 - W0 - i*dVext2*u2 + i*c2*dt*dot(Vext0)*u0) +=> u3

  // stage 3
  Field.get_value(time3, lfield); // Vext3
  HPW.htot(Proc, Bas, ONE, lfield, Wfn3, hWfn3); // W3
  if (rosenbrock) {
    Field.get_value(time3, timeh, ONE, -ONE, lfield); // dVext3 = Vext3 - Vext0
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn3, hWfn3); // W3 - i*dVext3*u3
    HPW.axpyz(Proc, Bas, -RUNIT, dWfn, hWfn3, hWfn3); // W3 - i*dVext3*u3 + i*c3*dt*dot(Vext0)*u0
  }
  HPW.xmyz(Proc, Bas, hWfn3, hWfn0, hWfn3); // W3 - W0 - i*dVext3*u3 + i*c3*dt*dot(Vext0)*u0

  // summation
  HPW.axpbyz(Proc, Bas, CTWO, hWfn1, CTWO, hWfn2, tWfn); // 2*W1 + 2*W2 - W3
  HPW.axpy(Proc, Bas, -RUNIT, hWfn3, tWfn);              // 
  prod2(Proc, Bas, HPW, timeh, Field, phi2dt1, phi2rk, tWfn, Wfn);

  HPW.axpbyz(Proc, Bas, -CFOUR, hWfn1, -CFOUR, hWfn2, tWfn); // -4*W1 - 4*W2 + 4*W3
  HPW.axpy(Proc, Bas, +CFOUR, hWfn3, tWfn);                  //
  prod2(Proc, Bas, HPW, timeh, Field, phi3dt1, phi3rk, tWfn, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cletdrb::prop4_orb(const clmpi& Proc, const clbas& Bas,
			const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  double lfield[9];
  double dEref;
  double timeh;
  double time0 = Field.time;
  double time1 = Field.time + etd_dt2;
  double time2 = Field.time + etd_dt2;
  double time3 = Field.time + etd_dt1;
  //  long NA2 = pow(Bas.ORMAS.nact, 2);
  //  long NA4 = pow(Bas.ORMAS.nact, 4);
  long NA2 = std::max(LONE, Bas.ORMAS.nact * Bas.ORMAS.nact);
  long NA4 = NA2 * NA2;

  if (time_of_h == 0) {
    timeh = Field.time;
  } else {
    timeh = Field.time + etd_dt2;
  }

  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  Eref = ZERO;
  zclear_omp_(&NA2, &Int1e[0]);
  zclear_omp_(&NA4, &Int2e[0]);
  zclear_omp_(&Bas.ORMAS.lcic, &hDiag[0]);

  // nodes
  dcomplex c1 = CHALF;
  dcomplex c3 = RUNIT;

  // initialization
  HPW.copy(Proc, Bas, Wfn, Wfn0); // u0
  HPW.clear(Proc, Bas, Wfn1);
  HPW.clear(Proc, Bas, Wfn3);

  // stage 0
  //20161127  HPW.clear(Proc, Bas, hWfn0);
  HPW.clear(Proc, Bas, tWfn);
  Field.get_value(timeh, lfield);
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn0, tWfn);
  Field.get_value(time0, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);  // hWfn0 = W0
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfn0, tWfn); // tWfn  = -i*h0*u0 + W0
  phi1dt2.prod_numer(Proc, Bas, timeh, c1, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn1);
  phi1dt1.prod_numer(Proc, Bas, timeh, c3, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn3);
  if (rosenbrock) {
    HPW.clear(Proc, Bas, dWfn);
    Field.get_der(timeh, lfield); // dot(Vext0)
    HPW.v1ext(Proc, Bas, RUNIT, lfield, Wfn0, dWfn); // dot(Vext0)*u0
    HPW.scal(Proc, Bas, tunit*etd_dt1, dWfn);        // -i*dt*dot(Vext0)*u0
    phid1dt2.prod_numer(Proc, Bas, timeh, c1*c1, Field, HPW, Eref, hDiag, Int1e, Int2e, dWfn, Wfn1);
    phid1dt1.prod_numer(Proc, Bas, timeh, c3*c3, Field, HPW, Eref, hDiag, Int1e, Int2e, dWfn, Wfn3);
  }

  // stage 1
  HPW.copy(Proc, Bas, Wfn1, tWfn);
  HPW.copy(Proc, Bas, Wfn0, Wfn1);
  phi1dt2.prod_denom(Proc, Bas, timeh, RUNIT, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn1);
  Field.get_value(time1, lfield); // Vext1
  HPW.htot(Proc, Bas, ONE, lfield, Wfn1, hWfn1); // W1
  HPW.axpy(Proc, Bas, -RUNIT, hWfn0, hWfn1);     // W1 - W0
  if (rosenbrock) {
    Field.get_value(time1, timeh, ONE, -ONE, lfield); // dVext1 = Vext1 - Vext0
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn1, hWfn1); // W1 - W0 - i*dVext1*u1 
    HPW.axpy(Proc, Bas, -CHALF, dWfn, hWfn1);         // W1 - W0 - i*dVext1*u1 + i*c1*dt*dot(Vext0)*u0
  }

  // stage 2
  HPW.copy(Proc, Bas, Wfn1, Wfn2);
  HPW.copy(Proc, Bas, hWfn1, tWfn);
  phi2dt2.prod2(Proc, Bas, timeh, RUNIT, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn2);
  Field.get_value(time2, lfield); // Vext2
  HPW.htot(Proc, Bas, ONE, lfield, Wfn2, hWfn2); // W2
  HPW.axpy(Proc, Bas, -RUNIT, hWfn0, hWfn2);     // W2 - W0
  if (rosenbrock) {
    Field.get_value(time2, timeh, ONE, -ONE, lfield); // dVext2 = Vext2 - Vext0
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn2, hWfn2); // W2 - W0 - i*dVext2*u2
    HPW.axpy(Proc, Bas, -CHALF, dWfn, hWfn2);         // W2 - W0 - i*dVext2*u2 + i*c2*dt*dot(Vext0)*u0
  }

  // stage 3
  HPW.copy(Proc, Bas, Wfn3, tWfn);
  HPW.copy(Proc, Bas, Wfn0, Wfn3);
  phi1dt1.prod_denom(Proc, Bas, timeh, RUNIT, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn3);
  HPW.copy(Proc, Bas, hWfn2, tWfn);
  phi2dt1.prod2(Proc, Bas, timeh, CTWO, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn3);
  Field.get_value(time3, lfield); // Vext3
  HPW.htot(Proc, Bas, ONE, lfield, Wfn3, hWfn3); // W3
  HPW.axpy(Proc, Bas, -RUNIT, hWfn0, hWfn3);     // W3 - W0
  if (rosenbrock) {
    Field.get_value(time3, timeh, ONE, -ONE, lfield); // dVext3 = Vext3 - Vext0
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn3, hWfn3); // W3 - W0 - i*dVext3*u3
    HPW.axpy(Proc, Bas, -RUNIT, dWfn, hWfn3);         // W3 - W0 - i*dVext3*u3 + i*c3*dt*dot(Vext0)*u0
  }

  // summation
  HPW.copy(Proc, Bas, Wfn3, Wfn);
  HPW.axpbyz(Proc, Bas, CTWO, hWfn1, -RUNIT, hWfn3, tWfn); // 2*W1 - W3 (2*W2 is included in hWfn3)
  phi2dt1.prod2(Proc, Bas, timeh, RUNIT, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn);
  HPW.xpyz(Proc, Bas, hWfn1, hWfn2, tWfn); 
  HPW.xmyz(Proc, Bas, tWfn, hWfn3, tWfn);  // W1 + W2 - W3 
  phi3dt1.prod2(Proc, Bas, timeh, -CFOUR, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn);
  //DEBUG
  //  std::cout << "ABORT for debug @ clerdrb::prop4_orb." << std::endl;
  //  abort();
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrb::prop4_orbci(const clmpi& Proc, const clbas& Bas,
			  const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  dcomplex tunit;
  double lfield[9];
  double dEref;
  double timeh;
  double time0 = Field.time;
  double time1 = Field.time + etd_dt2;
  double time2 = Field.time + etd_dt2;
  double time3 = Field.time + etd_dt1;
  //  long NA2 = pow(Bas.ORMAS.nact, 2);
  //  long NA4 = pow(Bas.ORMAS.nact, 4);
  long NA2 = std::max(LONE, Bas.ORMAS.nact * Bas.ORMAS.nact);
  long NA4 = NA2 * NA2;

  double Erefx;
  std::vector<dcomplex> Int1x(NA2), Int1y(NA2);
  std::vector<dcomplex> Int2x(NA4), Int2y(NA4);

  if (time_of_h == 0) {
    timeh = Field.time;
  } else {
    timeh = Field.time + etd_dt2;
  }

  if (clcontrol::icomp == 0) {
    tunit = -RUNIT;
  } else {
    tunit = -IUNIT;
  }

  Eref = ZERO;
  zclear_omp_(&NA2, &Int1e[0]);
  zclear_omp_(&NA4, &Int2e[0]);
  //  zclear_omp_(&Bas.ORMAS.lcic, &hDiag[0]);

  // nodes
  dcomplex c1 = CHALF;
  dcomplex c3 = RUNIT;

  // initialization
  HPW.copy(Proc, Bas, Wfn, Wfn0); // u0
  HPW.clear(Proc, Bas, Wfn1);
  HPW.clear(Proc, Bas, Wfn3);
  // stage 0
  HPW.clear(Proc, Bas, hWfn0);
  Field.get_value(timeh, lfield);
  HPW.h1tot(Proc, Bas, ONE, lfield, Wfn0, tWfn);
  Field.get_value(time0, lfield);
  HPW.htot(Proc, Bas, ONE, lfield, Wfn0, hWfn0);  // hWfn0 = W0
  HPW.axpyz(Proc, Bas, tunit, tWfn, hWfn0, tWfn); // tWfn  = -i*h0*u0 + W0
  if (cisplit != -1) {
    Eref = HPW.ene_act;
    //old    zcopy_omp_(&NA2, &HPW.den1[0], &Den1[0]);
    //old    zcopy_omp_(&NA2, &HPW.int1e[0], &Int1e[0]);
    //old    zcopy_omp_(&NA4, &HPW.int2e[0], &Int2e[0]);
    HPW.getden1(Den1);
    HPW.getint1(Int1e);
    HPW.getint2(Int2e);
    //    ormas_hdiag_old_(&Int1e[0], &Int2e[0], &hDiag[0]);
  }
  phi1dt2.prod_numer(Proc, Bas, timeh, c1, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn1);
  phi1dt1.prod_numer(Proc, Bas, timeh, c3, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn3);
  if (rosenbrock) {
    HPW.clear(Proc, Bas, dWfn);
    Field.get_der(timeh, lfield); // dot(Vext0)
    HPW.v1ext(Proc, Bas, RUNIT, lfield, Wfn0, dWfn); // dot(Vext0)*u0
    if (cisplit != -1 && clcontrol::oorot_type != 1) {
      //old      hprod_v1tot_mkint1_(&Wfn0.wfn[0], &dWfn.wfn[0], &HPW.int1e[0]);
      //old      dEref = hprod_ene_act1_(&HPW.int1e[0], &Den1[0]);
      //old      if (cisplit == 2) ormas_hcic1p_    (&dEref, &HPW.int1e[0], &Wfn0.wfn[Wfn0.size1], &dWfn.wfn[dWfn.size1]);
      hprod_v1tot_mkint1_(&Wfn0.wfn[0], &dWfn.wfn[0], &Int1e[0]);
      dEref = hprod_ene_act1_(&Int1e[0], &Den1[0]);
      if (cisplit == 2) ormas_hcic1p_    (&dEref, &Int1e[0], &Wfn0.wfn[Wfn0.size1], &dWfn.wfn[dWfn.size1]);
    }
    HPW.scal(Proc, Bas, tunit*etd_dt1, dWfn); // -i*dt*dot(Vext0)*u0
    phid1dt2.prod_numer(Proc, Bas, timeh, c1*c1, Field, HPW, Eref, hDiag, Int1e, Int2e, dWfn, Wfn1);
    phid1dt1.prod_numer(Proc, Bas, timeh, c3*c3, Field, HPW, Eref, hDiag, Int1e, Int2e, dWfn, Wfn3);
  }

  // stage 1
  HPW.clear(Proc, Bas, hWfn1);
  HPW.copy(Proc, Bas, Wfn1, tWfn);
  HPW.copy(Proc, Bas, Wfn0, Wfn1);
  phi1dt2.prod_denom(Proc, Bas, timeh, RUNIT, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn1);
  Field.get_value(time1, lfield); // Vext1
  HPW.htoto(Proc, Bas, ONE, lfield, Wfn1, hWfn1); // W1o
  HPW.axpyo(Proc, Bas, -RUNIT, hWfn0, hWfn1);     // W1o - W0o
  if (cisplit != -1) {
    Erefx = HPW.ene_act - Eref;
    //old    zxmyz_omp_(&NA2, &HPW.int1e[0], &Int1e[0], &Int1x[0]);
    //old    zxmyz_omp_(&NA4, &HPW.int2e[0], &Int2e[0], &Int2x[0]);
    HPW.getint1(Int1y);
    HPW.getint2(Int2y);
    zxmyz_omp_(&NA2, &Int1y[0], &Int1e[0], &Int1x[0]);
    zxmyz_omp_(&NA4, &Int2y[0], &Int2e[0], &Int2x[0]);

    ormas_hcic_(&Int1x[0], &Int2x[0], &Wfn1.wfn[Wfn1.size1], &hWfn1.wfn[hWfn1.size1], &Erefx);
    HPW.scalc(Proc, Bas, tunit, hWfn1); // -i*(H1 - H0)*C1
  }
  if (rosenbrock) {
    Field.get_value(time1, timeh, ONE, -ONE, lfield); // dVext1 = Vext1 - Vext0
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn1, hWfn1); // W1 - W0 - i*dVext1*u1
    HPW.axpy(Proc, Bas, -CHALF, dWfn, hWfn1);         // W1 - W0 - i*dVext1*u1 + i*c1*dt*dot(Vext0)*u0
  }

  // stage 2
  HPW.clear(Proc, Bas, hWfn2);
  HPW.copy(Proc, Bas, Wfn1, Wfn2);
  HPW.copy(Proc, Bas, hWfn1, tWfn);
  phi2dt2.prod2(Proc, Bas, timeh, RUNIT, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn2);
  Field.get_value(time2, lfield); // Vext2
  HPW.htoto(Proc, Bas, ONE, lfield, Wfn2, hWfn2); // W2o
  HPW.axpyo(Proc, Bas, -RUNIT, hWfn0, hWfn2);     // W2o - W0o
  if (cisplit != -1) {
    Erefx = HPW.ene_act - Eref;
    //old    zxmyz_omp_(&NA2, &HPW.int1e[0], &Int1e[0], &Int1x[0]);
    //old    zxmyz_omp_(&NA4, &HPW.int2e[0], &Int2e[0], &Int2x[0]);
    HPW.getint1(Int1y);
    HPW.getint2(Int2y);
    zxmyz_omp_(&NA2, &Int1y[0], &Int1e[0], &Int1x[0]);
    zxmyz_omp_(&NA4, &Int2y[0], &Int2e[0], &Int2x[0]);
    ormas_hcic_(&Int1x[0], &Int2x[0], &Wfn2.wfn[Wfn2.size1], &hWfn2.wfn[hWfn2.size1], &Erefx);
    HPW.scalc(Proc, Bas, tunit, hWfn2); // -i*(H2 - H0)*C2
  }
  if (rosenbrock) {
    Field.get_value(time2, timeh, ONE, -ONE, lfield); // dVext2 = Vext2 - Vext0
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn2, hWfn2); // W2 - W0 - i*dVext2*u2
    HPW.axpy(Proc, Bas, -CHALF, dWfn, hWfn2);         // W2 - W0 - i*dVext2*u2 + i*c2*dt*dot(Vext0)*u0
  }

  // stage 3
  HPW.clear(Proc, Bas, hWfn3);
  HPW.copy(Proc, Bas, Wfn3, tWfn);
  HPW.copy(Proc, Bas, Wfn0, Wfn3);
  phi1dt1.prod_denom(Proc, Bas, timeh, RUNIT, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn3);
  HPW.copy(Proc, Bas, hWfn2, tWfn);
  phi2dt1.prod2(Proc, Bas, timeh, CTWO, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn3);
  Field.get_value(time3, lfield); // Vext3
  HPW.htoto(Proc, Bas, ONE, lfield, Wfn3, hWfn3); // W3o
  HPW.axpyo(Proc, Bas, -RUNIT, hWfn0, hWfn3);     // W3o - W0o
  if (cisplit != -1) {
    Erefx = HPW.ene_act - Eref;
    //old    zxmyz_omp_(&NA2, &HPW.int1e[0], &Int1e[0], &Int1x[0]);
    //old    zxmyz_omp_(&NA4, &HPW.int2e[0], &Int2e[0], &Int2x[0]);
    HPW.getint1(Int1y);
    HPW.getint2(Int2y);
    zxmyz_omp_(&NA2, &Int1y[0], &Int1e[0], &Int1x[0]);
    zxmyz_omp_(&NA4, &Int2y[0], &Int2e[0], &Int2x[0]);
    ormas_hcic_(&Int1x[0], &Int2x[0], &Wfn3.wfn[Wfn3.size1], &hWfn3.wfn[hWfn3.size1], &Erefx);
    HPW.scalc(Proc, Bas, tunit, hWfn3); // -i*(H3 - H0)*C3
  }
  if (rosenbrock) {
    Field.get_value(time3, timeh, ONE, -ONE, lfield); // dVext3 = Vext3 - Vext0
    HPW.v1ext(Proc, Bas, tunit, lfield, Wfn3, hWfn3); // W3 - W0 - i*dVext3*u3
    HPW.axpy(Proc, Bas, -RUNIT, dWfn, hWfn3);         // W3 - W0 - i*dVext3*u3 + i*c3*dt*dot(Vext0)*u0
  }

  // summation
  HPW.copy(Proc, Bas, Wfn3, Wfn);
  HPW.axpbyz(Proc, Bas, CTWO, hWfn1, -RUNIT, hWfn3, tWfn); // 2*W1 - W3 (2*W2 is included in hWfn3)
  phi2dt1.prod2(Proc, Bas, timeh, RUNIT, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn);
  HPW.xpyz(Proc, Bas, hWfn1, hWfn2, tWfn); 
  HPW.xmyz(Proc, Bas, tWfn, hWfn3, tWfn);  // W1 + W2 - W3 
  phi3dt1.prod2(Proc, Bas, timeh, -CFOUR, Field, HPW, Eref, hDiag, Int1e, Int2e, tWfn, Wfn);
}
////////////////////////////////////////////////////////////////////////
void cletdrb::prod1(const clmpi& Proc, const clbas& Bas, clhprod& HPW, 
		    double time, const clfield& Field, const clh1rat& FacOrb, 
		    dcomplex FacCI, const clwfn& WIn, clwfn& WOut)
{
  HPW.copy(Proc, Bas, WIn, WOut);
  HPW.scalc(Proc, Bas, FacCI, WOut);
  FacOrb.prod(Proc, Bas, time, Field, HPW, WOut);
}
////////////////////////////////////////////////////////////////////////
void cletdrb::prod2(const clmpi& Proc, const clbas& Bas, clhprod& HPW, 
		     double time, const clfield& Field, const clh1rat& FacOrb, 
		     dcomplex FacCI, clwfn& WIn, clwfn& WOut)
{
  HPW.scalc(Proc, Bas, FacCI, WIn);
  FacOrb.prod(Proc, Bas, time, Field, HPW, WIn);
  HPW.xpyz(Proc, Bas, WIn, WOut, WOut);
}
//////////////////////////////////////////////////////////////////////////
//void cletdrb::prod_numer(const clmpi& Proc, const clbas& Bas, clhprod& HPW, 
//			 double time, const clfield& Field, dcomplex zfac, 
//			 const clh1rat& Rat, const clwfn& WIn, clwfn& WOut)
//{
//  Rat.prod_numer(Proc, Bas, time, Field, HPW, zfac, WIn, WOut);
//}
//////////////////////////////////////////////////////////////////////////
//void cletdrb::prod_denom(const clmpi& Proc, const clbas& Bas, clhprod& HPW, 
//			 double time, const clfield& Field, dcomplex zfac, 
//			 const clh1rat& Rat, clwfn& WIn, clwfn& WOut)
//{
//  Rat.prod_denom(Proc, Bas, time, Field, HPW, zfac, WIn, WOut);
//}
//////////////////////////////////////////////////////////////////////////
void cletdrb::gen_phi1(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  phi1rk = etd_dt1;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi1_n(max_dim_numer + 1);
  std::vector<dcomplex> phi1_d0(max_dim_denom);
  std::vector<dcomplex> phi1_d1(max_dim_denom);

  if (dim_numer1 == 1 && dim_denom1 == 1) {
    //               x + 6  
    // phi1(x) = [- -------]
    //              2 x - 6 
    // x0 = 3
    dfac = -1.0;
    phi1_n[0] = 6.0/2.0;
    phi1_n[1] = 1.0/2.0;
    roots[0] = 3.0;
  } else if (dim_numer1 == 2 && dim_denom1 == 2) {
    //               2              
    //              x  + 6 x + 60    
    // phi1(x) =  [----------------] 
    //                2              
    //             3 x  - 24 x + 60 
    dfac = 1.0;
    phi1_n[0] = 60.0/3.0;
    phi1_n[1] = 6.0/3.0;
    phi1_n[2] = 1.0/3.0;
    roots[0] =  4.0 - 2.0 * IUNIT;
    roots[1] =  4.0 + 2.0 * IUNIT;
  } else if (dim_numer1 == 2 && dim_denom1 == 3) {
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
  } else if (dim_numer1 == 3 && dim_denom1 == 3) {
    //                3       2                
    //               x  + 20 x  + 60 x + 840   
    // phi1(x) = [- --------------------------]
    //                 3       2               
    //              4 x  - 60 x  + 360 x - 840 
    dfac = -1.0;
    phi1_n[0] = 840.0/4.0;
    phi1_n[1] = 60.0/4.0;
    phi1_n[2] = 20.0/4.0;
    phi1_n[3] = 1.0/4.0;
    roots[0] = +5.648485971016884655 +0.000000000000000000 * IUNIT;
    roots[1] = +4.675757014491553676 +3.913489560603714335 * IUNIT;
    roots[2] = +4.675757014491553676 -3.913489560603714335 * IUNIT;
  } else {
    std::cout << "cletdrb:gen_exp. bad dim_numer1 or dim_denom1." << std::endl;
    abort();
  }

  phi1_d1[0] = +dfac;
  phi1_d0[0] = -dfac * roots[0];
  for (long idim = 1; idim < dim_denom1; idim ++) {
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
  for (long inumer = 0; inumer <= dim_numer1; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (long idenom = 0; idenom < dim_denom1; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi1dt1.gen(MPIP, IO, Bas, etd_dt1, cisplit, dim_numer1, dim_denom1, cn, cd0, cd1);

  // ##### dt * phi1(-i*h*dt/2) #####
  cn = phi1_n;
  cd0 = phi1_d0;
  cd1 = phi1_d1;
  tfac = tunit * etd_dt2; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (long inumer = 0; inumer <= dim_numer1; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (long idenom = 0; idenom < dim_denom1; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi1dt2.gen(MPIP, IO, Bas, etd_dt2, cisplit, dim_numer1, dim_denom1, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# phi1: real\n");
  // for (long ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi1_n, phi1_d0, phi1_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# phi1: imag\n");
  // for (long ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi1_n, phi1_d0, phi1_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrb::gen_phid1(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  phid1rk = etd_dt1 / TWO;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phid1_n(max_dim_numer + 1);
  std::vector<dcomplex> phid1_d0(max_dim_denom);
  std::vector<dcomplex> phid1_d1(max_dim_denom);

  if (dim_numer1 == 2 && dim_denom1 == 3) {
    //                   2		    
    //                  x  - 10 x + 60	    
    // phid1(x) = - -----------------------
    //               3       2		    
    //              x  - 12 x  + 60 x - 120
    dfac = -1.0;
    phid1_n[0] = 60.0;
    phid1_n[1] = -10.0;
    phid1_n[2] = 1.0;
    roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
    roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
    roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
  } else {
    if (rosenbrock) {
      std::cout << "cletdrb:gen_exp. bad dim_numer1 or dim_denom1." << std::endl;
      abort();
    } else {
      std::cout << "cletdrb:gen_exp. WARNING: bad dim_numer1 or dim_denom1." << std::endl;
      // taken from gen_phi3
      dfac = -1.0;
      phid1_n[0] = 10080.0/120.0;
      phid1_n[1] = -840.0/120.0;
      phid1_n[2] = 84.0/120.0;
      phid1_n[3] = 1.0/120.0;
      roots[0] = +7.653383973560514697 +0.000000000000000000 * IUNIT;
      roots[1] = +6.673308013219741319 +4.617378364686957504 * IUNIT;
      roots[2] = +6.673308013219741319 -4.617378364686957504 * IUNIT;
    }
  }

  phid1_d1[0] = +dfac;
  phid1_d0[0] = -dfac * roots[0];
  for (long idim = 1; idim < dim_denom1; idim ++) {
    phid1_d1[idim] = +RUNIT;
    phid1_d0[idim] = -roots[idim];
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

  // ##### dt * phid1(-i*h*dt) #####
  cn = phid1_n;
  cd0 = phid1_d0;
  cd1 = phid1_d1;
  tfac = tunit * etd_dt1; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (long inumer = 0; inumer <= dim_numer1; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (long idenom = 0; idenom < dim_denom1; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phid1dt1.gen(MPIP, IO, Bas, etd_dt1, cisplit, dim_numer1, dim_denom1, cn, cd0, cd1);

  // ##### dt * phid1(-i*h*dt/2) #####
  cn = phid1_n;
  cd0 = phid1_d0;
  cd1 = phid1_d1;
  tfac = tunit * etd_dt2; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (long inumer = 0; inumer <= dim_numer1; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (long idenom = 0; idenom < dim_denom1; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phid1dt2.gen(MPIP, IO, Bas, etd_dt2, cisplit, dim_numer1, dim_denom1, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# phid1: real\n");
  // for (long ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phid1_n, phid1_d0, phid1_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# phid1: imag\n");
  // for (long ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phid1_n, phid1_d0, phid1_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrb::gen_phi2(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  phi2rk = etd_dt1 / TWO;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi2_n(max_dim_numer + 1);
  std::vector<dcomplex> phi2_d0(max_dim_denom);
  std::vector<dcomplex> phi2_d1(max_dim_denom);

  if (dim_numer2 == 1 && dim_denom2 == 1) {
    //               x + 12  
    // phi2(x) = [- --------]
    //              6 x - 24 
    dfac = -1.0;
    phi2_n[0] = 12.0/6.0;
    phi2_n[1] = 1.0/6.0;
    roots[0] = 4.0;
  } else if (dim_numer2 == 2 && dim_denom2 == 2) {
    //                  2	       
    //                 x  + 180           
    // phi2(x) = [-------------------] 
    //                2	           
    //            12 x  - 120 x + 360 
    dfac = 1.0;
    phi2_n[0] = 180.0/12.0;
    phi2_n[1] = 0.0;
    phi2_n[2] = 1.0/12.0;
    roots[0] = +5.000000000000000000 +2.236067977499789805 * IUNIT;
    roots[1] = +5.000000000000000000 -2.236067977499789805 * IUNIT;
  } else if (dim_numer2 == 2 && dim_denom2 == 3) {
    //                    2		     
    //                 5 x  - 40 x + 420     	    
    // phi2(x) = - --------------------------
    //                3       2		     	    
    //             4 x  - 60 x  + 360 x - 840
    dfac = -1.0;
    phi2_n[0] = 420.0/4.0;
    phi2_n[1] = -40.0/4.0;
    phi2_n[2] = 5.0/4.0;
    roots[0] = +5.648485971016884655 +0.000000000000000000 * IUNIT;
    roots[1] = +4.675757014491553676 +3.913489560603714335 * IUNIT;
    roots[2] = +4.675757014491553676 -3.913489560603714335 * IUNIT;
  } else if (dim_numer2 == 2 && dim_denom2 == 3 && pade_equiv) {
    //                  2		    
    //                 x  - 10 x + 60	    
    // phi2(x) = - -----------------------
    //              3       2		    
    //             x  - 12 x  + 60 x - 120
    dfac = -1.0;
    phi2_n[0] = 60.0;
    phi2_n[1] = -10.0;
    phi2_n[2] = 1.0;
    roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
    roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
    roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
  } else if (dim_numer2 == 3 && dim_denom2 == 3) {
    //                 3       2		   
    //                x  + 40 x  - 140 x + 3360	   
    // phi2(x) = [- ------------------------------]
    //                  3        2		   
    //              20 x  - 360 x  + 2520 x - 6720 
    dfac = -1.0;
    phi2_n[0] = 3360.0/20.0;
    phi2_n[1] = -140.0/20.0;
    phi2_n[2] = 40.0/20.0;
    phi2_n[3] = 1.0/20.0;
    roots[0] = +6.651316808950019421 +0.000000000000000000 * IUNIT;
    roots[1] = +5.674341595524996507 +4.279971984629763249 * IUNIT;
    roots[2] = +5.674341595524996507 -4.279971984629763249 * IUNIT;
  } else {
    std::cout << "cletdrb:gen_exp. bad dim_numer2 or dim_denom2." << std::endl;
    abort();
  }

  phi2_d1[0] = +dfac;
  phi2_d0[0] = -dfac * roots[0];
  for (long idim = 1; idim < dim_denom2; idim ++) {
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
  for (long inumer = 0; inumer <= dim_numer2; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (long idenom = 0; idenom < dim_denom2; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi2dt1.gen(MPIP, IO, Bas, etd_dt1, cisplit, dim_numer2, dim_denom2, cn, cd0, cd1);

  // ##### dt * phi2(-i*h*dt/2) #####
  cn = phi2_n;
  cd0 = phi2_d0;
  cd1 = phi2_d1;
  tfac = tunit * etd_dt2; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (long inumer = 0; inumer <= dim_numer2; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (long idenom = 0; idenom < dim_denom2; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi2dt2.gen(MPIP, IO, Bas, etd_dt2, cisplit, dim_numer2, dim_denom2, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# phi2: real\n");
  // for (long ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi2_n, phi2_d0, phi2_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# phi2: imag\n");
  // for (long ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi2_n, phi2_d0, phi2_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void cletdrb::gen_phi3(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  phi3rk = etd_dt1 / SIX;

  dcomplex dfac;
  dcomplex roots[max_dim_denom + 1];

  std::vector<dcomplex> phi3_n(max_dim_numer + 1);
  std::vector<dcomplex> phi3_d0(max_dim_denom);
  std::vector<dcomplex> phi3_d1(max_dim_denom);

  if (dim_numer3 == 1 && dim_denom3 == 1) {
    //                x + 20   
    // phi3(x) = [- ----------]
    //              24 x - 120 
    dfac = -1.0;
    phi3_n[0] = 20.0/24.0;
    phi3_n[1] = 1.0/24.0;
    roots[0] = 5.0;
  } else if (dim_numer3 == 2 && dim_denom3 == 2) {
    //               2		       
    //              x  - 15 x + 420           
    // phi3(x) = [--------------------] 
    //                2		       
    //            60 x  - 720 x + 2520 
    dfac = 1.0;
    phi3_n[0] = 420.0/60.0;
    phi3_n[1] = -15.0/60.0;
    phi3_n[2] = 1.0/60.0;
    roots[0] = +6.000000000000000000 +2.449489742783177881 * IUNIT;
    roots[1] = +6.000000000000000000 -2.449489742783177881 * IUNIT;
  } else if (dim_numer3 == 2 && dim_denom3 == 3) {
    //                      2		         
    //                  11 x  - 140 x + 1120     
    // phi3(x) = - ------------------------------
    //                 3        2		 	     
    //             20 x  - 360 x  + 2520 x - 6720
    dfac = -1.0;
    phi3_n[0] = 1120.0/20.0;
    phi3_n[1] = -140.0/20.0;
    phi3_n[2] = 11.0/20.0;
    roots[0] = +6.651316808950019421 +0.000000000000000000 * IUNIT;
    roots[1] = +5.674341595524996507 +4.279971984629763249 * IUNIT;
    roots[2] = +5.674341595524996507 -4.279971984629763249 * IUNIT;
  } else if (dim_numer3 == 2 && dim_denom3 == 3 && pade_equiv) {
    //                    2		     
    //                   x  - 10 x + 40	     
    // phi3(x) = - --------------------------
    //                3       2		     
    //             2 x  - 24 x  + 120 x - 240
    dfac = -1.0;
    phi3_n[0] = 40.0/2.0;
    phi3_n[1] = -10.0/2.0;
    phi3_n[2] = 1.0/2.0;
    roots[0] = +4.644370709252173768 +0.000000000000000000 * IUNIT;
    roots[1] = +3.677814645373918445 +3.508761919567442433 * IUNIT;
    roots[2] = +3.677814645373918445 -3.508761919567442433 * IUNIT;
  } else if (dim_numer3 == 3 && dim_denom3 == 3) {
    //                   3       2		       
    //                  x  + 84 x  - 840 x + 10080     	   
    // phi3(x) = [- ----------------------------------]
    //                   3         2		       
    //              120 x  - 2520 x  + 20160 x - 60480 
    dfac = -1.0;
    phi3_n[0] = 10080.0/120.0;
    phi3_n[1] = -840.0/120.0;
    phi3_n[2] = 84.0/120.0;
    phi3_n[3] = 1.0/120.0;
    roots[0] = +7.653383973560514697 +0.000000000000000000 * IUNIT;
    roots[1] = +6.673308013219741319 +4.617378364686957504 * IUNIT;
    roots[2] = +6.673308013219741319 -4.617378364686957504 * IUNIT;
  } else {
    std::cout << "cletdrb:gen_exp. bad dim_numer3 or dim_denom3." << std::endl;
    abort();
  }

  phi3_d1[0] = +dfac;
  phi3_d0[0] = -dfac * roots[0];
  for (long idim = 1; idim < dim_denom3; idim ++) {
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
  for (long inumer = 0; inumer <= dim_numer3; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (long idenom = 0; idenom < dim_denom3; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi3dt1.gen(MPIP, IO, Bas, etd_dt1, cisplit, dim_numer3, dim_denom3, cn, cd0, cd1);

  // ##### dt * phi3(-i*h*dt/2) #####
  cn = phi3_n;
  cd0 = phi3_d0;
  cd1 = phi3_d1;
  tfac = tunit * etd_dt2; // BE CAREFUL!
  ttmp = RUNIT * etd_dt1; // BE CAREFUL!
  for (long inumer = 0; inumer <= dim_numer3; inumer ++) {
    cn[inumer] *= ttmp;
    ttmp *= tunit;
  }
  for (long idenom = 0; idenom < dim_denom3; idenom ++) {
    cd1[idenom] *= tfac;
  }
  phi3dt2.gen(MPIP, IO, Bas, etd_dt2, cisplit, dim_numer3, dim_denom3, cn, cd0, cd1);

  //DEBUG
  // dcomplex vala;
  // double dx = 0.001;
  // printf("\n");
  // printf("\n");
  // printf("# phi3: real\n");
  // for (long ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi3_n, phi3_d0, phi3_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, real(valx), real(vala));
  // }
  // printf("\n");
  // printf("\n");
  // printf("# phi3: imag\n");
  // for (long ix = 10; ix <= 500; ix ++) {
  //   double xr = dx * ix;
  //   dcomplex x = -xr * IUNIT;
  //   dcomplex valx = (exp(x)-1.0)/x;
  //   vala = test_pade(phi3_n, phi3_d0, phi3_d1, x);
  //   printf("%15.8f%15.8f%15.8f\n", xr, imag(valx), imag(vala));
  // }
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
