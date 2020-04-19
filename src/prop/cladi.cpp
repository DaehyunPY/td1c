////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
cladi::cladi()
{
}
////////////////////////////////////////////////////////////////////////
cladi::~cladi()
{
}
////////////////////////////////////////////////////////////////////////
cladi::cladi(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
	     double dt, int td_type)
	     
{
  gen(MPIP, IO, Bas, dt, td_type, LONE);
}
////////////////////////////////////////////////////////////////////////
cladi::cladi(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
	     double dt, int td_type, int split_type)
	     
{
  gen(MPIP, IO, Bas, dt, td_type, split_type);
}
////////////////////////////////////////////////////////////////////////
void cladi::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		double dt, int td_type, int split_type)
{
  printf("cladi: instance created.\n");

  IO.read_info("adi_nddt", LONE, adi_nddt);
//  if (adi_nddt != 1) {
//    std::cout << "cladi::gen. adi_nddt must be 1." << std::endl;
//    abort();
//  }

  dtime = dt / adi_nddt;
  icomp = td_type;

  if (split_type == 1) {
    dokin = 1;
    donuc = 1;
    doext = 1;
  } else if (split_type == 2) {
    dokin = 1;
    donuc = 0;
    doext = 0;
  } else if (split_type == 3) {
    dokin = 1;
    donuc = 1;
    doext = 0;
  } else if (split_type == 4) {
    dokin = 1;
    donuc = 0;
    doext = 1;
    std::cout << "cladi::gen: split_type = 4 is nyi..." << std::endl;
    abort();
  } else {
    std::cout << "cladi::gen: bad split_type" << std::endl;
    abort();
  }

  int size;
  size = (2 * Bas.GRad.nmax + 1) * (Bas.GRad.nrad - 1) * (Bas.Sph1.lmax + 1);
  tadi1.resize(size);
  size = (3 * Bas.GRad.nmax + 1) * (Bas.GRad.nrad - 1) * (Bas.Sph1.lmax + 1);
  tadi2.resize(size);
  size = (Bas.GRad.nrad - 1) * (Bas.Sph1.lmax + 1);
  tpiv2.resize(size);
  adi_gen_tadi_(&icomp, &dokin, &donuc, &dtime, &tadi1[0], &tadi2[0], &tpiv2[0]);
}
////////////////////////////////////////////////////////////////////////
void cladi::prop(const clmpi& Proc, const clbas& Bas, 
		 double time, const clfield& Field, clwfn& Wfn)
{
  //old  prop1(Proc, Bas, time, Field, Wfn);
  //old  prop2(Proc, Bas, time, Field, Wfn);
  if (Field.gauge.compare("length1") == 0) {
    propl(Proc, Bas, time, Field, Wfn);
  } else if (Field.gauge.compare("length2") == 0) {
    propl(Proc, Bas, time, Field, Wfn);
  } else if (Field.gauge.compare("velocity1") == 0) {
    propv(Proc, Bas, time, Field, Wfn);
  } else if (Field.gauge.compare("velocity2") == 0) {
    propv(Proc, Bas, time, Field, Wfn);
  }
}
////////////////////////////////////////////////////////////////////////
void cladi::propl(const clmpi& Proc, const clbas& Bas, 
		  double time, const clfield& Field, clwfn& Wfn)
{
  //  std::cout << "cladi::propl" << std::endl;
  double tnow;
  double dt2 = dtime * HALF;
  for (int istep = 0; istep < adi_nddt; istep ++) {
    if (doext == 1) {
      tnow = time + dtime * istep;
      laser_lgauge(Proc, Bas, 0, tnow, dt2, Field, Wfn);
    }
    t_explicit(Proc, Bas, Wfn);
    t_implicit(Proc, Bas, Wfn);
    if (doext == 1) {
      tnow = time + dtime * (istep + 1);
      laser_lgauge(Proc, Bas, LONE, tnow, dt2, Field, Wfn);
    }
  }
}
////////////////////////////////////////////////////////////////////////
void cladi::propv(const clmpi& Proc, const clbas& Bas, 
		  double time, const clfield& Field, clwfn& Wfn)
{
  //  std::cout << "cladi::propv" << std::endl;
  double tnow;
  double dt2 = dtime * HALF;
  for (int istep = 0; istep < adi_nddt; istep ++) {
    if (doext == 1) {
      tnow = time + dtime * istep;
      laser_vgauge(Proc, Bas, 0, tnow, dt2, Field, Wfn);
    }
    t_explicit(Proc, Bas, Wfn);
    t_implicit(Proc, Bas, Wfn);
    if (doext == 1) {
      tnow = time + dtime * (istep + 1);
      laser_vgauge(Proc, Bas, LONE, tnow, dt2, Field, Wfn);
    }
  }
}
//////////////////////////////////////////////////////////////////////////
void cladi::t_explicit(const clmpi& MPIP, const clbas& Bas, clwfn& Wfn)
{
  adi_t_explicit_(&tadi1[0], &Wfn.wfn[0], &clhprod::orb[0]);
}
//////////////////////////////////////////////////////////////////////////
void cladi::t_implicit(const clmpi& MPIP, const clbas& Bas, clwfn& Wfn)
{
  adi_t_implicit_(&tadi2[0], &tpiv2[0], &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void cladi::laser_lgauge(const clmpi& MPIP, const clbas& Bas, int istag, 
			 double time, double dt, const clfield& Field, clwfn& Wfn)
{
  double lfield[9];
  Field.get_value(time, lfield);
  adi_laser_lgauge_(&icomp, &istag, &dt, lfield, &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void cladi::laser_vgauge(const clmpi& MPIP, const clbas& Bas, int istag, 
			 double time, double dt, const clfield& Field, clwfn& Wfn)
{
  double lfield[9];
  Field.get_value(time, lfield);
  adi_laser_vgauge_(&icomp, &istag, &dt, lfield, &Wfn.wfn[0]);
}
////////////////////////////////////////////////////////////////////////
//old void cladi::prop1(const clmpi& Proc, const clbas& Bas,
//old 		  double time, const clfield& Field, clwfn& Wfn)
//old {
//old   // DEBUG
//old   //  std::cout << "Do adi_gen_tadi @ every time step!" << std::endl;
//old   //  adi_gen_tadi_(&icomp, &dokin, &donuc, 
//old   //		&Bas.GRad.nrad, &Bas.Sph1.nsph, &Bas.Sph1.lmax, &Bas.Sph1.mmax,
//old   //		&Bas.GRad.nmax, &Bas.znuc, &dtime, &Bas.GRad.xrad[0],
//old   //		&Bas.GRad.kinetic(0,0), &tadi1[0], &tadi2[0], &tpiv2[0]);
//old   // DEBUG
//old   std::cout << "cladi::prop1" << std::endl;
//old 
//old   double tnow;
//old   for (int istep = 0; istep < adi_nddt; istep ++) {
//old     tnow = time + dtime * (istep + HALF);
//old     t_explicit(Proc, Bas, Wfn);
//old     if (doext == 1) v_cayley(Proc, Bas, tnow, dtime, Field, Wfn);
//old     t_implicit(Proc, Bas, Wfn);
//old   }
//old }
//old ////////////////////////////////////////////////////////////////////////
//old void cladi::prop2(const clmpi& Proc, const clbas& Bas, 
//old 		  double time, const clfield& Field, clwfn& Wfn)
//old {
//old   // DEBUG
//old   //  std::cout << "Do adi_gen_tadi @ every time step!" << std::endl;
//old   //  adi_gen_tadi_(&icomp, &dokin, &donuc, 
//old   //		&Bas.GRad.nrad, &Bas.Sph1.nsph, &Bas.Sph1.lmax, &Bas.Sph1.mmax,
//old   //		&Bas.GRad.nmax, &Bas.znuc, &dtime, &Bas.GRad.xrad[0],
//old   //		&Bas.GRad.kinetic(0,0), &tadi1[0], &tadi2[0], &tpiv2[0]);
//old   // DEBUG
//old   std::cout << "cladi::prop2" << std::endl;
//old 
//old   double tnow;
//old   double dt2 = dtime * HALF;
//old   for (int istep = 0; istep < adi_nddt; istep ++) {
//old     if (doext == 1) {
//old       tnow = time + dtime * istep;
//old       v_cayley(Proc, Bas, tnow, dt2, Field, Wfn);
//old     }
//old     t_explicit(Proc, Bas, Wfn);
//old     t_implicit(Proc, Bas, Wfn);
//old     if (doext == 1) {
//old       tnow = time + dtime * (istep + 1);
//old       v_cayley(Proc, Bas, tnow, dt2, Field, Wfn);
//old     }
//old   }
//old }
//old //////////////////////////////////////////////////////////////////////////
//old void cladi::v_cayley(const clmpi& MPIP, const clbas& Bas, 
//old 		     double time, double dt, const clfield& Field, clwfn& Wfn)
//old {
//old   Bas.sph2ang1(MPIP, Wfn.wfn, Wfn.wfng);
//old 
//old   Field.get_value(time, lfield);
//old   adi_v_cayley_(&icomp, &Bas.GRad.nrad, &Bas.Sph1.nang, &Bas.ORMAS.nfun, &dt, lfield,
//old 		&Bas.grid(0,0), &Wfn.wfng[0]);
//old 
//old   Bas.ang2sph1(MPIP, Wfn.wfng, Wfn.wfn);
//old }
//old //////////////////////////////////////////////////////////////////////////
