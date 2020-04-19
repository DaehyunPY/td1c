////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
void init(const clmpi& Proc, const clio& IO, const clbas& Bas, clwfn& Wfn)
{
  if (IO.job_type.compare("init") != 0 &&
      IO.job_type.compare("both") != 0 ) return;
  clock_t time0 = clock();

  clfield Field(Proc, IO);
  clhprod HPW(Proc, IO, Bas, Field);

  //DEBUG
  bool readno;
  IO.read_info("readno", false, readno);
  if (readno) {
    clwfn tWfn(Proc, IO, Bas);
    HPW.getno(Proc, IO, Bas, Wfn, tWfn);
    printf("lm components of natural orbitals:\n");
    tWfn.print_wang(Bas);
    HPW.copyo(Proc, Bas, tWfn, Wfn);
    ormas_cic0_(&Wfn.wfn[Wfn.size1]);
  }
  //DEBUG

  HPW.set_wfn0(Proc, IO, Bas, Wfn);

  std::string prop_type;
  IO.read_info("prop_type", "etdrb", prop_type);

  // initial CI diagonalization
  bool cidiag;
  IO.read_info("init_cidiag", false, cidiag);
  if (cidiag && Bas.ORMAS.ndetx > 1) {
    HPW.cidiag(Proc, Bas, Wfn);
  }

  if (prop_type.compare("split4") == 0) {
    std::cout << "init_split4 nyi." << std::endl;
    abort();
  } else if (prop_type.compare("split2") == 0 ||
	     prop_type.compare("split2_itr") == 0 ||
	     prop_type.compare("split2_itr_crk") == 0) {
    clh1itr H1P;
    clcrk H2P;
    init_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_itr_irk") == 0) {
    clh1itr H1P;
    clirk H2P;
    init_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_adi") == 0 ||
	     prop_type.compare("split2_adi_crk") == 0 ) {
    clh1adi H1P;
    clcrk H2P;
    init_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_adi_irk") == 0 ) {
    clh1adi H1P;
    clirk H2P;
    init_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_exp") == 0 ||
	     prop_type.compare("split2_exp_crk") == 0) {
    clh1exp H1P;
    clcrk H2P;
    init_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_exp_irk") == 0) {
    clh1exp H1P;
    clirk H2P;
    init_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("etdrb") == 0) {
    cletdrb H12P;
    init_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("etdrk") == 0) {
    cletdrk H12P;
    init_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("etdrk2") == 0) {
    cletdrk2 H12P;
    init_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("etdrk4") == 0) {
    cletdrk4 H12P;
    init_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("lawrk") == 0) {
    cllawrk H12P;
    init_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("lawab") == 0) {
    cllawab H12P;
    init_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else {
    std::cout << "init: bad prop_type" << std::endl;
    abort();
  }

  bool fock_diag;
  double lfield[9];
  std::vector<double> Eig(Bas.ORMAS.nfun);
  IO.read_info("fock_diag", false, fock_diag);
  Field.get_value(lfield);

  if (fock_diag) {
    HPW.fockdiag(Proc, Bas, lfield, Wfn, Eig);
  } else {
    HPW.orbene(Proc, Bas, lfield, Wfn, Eig);
  }

  printf("orbital energies:\n");
  for (int ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
    printf("%5d %25.15e\n", ifun, Eig[ifun]);
  }

  IO.read_info("post_cidiag", false, cidiag);
  if (cidiag && Bas.ORMAS.ndetx > 1) {
    HPW.cidiag(Proc, Bas, Wfn);
  }

  Wfn.write(Proc, IO, Bas);

  typhys Phys;
  HPW.spin(Proc, Bas, Wfn, Phys);
  printf("init: <S_z> = %25.15e\n", Phys.sz);
  printf("init: <S^2> = %25.15e\n", Phys.s2);
  HPW.oang(Proc, Bas, Wfn, Phys);
  printf("init: <L_z> = %25.15e\n", Phys.lz);
  printf("init: <L^2> = %25.15e\n", Phys.l2);

  //  if (clcontrol::psp && clcontrol::psp_type==2) {
  //    Bas.ppgenkb(Wfn.wfn);
  //  }

//  HPW.ladapt(Proc, Bas, Wfn);
//  HPW.oang(Proc, Bas, Wfn, Phys);
//  printf("init: <L_z> after projection = %20.10f\n", Phys.lz);
//  printf("init: <L^2> after projection = %20.10f\n", Phys.l2);
//  test = Phys.ene;
//  HPW.dipole3d(Proc, Bas, Wfn, Phys);
//  HPW.energy(Proc, Bas, lfield, Wfn, Phys);
//  test = Phys.ene - test;
//  printf("scf proj: %6ld %20.10f%20.10f%20.10f%20.10f%20.10f%20.10f%20.10f%20.5e\n",
//	 ncyc, Phys.dip[0], Phys.dip[1], Phys.dip[2], 
//	 Phys.acc[0], Phys.acc[1], Phys.acc[2], Phys.ene, test);

  printf("lm components of initial orbitals:\n");
  Wfn.print_wang(Bas);

  if (IO.iprint > 0) {
    clwfn tWfn(Proc, IO, Bas);
    HPW.getno(Proc, IO, Bas, Wfn, tWfn);
    printf("lm components of natural orbitals:\n");
    tWfn.print_wang(Bas);
  }
  //nyi  DEBUG: full-diagonalization
  //nyi  Ham.fock1(Bas);
  //nyi  Ham.fock2(Bas, Wfn);

  bool print_orb, print_cic;
  int print_orb_nstep;
  std::string fname;
  IO.read_info("print_orb", true, print_orb);
  if (print_orb) {
    fname = IO.name + ".orbital";
    IO.read_info("print_orb_nstep", LONE*10, print_orb_nstep);
    Wfn.print_orb(print_orb_nstep, fname, Bas);
//    printf("initial density interpolated:\n");
//    std::vector<dcomplex> rrad(Bas.GRad.nrad);
//    HPW.mkrrad(Wfn, rrad);
//    Wfn.print_orb(print_orb_nstep, IO, Bas, rrad);
  }

  IO.read_info("print_cic", true, print_cic);
  if (print_cic) {
    fname = IO.name + ".cicoeff";
    Wfn.print_cic(fname, Bas);
  }

  clock_t time1 = clock();
  std::cout << "# init: " << (double)(time1 - time0) / CLOCKS_PER_SEC << std::endl;
}
////////////////////////////////////////////////////////////////////////
void init_split2(const clmpi& Proc, const clio& IO, const clbas& Bas, clfield& Field, 
		 clh1prop& H1P, clh2prop& H2P, clhprod& HPW, clwfn& Wfn)
{
  int STEP1 = LONE;
  int STEP2 = LTWO;

  bool doh1 = clcontrol::split_type != 0;
  bool doh2 = Bas.ORMAS.neltot[2] > 1 || clcontrol::split_type == 0;  

  int init_maxcyc;
  double init_dstep;
  double init_tolene;
  IO.read_info("init_dstep", 0.01, init_dstep);
  IO.read_info("init_tolene", 1.E-13, init_tolene);
  IO.read_info("init_maxcyc", LONE*1000000, init_maxcyc);

  typhys Phys;
  int ncyc = 0;
  double test = ZERO;
  double ene1 = -1.0E+10;

  Field.step = ZERO;
  Field.time = ZERO;
  if (init_dstep > ZERO) Field.dtime = init_dstep;

  double lfield[9];
  H1P.gen(Proc, IO, Bas, Field, HPW);
  H2P.gen(Proc, IO, Bas, Field, HPW);

  for (int icyc = 0; icyc < init_maxcyc; icyc ++) {
    Field.get_value(lfield);
    HPW.dipole(Proc, Bas, lfield, Wfn, Phys);
    HPW.energy(Proc, Bas, lfield, Wfn, Phys);

    // convergence check
    ncyc = icyc;
    test = Phys.ene - ene1;
    if (get_abs(test) < init_tolene) {
      break;
    } else {
      ene1 = Phys.ene;
    }
 
    if (icyc == 0) {
      printf("cycle %10d %20.10f%20.10f%20.10f%20.10f\n", 
 	     icyc, Phys.dip[2], Phys.vel[2], Phys.acc[2], Phys.ene);
    } else {
      printf("cycle %10d %20.10f%20.10f%20.10f%20.10f%20.5e\n",
 	     icyc, Phys.dip[2], Phys.vel[2], Phys.acc[2], Phys.ene, test);
    }

    if (doh1) H1P.prop(Proc, Bas, Field, STEP1, HPW, Wfn);
    if (doh2) H2P.prop(Proc, Bas, Field, HPW, Wfn);
    if (doh1) H1P.prop(Proc, Bas, Field, STEP2, HPW, Wfn);

//  HPW.L_adapt(Proc, Bas, LZERO, Wfn);
    Wfn.madapt(Proc, Bas);
    Wfn.ladapt(Proc, Bas);
    Wfn.orth(Bas);
  }
  printf("scf done: %6ld %20.10f%20.10f%20.10f%20.10f%20.5e\n",
 	 ncyc, Phys.dip[2], Phys.vel[2], Phys.acc[2], Phys.ene, test);
}
////////////////////////////////////////////////////////////////////////
void init_prop12(const clmpi& Proc, const clio& IO, const clbas& Bas, clfield& Field, 
		 clh2prop& H2P, clhprod& HPW, clwfn& Wfn)
{
  int init_maxcyc;
  double init_dstep;
  double init_tolene;
  double init_tolwfn;
  IO.read_info("init_dstep", 0.01, init_dstep);
  IO.read_info("init_tolene", 1.E-12, init_tolene);
  IO.read_info("init_tolwfn", 1.E-6, init_tolwfn);
  IO.read_info("init_maxcyc", LONE*1000000, init_maxcyc);

  typhys Phys;
  int ncyc = 0;
  double test = ZERO;
  double ene1 = -1.0E+10;

  Field.step = ZERO;
  Field.time = ZERO;
  if (init_dstep > ZERO) Field.dtime = init_dstep;

  double lfield[9];
  H2P.gen(Proc, IO, Bas, Field, HPW);

  for (int icyc = 0; icyc < init_maxcyc; icyc ++) {
    Field.get_value(lfield);

    HPW.dipole(Proc, Bas, lfield, Wfn, Phys);

    //HPW.energy(Proc, Bas, lfield, Wfn, Phys);
    HPW.chkconv(lfield, Wfn, Phys);

    // convergence check
    ncyc = icyc;
    test = Phys.ene - ene1;
    if (get_abs(test) < init_tolene && Phys.ipx[2] < init_tolwfn) {
      break;
    } else {
      ene1 = Phys.ene;
    }
 
    if (icyc == 0) {
      printf("cycle %10d %20.10f%20.10f%20.10f%20.10f%20.5e%20.5e%20.5e\n", 
 	     icyc, Phys.dip[2], Phys.vel[2], Phys.acc[2], Phys.ene, 
	     Phys.ipx[0], Phys.ipx[1], Phys.ipx[2]);
    } else {
      printf("cycle %10d %20.10f%20.10f%20.10f%20.10f%20.5e%20.5e%20.5e%20.5e\n",
 	     icyc, Phys.dip[2], Phys.vel[2], Phys.acc[2], Phys.ene, test, 
	     Phys.ipx[0], Phys.ipx[1], Phys.ipx[2]);
    }

    H2P.prop(Proc, Bas, Field, HPW, Wfn);

    //NEW
    HPW.ProjHigh(Wfn);
    //NEW

//  HPW.L_adapt(Proc, Bas, LZERO, Wfn);
    Wfn.madapt(Proc, Bas);
    Wfn.ladapt(Proc, Bas);
    Wfn.orth(Bas);

//DEBUG CI diagonalization
//    HPW.cidiag(Proc, Bas, Wfn);
//    Wfn.madapt(Proc, Bas);
//    Wfn.ladapt(Proc, Bas);
//    Wfn.orth(Bas);
//DEBUG 
  }

  Field.get_value(lfield);
  HPW.dipole(Proc, Bas, lfield, Wfn, Phys);
  HPW.energy(Proc, Bas, lfield, Wfn, Phys);
  test = Phys.ene - ene1;
  printf("scf done: %6ld %20.10f%20.10f%20.10f%20.10f%20.5e\n",
 	 ncyc + 1, Phys.dip[2], Phys.vel[2], Phys.acc[2], Phys.ene, test);
}
////////////////////////////////////////////////////////////////////////
