////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
// Orimo_ECS
#include "surff.hpp"
// Orimo_ECS
////////////////////////////////////////////////////////////////////////
void tdse(const clmpi& Proc, const clio& IO, const clbas& Bas, clwfn& Wfn)
{
  if (IO.job_type.compare("td") != 0 &&
      IO.job_type.compare("both") != 0 ) return;
  clock_t time0 = clock();

  clfield Field(Proc, IO);
  Field.init();

  clhprod HPW(Proc, IO, Bas, Field);
  HPW.set_wfn0(Proc, IO, Bas, Wfn);

  bool kick_orb;
  IO.read_info("kick_orb", false, kick_orb);
  if (kick_orb) {
    double knorm;
    IO.read_info("knorm", 0.01, knorm);
    knorm *= 0.5291772086;
    HPW.kickk(Proc, IO, Bas, Wfn, knorm);
  }

  std::string prop_type;
  IO.read_info("prop_type", "etdrb", prop_type);

  if (prop_type.compare("split4") == 0 ||
      prop_type.compare("split4_irk") == 0 ||
      prop_type.compare("split4_crk") == 0) {
    tdse_split4(Proc, IO, Bas, Field, HPW, Wfn);
  } else if (prop_type.compare("split2") == 0 ||
	     prop_type.compare("split2_itr") == 0 ||
	     prop_type.compare("split2_itr_crk") == 0) {
    clh1itr H1P;
    clcrk H2P;
    tdse_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_itr_irk") == 0) {
    clh1itr H1P;
    clirk H2P;
    tdse_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_adi") == 0 ||
	     prop_type.compare("split2_adi_crk") == 0 ) {
    clh1adi H1P;
    clcrk H2P;
    tdse_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_adi_irk") == 0 ) {
    clh1adi H1P;
    clirk H2P;
    tdse_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_exp") == 0 ||
	     prop_type.compare("split2_exp_crk") == 0) {
    clh1exp H1P;
    clcrk H2P;
    tdse_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("split2_exp_irk") == 0) {
    clh1exp H1P;
    clirk H2P;
    tdse_split2(Proc, IO, Bas, Field, H1P, H2P, HPW, Wfn);
  } else if (prop_type.compare("etdrb") == 0) {
    cletdrb H12P;
    tdse_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("etdrk") == 0) {
    cletdrk H12P;
    tdse_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("etdrk2") == 0) {
    cletdrk2 H12P;
    tdse_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("etdrk4") == 0) {
    cletdrk4 H12P;
    tdse_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("lawrk") == 0) {
    cllawrk H12P;
    tdse_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else if (prop_type.compare("crnic") == 0) {
    clcrnic H12P;
    tdse_prop12(Proc, IO, Bas, Field, H12P, HPW, Wfn);
  } else {
    std::cout << "tdse: bad prop_type" << std::endl;
    abort();
  }

  Wfn.write(Proc, IO, Bas);
// Sato_tSURFF
// Orimo_ECS
  if (clcontrol::tsurff) Wfn.surff->print_amplitude(IO);
// Orimo_ECS
// Sato_tSURFF
  clock_t time1 = clock();
  std::cout << "# tdse: " << (double)(time1 - time0) / CLOCKS_PER_SEC << std::endl;

  typhys Phys;
  HPW.spin(Proc, Bas, Wfn, Phys);
  printf("tdse: <S_z> = %25.15e\n", Phys.sz);
  printf("tdse: <S^2> = %25.15e\n", Phys.s2);
  HPW.oang(Proc, Bas, Wfn, Phys);
  printf("tdse: <L_z> = %25.15e\n", Phys.lz);
  printf("tdse: <L^2> = %25.15e\n", Phys.l2);
  printf("lm-components of final orbitals:\n");
  Wfn.print_wang(Bas);

  bool print_orb, print_cic;
  long print_orb_nstep;
  std::string fname;
  IO.read_info("print_orb", false, print_orb);
  if (print_orb) {
    fname = IO.name + ".orbital";
    IO.read_info("print_orb_nstep", LONE*10, print_orb_nstep);
    Wfn.print_orb(print_orb_nstep, fname, Bas);
  }

  IO.read_info("print_cic", false, print_cic);
  if (print_cic) {
    fname = IO.name + ".cicoeff";
    Wfn.print_cic(fname, Bas);
  }

}
////////////////////////////////////////////////////////////////////////
void tdse_split4(const clmpi& Proc, const clio& IO, const clbas& Bas, clfield& Field, 
		 clhprod& HPW, clwfn& Wfn)
{
  const double w1 = ONE / (TWO - pow(TWO, THIRD));
  const double w2 = ONE - TWO * w1;
  const double dt1 = w1 * Field.dtime;
  const double dt2 = w2 * Field.dtime;

  std::string prop_type;
  IO.read_info("prop_type", "split4", prop_type);
  if (prop_type.compare("split4") == 0 ||
      prop_type.compare("split4_irk") == 0) {
    clh1itr H1P1(Proc, IO, Bas, Field, dt1, HPW);
    clh1itr H1P2(Proc, IO, Bas, Field, dt2, HPW);
    clirk H2P(Proc, IO, Bas, Field, HPW);
    tdse_split4x(Proc, IO, Bas, Field, H1P1, H1P2, H2P, dt1, dt2, HPW, Wfn);
  } else if (prop_type.compare("split4_crk") == 0) {
    clh1itr H1P1(Proc, IO, Bas, Field, dt1, HPW);
    clh1itr H1P2(Proc, IO, Bas, Field, dt2, HPW);
    clcrk H2P(Proc, IO, Bas, Field, HPW);
    tdse_split4x(Proc, IO, Bas, Field, H1P1, H1P2, H2P, dt1, dt2, HPW, Wfn);
  } else {
    std::cout << "tdse_split4: bad prop_type" << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void tdse_split4x(const clmpi& Proc, const clio& IO, const clbas& Bas, clfield& Field, 
		  clh1prop& H1P1, clh1prop& H1P2, clh2prop& H2P, double dt1, double dt2,
		  clhprod& HPW, clwfn& Wfn)
{
  double time;
  H1P1.gen(Proc, IO, Bas, Field, dt1, HPW);
  H1P2.gen(Proc, IO, Bas, Field, dt2, HPW);
  H2P.gen(Proc, IO, Bas, Field, HPW);

  bool doh1 = clcontrol::split_type != 0;
  bool doh2 = Bas.ORMAS.neltot[2] > 1 || clcontrol::split_type == 0;

  while (! Field.finished()) {
    tdse_print(Proc, IO, Bas, Field, HPW, Wfn);

    time = Field.time;
    if (doh1) H1P1.prop(Proc, Bas, Field, time, dt1, HPW, Wfn);
    if (doh2) H2P.prop(Proc, Bas, Field, time, dt1, HPW, Wfn);
    if (doh1) H1P1.prop(Proc, Bas, Field, time + dt1 * HALF, dt1, HPW, Wfn);

    time += dt1;
    if (doh1) H1P2.prop(Proc, Bas, Field, time, dt2, HPW, Wfn);
    if (doh2) H2P.prop(Proc, Bas, Field, time, dt2, HPW, Wfn);
    if (doh1) H1P2.prop(Proc, Bas, Field, time + dt2 * HALF, dt2, HPW, Wfn);

    time += dt2;
    if (doh1) H1P1.prop(Proc, Bas, Field, time, dt1, HPW, Wfn);
    if (doh2) H2P.prop(Proc, Bas, Field, time, dt1, HPW, Wfn);
    if (doh1) H1P1.prop(Proc, Bas, Field, time + dt1 * HALF, dt1, HPW, Wfn);

    Wfn.mask(Bas);

    Field.step ++;
    Field.time += Field.dtime;
  }
}
////////////////////////////////////////////////////////////////////////
void tdse_split2(const clmpi& Proc, const clio& IO, const clbas& Bas, clfield& Field, 
		 clh1prop& H1P, clh2prop& H2P, clhprod& HPW, clwfn& Wfn)
{
  H1P.gen(Proc, IO, Bas, Field, HPW);
  H2P.gen(Proc, IO, Bas, Field, HPW);

  long STEP1 = LONE;
  long STEP2 = LTWO;

  bool doh1 = clcontrol::split_type != 0;
  bool doh2 = Bas.ORMAS.neltot[2] > 1 || clcontrol::split_type != 1;

// Sato_tSURFF
//// Orimo_ECS
//  clwfn swfn(Proc, IO, Bas);
//// Orimo_ECS
// Sato_tSURFF

  while (! Field.finished()) {
    tdse_print(Proc, IO, Bas, Field, HPW, Wfn);

// Sato_tSURFF
//// Orimo_ECS
//    if(Wfn.surff->nprint > 0){
//      swfn.copy(Bas, Wfn);
//      if (doh2) H2P.prop(Proc, Bas, Field, HPW, swfn); // calculate v2xmat of swfn
//      Wfn.surff->v2xmat = swfn.surff->v2xmat; // insert correct v2xmat from swfn.surff->v2xmat
//      Wfn.surff->prop(Proc, Wfn, Bas, Field, HPW);
//    }
//// Orimo_ECS
// Sato_tSURFF

    if (doh1) H1P.prop(Proc, Bas, Field, STEP1, HPW, Wfn);
    if (doh2) H2P.prop(Proc, Bas, Field, HPW, Wfn);
    if (doh1) H1P.prop(Proc, Bas, Field, STEP2, HPW, Wfn);

    Wfn.mask(Bas);

    Field.step ++;
    Field.time += Field.dtime;
  }
}
////////////////////////////////////////////////////////////////////////
void tdse_prop12(const clmpi& Proc, const clio& IO, const clbas& Bas, clfield& Field, 
		 clh2prop& H2P, clhprod& HPW, clwfn& Wfn)
{
  H2P.gen(Proc, IO, Bas, Field, HPW);
//std::cout << "after H2P.gen" << std::endl;
//abort();

  while (! Field.finished()) {
    tdse_print(Proc, IO, Bas, Field, HPW, Wfn);
    //DEBUG
    //    std::cout << "ABORT for debug @ tdse_prop12." << std::endl;
    //    abort();
    //DEBUG
    H2P.prop(Proc, Bas, Field, HPW, Wfn);

    //NEW
    HPW.ProjHigh(Wfn);
    //NEW

    Wfn.mask(Bas);
    Field.step ++;
    Field.time += Field.dtime;
    //DEBUG
    //    std::cout << "ABORT for debug @ tdse_prop12." << std::endl;
    //    abort();
    //DEBUG
  }
}
////////////////////////////////////////////////////////////////////////
void tdse_print(const clmpi& Proc, const clio& IO, const clbas& Bas, 
		const clfield& Field, clhprod& HPW, const clwfn& Wfn)
{
  //RTDISP
  //  rtdisp_print(Proc, IO, Bas, Field, HPW, Wfn);
  //  return;
  //RTDISP

  typhys Phys;
  double lfield[9];
  Field.get_value(lfield);

  if (IO.nprint_ene > 0 && Field.step % IO.nprint_ene == 0) {
    HPW.energy(Proc, Bas, lfield, Wfn, Phys);
    fprintf(IO.fp_ene, "%10ld %20.10e %20.10e %20.10e %20.10e %20.10e\n",
	    Field.step, Field.time, Field.ncyc(), lfield[5], lfield[8], Phys.ene);
  }

  if (IO.nprint_dip > 0 && Field.step % IO.nprint_dip == 0) {
    HPW.dipole(Proc, Bas, lfield, Wfn, Phys);
    if (IO.nprint_opx > 0) {
      fprintf(IO.fp_op0, "%10ld %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
	      Field.step, Field.time, Field.ncyc(), lfield[5], lfield[8], Phys.dip[0], Phys.vel[0], Phys.acc[0]);
      fprintf(IO.fp_opd, "%10ld %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
	      Field.step, Field.time, Field.ncyc(), lfield[5], lfield[8], Phys.dip[1], Phys.vel[1], Phys.acc[1]);
    }
    fprintf(IO.fp_op1, "%10ld %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
	    Field.step, Field.time, Field.ncyc(), lfield[5], lfield[8], Phys.dip[2], Phys.vel[2], Phys.acc[2]);
  }

//  if (IO.nprint_opx > 0 && Field.step % IO.nprint_opx == 0) {
//    HPW.dipole_ion(Proc, Bas, lfield, Wfn, Phys);
//    fprintf(IO.fp_opx, "%10ld %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
//	    Field.step, Field.time, Field.ncyc(), lfield[5], lfield[8], Phys.dip[2], Phys.vel[2], Phys.acc[2]);
//  }

  if (IO.nprint_ion > 0 && Field.step % IO.nprint_ion == 0) {
    HPW.normx(Proc, Bas, lfield, Wfn, Phys);
    fprintf(IO.fp_ipx, "%10ld %20.10e %20.10e %20.10e %20.10e", 
	    Field.step, Field.time, Field.ncyc(), lfield[5], lfield[8]);
    fprintf(IO.fp_ipx, " %20.10e", Phys.ipx[0]);

    HPW.ionpx(Proc, Bas, lfield, Wfn, Phys);
    for (long iipx = 0; iipx <= clcontrol::maxipx; iipx ++) {
      fprintf(IO.fp_ipx, " %20.10e", Phys.ipx[iipx]);
    }
    fprintf(IO.fp_ipx, "\n");
  }

  if (IO.nprint_ipd > 0 && Field.step % IO.nprint_ipd == 0) {
    fprintf(IO.fp_ipd, "%10ld %20.10e %20.10e %20.10e %20.10e", 
	    Field.step, Field.time, Field.ncyc(), lfield[5], lfield[8]);
    fprintf(IO.fp_ipd, " %20.10e", Phys.ipx[0]);

    HPW.ionpd(Proc, Bas, lfield, Wfn, Phys);
    for (long iipx = 0; iipx <= clcontrol::maxipx; iipx ++) {
      fprintf(IO.fp_ipd, " %20.10e", Phys.ipx[iipx]);
    }
    fprintf(IO.fp_ipd, "\n");
  }

  if (IO.nprint_rrad > 0 && Field.step % IO.nprint_rrad == 0) {
    double drrad, rradp, rradm;
    std::vector<dcomplex> rrad(Bas.GRad.nrad + 1);
    std::vector<dcomplex> rradpw((Bas.GRad.nrad + 1)*(Bas.GAng.lmax1 + 1));
    HPW.mkrradx(Wfn, rrad, rradpw);
    //    HPW.mkrrad(Wfn, rrad, rradpw);
    //    HPW.mkrrad1(Wfn, rrad, rradpw);
    fprintf(IO.fp_rrad, "###\n");
    for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
      drrad = rrad[irad].real();
// Orimo_ECS
      if(Bas.GRad.inf_range  ){ 
	long infx = Bas.GRad.nrad - Bas.GRad.ndvr[Bas.GRad.nfe-1];
	if(irad >= infx) drrad = drrad * exp(-Bas.GRad.exp_factor*2*(Bas.GRad.xrad[irad] - Bas.GRad.xrad[infx]));
      }
// Orimo_ECS
      fprintf(IO.fp_rrad, " %20.10e %20.10e", Bas.GRad.xrad[irad], drrad);

      for (long l = 0; l <= Bas.GAng.lmax1; l++) {
	drrad = rradpw[(Bas.GRad.nrad+1)*l+irad].real();
	fprintf(IO.fp_rrad, "%20.10e", drrad);
      }
      fprintf(IO.fp_rrad, "\n");
    }
    fprintf(IO.fp_rrad, "\n");
    fprintf(IO.fp_rrad, "\n");
  }

// Sato_tSURFF
// Orimo_ECS
  if (clcontrol::tsurff) {
    if(Wfn.surff->nprint > 0 && Field.step % Wfn.surff->nprint == 0){
      Wfn.surff->calc_me_spectrum(Proc, Wfn, Bas, HPW);
      Wfn.surff->print_file(IO);
      if(Wfn.surff->is_m_resolved){
	Wfn.surff->print_m_resolved_ang_spec(IO);
      }
    }
    HPW.htotx(Proc, Bas, ONE, lfield, Wfn);
    Wfn.surff->rec_v2xmat(ONE);
    Wfn.surff->prop(Proc, Wfn, Bas, Field, HPW);
  }
// Orimo_ECS
// Sato_tSURFF

  if (IO.nprint_full > 0 && Field.step % IO.nprint_full == 0) {
    double wfnr, wfni;
    // write out orbitals
    fprintf(IO.fp_torb, "###\n");
    fprintf(IO.fp_torb, "%10ld\n", Bas.GRad.nrad);
    fprintf(IO.fp_torb, "%10ld\n", Bas.GAng.lmax1);
    fprintf(IO.fp_torb, "%10ld\n", Bas.GAng.mmax1);
    fprintf(IO.fp_torb, "%10ld\n", Bas.ORMAS.nfun);
    long ind;
    for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
      for (long l = 0; l <= Bas.GAng.lmax1; l ++) {
	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
             	                  + l * (Bas.GRad.nrad - 1) + irad - 1;
	  wfnr = real(Wfn.wfn[ind]);
	  wfni = imag(Wfn.wfn[ind]);
	  fprintf(IO.fp_torb, "%25.15le%25.15le\n", wfnr, wfni);
	}
      }
    }
    fprintf(IO.fp_torb, "\n");
    fprintf(IO.fp_torb, "\n");

    fprintf(IO.fp_tcic, "###\n");
    fprintf(IO.fp_tcic, "%10ld\n", Bas.ORMAS.lcic);
    for (long idet = 0; idet < Bas.ORMAS.lcic; idet ++) {
      wfnr = real(Wfn.wfn[Wfn.size1+idet]);
      wfni = imag(Wfn.wfn[Wfn.size1+idet]);
      fprintf(IO.fp_tcic, "%25.15le%25.15le\n", wfnr, wfni); 
    }
    fprintf(IO.fp_tcic, "\n");
    fprintf(IO.fp_tcic, "\n");
  }
}
////////////////////////////////////////////////////////////////////////
//void tdse_print_rho(const clmpi& Proc, const clio& IO, const clbas& Bas, const std::vector<dcomplex>& vec)
//{
//  FILE *fpo;
//  fpo = fopen(FNAME.c_str(), "w");
//
//  long ind, irad, num_dvr;
//  double x_ll, x_ul, dpt;
//  double point, orbr, orbi;
//  dcomplex vec_val;
//  std::vector<dcomplex> bas_val(Bas.GRad.nmax + 1);
//
//  for (long ife = 0; ife < femax; ife ++) {
//    x_ll = Bas.GRad.get_x0(ife);
//    x_ul = Bas.GRad.get_x1(ife);
//    dpt = (x_ul - x_ll) / npt;
//    num_dvr = Bas.GRad.get_ndvr(ife);
//
//    point = x_ll;
//    while (point + dpt * HALF < x_ul) {
//
//      for(long idvr = 0; idvr <= num_dvr; idvr ++) {
//	irad = Bas.GRad.mapf[ife] + idvr;
//	bas_val[idvr] = Bas.GRad.get_val(ife, idvr, point);
//      }
//
//      fprintf(fpo, "%20.10f", point);
//      for (long l = 0; l < ndata; l ++) {
//	vec_val = CZERO;
//	for(long idvr = 0; idvr <= num_dvr; idvr ++) {
//	  irad = Bas.GRad.mapf[ife] + idvr;
//	  if (irad > 0 && irad < Bas.GRad.nrad) {
//	    ind = l * (Bas.GRad.nrad - 1) + irad - 1;
//	    if (nwrad == 1) {
//	      vec_val += vec[ind] * bas_val[idvr] / sqrt(Bas.GRad.wrad[irad]);
//	    } else if (nwrad == 2) {
//	      vec_val += vec[ind] * bas_val[idvr] / Bas.GRad.wrad[irad];
//	    }
//	  }
//	}
//	orbr = real(vec_val);
//	orbi = imag(vec_val);
//	fprintf(fpo, "%15.5e%15.5e", orbr, orbi);
//      }
//      fprintf(fpo, "\n");
//      point += dpt;
//    }
//  }
//
//  fclose(fpo);
//}
////////////////////////////////////////////////////////////////////////
//void rtdisp_print(const clmpi& Proc, const clio& IO, const clbas& Bas, 
//		  const clfield& Field, clhprod& HPW, const clwfn& Wfn)
//{
//  typhys Phys;
//  double lfield[9];
//  if (Field.step > 0) {
//    Field.get_value(lfield);
//    HPW.dipole(Proc, Bas, lfield, Wfn, Phys);
//    fprintf(IO.fp_op1, "%10ld %20.10e %20.10e %20.10e %20.10e", Field.step, Field.time, Field.ncyc(), lfield[5], Phys.dip[2]);
//
//    double resp, ufreq;
//    double maxresp = ZERO;
//    for (int i = 0; i <= 20; i ++) {
//      ufreq = 0.5 * i;
//      //      resp = exp(-ufreq*(Field.step-1)*Field.dtime) * Phys.dip[2] / Field.Fields[0].famp;
//      resp = exp(-ufreq*(Field.step-1)*Field.dtime) * Phys.dip[2];
//      fprintf(IO.fp_op1, "%20.10e", resp);
//    }
//    fprintf(IO.fp_op1, "\n");
//  }
//}
////////////////////////////////////////////////////////////////////////
