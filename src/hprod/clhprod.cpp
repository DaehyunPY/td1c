////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
// static members
////////////////////////////////////////////////////////////////////////
int clhprod::num_hprod = 0;
////////////////////////////////////////////////////////////////////////
clhprod::clhprod()
{
}
////////////////////////////////////////////////////////////////////////
clhprod::~clhprod()
{
  num_hprod --;
  std::cout << "clhprod::~clhprod: num_hprod = " << num_hprod << std::endl;
  if (num_hprod != 0) abort();
  hprod_final_();
}
////////////////////////////////////////////////////////////////////////
clhprod::clhprod(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		 const clfield& Field)
{
  gen(MPIP, IO, Bas, Field);
}
////////////////////////////////////////////////////////////////////////
void clhprod::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		  const clfield& Field)
{
  num_hprod ++;
  std::cout << "clhprod::gen: num_hprod = " << num_hprod << std::endl;
  if (num_hprod != 1) abort();

  IO.read_info("projhigh", false, projhigh);
  if (projhigh) {
    IO.read_info("projhigh_cutoff", projhigh_cutoff);
  } else {
    projhigh_cutoff = ZERO;
  }

  // c++/fortran bindings
  hprod_bind_(&ene_fcore, &ene_dcore, &ene_core, &ene_act, &ene_tot,
	      &dip_exp, &vel_exp, &acc_exp, &projhigh, &projhigh_cutoff);
  // allocate static arrays
  hprod_init_();
}
//////////////////////////////////////////////////////////////////////////
void clhprod::set_wfn0(const clmpi& MPIP, const clio& IO, const clbas& Bas, const clwfn& Wfn) const
{
  std::string v2xfc_type;
  IO.read_info("v2xfc_type", "kli", v2xfc_type);
  hprod_set_wfn0_(&Wfn.wfn[0], &Wfn.wfn[Wfn.size1], v2xfc_type.c_str(), v2xfc_type.length());
}
//////////////////////////////////////////////////////////////////////////
void clhprod::ProjHigh(clwfn& Wfn) const
{
  hprod_projhigh_(&Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::htot(const clmpi& MPIP, const clbas& Bas, double dt,
		   const double* lfield, const clwfn& Wfn, clwfn& hWfn)
{
  hprod_htot_(&dt, lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &hWfn.wfn[0], &hWfn.wfn[hWfn.size1]);
  hWfn.madapt(MPIP, Bas);
  hWfn.ladapt(MPIP, Bas);
  //  ProjHigh(hWfn);

  //debug static int ncalled = 0;
  //debug std::cout << "ncalled = " << ncalled << std::endl;
  //debug ncalled ++;
  //debug if (ncalled == 10) {
  //debug   int npt = 5;
  //debug   std::string fname;
  //debug   fname = "test.orb";
  //debug   Wfn.print_orb(npt, fname, Bas);
  //debug   fname = "test.horb";
  //debug   hWfn.print_orb(npt, fname, Bas);
  //debug   exit(1);
  //debug }
}
//////////////////////////////////////////////////////////////////////////
void clhprod::htotx(const clmpi& MPIP, const clbas& Bas, double dt,
		   const double* lfield, const clwfn& Wfn)
{
  hprod_htotx_(&dt, lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::htoto(const clmpi& MPIP, const clbas& Bas, double dt,
		    const double* lfield, const clwfn& Wfn, clwfn& hWfn)
{
  hprod_htoto_(&dt, lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &hWfn.wfn[0]);
  hWfn.madapt(MPIP, Bas);
  hWfn.ladapt(MPIP, Bas);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::energy(const clmpi& MPIP, const clbas& Bas, const double* lfield,
		     const clwfn& Wfn, typhys& Phys)
{
  hprod_energy_(lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1]);
  Phys.ene = ene_tot;
}
//////////////////////////////////////////////////////////////////////////
void clhprod::chkconv(const double* lfield, const clwfn& Wfn, typhys& Phys)
{
  hprod_chkconv_(lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &Phys.ene, &Phys.ipx[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::h1tot(const clmpi& MPIP, const clbas& Bas, double dt,
		    const double* lfield, const clwfn& Wfn, clwfn& hWfn)
{
  bool dofc = false;
  clear(MPIP, Bas, hWfn);
  hprod_h1tot_(&dofc, lfield, &Wfn.wfn[0], &hWfn.wfn[0]);
  hWfn.madapt(MPIP, Bas);
  hWfn.ladapt(MPIP, Bas);
  //  ProjHigh(hWfn);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::v1tot(const clmpi& MPIP, const clbas& Bas, double dt,
		    const double* lfield, const clwfn& Wfn, clwfn& hWfn)
{
  hprod_v1tot_(&dt, lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &hWfn.wfn[0], &hWfn.wfn[hWfn.size1]);
  hWfn.madapt(MPIP, Bas);
  hWfn.ladapt(MPIP, Bas);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::h1add(const clmpi& MPIP, const clbas& Bas, dcomplex zfac,
		    const double* lfield, const clwfn& Wfn, clwfn& hWfn)
{
  bool dofc = false;
  hprod_h1add_(&dofc, &zfac, lfield, &Wfn.wfn[0], &hWfn.wfn[0]);
  //  hWfn.madapt(MPIP, Bas);
  //  hWfn.ladapt(MPIP, Bas);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::v1ext(const clmpi& MPIP, const clbas& Bas, dcomplex zfac,
		    const double* lfield, const clwfn& Wfn, clwfn& hWfn)
{
  bool dofc = false;
  hprod_v1ext_(&dofc, &zfac, lfield, &Wfn.wfn[0], &hWfn.wfn[0]);
  //  hWfn.madapt(MPIP, Bas);
  //  hWfn.ladapt(MPIP, Bas);
  //  ProjHigh(hWfn);
}
////////////////////////////////////////////////////////////////////////
void clhprod::fulldiag(const clmpi& MPIP, const clbas& Bas, clwfn& Wfn)
{
  hprod_fulldiag_(&Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::fockdiag(const clmpi& MPIP, const clbas& Bas, const double* lfield,
		       clwfn& Wfn, std::vector<double>& Eig)
{
  hprod_fockdiag_(lfield, &Wfn.wfn[0],  &Wfn.wfn[Wfn.size1], &Eig[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::mkrrad(const clwfn& Wfn, std::vector<dcomplex>& rrad, std::vector<dcomplex>& rradpw)
{
  hprod_mkrrad_(&LZERO, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &rrad[0], &rradpw[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::mkrradx(const clwfn& Wfn, std::vector<dcomplex>& rrad, std::vector<dcomplex>& rradpw)
{
  if (clcontrol::rrad_type == -1) {
    hprod_mkrrad_(&LZERO, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &rrad[0], &rradpw[0]);
  } else if (clcontrol::rrad_type == 0) {
    hprod_mkrrad0_(&LZERO, &clcontrol::radipx, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &rrad[0], &rradpw[0]);
  } else if (clcontrol::rrad_type == 1) {
    hprod_mkrrad1_(&LZERO, &clcontrol::radipx, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &rrad[0], &rradpw[0]);
  } else if (clcontrol::rrad_type == 2) {
    hprod_mkrrad2_(&LZERO, &clcontrol::rrad_rion, &clcontrol::rrad_r2in, &clcontrol::rrad_r2out, 
		   &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &rrad[0], &rradpw[0]);
  } else if (clcontrol::rrad_type == 3) {
    hprod_mkrrad3_(&LZERO, &clcontrol::rrad_rion, &clcontrol::rrad_r2in, &clcontrol::rrad_r2out, 
		   &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &rrad[0], &rradpw[0]);
  }
}
//////////////////////////////////////////////////////////////////////////
void clhprod::mkrrad0(const clwfn& Wfn, std::vector<dcomplex>& rrad, std::vector<dcomplex>& rradpw)
{
  hprod_mkrrad0_(&LZERO, &clcontrol::radipx, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &rrad[0], &rradpw[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::mkrrad1(const clwfn& Wfn, std::vector<dcomplex>& rrad, std::vector<dcomplex>& rradpw)
{
  hprod_mkrrad1_(&LZERO, &clcontrol::radipx, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &rrad[0], &rradpw[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::spin(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn, typhys& Phys)
{
  double tsz, ts2;
  ormas_spin_(&Wfn.wfn[Wfn.size1], &tsz, &ts2);
  Phys.sz = tsz;
  Phys.s2 = ts2;
}
//////////////////////////////////////////////////////////////////////////
void clhprod::oang(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn, typhys& Phys)
{
  double tlz, tl2;
  hprod_oang_(&Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &tlz, &tl2);
  Phys.lz = tlz;
  Phys.l2 = tl2;
}
//////////////////////////////////////////////////////////////////////////
void clhprod::orbene(const clmpi& MPIP, const clbas& Bas, const double* lfield,
		     const clwfn& Wfn, std::vector<double>& Eig)
{
  hprod_orbene_(lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &Eig[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::cidiag(const clmpi& MPIP, const clbas& Bas, clwfn& Wfn)
{
  hprod_cidiag_(&Wfn.wfn[0], &Wfn.wfn[Wfn.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::dipole(const clmpi& MPIP, const clbas& Bas, const double* lfield,
		     const clwfn& Wfn, typhys& Phys)
{
  hprod_dipole_(lfield, &Wfn.wfn[0],  &Wfn.wfn[Wfn.size1], Phys.dip, Phys.vel, Phys.acc);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::dipole_ipx(const clmpi& MPIP, const clbas& Bas, const double* lfield,
			 const clwfn& Wfn, typhys* OpIpx)
{
  double *dip, *vel, *acc;
  dip = new double [clcontrol::maxipx+1];
  vel = new double [clcontrol::maxipx+1];
  acc = new double [clcontrol::maxipx+1];
  hprod_dipole_ipx_(&clcontrol::maxipx, &clcontrol::radipx, lfield, &Wfn.wfn[0],
		    &Wfn.wfn[Wfn.size1], dip, vel, acc);
  for (int iipx = 0; iipx <= clcontrol::maxipx; iipx ++) {
    OpIpx[iipx].dip[2] = dip[iipx];
    OpIpx[iipx].vel[2] = vel[iipx];
    OpIpx[iipx].acc[2] = acc[iipx];
  }
  delete [] acc;
  delete [] vel;
  delete [] dip;
}
//////////////////////////////////////////////////////////////////////////
void clhprod::dipole_ipd(const clmpi& MPIP, const clbas& Bas, const double* lfield,
			 const clwfn& Wfn, typhys* OpIpx)
{
  double *dip, *vel, *acc;
  dip = new double [clcontrol::maxipd+1];
  vel = new double [clcontrol::maxipd+1];
  acc = new double [clcontrol::maxipd+1];
  hprod_dipole_ipd_(&clcontrol::maxipd, &clcontrol::radipx, lfield, &Wfn.wfn[0],
		    &Wfn.wfn[Wfn.size1], dip, vel, acc);
  for (int iipx = 0; iipx <= clcontrol::maxipd; iipx ++) {
    OpIpx[iipx].dip[2] = dip[iipx];
    OpIpx[iipx].vel[2] = vel[iipx];
    OpIpx[iipx].acc[2] = acc[iipx];
  }
  delete [] acc;
  delete [] vel;
  delete [] dip;
}
//////////////////////////////////////////////////////////////////////////
void clhprod::enepole(const clmpi& MPIP, const clbas& Bas, const double* lfield, 
		      const clwfn& Wfn, typhys& Phys)
{
  hprod_enepole_(lfield, &Wfn.wfn[0],  &Wfn.wfn[Wfn.size1]);
  Phys.dip[2] = dip_exp;
  Phys.vel[2] = vel_exp;
  Phys.acc[2] = acc_exp;
}
//////////////////////////////////////////////////////////////////////////
void clhprod::normx(const clmpi& MPIP, const clbas& Bas, const double* lfield,
		    const clwfn& Wfn, typhys& Phys)
{
  hprod_ipx_(&LZERO, &Bas.GRad.xrad[Bas.GRad.nrad], lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &Phys.ipx[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::ionpx(const clmpi& MPIP, const clbas& Bas, const double* lfield,
		    const clwfn& Wfn, typhys& Phys)
{
  hprod_ipx_(&clcontrol::maxipx, &clcontrol::radipx, lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &Phys.ipx[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::ionpd(const clmpi& MPIP, const clbas& Bas, const double* lfield,
		    const clwfn& Wfn, typhys& Phys)
{
  hprod_ipd_(&clcontrol::maxipd, &clcontrol::radipx, lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &Phys.ipx[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::getno(const clmpi& MPIP, const clio& IO, const clbas& Bas,
		    const clwfn& Wfn, clwfn& tWfn)
{
  hprod_norb_(&Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &tWfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::kickk(const clmpi& MPIP, const clio& IO, const clbas& Bas, clwfn& Wfn, double knorm) const
{
  hprod_kickk_(&knorm, &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::getden1(std::vector<dcomplex>& den1) const
{
  hprod_getden1_(&den1[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::getint1(std::vector<dcomplex>& int1) const
{
  hprod_getint1_(&int1[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::getint2(std::vector<dcomplex>& int2) const
{
  hprod_getint2_(&int2[0]);
}
////////////////////////////////////////////////////////////////////////
void clhprod::op1tr_init(const clio& IO, const clbas& Bas, const clwfn& Wfn)
{
  IO.read_info("op1tr_nfun", Bas.ORMAS.nfun, op1tr_nfun);

  dQ.resize(op1tr_nfun);
  vQ.resize(op1tr_nfun);
  aQ.resize(op1tr_nfun);
  dP.resize(op1tr_nfun);
  vP.resize(op1tr_nfun);
  aP.resize(op1tr_nfun);
  vdP.resize(op1tr_nfun*op1tr_nfun);
  vvP.resize(op1tr_nfun*op1tr_nfun);
  vaP.resize(op1tr_nfun*op1tr_nfun);
  for (int ifun = 0; ifun < op1tr_nfun; ifun ++) {
    dP[ifun] = &vdP[op1tr_nfun*ifun];
    vP[ifun] = &vvP[op1tr_nfun*ifun];
    aP[ifun] = &vaP[op1tr_nfun*ifun];
  }
  hprod_op1tr_init_(&op1tr_nfun,&Wfn.wfn[0],&vdP[0],&vvP[0],&vaP[0],&dQ[0],&vQ[0],&aQ[0]);
}
////////////////////////////////////////////////////////////////////////
void clhprod::op1tr_final()
{
  hprod_op1tr_final_();
}
////////////////////////////////////////////////////////////////////////
void clhprod::op1tr(const double* lfield, typhys& Phys) const
{
  hprod_op1tr_(lfield, Phys.dip, Phys.vel, Phys.acc);
}
////////////////////////////////////////////////////////////////////////
void clhprod::op1tr_print(const clio& IO, const clbas& Bas, const clfield& Field, 
			  const double* lfield, const typhys& Phys) const
{
  fprintf(IO.fp_atrQ,"%10d %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e",
	  Field.step,Field.time,Field.ncyc(),lfield[5],lfield[8],Phys.acc[0],Phys.acc[1],Phys.acc[2]);
  for (int ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++){
    fprintf(IO.fp_atrQ,"%20.10e",aQ[ifun].real());
  }
  fprintf(IO.fp_atrQ,"\n");

  fprintf(IO.fp_atrP,"%10d %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e",
	  Field.step,Field.time,Field.ncyc(),lfield[5],lfield[8],Phys.acc[0],Phys.acc[1],Phys.acc[2]);
  for (int ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++){
    for (int jfun = ifun+1; jfun < Bas.ORMAS.nfun; jfun ++){
      if (Bas.lval[ifun] != Bas.lval[jfun] &&
	  Bas.mval[ifun] == Bas.mval[jfun]) {
	fprintf(IO.fp_atrP,"%20.10e",2*aP[ifun][jfun].real());
      }
    }
  }
  fprintf(IO.fp_atrP,"\n");
}
////////////////////////////////////////////////////////////////////////
void clhprod::op1tr_printp_tag(const clio& IO, const clbas& Bas) const
{
  fprintf(IO.fp_atrP,"# 1: istep\n");
  fprintf(IO.fp_atrP,"# 2: time/au\n");
  fprintf(IO.fp_atrP,"# 3: time/cyc\n");
  fprintf(IO.fp_atrP,"# 4: electric field\n");
  fprintf(IO.fp_atrP,"# 5: vector potential\n");
  fprintf(IO.fp_atrP,"# 6: nuclear force\n");
  fprintf(IO.fp_atrP,"# 7: correction\n");
  fprintf(IO.fp_atrP,"# 8: total acceleration\n");

  int col = 9;
  for (int ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++){
    for (int jfun = ifun+1; jfun < Bas.ORMAS.nfun; jfun ++){
      if (Bas.lval[ifun] != Bas.lval[jfun] &&
	  Bas.mval[ifun] == Bas.mval[jfun]) {
	fprintf(IO.fp_atrP,"# %5d: %5d <--> %5d\n",col,ifun,jfun);
	col ++;
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////
void clhprod::op1tr_printq_tag(const clio& IO, const clbas& Bas) const
{
  fprintf(IO.fp_atrQ,"#    1: istep\n");
  fprintf(IO.fp_atrQ,"#    2: time/au\n");
  fprintf(IO.fp_atrQ,"#    3: time/cyc\n");
  fprintf(IO.fp_atrQ,"#    4: electric field\n");
  fprintf(IO.fp_atrQ,"#    5: vector potential\n");
  fprintf(IO.fp_atrQ,"#    6: nuclear force\n");
  fprintf(IO.fp_atrQ,"#    7: correction\n");
  fprintf(IO.fp_atrQ,"#    8: total acceleration\n");

  int col = 9;
  for (int ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++){
    fprintf(IO.fp_atrQ,"# %5d: %5d\n",col,ifun);
    col ++;
  }
}
// tdcis-teramura
////////////////////////////////////////////////////////////////////////
void clhprod::phys_tdcis(const clmpi& MPIP, const clbas& Bas, const double* lfield,
			 const clwfn& Wfn, typhys& Phys, const double dt)
{
  hprod_phys_tdcis_(&dt, &clcontrol::radipx, lfield, &Wfn.wfn[0], 
		    &Wfn.wfn[Wfn.size1], Phys.dip, Phys.vel, Phys.acc);
}
////////////////////////////////////////////////////////////////////////
void clhprod::orbin_tdcis_init()
{
  hprod_orbin_tdcis_init_();
}
////////////////////////////////////////////////////////////////////////
// tdcis-teramura
