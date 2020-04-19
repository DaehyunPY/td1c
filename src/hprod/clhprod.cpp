////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
// static members
////////////////////////////////////////////////////////////////////////
long clhprod::num_hprod = 0;
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
void clhprod::dipole_ion(const clmpi& MPIP, const clbas& Bas, const double* lfield,
			 const clwfn& Wfn, typhys& Phys)
{
  hprod_dipole_ion_(&clcontrol::radipx, lfield, &Wfn.wfn[0],
		    &Wfn.wfn[Wfn.size1], Phys.dip, Phys.vel, Phys.acc);
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
  hprod_ipd_(&clcontrol::maxipx, &clcontrol::radipx, lfield, &Wfn.wfn[0], &Wfn.wfn[Wfn.size1], &Phys.ipx[0]);
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
