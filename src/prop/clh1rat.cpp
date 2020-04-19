////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
clh1rat::clh1rat()
{
}
////////////////////////////////////////////////////////////////////////
clh1rat::~clh1rat()
{
}
////////////////////////////////////////////////////////////////////////
void clh1rat::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		  double dt,
		  int ci_type,
		  int dimn,
		  int dimd,
		  std::vector<dcomplex>& cn,
		  std::vector<dcomplex>& cd0,
		  std::vector<dcomplex>& cd1)
{
  dtime = dt;
  h1ci_type = ci_type;
  dim_numer = dimn;
  dim_denom = dimd;
  ncoeff = cn;
  dcoeff0 = cd0;
  dcoeff1 = cd1;

  IO.read_info("h1rat_prodci_type", 2, prodci_type);
  if (prodci_type < 0 || prodci_type > 2) abort();

  int len_piv = (Bas.GRad.nrad - 1) * (Bas.GAng.lmax1 + 1);
  int len_inv = (3 * Bas.GRad.nmax + 1) * (Bas.GRad.nrad - 1) * (Bas.GAng.lmax1 + 1);
  tpiv.resize(dim_denom * len_piv);
  tinv.resize(dim_denom * len_inv);
  h1rat_gen_(&dim_denom, &dcoeff0[0], &dcoeff1[0], &tinv[0], &tpiv[0]);

  numer_limit = ncoeff[0];
  denom_limit = RUNIT;
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    denom_limit *= dcoeff0[idenom];
  }
  coeff_limit = numer_limit / denom_limit;
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod(const clmpi& Proc, const clbas& Bas, double time,
		   const clfield& Field, clhprod& HPW, clwfn& Wfn) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  if (clcontrol::split_type == 0) {
    HPW.scalo(Proc, Bas, coeff_limit, Wfn);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_(&CTRUE, &CTRUE, lfield, &dtime, &dim_numer, &dim_denom, 
		&ncoeff[0], &dcoeff1[0], &tpiv[0], &tinv[0], &Wfn.wfn[0]);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac,
		    const clfield& Field, clhprod& HPW, clwfn& Wfn) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  if (h1ci_type == -1) {
    HPW.scalc(Proc, Bas, zfac*coeff_limit, Wfn);
  } else {
    std::cout << "clh1rat::prod: bad h1ci_type." << std::endl;
    abort();
  }

  if (clcontrol::split_type == 0) {
    HPW.scalo(Proc, Bas, zfac*coeff_limit, Wfn);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_(&CTRUE, &CTRUE, lfield, &dtime, &dim_numer, &dim_denom, 
		&ncoeff[0], &dcoeff1[0], &tpiv[0], &tinv[0], &Wfn.wfn[0]);
    HPW.scalo(Proc, Bas, zfac, Wfn);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod2(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac,
		    const clfield& Field, clhprod& HPW, clwfn& WIn, clwfn& WOut) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  if (h1ci_type == -1) {
    HPW.axpyc(Proc, Bas, zfac*coeff_limit, WIn, WOut);
  } else {
    std::cout << "clh1rat::prod2 bad h1ci_type." << std::endl;
    abort();
  }

  if (clcontrol::split_type == 0) {
    HPW.axpyo(Proc, Bas, zfac*coeff_limit, WIn, WOut);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_(&CTRUE, &CTRUE, lfield, &dtime, &dim_numer, &dim_denom, 
		&ncoeff[0], &dcoeff1[0], &tpiv[0], &tinv[0], &WIn.wfn[0]);
    HPW.axpyo(Proc, Bas, zfac, WIn, WOut);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod2(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac,
		    const clfield& Field, clhprod& HPW, const std::vector<dcomplex>& hDiag,
		    const std::vector<dcomplex>& Int1e, const std::vector<dcomplex>& Int2e, 
		    clwfn& WIn, clwfn& WOut) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  double Eref = ZERO;
  if (h1ci_type == -1) {
    HPW.axpyc(Proc, Bas, zfac*coeff_limit, WIn, WOut);
  } else if (h1ci_type == 0) {
//old    h1rat_prodcid_(&CTRUE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
//old		       &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &Eref, &hDiag[0], 
//old		       &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
  } else {
    if (prodci_type == 0) {
      h1rat_prodci_(&CTRUE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		    &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		    &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 1) {
      h1rat_prodci1_(&CTRUE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 2) {
      h1rat_prodci2_(&CTRUE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    }
  }

  if (clcontrol::split_type == 0) {
    HPW.axpyo(Proc, Bas, zfac*coeff_limit, WIn, WOut);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_(&CTRUE, &CTRUE, lfield, &dtime, &dim_numer, &dim_denom, &ncoeff[0], &dcoeff1[0], &tpiv[0], &tinv[0], &WIn.wfn[0]);
    HPW.axpyo(Proc, Bas, zfac, WIn, WOut);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod2(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac,
		    const clfield& Field, clhprod& HPW, double Eref, const std::vector<dcomplex>& hDiag,
		    const std::vector<dcomplex>& Int1e, const std::vector<dcomplex>& Int2e, 
		    clwfn& WIn, clwfn& WOut) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  if (h1ci_type == -1) {
    HPW.axpyc(Proc, Bas, zfac*coeff_limit, WIn, WOut);
  } else if (h1ci_type == 0) {
//old    h1rat_prodcid_(&CTRUE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
//old		       &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &Eref, &hDiag[0], 
//old		       &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
  } else {
    if (prodci_type == 0) {
      h1rat_prodci_(&CTRUE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		    &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		    &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 1) {
      h1rat_prodci1_(&CTRUE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 2) {
      h1rat_prodci2_(&CTRUE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    }
  }

  if (clcontrol::split_type == 0) {
    HPW.axpyo(Proc, Bas, zfac*coeff_limit, WIn, WOut);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_(&CTRUE, &CTRUE, lfield, &dtime, &dim_numer, &dim_denom, &ncoeff[0], &dcoeff1[0], &tpiv[0], &tinv[0], &WIn.wfn[0]);
    HPW.axpyo(Proc, Bas, zfac, WIn, WOut);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod_numer(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac, 
			 const clfield& Field, clhprod& HPW, const clwfn& WIn, clwfn& WOut) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  if (h1ci_type == -1) {
    HPW.axpyc(Proc, Bas, zfac*numer_limit, WIn, WOut);
  } else {
    std::cout << "clh1rat::prod2 bad h1ci_type." << std::endl;
    abort();
  }

  if (clcontrol::split_type == 0) {
    HPW.axpyo(Proc, Bas, zfac*numer_limit, WIn, WOut);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_numer_(lfield, &dtime, &dim_numer, &zfac, &ncoeff[0], &WIn.wfn[0], &WOut.wfn[0]);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod_numer(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac, 
			 const clfield& Field, clhprod& HPW, const std::vector<dcomplex>& hDiag,
			 const std::vector<dcomplex>& Int1e, const std::vector<dcomplex>& Int2e, 
			 const clwfn& WIn, clwfn& WOut) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  double Eref = ZERO;
  if (h1ci_type == -1) {
    HPW.axpyc(Proc, Bas, zfac*numer_limit, WIn, WOut);
  } else if (h1ci_type == 0) {
//old    h1rat_prodcid_(&CTRUE, &CFALSE, &zfac, &dtime, &dim_numer, &dim_denom, 
//old		       &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &Eref, &hDiag[0], 
//old		       &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
  } else {
    if (prodci_type == 0) {
      h1rat_prodci_(&CTRUE, &CFALSE, &zfac, &dtime, &dim_numer, &dim_denom, 
		    &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		    &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 1) {
      h1rat_prodci1_(&CTRUE, &CFALSE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 2) {
      h1rat_prodci2_(&CTRUE, &CFALSE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    }
  }

  if (clcontrol::split_type == 0) {
    HPW.axpyo(Proc, Bas, zfac*numer_limit, WIn, WOut);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_numer_(lfield, &dtime, &dim_numer, &zfac, &ncoeff[0], &WIn.wfn[0], &WOut.wfn[0]);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod_numer(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac, 
			 const clfield& Field, clhprod& HPW, double Eref, const std::vector<dcomplex>& hDiag,
			 const std::vector<dcomplex>& Int1e, const std::vector<dcomplex>& Int2e, 
			 const clwfn& WIn, clwfn& WOut) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  if (h1ci_type == -1) {
    HPW.axpyc(Proc, Bas, zfac*numer_limit, WIn, WOut);
  } else if (h1ci_type == 0) {
//old    h1rat_prodcid_(&CTRUE, &CFALSE, &zfac, &dtime, &dim_numer, &dim_denom, 
//old		       &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &Eref, &hDiag[0], 
//old		       &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
  } else {
    if (prodci_type == 0) {
      h1rat_prodci_(&CTRUE, &CFALSE, &zfac, &dtime, &dim_numer, &dim_denom, 
		    &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		    &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 1) {
      h1rat_prodci1_(&CTRUE, &CFALSE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 2) {
      h1rat_prodci2_(&CTRUE, &CFALSE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    }
  }

  if (clcontrol::split_type == 0) {
    HPW.axpyo(Proc, Bas, zfac*numer_limit, WIn, WOut);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_numer_(lfield, &dtime, &dim_numer, &zfac, &ncoeff[0], &WIn.wfn[0], &WOut.wfn[0]);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod_denom(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac, 
			 const clfield& Field, clhprod& HPW, clwfn& WIn, clwfn& WOut) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  if (h1ci_type == -1) {  
    HPW.axpyc(Proc, Bas, zfac/denom_limit, WIn, WOut);
  } else {
    std::cout << "clh1rat::prod2 bad h1ci_type." << std::endl;
    abort();
  }

  if (clcontrol::split_type == 0) {
    HPW.axpyo(Proc, Bas, zfac/denom_limit, WIn, WOut);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_(&CFALSE, &CTRUE, lfield, &dtime, &dim_numer, &dim_denom, &ncoeff[0], &dcoeff1[0], &tpiv[0], &tinv[0], &WIn.wfn[0]);
    HPW.axpyo(Proc, Bas, zfac, WIn, WOut);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod_denom(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac, 
			 const clfield& Field, clhprod& HPW, const std::vector<dcomplex>& hDiag,
			 const std::vector<dcomplex>& Int1e, const std::vector<dcomplex>& Int2e, 
			 clwfn& WIn, clwfn& WOut) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  double Eref = ZERO;
  if (h1ci_type == -1) {  
    HPW.axpyc(Proc, Bas, zfac/denom_limit, WIn, WOut);
  } else if (h1ci_type == 0) {
//old    h1rat_prodcid_(&CFALSE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
//old		       &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &Eref, &hDiag[0], 
//old		       &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
  } else {
    if (prodci_type == 0) {
      h1rat_prodci_(&CFALSE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		    &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		    &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 1) {
      h1rat_prodci1_(&CFALSE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 2) {
      h1rat_prodci2_(&CFALSE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    }
  }

  if (clcontrol::split_type == 0) {
    HPW.axpyo(Proc, Bas, zfac/denom_limit, WIn, WOut);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_(&CFALSE, &CTRUE, lfield, &dtime, &dim_numer, &dim_denom, &ncoeff[0], &dcoeff1[0], &tpiv[0], &tinv[0], &WIn.wfn[0]);
    HPW.axpyo(Proc, Bas, zfac, WIn, WOut);
  }
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::prod_denom(const clmpi& Proc, const clbas& Bas, double time, dcomplex zfac, 
			 const clfield& Field, clhprod& HPW, double Eref, const std::vector<dcomplex>& hDiag,
			 const std::vector<dcomplex>& Int1e, const std::vector<dcomplex>& Int2e, 
			 clwfn& WIn, clwfn& WOut) const
{
  if (clcontrol::iprojfc != 2) {
    std::cout << "clh1rat::prod: projfc is nyi." << std::endl;
    abort();
  }

  if (h1ci_type == -1) {  
    HPW.axpyc(Proc, Bas, zfac/denom_limit, WIn, WOut);
  } else if (h1ci_type == 0) {
//old    h1rat_prodcid_(&CFALSE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
//old		       &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &Eref, &hDiag[0], 
//old		       &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
  } else {
    if (prodci_type == 0) {
      h1rat_prodci_(&CFALSE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		    &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		    &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 1) {
      h1rat_prodci1_(&CFALSE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    } else if (prodci_type == 2) {
      h1rat_prodci2_(&CFALSE, &CTRUE, &zfac, &dtime, &dim_numer, &dim_denom, 
		     &ncoeff[0], &dcoeff0[0], &dcoeff1[0], &h1ci_type, &Eref, &hDiag[0], 
		     &Int1e[0], &Int2e[0], &WIn.wfn[WIn.size1], &WOut.wfn[WIn.size1]);
    }
  }

  if (clcontrol::split_type == 0) {
    HPW.axpyo(Proc, Bas, zfac/denom_limit, WIn, WOut);
  } else {
    double lfield[9];
    Field.get_value(time, lfield);
    h1rat_prod_(&CFALSE, &CTRUE, lfield, &dtime, &dim_numer, &dim_denom, &ncoeff[0], &dcoeff1[0], &tpiv[0], &tinv[0], &WIn.wfn[0]);
    HPW.axpyo(Proc, Bas, zfac, WIn, WOut);
  }
}
//////////////////////////////////////////////////////////////////////////
dcomplex clh1rat::get_numer_limit() const
{
  return numer_limit;
}
//////////////////////////////////////////////////////////////////////////
dcomplex clh1rat::get_denom_limit() const
{
  return denom_limit;
}
//////////////////////////////////////////////////////////////////////////
dcomplex clh1rat::get_coeff_limit() const
{
  return coeff_limit;
}
//////////////////////////////////////////////////////////////////////////
void clh1rat::print() const
{
  printf("# clh1rat::ncoeff %5d\n", dim_numer);
  for (int inumer = 0; inumer <= dim_numer; inumer ++) {
    printf("%5d%20.10e%20.10e\n", inumer, 
  	   real(ncoeff[inumer]), imag(ncoeff[inumer]));
  }
  printf("# clh1rat::dcoeff %5d\n", dim_denom);
  for (int idenom = 0; idenom < dim_denom; idenom ++) {
    printf("%5d%20.10e%20.10e%20.10e%20.10e\n", idenom, 
  	   real(dcoeff0[idenom]), imag(dcoeff0[idenom]),
  	   real(dcoeff1[idenom]), imag(dcoeff1[idenom]));
  }
  printf("# clh1rat. limit of dt-->0 (numer) = %20.10f%20.10f\n", 
	 real(numer_limit), imag(numer_limit));
  printf("# clh1rat. limit of dt-->0 (denom) = %20.10f%20.10f\n", 
	 real(denom_limit), imag(denom_limit));
  printf("# clh1rat. limit of dt-->0 = %20.10f%20.10f\n", 
	 real(coeff_limit), imag(coeff_limit));
}
//////////////////////////////////////////////////////////////////////////
