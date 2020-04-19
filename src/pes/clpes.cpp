////////////////////////////////////////////////////////////////////////
// Photoelectron spectra
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
#include "gsl/gsl_specfunc.h"
////////////////////////////////////////////////////////////////////////
clpes::clpes()
{
  std::cout << "clpes" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clpes::~clpes()
{
  std::cout << "~clpes" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clpes::clpes(const clmpi& MPIP, const clio& IO, 
	     const clbas& Bas, const clhprod& HP)
{
  std::cout << "clpes" << std::endl;
  gen(MPIP, IO, Bas, HP);
}
////////////////////////////////////////////////////////////////////////
void clpes::gen(const clmpi& MPIP, const clio& IO, 
		const clbas& Bas, const clhprod& HP)
{
  IO.read_info("pes_r_min", 0e0, pes_r_min);
  IO.read_info("pes_k_min", 0e0, pes_k_min);
  IO.read_info("pes_k_max", 10e0, pes_k_max);
  IO.read_info("pes_numk", 1000, pes_numk);
  pes_k_step = (pes_k_max - pes_k_min) / pes_numk;

  pes_llr = Bas.GRad.get_irad(pes_r_min);
  pes_ulr = Bas.GRad.nrad - 1;

  pes_psik.resize((pes_numk + 1) * (Bas.GAng.lmax1 + 1) * Bas.ORMAS.nfun);
  pes_rhok.resize((pes_numk + 1) * Bas.ORMAS.nfun * Bas.ORMAS.nfun);
  pes_bess.resize((Bas.GAng.lmax1 + 1) * (pes_numk + 1) * (Bas.GRad.nrad - 1));
  gen_bess(MPIP, Bas, HP);
  pes_bind_(&pes_numk, &pes_llr, &pes_ulr, &pes_k_min, &pes_k_max, &pes_k_step, 
	    &pes_bess[0], &pes_psik[0], &pes_rhok[0]);
}
////////////////////////////////////////////////////////////////////////
void clpes::gen_bess(const clmpi& MPIP, const clbas& Bas, const clhprod& HP)
{
  //
  // pess_bess[l][ik][irad] = bessel_function[l][ik][irad] * rval[irad] * sqrt(wval[irad])
  //
  int lkr;
  int nrad = Bas.GRad.nrad;
  int lmax = Bas.GAng.lmax1;
  //  std::vector<double> bessel(lmax + 1);
  double bessel[MPIP.nthr][lmax + 1];
  double rval;  // radial coordinate at irad
  double wval;  // radial weight at irad
  double kval;  // momentum value at ik 
  double rwval; // rval * wval
  double krval; // rval * kval

  int ithr, ik, irad, l;
  #pragma omp parallel for private(ithr, ik, irad, l, kval, rval, wval, rwval, krval, lkr)
  for (ik = 0; ik <= pes_numk; ik ++) {
    ithr = omp_get_thread_num();
    kval = pes_k_min + ik * pes_k_step;
    for (irad = pes_llr; irad < pes_ulr; irad ++) {
      rval = Bas.GRad.xrad[irad];
      wval = Bas.GRad.wrad[irad];
      rwval = rval * sqrt(wval);
      krval = kval * rval;

      gsl_sf_bessel_jl_array(lmax, krval, &bessel[ithr][0]);
      for (l = 0; l <= lmax; l ++) {
	lkr = l * (pes_numk + 1) * (nrad - 1) 
  	                    + ik * (nrad - 1)
     	                          + irad - 1;
	pes_bess[lkr] = bessel[ithr][l] * rwval;
      }
    }
  }

  //DEBUG
  //int ik = 0;
  //kval = pes_k_min + ik * pes_k_step;
  //printf("clpes::gen_bess: bessel functions for k = %20.10f\n", kval);
  //for (int irad = pes_llr; irad <= pes_ulr; irad ++) {
  //  rval = Bas.GRad.xrad[irad];
  //  wval = Bas.GRad.wrad[irad];
  //  rwval = rval * sqrt(wval);
  //  printf("%20.10f", rval);
  //  for (int l = 0; l <= lmax; l ++) {
  //    lkr = l * (pes_numk + 1) * (nrad - 1) 
  //	                  + ik * (nrad - 1)
  //     	                        + irad - 1;
  //    bessel[l] = pes_bess[lkr] / rwval;
  //    printf("%20.10f", bessel[l]);
  //  }
  //  printf("\n");
  //}
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void clpes::spec1_k(const clmpi& MPIP, const clio& IO, 
		    const clbas& Bas, const clhprod& HP, const clwfn& Wfn) const
{
  pes_spec1_k_(&Wfn.wfn[0], &Wfn.wfn[Wfn.size1]);
}
////////////////////////////////////////////////////////////////////////
void clpes::spec1_kz(const clmpi& MPIP, const clio& IO, 
		    const clbas& Bas, const clhprod& HP, const clwfn& Wfn) const
{
  pes_spec1_kz_(&Wfn.wfn[0], &Wfn.wfn[Wfn.size1]);
}
////////////////////////////////////////////////////////////////////////
