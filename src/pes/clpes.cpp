////////////////////////////////////////////////////////////////////////
// Photoelectron spectra
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
#include "gsl/gsl_specfunc.h"
#include "radau_alglib.hpp"
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

  krad.resize(pes_numk + 1);
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

  long ithr, ik, irad, l;
  #pragma omp parallel for private(ithr, ik, irad, l, kval, rval, wval, rwval, krval, lkr)
  for (ik = 0; ik <= pes_numk; ik ++) {
    ithr = (long) omp_get_thread_num();
    kval = pes_k_min + ik * pes_k_step;
    krad[ik] = kval;    
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
////////////////////////////////////////////////////////////////////////
void clpes::spec1_k_angres(const clmpi& MPIP, const clio& IO, 
		    const clbas& Bas, const clhprod& HP, const clwfn& Wfn) const
{


  // calculate < j_kl | phi_n>
  int lnum = Bas.GAng.lmax1 + 1;

  std::vector<dcomplex> psikl((pes_numk+1) * lnum * Bas.ORMAS.nfun);

  pes_spec1_k_psik_(&Wfn.wfn[0], &psikl[0]);


  // get angular parts of spehrical harmonics
  long pes_num_kang;
  IO.read_info("pes_num_kang", pes_num_kang);
  std::vector<double> kang(pes_num_kang);
  std::vector<double> wang(pes_num_kang);
  //double kaunit = PI / (pes_num_kang-1);
  //for(int i = 0; i < pes_num_kang; i++){
  //  kang[i] = i * kaunit;
  //  wang[i] = sin(kang[ikang]) * kaunit;
  //}



  {  
    //function : init_kang
    //--------------------------------------------
    // taken from surff class
    //--------------------------------------------

    int num_kang = pes_num_kang;
    kang.resize(num_kang);
    wang.resize(num_kang);
    
    int N = num_kang; 
    
    double left = -1.0; // left boundary
    double right = 1.0;
    double range = right - left; // range

    double mu0; //0th moment of weighting function in the integral region
    mu0 = range;

    alglib::real_1d_array alpha;
    alglib::real_1d_array beta;
    alpha.setlength(N-1);
    beta.setlength(N-1);   
    double di;
    for(int i = 0; i < N-1; ++i){
      di = static_cast<double>(i);
      alpha[i] = 0.5 * (2 * left + range);
      beta[i] = di*di / (4*di*di - 1) * std::pow( 0.5 * range, 2);
    }
    
    long info;
    alglib::real_1d_array x;
    alglib::real_1d_array w;
    x.setlength(N);
    w.setlength(N);
    alglib::gqgenerategausslobattorec(alpha, beta, mu0, left, right, N, info, x, w);
  
    std::cout << "# gauss_lobatto quadrature:" << std::endl;
    std::cout << "# cos(kang)     kang     wang" << std::endl;
    int i = 0;
    kang[i] = PI;
    wang[i] = w[i];
    std::cout << "      " << x[i] << " " << kang[i] << " " << wang[i] << std::endl;
    for(int i = 1; i < N-1; ++i){
      kang[i] = std::acos(x[i]);
      wang[i] = w[i];
      std::cout << "      " << x[i] << " " << kang[i] << " " << wang[i] << std::endl;
    }
    i = N-1;
    kang[i] = 0;
    wang[i] = w[i];
    std::cout << "      " << x[i] << " " << kang[i] << " " << wang[i] << std::endl;
  
  }


  long mnum = Bas.GAng.mmax1 + 1;
  std::vector<double> sph_ang(mnum * lnum * pes_num_kang); // P_lm(cos(theta_k))

  for(int mval = 0; mval < mnum; mval++){ 
    for(int l = mval; l < lnum; l++){
      for(int ikang = 0; ikang < pes_num_kang; ikang++){ 
	sph_ang[mval * lnum * pes_num_kang + l * pes_num_kang + ikang] =  
	  std::sqrt(TWO * PI) * gsl_sf_legendre_sphPlm(l, mval, std::cos(kang[ikang]));
      }
    }
  }  
  
  // calculate angular resolved pes
  std::vector<dcomplex> orbpes((pes_numk+1) *  pes_num_kang * Bas.ORMAS.nfun);
  std::vector<dcomplex> ipow(lnum);
  for(int l = 0; l < lnum; ++l){
    ipow[l] = pow(-IUNIT, l);
  }
  
  double invpi = 1.0 / std::pow(2*PI, 1.5);
  int mval = 0;
  for(int ifun = 0; ifun < Bas.ORMAS.nfun; ++ifun){
    mval = std::abs(Bas.mval[ifun]);
    for(int l = 0; l < lnum; ++l){
      for(int ikang = 0; ikang < pes_num_kang; ++ikang){
	for(int ik = 0; ik < (pes_numk+1); ++ik){
	  orbpes[ifun * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ik] +=
	    invpi * FOUR * PI * ipow[l] * sph_ang[mval * lnum * pes_num_kang + l * pes_num_kang + ikang]
	    * psikl[ifun * lnum * (pes_numk+1) + l * (pes_numk+1) + ik];
	}
      }
    }
  }



  //////////////////////////////////////////////////////////////////
  // physical value
  //////////////////////////////////////////////////////////////////

  int norb = Bas.ORMAS.nfun;

  std::vector<dcomplex> rhok(norb * norb * pes_num_kang * (pes_numk+1));

#pragma omp parallel default(shared)
  {
    long ithr = omp_get_thread_num();
    long llkang, ulkang;
    MPIP.omp_divide(ithr, 0, pes_num_kang, llkang, ulkang);

    for(int iorb = 0; iorb < norb; ++iorb){
      for(int jorb = 0; jorb < norb; ++jorb){
	if(Bas.mval[iorb] == Bas.mval[jorb]){
	  ////////////////////////
	  // omp parallelized
	  ////////////////////////
	  for(int ikang = llkang; ikang < ulkang; ++ikang){
	    for(int ikrad = 0; ikrad < pes_numk + 1; ++ikrad){
	      rhok[iorb * norb * pes_num_kang * (pes_numk+1) + jorb * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ikrad]
		= conj(orbpes[iorb * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ikrad]) 
		* orbpes[jorb * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ikrad];
	      //std::cout << "rhok_first " << rhok[iorb * norb * pes_num_kang * (pes_numk+1) + jorb * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ikrad]
	      //	<< std::endl;
	    }
	  }
	}
      }
    }
  }

  

  std::vector<dcomplex> den1mat(norb * norb, 0.0);
  ormas_mkden1_(&Wfn.wfn[Wfn.size1], &den1mat[0]);
  
  std::vector<dcomplex> ang_spec(pes_num_kang * (pes_numk+1));
#pragma omp parallel default(shared)
  {
    dcomplex den1;
    int ncore = Bas.ORMAS.ncore;
    
    long ithr = omp_get_thread_num();
    long llkang, ulkang;
    MPIP.omp_divide(ithr, 0, pes_num_kang, llkang, ulkang);
    
    for(int iorb = 0; iorb < ncore; ++iorb){
      ////////////////////////
      // omp parallelized
      ////////////////////////
      for(int ikang = llkang; ikang < ulkang; ++ikang){
	for(int ikrad = 0; ikrad < pes_numk+1; ++ikrad){
	  ang_spec[ikang * (pes_numk+1) + ikrad] += 
	    rhok[iorb * norb * pes_num_kang * (pes_numk+1) + iorb * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ikrad]
	    * 2.0;
	}
      }
    }
  }
  

  std::vector<double> mom_spec(pes_numk + 1);

#pragma omp parallel default(shared)
  {
    dcomplex den1;
    int ncore = Bas.ORMAS.ncore;
    
    long ithr = omp_get_thread_num();
    long llkang, ulkang;
    MPIP.omp_divide(ithr, 0, pes_num_kang, llkang, ulkang);
    
    for(int iorb = ncore; iorb < norb; ++iorb){
      for(int jorb = ncore; jorb < norb; ++jorb){

	den1 = den1mat[(iorb - ncore) * (norb - ncore) + (jorb - ncore)];
	////////////////////////
	// omp parallelized
	////////////////////////
	for(int ikang = llkang; ikang < ulkang; ++ikang){
	  for(int ikrad = 0; ikrad < pes_numk+1; ++ikrad){
	    ang_spec[ikang * (pes_numk+1) + ikrad] += 
	      rhok[iorb * norb * pes_num_kang * (pes_numk+1) + jorb * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ikrad]
	      * den1;
	  }
	}
      }
    }
  }

  for(int ikang = 0; ikang < pes_num_kang; ++ikang){
    for(int ikrad = 0; ikrad < pes_numk+1; ++ikrad){
      mom_spec[ikrad] += real(ang_spec[ikang * (pes_numk+1) + ikrad]) * wang[ikang];
    }
  }    


//   //calculate momentum and energy distribution for 1 electron
//   std::vector<double> mom_spec(pes_numk + 1);
//   double del_ang = kang[1] - kang[0];
//   for(int ikang = 0; ikang < pes_num_kang; ++ikang){
//     for(int ik = 0; ik < (pes_numk+1); ++ik){
//       mom_spec[ik] += (std::abs(orbpes[0 * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ik]
// 			       * orbpes[0 * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ik]) 
// 		       * sin(kang[ikang]) * del_ang);
//     }
//   }

  std::vector<double> ene_spec(pes_numk + 1);
  double small = 10E-10;

  for(int ikrad = 0; ikrad < pes_numk+1; ++ikrad){
    if(krad[ikrad] > small){
      ene_spec[ikrad] = mom_spec[ikrad] / krad[ikrad];
    }
    else {
      ene_spec[ikrad] = 0.0;
    }
  }


  // print angular distribution momentum and energy distribution
  std::string fname;
  fname = IO.name + ".d_ang_spec";
  FILE *fpo;
  fpo = fopen(fname.c_str(), "w");
  
  fprintf(fpo, "# krad kang real[spectrum amplitude (iorb)] imag[spectrum amplitude (iorb)] ... \n" );
  fprintf(fpo, "%ld\n%ld\n", pes_numk+1, pes_num_kang);
  double tmp;
  for(int ikrad = 0; ikrad  < (pes_numk+1); ikrad++){
    for(int ikang = 0; ikang < pes_num_kang; ikang++){
      fprintf(fpo, " %20.10E    %20.10E   ",  krad[ikrad], kang[ikang]);
      tmp = real(ang_spec[ikang * (pes_numk+1) + ikrad]);
      fprintf(fpo, " %20.10E ",  tmp);
      fprintf(fpo, "\n");
    }
  }
  fclose(fpo);



  // calculate norm of orbpes
  double norm = 0.0;
  double delk = pes_k_step;
  for(int ikrad = 0; ikrad < pes_numk+1; ++ikrad){
    norm += mom_spec[ikrad] * krad[ikrad] * krad[ikrad] * delk; 
  }
  std::cout << "# norm of orbpes = " << norm << std::endl;


  // print momentum and energy distribution
  fname = IO.name + ".d_me_spec";
  fpo = fopen(fname.c_str(), "w");
  
  fprintf(fpo, "# norm = %20.10E, k rho(k) E rho(E) \n", norm);
  double ene;
  double rho;
  for(int ikrad = 0; ikrad  < pes_numk+1; ikrad++){
    ene = HALF * krad[ikrad] * krad[ikrad];
    fprintf(fpo, " %20.10E    %20.10E   %20.10E    %20.10E", krad[ikrad], mom_spec[ikrad], ene, ene_spec[ikrad]);
    fprintf(fpo, "\n");
  }
  fclose(fpo);


  // print surff
  fname = IO.name + ".d_surff";
  fpo = fopen(fname.c_str(), "w");
  
  fprintf(fpo, "# krad kang real[spectrum amplitude (iorb)] imag[spectrum amplitude (iorb)] ... \n" );
  fprintf(fpo, "%ld\n%ld\n", pes_numk+1, pes_num_kang);
  dcomplex dtmp;
  for(int ikrad = 0; ikrad  < (pes_numk+1); ikrad++){
    for(int ikang = 0; ikang < pes_num_kang; ikang++){
      fprintf(fpo, " %20.10E    %20.10E   ",  krad[ikrad], kang[ikang]);
      for(long iorb = 0; iorb < norb; ++iorb){
	dtmp = orbpes[iorb * pes_num_kang * (pes_numk+1) + ikang * (pes_numk+1) + ikrad];
	fprintf(fpo, " %20.10E    %20.10E   ",  dtmp.real(), dtmp.imag());
      }
      fprintf(fpo, "\n");
    }
  }
  fclose(fpo);

  
  return;
  
  

}
////////////////////////////////////////////////////////////////////////




