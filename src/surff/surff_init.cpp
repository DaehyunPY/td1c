#include "td1c.hpp"
#include "wrapper.hpp"
#include "surff.hpp"

#include "gsl/gsl_specfunc.h"
#include <dirent.h>
#include <sys/stat.h>
#include "radau_alglib.hpp"

//////////////////////////////////////////////////////////////////////////
Surff::Surff(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
// Sato_tSURFF
  //  IO.read_info("nprint_surff", 0, nprint);
  //  if(nprint <= 0){return;}
  if (! clcontrol::tsurff) {
    std::cout << "don't instantiate Surff when tsurff is false." << std::endl;
    abort();
  }
  IO.read_info("nprint_surff", 0, nprint);
// Sato_tSURFF

  init(MPIP, IO, Bas);
  return;
}
//////////////////////////////////////////////////////////////////////////
Surff::~Surff(){
// Sato_tSURFF
//DEBUG  if(nprint <= 0){return;}
// Sato_tSURFF
  delete integ;
}
//////////////////////////////////////////////////////////////////////////
void Surff::init(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  struct stat buf;
  int ret;
  char dir[256];
// Sato_TSURFF
//  char mkdir[512];
  char command[256];
// Sato_TSURFF


// Sato_TSURFF
  //  snprintf(dir,256,"./SPECTRUM");
  sprintf(dir,"./%s.SPECTRUM",IO.name.c_str());
  ret=stat(dir, &buf);
  if(ret != 0) {
    sprintf(command,"mkdir %s",dir);
    system(command); 
  }
  else if(ret == 0) {
    sprintf(command,"rm -rf %s/*",dir);
    system(command); 
  }
  std::cout << command << std::endl;
// Sato_TSURFF

  IO.read_info("num_krad", 1, num_krad);
  krad.resize(num_krad);
  double emin, emax;
// Sato_tSURFF
  double emin_ev, emax_ev;
  IO.read_info("emin", emin_ev);
  IO.read_info("emax", emax_ev);
  emin = emin_ev/27.2113845;
  emax = emax_ev/27.2113845;
// Sato_tSURFF

  if(emin < 10E-7 ){
// Sato_tSURFF
//  double eunit = (emax - emin) / 27.2 / (num_krad-1);
    double eunit = (emax - emin)/(num_krad-1);
// Sato_tSURFF
    emin += eunit;
  }
// Sato_tSURFF
//double eunit = (emax - emin) / 27.2 / (num_krad-1);
  double eunit = (emax - emin)/(num_krad-1);
// Sato_tSURFF
  for(int i = 0; i < num_krad; i++){
    krad[i] = std::sqrt((emin + i * eunit) * 2);
  }
  mom_spec.resize(num_krad, 0.0);
  ene_spec.resize(num_krad, 0.0);

  IO.read_info("num_kang", 1, num_kang);  
//   kang.resize(num_kang);
//   wang.resize(num_kang);
//   double kaunit = PI / (num_kang-1);
//   for(int i = 0; i < num_kang; i++){
//     kang[i] = i * kaunit;
//     wang[i] = sin(kang[i]) * kaunit;
//   }
  init_kang(MPIP, IO);

  norb = Bas.ORMAS.nfun;
  lnum = Bas.GAng.lmax1+1;
  orb_surf.resize(norb * lnum);
  dorb_surf.resize(norb * lnum);

  ksize = num_krad * num_kang;
  orbpes.resize(ksize * norb, 0.0); // SPCETRUM : krad and kang

  //TEST  //TEST  //TEST  //TEST
  //orbpes.resize(ksize * norb, 1.0); // SPCETRUM : krad and kang
  //TEST  //TEST  //TEST  //TEST

  orbpes_dt.resize(ksize * norb, 0.0); // SPCETRUM : krad and kang
  torbpes.resize(ksize * norb, 0.0); // SPCETRUM : krad and kang

  ang_spec.resize(num_kang * num_krad);
  rhok.resize(norb * norb * num_kang * num_krad);


  mnum = Bas.GAng.mmax1 + 1;
  sph_ang.resize(mnum * lnum * num_kang); // P_lm(cos(theta_k))
  sbessel.resize(lnum * num_krad);
  sbessel_dr.resize(lnum * num_krad);
  alph_lm.resize(lnum * (2 * (mnum - 1) + 1));


  for(int mval = 0; mval < mnum; mval++){ 
    for(int l = mval; l < lnum; l++){
      for(int ikang = 0; ikang < num_kang; ikang++){ 
	sph_ang[mval * lnum * num_kang + l * num_kang + ikang] =  
	  std::sqrt(TWO * PI) * gsl_sf_legendre_sphPlm(l, mval, std::cos(kang[ikang]));
	//gsl_sf_legendre_Plm_array(lnum-1, mmin, std::cos(kang[ikang]), &ass_legendre[0]);
      }
    }
  }

  //DEBUG
  //old double delang = kang[1] - kang[0];
  double norm;
//   std::cout << std::endl;   
//   for(int mval = 0; mval < mnum; mval++){ 
//     for(int l = 0; l < lnum; l++){
//       for(int l2 = 0; l2 < lnum; l2++){
// 	norm = 0.0;
// 	for(int ikang = 0; ikang < num_kang; ikang++){ 	  
// 	  //old norm += sph_ang[mval * lnum * num_kang + l2 * num_kang + ikang] * sph_ang[mval * lnum * num_kang + l * num_kang + ikang] 
// 	  //old   * sin(kang[ikang]) * delang;
// 	  norm += sph_ang[mval * lnum * num_kang + l2 * num_kang + ikang] * sph_ang[mval * lnum * num_kang + l * num_kang + ikang] 
// 	    * wang[ikang];
// 	    }
// 	std::cout << l << " vs " << l2 <<" : norm = " << norm << std::endl;
//       }
//     }
//   }
//  abort();
  
//   std::cout << std::endl;    
//   for(int mval = 0; mval < mnum; mval++){ 
//     for(int ikang = 0; ikang < num_kang; ikang++){ 
//       std::cout << kang[ikang];
//       for(int l = mval; l < lnum; l++){
// 	std::cout << " " << sph_ang[mval * lnum * num_kang + l * num_kang + ikang];
//       }
//     std::cout << std::endl;
//     }
//     std::cout << std::endl;
//   }
//   abort();  
  //DEBUG
  
  //DEBUG
  //   int mval = 0;
  //   int l = 2;
  //   for(int ikang = 0; ikang < num_kang; ikang++){ 
  //     std::cout << kang[ikang] << "  " << sph_ang[mval * lnum * num_kang + l * num_kang + ikang] << std::endl;
  //   }
  //   abort();
  //DEBUG


  irad_surf = -1;
  double srad;
  IO.read_info("srad", srad);
  for(int irad = 0; irad <  Bas.GRad.nrad + 1; ++irad){
    if(std::abs(Bas.GRad.xrad[irad] - srad ) < 10E-10){
      irad_surf = irad;
      rrad_surf = Bas.GRad.xrad[irad];
      std::cout << "# irad_surf = " << irad_surf << std::endl;
      std::cout << "# rrad_surf = " << rrad_surf << std::endl;
      break;
    }
  }
  if(irad_surf < 0){
    std::cout << "bad srad" << std::endl;
    abort();
  }
  
  double sph_bess[lnum + 1];
  double sph_bess_minus1;
  double sph_bess_dr[lnum + 1];
  double krval;
  for(int ikrad = 0; ikrad < num_krad; ++ikrad){
    krval = krad[ikrad] * rrad_surf;
    gsl_sf_bessel_jl_array(lnum, krval, &sph_bess[0]);
    sph_bess_minus1 = -1.0 * sph_bess[1];
    
    //sph_bess_dr[0] = krad[ikrad] * HALF * (sph_bess_minus1 - (sph_bess[0] + krval * sph_bess[1]) / krval); //calc differenciation of j_0
    sph_bess_dr[0] = krad[ikrad] * ( -sph_bess[1]); //calc differenciation of j_0
    sbessel_dr[0 * num_krad + ikrad] = sph_bess_dr[0]; 
    sbessel[0 * num_krad + ikrad] = sph_bess[0];
    for(int l = 1; l < lnum; l++){
      //sph_bess_dr[l] = krad[ikrad] * HALF * (sph_bess[l-1] - (sph_bess[l] + krval * sph_bess[l+1]) / krval); //calc differenciation
      sph_bess_dr[l] = krad[ikrad] * (-sph_bess[l+1] + double(l) / krval * sph_bess[l]); //calc differenciation
      sbessel_dr[l * num_krad + ikrad] = sph_bess_dr[l];
      sbessel[l * num_krad + ikrad] = sph_bess[l];
    }
  }

  // DEBUG
//   std::cout << "debug "<<std::endl;
//   for(int ikrad = 0; ikrad < num_krad; ++ikrad){
//     std::cout << krad[ikrad] ;
//     for(int l = 0; l < lnum; ++l){
//       std::cout << " " << sbessel[l * num_krad + ikrad];
//     }  
//     std::cout << std::endl;
//   }
//   std::cout << std::endl;
  
//   for(int ikrad = 0; ikrad < num_krad; ++ikrad){
//     std::cout << krad[ikrad] ;
//     for(int l = 0; l < lnum; ++l){
//       std::cout << " " << sbessel_dr[l * num_krad + ikrad];
//     }  
//     std::cout << std::endl;
//   }
//   std::cout << std::endl;
//   abort();
// DEBUG   


  alph_lm = Bas.alph_lm;
  //DEBUG
//   for(int l = 0; l < lnum; l++){
//     printf("%5d", l);
//   }
//   std::cout << std::endl;
//   for(int mval = -mnum + 1; mval < mnum; mval++){
//     printf("%5d", mval);
//     int midx = mval + (mnum - 1);
//     for(int l = 0; l < lnum; l++){
//       printf("%20.10E", alph_lm[midx * lnum + l]);
//     }
//     std::cout << std::endl;
//   }
//   std::cout << "done" << std::endl;
//   abort();
  //DEBUG

  
  vphase.resize(ksize, 0.0);
  tdfac.resize(ksize, 0.0);

  v2xmat.resize(norb * norb, 0.0);

  //integral
  IO.read_info("integral_type", "EI", integ_type);

  if(integ_type == "trap"){
    integ = new Integ_trap(ksize, norb);
  }
  else if(integ_type == "simpson"){
    integ = new Integ_simpson;
  } 
  else if(integ_type == "CN"){
    integ = new Integ_CN(ksize, norb);
  } 
  else if(integ_type == "EI"){
    integ = new Integ_EI(ksize, norb);
    int iren;
    IO.read_info("renormarize", 0, iren);
    integ->set_renormarize(iren);
  } 
  else {
    std::cout << "bad integral_type. abort" << std::endl;
    abort();
  }
  
  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::init_kang(const clmpi&, const clio&)
{
  
  kang.resize(num_kang);
  wang.resize(num_kang);
  
  //The number of nodes including both boundary (>= ?)
  int N = num_kang; 
  
  double left = -1.0; // left boundary
  double right = 1.0;
  double range = right - left; // range
  
  // http://people.math.gatech.edu/~jeanbel/6580/orthogPol13.pdf
  // monic !!
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
  
  // sato:2019/04/22
  // see the definition of ae_int_t in alglib.
  // int info;
  long info;
  // sato:2019/04/22
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
  
  return;
}
//////////////////////////////////////////////////////////////////////////

