#include "td1c.hpp"
#include "wrapper.hpp"
#include "surff.hpp"


//////////////////////////////////////////////////////////////////////////
void Surff::rec_v2xmat(double dtime)
{
//Sato_tSURFF
//  if(nprint <= 0){return;}
//Sato_tSURFF
  
  surff_rec_v2xmat_(&v2xmat[0]);
  for(int iorb = 0; iorb < norb; ++iorb){
    for(int jorb = 0; jorb < norb; ++jorb){
      v2xmat[iorb * norb + jorb] = v2xmat[iorb * norb + jorb] / dtime;
    }
  }
  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::rec_v2xmat(double dtime, const std::vector<dcomplex>& xmat)
{
//Sato_tSURFF
//  if(nprint <= 0){return;}
//Sato_tSURFF
  
  for(int iorb = 0; iorb < norb; ++iorb){
    for(int jorb = 0; jorb < norb; ++jorb){
      v2xmat[iorb * norb + jorb] = xmat[iorb * norb + jorb] / dtime;
    }
  }

  //DEBUG
  //printf("v11.real %20.15E \n", v2xmat[0].real());

  //double norm = 0;
  //for(int iorb = 0; iorb < norb; ++iorb){
  //norm += std::pow(std::abs(orbpes[iorb * ksize]), 2.0);
  //}
  //std::cout <<  "kconst "<< norm << std::endl;

  //std::cout <<  "v21.real "<< v2xmat[1 * norb + 0].real() << std::endl;
  //std::cout <<  "v21.imag "<< v2xmat[1 * norb + 0].imag() << std::endl;
  //DEBUG

}
//////////////////////////////////////////////////////////////////////////
void Surff::prop(const clmpi& MPIP, const clwfn& Wfn, const clbas& Bas, const clfield& Field, const clhprod& HP)
{
//Sato_tSURFF
//  if(nprint <= 0){return;}
//Sato_tSURFF

  calc_dt(MPIP, Wfn, Bas, Field, orbpes_dt);

  //prop_trap(MPIP, Bas, Field, orbpes_dt, orbpes);

  ///////test
  //std::fill(orbpes_dt.begin(), orbpes_dt.end(), 0.0);
  ///////test

  std::vector<dcomplex> den1mat(norb * norb, 0.0);
  ormas_mkden1_(&Wfn.wfn[Wfn.size1], &den1mat[0]);


  //default integ->prop(MPIP, Bas, Field, v2xmat, orbpes_dt, orbpes);
  if(integ_type == "EI"){
    integ->prop(MPIP, Bas, Field, v2xmat, den1mat, orbpes_dt, orbpes);
  }
  else{
    integ->prop(MPIP, Bas, Field, v2xmat, orbpes_dt, orbpes);
  }
  
  // for remove_Q
  //prop_trap_remove_Q(MPIP, Bas, Field, orbpes_dt, orbpes);
  
  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::prop_trap(const clmpi& MPIP, const clbas& Bas, const clfield& Field,
		      std::vector<dcomplex> &tmppes_dt, std::vector<dcomplex> &tmppes)
{

  //propagate
  double dt = Field.dtime;
  
  torbpes = tmppes;
  
#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);

    for(int iorb = 0; iorb < norb; ++iorb){
      ////////////////////////
      // omp parallelized
      ////////////////////////
      for(int ikang = llkang; ikang < ulkang; ikang++){
	for(int ikrad = 0; ikrad  < num_krad; ikrad++){
	  tmppes[iorb * num_kang * num_krad + ikang * num_krad + ikrad] 
	    += dt * tmppes_dt[iorb * num_kang * num_krad + ikang * num_krad + ikrad];
	}
      }
    }
  }


#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);

    for(int iorb = 0; iorb < norb; ++iorb){
      for(int jorb = 0; jorb < norb; ++jorb){
	////////////////////////
	// omp parallelized
	////////////////////////
	for(int ikang = llkang; ikang < ulkang; ikang++){
	  for(int ikrad = 0; ikrad  < num_krad; ikrad++){
	    tmppes[iorb * num_kang * num_krad + ikang * num_krad + ikrad] 
	      += dt * torbpes[jorb * num_kang * num_krad + ikang * num_krad + ikrad] * v2xmat[iorb * norb + jorb];
	  }
	}
      }
    }
  }

  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::prop_trap_remove_Q(const clmpi& MPIP, const clbas& Bas, const clfield& Field,
		      std::vector<dcomplex> &tmppes_dt, std::vector<dcomplex> &tmppes)
{

  //propagate
  double dt = Field.dtime;
  
  torbpes = tmppes;
  
#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);

    for(int iorb = 0; iorb < norb; ++iorb){
      ////////////////////////
      // omp parallelized
      ////////////////////////
      for(int ikang = llkang; ikang < ulkang; ikang++){
	for(int ikrad = 0; ikrad  < num_krad; ikrad++){
	  tmppes[iorb * num_kang * num_krad + ikang * num_krad + ikrad] 
	    += dt * tmppes_dt[iorb * num_kang * num_krad + ikang * num_krad + ikrad];
	}
      }
    }
  }


  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::calc_dt(const clmpi& MPIP, const clwfn& Wfn, const clbas& Bas, const clfield& Field, std::vector<dcomplex> &tmppes_dt)
{
  double field[9];
  Field.get_value(field);
  double afield = -field[8];
  dcomplex ivec_pot = IUNIT * afield;
  double fpi = FOUR * PI;

  prop_trap_vphase(MPIP, afield, Field.time, Field.dtime);
  calc_orb(Bas, Wfn, irad_surf, orb_surf, dorb_surf);
  std::vector<dcomplex> ipow(lnum);
  for(int l = 0; l < lnum; ++l){
    ipow[l] = pow(-IUNIT, l);
  }
  std::fill(tmppes_dt.begin(), tmppes_dt.end(), 0.0);

#pragma omp parallel default(shared)
  {
    int mval;
    dcomplex tmp = 0.0;
    dcomplex factl = 0.0;
    double legendre_ang = 0.0;
    dcomplex omega_rad;

    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);
    
    for(int iorb = 0; iorb < norb; ++iorb){
      mval = std::abs(Bas.mval[iorb]);
      for(int orbl = mval; orbl < lnum; orbl++){
	factl = fpi * ipow[orbl];
	////////////////////////
	// omp parallelized
	////////////////////////
	for(int ikang = llkang; ikang < ulkang; ikang++){ //kang l loop
	  //std::cout << sph_ang.size() <<  " " << mval * lnum * num_kang + orbl * num_kang + ikang << std::endl;
	  legendre_ang = sph_ang[mval * lnum * num_kang + orbl * num_kang + ikang];
	  for(int ikrad = 0; ikrad < num_krad; ikrad++){
	    //omega_rad = <j_kl'| Omega(r) |phi_nl'>
	    omega_rad = HALF * rrad_surf * rrad_surf 
	      * ( conj(sbessel_dr[orbl * num_krad + ikrad]) * orb_surf[iorb * lnum + orbl]
		  - conj(sbessel[orbl * num_krad + ikrad]) * dorb_surf[iorb * lnum + orbl] );
	    tmp = factl * legendre_ang * omega_rad;
	    tmppes_dt[iorb * num_kang * num_krad + ikang * num_krad + ikrad] += tmp;
	  }
	}
      }
    }
  }
  //std::cout << " 1st "<< tmppes_dt[0] << std::endl;
  

  //DEBUG
//   int iorb = 0;
//   for(int ikrad = 0; ikrad  < num_krad; ikrad++){    
//     for(int ikang = 0; ikang < num_kang; ikang++){
//       std::cout << krad[ikrad] << " " << kang[ikang] << " " << tmppes_dt[iorb * num_kang * num_krad + ikang * num_krad + ikrad] << std::endl; 
//     }
//   }
    

#pragma omp parallel default(shared)
  {
    dcomplex tmp0, tmp1;
    dcomplex factl0, factl1;
    double legendre_ang0;
    double legendre_ang1;
    dcomplex delta_rad0, delta_rad1;
    int smval, mval, midx;
    int mmax = Bas.GAng.mmax1;

    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);

    for(int iorb = 0; iorb < norb; ++iorb){
      smval = Bas.mval[iorb];
      mval = std::abs(Bas.mval[iorb]);
      midx = mmax + smval;
      for(int orbl = mval; orbl < lnum - 1; orbl++){
	//old factl0 = - ivec_pot * fpi * ipow[orbl] * alph_lm[midx * lnum + orbl];
	//old factl1 = - ivec_pot * fpi * ipow[orbl + 1] * alph_lm[midx * lnum + orbl];
	factl0 = ivec_pot * fpi * ipow[orbl] * alph_lm[midx * lnum + orbl];
	factl1 = ivec_pot * fpi * ipow[orbl + 1] * alph_lm[midx * lnum + orbl];

	////////////////////////
	// omp parallelized
	////////////////////////
	for(int ikang = llkang; ikang < ulkang; ikang++){ //kang l loop
	  legendre_ang0 = sph_ang[mval * lnum * num_kang + orbl * num_kang + ikang];
	  legendre_ang1 = sph_ang[mval * lnum * num_kang + (orbl+1) * num_kang + ikang];
	  for(int ikrad = 0; ikrad  < num_krad; ikrad++){
	    //delta_rad0 = <j_k(orbl)| delta(r) |phi_n(orbl+1)>
	    //delta_rad1 = <j_k(orbl+1)| delta(r) |phi_n(orbl)>
	    delta_rad0 = rrad_surf * rrad_surf * conj(sbessel[orbl * num_krad + ikrad]) * orb_surf[iorb * lnum + orbl + 1];
	    delta_rad1 = rrad_surf * rrad_surf * conj(sbessel[(orbl + 1) * num_krad + ikrad]) * orb_surf[iorb * lnum + orbl];
	    tmp0  = factl0 * legendre_ang0 * delta_rad0;
	    tmp1  = factl1 * legendre_ang1 * delta_rad1;
	    tmppes_dt[iorb * num_kang * num_krad + ikang * num_krad + ikrad] += tmp0 + tmp1;
	  }
	}
      }
    } 
  }  
  //std::cout << " 2rd "<< tmppes_dt[0] << std::endl;

  
  // calculate orbpes_dt
#pragma omp parallel default(shared)
  {
    dcomplex dfact;
    double invpi = 1.0 / std::pow(2*PI, 1.5);
    
    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);

    for(int iorb = 0; iorb < norb; ++iorb){
      ////////////////////////
      // omp parallelized
      ////////////////////////
      for(int ikang = llkang; ikang < ulkang; ikang++){ //kang l loop
	for(int ikrad = 0; ikrad  < num_krad; ikrad++){
	  dfact = invpi * IUNIT * exp(IUNIT * vphase[ikang * num_krad + ikrad]);
	  tmppes_dt[iorb * num_kang * num_krad + ikang * num_krad + ikrad] 
	    = dfact * tmppes_dt[iorb * num_kang * num_krad + ikang * num_krad + ikrad];
	}
      }
    }
  }
  //std::cout << " 3rd "<< tmppes_dt[0] << std::endl;


  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::prop_trap_vphase(const clmpi& MPIP, double afield, double time, double dtime)
{

#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);

    ////////////////////////
    // omp parallelized
    ////////////////////////
    for(int ikang = llkang; ikang < ulkang; ikang++){ 
      for(int ikrad = 0; ikrad < num_krad; ++ikrad){
	tdfac[ikang * num_krad + ikrad] += 
	  dtime * ( - afield * krad[ikrad] * std::cos(kang[ikang]));
	//old tdfac[ikang * num_krad + ikrad] += 
	//old dtime * (HALF * afield * afield - afield * krad[ikrad] * std::cos(kang[ikang]));
	
	vphase[ikang * num_krad + ikrad] = 
	  ( HALF * krad[ikrad] * krad[ikrad] * time + tdfac[ikang * num_krad + ikrad]);
      }
    }
  }
  
  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::calc_orb(const clbas &Bas, const clwfn &Wfn, 
		     int sirad, std::vector<dcomplex> &norb, std::vector<dcomplex> &dorb)
{
  calc_orb_(&sirad, &Wfn.wfn[0], &norb[0], &dorb[0]);
  return;
}
//////////////////////////////////////////////////////////////////////////

