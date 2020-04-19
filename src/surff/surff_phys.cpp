#include "td1c.hpp"
#include "wrapper.hpp"
#include "surff.hpp"

//////////////////////////////////////////////////////////////////////////
void Surff::calc_me_spectrum(const clmpi& MPIP, const clwfn &Wfn, const clbas &Bas, clhprod &HPW)
{ 
//Sato_tSURFF
  //  if(nprint <= 0){return;}
//Sato_tSURFF
  
  calc_momentum_spectrum(MPIP, &Wfn.wfn[Wfn.size1], Bas, HPW);
  //calc_momentum_spectrum_1e(MPIP, &Wfn.wfn[Wfn.size1], Bas, HPW);
  convert_energy_spectrum(); 
  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::calc_momentum_spectrum_1e(const clmpi &MPIP, const dcomplex* cic, const clbas &Bas, clhprod &HPW)
{
////////////////////////
// for 1 electron
////////////////////////
  std::fill(mom_spec.begin(), mom_spec.end(), 0.0);  
  for(int ikang = 0; ikang < num_kang; ++ikang){
    for(int ikrad = 0; ikrad < num_krad; ++ikrad){
      mom_spec[ikrad] += std::abs(orbpes[0 * num_kang * num_krad + ikang * num_krad + ikrad]
				  * orbpes[0 * num_kang * num_krad + ikang * num_krad + ikrad]) * wang[ikang];
    }
  }

  
  norm = 0.0;
  double delk = 0.0;
  for(int ikrad = 0; ikrad < num_krad-1; ++ikrad){
    delk = krad[ikrad + 1] - krad[ikrad];
    norm += mom_spec[ikrad] * krad[ikrad] * krad[ikrad] * delk; 
  }
  std::cout << "# norm of orbpes = " << norm << std::endl;
  
  return;
  
}
//////////////////////////////////////////////////////////////////////////
 void Surff::calc_momentum_spectrum(const clmpi &MPIP, const dcomplex* cic, const clbas &Bas, clhprod &HPW){
  
   std::vector<dcomplex> den1mat(norb * norb, 0.0);
   ormas_mkden1_(&cic[0], &den1mat[0]);
   
   
   //not required. std::fill(rhok.begin(), rhok.end(), 0.0);
#pragma omp parallel default(shared)
   {
    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);
    
    for(int iorb = 0; iorb < norb; ++iorb){
      for(int jorb = 0; jorb < norb; ++jorb){
	if(Bas.mval[iorb] == Bas.mval[jorb]){

	  ////////////////////////
	  // omp parallelized
	  ////////////////////////
	  for(int ikang = llkang; ikang < ulkang; ++ikang){
	    for(int ikrad = 0; ikrad < num_krad; ++ikrad){
	      rhok[iorb * norb * num_kang * num_krad + jorb * num_kang * num_krad + ikang * num_krad + ikrad]
		= conj(orbpes[iorb * num_kang * num_krad + ikang * num_krad + ikrad]) 
		* orbpes[jorb * num_kang * num_krad + ikang * num_krad + ikrad];
	    }
	  }	
	}
      }
    }
  }
  
  std::fill(ang_spec.begin(), ang_spec.end(), 0.0);
#pragma omp parallel default(shared)
  {
    dcomplex den1;
    int ncore = Bas.ORMAS.ncore;

    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);
    
    for(int iorb = 0; iorb < ncore; ++iorb){
      ////////////////////////
      // omp parallelized
      ////////////////////////
      for(int ikang = llkang; ikang < ulkang; ++ikang){
	for(int ikrad = 0; ikrad < num_krad; ++ikrad){
	  ang_spec[ikang * num_krad + ikrad] += 
	    rhok[iorb * norb * num_kang * num_krad + iorb * num_kang * num_krad + ikang * num_krad + ikrad]
	    * 2.0;
	}
      }
    }
  }
   
#pragma omp parallel default(shared)
  {
    dcomplex den1;
    int ncore = Bas.ORMAS.ncore;

    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);
 
    for(int iorb = ncore; iorb < norb; ++iorb){
      for(int jorb = ncore; jorb < norb; ++jorb){
	//den1 = den1mat[(iorb - nfcore) * norb + (jorb - nfcore)];
	den1 = den1mat[(iorb - ncore) * (norb - ncore) + (jorb - ncore)];

	////////////////////////
	// omp parallelized
	////////////////////////
	for(int ikang = llkang; ikang < ulkang; ++ikang){
	  for(int ikrad = 0; ikrad < num_krad; ++ikrad){
	    ang_spec[ikang * num_krad + ikrad] += 
	      rhok[iorb * norb * num_kang * num_krad + jorb * num_kang * num_krad + ikang * num_krad + ikrad]
	      * den1;
	  }
	}	
      }
    }
  }
  

  static std::vector<double> tmp_mom(MPIP.nthr * num_krad);
  std::fill(tmp_mom.begin(), tmp_mom.end(), 0.0);

#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llkang, ulkang;
    MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);

    ////////////////////////
    // omp parallelized
    ////////////////////////
    for(int ikang = llkang; ikang < ulkang; ++ikang){
      for(int ikrad = 0; ikrad < num_krad; ++ikrad){
	tmp_mom[ithr * num_krad + ikrad] += real(ang_spec[ikang * num_krad + ikrad]) * wang[ikang];
      }
    }
  }  
 
  std::fill(mom_spec.begin(), mom_spec.end(), 0.0);  
#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llkrad, ulkrad;
    MPIP.omp_divide(ithr, 0, num_krad, llkrad, ulkrad);
    
    for(int jthr = 0; jthr < MPIP.nthr; ++jthr){
      
      ////////////////////////
      // omp parallelized
      ////////////////////////
      for(int ikrad = llkrad; ikrad < ulkrad; ikrad++){
	mom_spec[ikrad] += tmp_mom[jthr * num_krad + ikrad];
      }
    }
  }
  

//bugs #pragma omp parallel default(shared)
//bugs   {
//bugs     int ithr = omp_get_thread_num();
//bugs     int llkang, ulkang;
//bugs     MPIP.omp_divide(ithr, 0, num_kang, llkang, ulkang);
//bugs
//bugs     ////////////////////////
//bugs     // omp parallelized
//bugs     ////////////////////////
//bugs     for(int ikang = llkang; ikang < ulkang; ++ikang){
//bugs       for(int ikrad = 0; ikrad < num_krad; ++ikrad){
//bugs 	mom_spec[ikrad] += real(ang_spec[ikang * num_krad + ikrad]) * wang[ikang];
//bugs       }
//bugs     }
//bugs   }  
  

//old double del_ang = kang[1] - kang[0];
////////////////////////
// for 1 electron
////////////////////////
//   std::fill(mom_spec.begin(), mom_spec.end(), 0.0);  
//   for(int ikang = 0; ikang < num_kang; ++ikang){
//     for(int ikrad = 0; ikrad < num_krad; ++ikrad){
//       //old mom_spec[ikrad] += std::abs(orbpes[0 * num_kang * num_krad + ikang * num_krad + ikrad]
//       //old				  * orbpes[0 * num_kang * num_krad + ikang * num_krad + ikrad]) * sin(kang[ikang]) * del_ang;
//       mom_spec[ikrad] += std::abs(orbpes[0 * num_kang * num_krad + ikang * num_krad + ikrad]
// 				  * orbpes[0 * num_kang * num_krad + ikang * num_krad + ikrad]) * wang[ikang];
//     }
//   }

  
  norm = 0.0;
  double delk = 0.0;
  for(int ikrad = 0; ikrad < num_krad-1; ++ikrad){
    delk = krad[ikrad + 1] - krad[ikrad];
    norm += mom_spec[ikrad] * krad[ikrad] * krad[ikrad] * delk; 
  }
  std::cout << "# norm of orbpes = " << norm << std::endl;

  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::convert_energy_spectrum(){

  double small = 10E-10;
  for(int ikrad = 0; ikrad < num_krad; ++ikrad){
    if(krad[ikrad] > small){
      ene_spec[ikrad] = mom_spec[ikrad] / krad[ikrad];
    }
    else {
      ene_spec[ikrad] = 0.0;
    }
  }  
  
  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::print(const clio &IO)
{
//Sato_tSURFF
//  if(nprint <= 0){return;}
//Sato_tSURFF

  printf("\n");
  printf("# spectrum amplitude \n");
  printf("# krad kang real[spectrum amplitude (iorb)] imag[spectrum amplitude (iorb)] ... \n" );
  dcomplex tmp;
  for(int ikrad = 0; ikrad  < num_krad; ikrad++){
    for(int ikang = 0; ikang < num_kang; ikang++){ //kang l loop
      printf(" %20.10E    %20.10E   ",  krad[ikrad], kang[ikang]);
      for(int iorb = 0; iorb < norb; ++iorb){
	tmp = orbpes[iorb * num_kang * num_krad + ikang * num_krad + ikrad];
	printf(" %20.10E    %20.10E   ",  tmp.real(), tmp.imag());
      }
      printf("\n");
    }
  }
  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::print_amplitude(const clio &IO)
{
//Sato_tSURFF
//  if (nprint <= 0) return;
//Sato_tSURFF

  std::string fname;
  fname = strsprintf("%s.surff", IO.name.c_str());
  
  //fname = IO.name + ".surff";
  FILE *fpo;
  fpo = fopen(fname.c_str(), "w");
  
  fprintf(fpo, "# krad kang real[spectrum amplitude (iorb)] imag[spectrum amplitude (iorb)] ... \n" );
  fprintf(fpo, "%d\n%d\n", num_krad, num_kang);
  dcomplex tmp;
  for(int ikrad = 0; ikrad  < num_krad; ikrad++){
    for(int ikang = 0; ikang < num_kang; ikang++){
      fprintf(fpo, " %20.10E    %20.10E   ",  krad[ikrad], kang[ikang]);
      for(int iorb = 0; iorb < norb; ++iorb){
	tmp = orbpes[iorb * num_kang * num_krad + ikang * num_krad + ikrad];
	fprintf(fpo, " %20.10E    %20.10E   ",  tmp.real(), tmp.imag());
      }
      fprintf(fpo, "\n");
    }
  }
  fclose(fpo);

}
//////////////////////////////////////////////////////////////////////////
void Surff::print_file(const clio &IO)
{
//Sato_tSURFF
//  if(nprint <= 0){return;}
//Sato_tSURFF

  static int iprint = 0;
  
  FILE *fpo;  
  
  std::string fname;
// Sato_TSURFF
  //  fname = strsprintf("./SPECTRUM/%s_%08d.ang_spec", IO.name.c_str(), iprint);
  fname = strsprintf("./%s.SPECTRUM/%s_%08d.ang_spec", IO.name.c_str(), IO.name.c_str(), iprint);
// Sato_TSURFF
  fpo = fopen(fname.c_str(), "w");
  
  //DEBUG
//   for(int i = 0; i < norb; ++i){
//     for(int j = 0; j < norb; ++j){
//       fprintf(fpo, "%+15.5E ", v2xmat[i * norb + j].real());
//     }
//     fprintf(fpo, "\n");
//   }
//   fprintf(fpo, "\n");
//   fprintf(fpo, "\n");
//   for(int i = 0; i < norb; ++i){
//     for(int j = 0; j < norb; ++j){
//       fprintf(fpo, "%+15.5E ", v2xmat[i * norb + j].imag());
//     }
//     fprintf(fpo, "\n");
//   }
//   fprintf(fpo, "\n");
  //Debug

  fprintf(fpo, "# krad kang momentum spectrum  \n" );
  fprintf(fpo, "%d\n%d\n", num_krad, num_kang);
  double dtmp;
  for(int ikrad = 0; ikrad  < num_krad; ikrad++){
    for(int ikang = 0; ikang < num_kang; ikang++){
      fprintf(fpo, " %20.10E    %20.10E   ",  krad[ikrad], kang[ikang]);
      dtmp = real(ang_spec[ikang * num_krad + ikrad]);
      fprintf(fpo, " %20.10E    ",  dtmp);
      fprintf(fpo, "\n");
    }
  }  
  fclose(fpo);

// Sato_TSURFF
//  fname = strsprintf("./SPECTRUM/%s_%08d.me_spec", IO.name.c_str(), iprint);
  fname = strsprintf("./%s.SPECTRUM/%s_%08d.me_spec", IO.name.c_str(), IO.name.c_str(), iprint);
// Sato_TSURFF
  //fname = IO.name + ".me_spec";
  fpo = fopen(fname.c_str(), "w");
  
  fprintf(fpo, "# norm = %20.10E, k rho(k) E rho(E) \n", norm);
  double ene;
  double rho;
  for(int ikrad = 0; ikrad  < num_krad; ikrad++){
    ene = HALF * krad[ikrad] * krad[ikrad];
    fprintf(fpo, " %20.10E    %20.10E   %20.10E    %20.10E", krad[ikrad], mom_spec[ikrad], ene, ene_spec[ikrad]);
    fprintf(fpo, "\n");
  }
  fclose(fpo);


  ++iprint;
  return;
}
//////////////////////////////////////////////////////////////////////////
void Surff::calc_t_surff(const clmpi& MPIP, const clwfn& Wfn, const clbas& Bas, const clfield& Field)
{


}
