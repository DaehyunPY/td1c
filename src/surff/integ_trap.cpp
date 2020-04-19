#include "td1c.hpp"
#include "wrapper.hpp"
#include "surff.hpp"
//////////////////////////////////////////////////////////////////////////
Integ_trap::Integ_trap(int vir_ksize, int vir_norb) : Integ(vir_ksize, vir_norb)
{
  std::cout << "# Integ_trap()" << std::endl;

  tmp_opes.resize(ksize * norb, 0.0);
  return;
};
//////////////////////////////////////////////////////////////////////////
Integ_trap::~Integ_trap()
{
};
//////////////////////////////////////////////////////////////////////////
void Integ_trap::prop(const clmpi& MPIP, const clbas& Bas, const clfield& Field, const std::vector<dcomplex> &v2xmat,
		    const std::vector<dcomplex> &den1, std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes)
{
}
//////////////////////////////////////////////////////////////////////////
void Integ_trap::prop(const clmpi& MPIP, const clbas& Bas, const clfield& Field, const std::vector<dcomplex> &v2xmat,
		      std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes)
{
  //propagate
  double dt = Field.dtime;
  
  tmp_opes = opes;
  
#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(int iorb = 0; iorb < norb; ++iorb){
      ////////////////////////
      // omp parallelized
      ////////////////////////
      for(int ik = llk; ik  < ulk; ik++){
	  opes[iorb * ksize + ik] += dt * opes_dt[iorb * ksize + ik];
      }
    }
  }


#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(int iorb = 0; iorb < norb; ++iorb){
      for(int jorb = 0; jorb < norb; ++jorb){
	////////////////////////
	// omp parallelized
	////////////////////////
	for(int ik = llk; ik  < ulk; ik++){
	  opes[iorb * ksize + ik]
	    += dt * tmp_opes[jorb * ksize + ik] * v2xmat[iorb * norb + jorb];
	}
      }
    }
  }


//old unparallelized
//   for(int iorb = 0; iorb < norb; ++iorb){
//     for(int ik = 0; ik  < ksize; ik++){
//       opes[iorb * ksize + ik] += dt * opes_dt[iorb * ksize + ik];
//     }
//   }

//   for(int iorb = 0; iorb < norb; ++iorb){
//     for(int jorb = 0; jorb < norb; ++jorb){
//       for(int ik = 0; ik  < ksize; ik++){
// 	opes[iorb * ksize + ik]
// 	  += dt * tmp_opes[jorb * ksize + ik] * v2xmat[iorb * norb + jorb];
//       }
//     }
//   }

  return;

};
//////////////////////////////////////////////////////////////////////////

