#include "td1c.hpp"
#include "wrapper.hpp"
#include "surff.hpp"
//////////////////////////////////////////////////////////////////////////
Integ_CN::Integ_CN(int vir_ksize, int vir_norb) : Integ(vir_ksize, vir_norb)
{
  std::cout << "# Integ_CN()" << std::endl;

  tmp_opes.resize(ksize * norb, 0.0);
  exp_vmat.resize(norb * norb, 0.0);
  imp_vmat.resize(norb * norb, 0.0);

  old_v2xmat.resize(norb * norb, 0.0);
  old_opes_dt.resize(ksize * norb, 0.0);
  
}
//////////////////////////////////////////////////////////////////////////
Integ_CN::~Integ_CN()
{
}
//////////////////////////////////////////////////////////////////////////
void Integ_CN::rec_old_state(const clmpi& MPIP, const std::vector<dcomplex> &v2xmat, std::vector<dcomplex> &opes_dt)
{

  old_v2xmat = v2xmat;
  old_opes_dt = opes_dt;

  return;
}
//////////////////////////////////////////////////////////////////////////
void Integ_CN::prop(const clmpi& MPIP, const clbas& Bas, const clfield& Field, const std::vector<dcomplex> &v2xmat,
		    const std::vector<dcomplex> &den1, std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes)
{
}
//////////////////////////////////////////////////////////////////////////
void Integ_CN::prop(const clmpi& MPIP, const clbas& Bas, const clfield& Field, const std::vector<dcomplex> &v2xmat,
		      std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes)
{

  double dt = Field.dtime;
  std::fill(tmp_opes.begin(), tmp_opes.end(), 0.0);


  for(int iorb = 0; iorb < norb; ++iorb){
    for(int jorb = 0; jorb < norb; ++jorb){
      exp_vmat[iorb * norb + jorb] = + HALF * dt * v2xmat[iorb * norb + jorb];
      imp_vmat[iorb * norb + jorb] = - HALF * dt * old_v2xmat[iorb * norb + jorb];
    }  
    exp_vmat[iorb * norb + iorb] += ONE;
    imp_vmat[iorb * norb + iorb] += ONE;
  }


  static long long dim = norb;
  static long long lda= norb;
  static std::vector<long long> ipiv(norb);
  static long long lwork = norb;
  static std::vector<dcomplex>  work(norb*norb);
  static long long info;

// Sato_TSURFF
  long norbl = norb;
  futil_gmatinv_(&norbl, &ZERO, &imp_vmat[0], &imp_vmat[0]);
//  zgetrf(&dim, &dim, &imp_vmat[0], &lda, &ipiv[0], &info);
//  //info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, dim, dim, &imp_vmat[0], lda, &ipiv[0]);
//  if(info != 0){
//    std::cout << "zgetrf_ failed at Integ_CN::prop" << std::endl;
//    abort();
//  }
//
//  zgetri(&dim, &imp_vmat[0], &lda, &ipiv[0], &work[0], &lwork, &info);
//  //info = LAPACKE_zgetri(LAPACK_COL_MAJOR, dim, &imp_vmat[0], lda, &ipiv[0]);
//  if(info != 0){
//    std::cout << "zgetri_ failed at Integ_CN::prop" << std::endl;
//    abort();
//  }
// Sato_TSURFF

//   for(long iorb = 0; iorb < norb; ++iorb){
//     for(long jorb = 0; jorb < norb; ++jorb){
//       std::cout << imp_vmat[iorb * norb + jorb];
//     }
//   }
  
#pragma omp parallel default(shared)
  {
    long ithr = omp_get_thread_num();
    long llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(long iorb = 0; iorb < norb; ++iorb){
      for(long jorb = 0; jorb < norb; ++jorb){
	if(Bas.mval[iorb] == Bas.mval[jorb]){
	  ////////////////////////
	  // omp parallelized
	  ////////////////////////
	  for(int ik = llk; ik < ulk; ik++){
	    tmp_opes[iorb * ksize + ik] += opes[jorb * ksize + ik] * exp_vmat[iorb * norb + jorb]; 
	    //bugs tmp_opes[iorb * ksize + ik] += opes[jorb * ksize + ik] * exp_vmat[jorb * norb + iorb]; 
	  }
	}
      }
      for(int ik = llk; ik < ulk; ik++){
	tmp_opes[iorb * ksize + ik] += HALF * dt * (opes_dt[iorb * ksize + ik] + old_opes_dt[iorb * ksize + ik]);
      }
    }    
  }


  std::fill(opes.begin(), opes.end(), 0.0);
#pragma omp parallel default(shared)
  {
    long ithr = omp_get_thread_num();
    long llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(long iorb = 0; iorb < norb; ++iorb){
      for(long jorb = 0; jorb < norb; ++jorb){
	if(Bas.mval[iorb] == Bas.mval[jorb]){
	  ////////////////////////
	  // omp parallelized
	  ////////////////////////
	  for(int ik = llk; ik < ulk; ik++){
	    opes[iorb * ksize + ik] += tmp_opes[jorb * ksize + ik] * imp_vmat[iorb * norb + jorb];
	    //bugs opes[iorb * ksize + ik] += tmp_opes[jorb * ksize + ik] * imp_vmat[jorb * norb + iorb];
	  }
	}
      }
    }    
  }

  rec_old_state(MPIP, v2xmat, opes_dt);
  return;
}
//////////////////////////////////////////////////////////////////////////
void Integ_CN::prop_old(const clmpi& MPIP, const clbas& Bas, const clfield& Field, const std::vector<dcomplex> &v2xmat,
		      std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes)
{

  double dt = Field.dtime;
  std::fill(tmp_opes.begin(), tmp_opes.end(), 0.0);


  for(int iorb = 0; iorb < norb; ++iorb){
    for(int jorb = 0; jorb < norb; ++jorb){
      exp_vmat[iorb * norb + jorb] = + HALF * dt * v2xmat[iorb * norb + jorb];
      imp_vmat[iorb * norb + jorb] = - HALF * dt * v2xmat[iorb * norb + jorb];
    }  
    exp_vmat[iorb * norb + iorb] += ONE;
    imp_vmat[iorb * norb + iorb] += ONE;
  }


  long long dim = norb;
  long long lda= norb;
  static std::vector<long long> ipiv(norb);
  long long lwork = norb;
  static std::vector<dcomplex>  work(norb);
  long long info;

// Sato_TSURFF
  long norbl = norb;
  futil_gmatinv_(&norbl, &ZERO, &imp_vmat[0], &imp_vmat[0]);
//  zgetrf(&dim, &dim, &imp_vmat[0], &lda, &ipiv[0], &info);
//  //info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, dim, dim, &imp_vmat[0], lda, &ipiv[0]);
//  if(info != 0){
//    std::cout << "zgetrf_ failed at Integ_CN::prop" << std::endl;
//    abort();
//  }
//
//  zgetri(&dim, &imp_vmat[0], &lda, &ipiv[0], &work[0], &lwork, &info);
//  //info = LAPACKE_zgetri(LAPACK_COL_MAJOR, dim, &imp_vmat[0], lda, &ipiv[0]);
//  if(info != 0){
//    std::cout << "zgetri_ failed at Integ_CN::prop" << std::endl;
//    abort();
//  }
// Sato_TSURFF

//   for(long iorb = 0; iorb < norb; ++iorb){
//     for(long jorb = 0; jorb < norb; ++jorb){
//       std::cout << imp_vmat[iorb * norb + jorb];
//     }
//   }
  
#pragma omp parallel default(shared)
  {
    long ithr = omp_get_thread_num();
    long llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(long iorb = 0; iorb < norb; ++iorb){
      for(long jorb = 0; jorb < norb; ++jorb){
	if(Bas.mval[iorb] == Bas.mval[jorb]){
	  ////////////////////////
	  // omp parallelized
	  ////////////////////////
	  for(int ik = llk; ik < ulk; ik++){
	    tmp_opes[iorb * ksize + ik] += opes[jorb * ksize + ik] * exp_vmat[iorb * norb + jorb]; 
	  }
	}
      }
      for(int ik = llk; ik < ulk; ik++){
	tmp_opes[iorb * ksize + ik] += dt * opes_dt[iorb * ksize + ik];
      }
    }    
  }


  std::fill(opes.begin(), opes.end(), 0.0);
#pragma omp parallel default(shared)
  {
    long ithr = omp_get_thread_num();
    long llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(long iorb = 0; iorb < norb; ++iorb){
      for(long jorb = 0; jorb < norb; ++jorb){
	if(Bas.mval[iorb] == Bas.mval[jorb]){
	  ////////////////////////
	  // omp parallelized
	  ////////////////////////
	  for(int ik = llk; ik < ulk; ik++){
	    opes[iorb * ksize + ik] += tmp_opes[jorb * ksize + ik] * imp_vmat[iorb * norb + jorb];
	  }
	}
      }
    }    
  }


////////////////////////
// unprallelized version
////////////////////////

//   int llk = 0;
//   int ulk = ksize;
  
//     for(long iorb = 0; iorb < norb; ++iorb){
//       for(long jorb = 0; jorb < norb; ++jorb){
// 	////////////////////////
// 	// omp parallelized
// 	////////////////////////
// 	for(int ik = llk; ik < ulk; ik++){
// 	  tmp_opes[iorb * ksize + ik] += opes[jorb * ksize + ik] * exp_vmat[jorb * norb + iorb];
// 	}
//       }
//       for(int ik = llk; ik < ulk; ik++){
// 	tmp_opes[iorb * ksize + ik] += dt * opes_dt[iorb * ksize + ik];
//       }
//     }    

//     std::fill(opes.begin(), opes.end(), 0.0);
//     for(long iorb = 0; iorb < norb; ++iorb){
//       for(long jorb = 0; jorb < norb; ++jorb){
// 	////////////////////////
// 	// omp parallelized
// 	////////////////////////
// 	for(int ik = llk; ik < ulk; ik++){
// 	  opes[iorb * ksize + ik] += tmp_opes[jorb * ksize + ik] * imp_vmat[jorb * norb + iorb];
// 	}
//       }
//     }    


  return;
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

