#include "td1c.hpp"
#include "wrapper.hpp"
#include "surff.hpp"
//////////////////////////////////////////////////////////////////////////
Integ_EI::Integ_EI(int vir_ksize, int vir_norb) : Integ(vir_ksize, vir_norb)
{
  std::cout << "# Integ_EI()" << std::endl;
 
  renormarize = 0;
  return;
}
//////////////////////////////////////////////////////////////////////////
Integ_EI::~Integ_EI()
{
}
//////////////////////////////////////////////////////////////////////////
void Integ_EI::set_renormarize(int iren)
{
  renormarize = iren;
  return;
} 
//////////////////////////////////////////////////////////////////////////
void Integ_EI::prop(const clmpi& MPIP, const clbas& Bas, const clfield& Field, const std::vector<dcomplex> &v2xmat,
		    std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes)
{
}
//////////////////////////////////////////////////////////////////////////
void Integ_EI::prop(const clmpi& MPIP, const clbas& Bas, const clfield& Field, const std::vector<dcomplex> &v2xmat,
		    const std::vector<dcomplex> &den1, std::vector<dcomplex> &opes_dt, std::vector<dcomplex> &opes)
{
  double dt = Field.dtime;

  static std::vector<dcomplex> eval(norb);
  static std::vector<dcomplex> evl(norb*norb);
  static std::vector<dcomplex> evr(norb*norb);
  eigen_decompose(dt, v2xmat, eval, evl, evr);

  static int nstep = 0;
  //  int pstep = 100;
  //  if(nstep % pstep == 0){
  //    printf("eig%d", pstep);
  //    for(int i = 0; i < norb; ++i){
  //      printf("%20.10E %20.10E ", eval[i].real(), eval[i].imag());
  //    }
  //    std::cout << std::endl;
  //  }
  ++ nstep;

  static std::vector<dcomplex> exp_eval(norb);
  calc_pade1(eval, exp_eval);
  //for(int iorb = 0; iorb < norb; ++iorb){
  //  exp_eval[iorb] = std::exp(eval[iorb]);
  //}

  static std::vector<dcomplex> inv_evr(norb*norb);  
  inv_evr = evr;
  inv_mat(inv_evr);


  ////////////////////////////
  // homogenious part
  ////////////////////////////
  static std::vector<dcomplex> homo_opes(norb * ksize);
  static std::vector<dcomplex> tmp_homo_opes(norb * ksize);
  
  std::fill(tmp_homo_opes.begin(), tmp_homo_opes.end(), 0.0);
#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(int iorb = 0; iorb < norb; ++iorb){
      for(int jorb = 0; jorb < norb; ++jorb){
	//	if(Bas.mval[iorb] == Bas.mval[jorb]){
	  ////////////////////////
	  // omp parallelized
	  ////////////////////////
	for(int ik = llk; ik < ulk; ik++){
	  tmp_homo_opes[iorb * ksize + ik] += evr[iorb * norb + jorb] * opes[jorb * ksize + ik];
	  //old2 tmp_homo_opes[iorb * ksize + ik] += inv_evr[iorb * norb + jorb] * opes[jorb * ksize + ik];
	  //old3 tmp_homo_opes[iorb * ksize + ik] += inv_evr[jorb * norb + iorb] * opes[jorb * ksize + ik];

	}
	//}
      }
      for(int ik = llk; ik < ulk; ik++){
	tmp_homo_opes[iorb * ksize + ik] *= exp_eval[iorb];
      }
    }
  }

  std::fill(homo_opes.begin(), homo_opes.end(), 0.0);
#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(int iorb = 0; iorb < norb; ++iorb){
      for(int jorb = 0; jorb < norb; ++jorb){
	//	if(Bas.mval[iorb] == Bas.mval[jorb]){
	////////////////////////
	// omp parallelized
	////////////////////////
	for(int ik = llk; ik < ulk; ik++){
	  homo_opes[iorb * ksize + ik] += inv_evr[iorb * norb + jorb] * tmp_homo_opes[jorb * ksize + ik];
	  //old2 homo_opes[iorb * ksize + ik] += evr[iorb * norb + jorb] * tmp_homo_opes[jorb * ksize + ik];
	  //old3 homo_opes[iorb * ksize + ik] += evr[jorb * norb + iorb] * tmp_homo_opes[jorb * ksize + ik];
	}
	//}
      }
    }
  }
  
  
  ////////////////////////////
  // inhomogenious part
  ////////////////////////////
  static std::vector<dcomplex> inhomo_opes(norb * ksize);
  static std::vector<dcomplex> tmp_inhomo_opes(norb * ksize);  
  static std::vector<dcomplex> pade_iexp(norb);
  
  calc_pade2(eval, pade_iexp);
  
  std::fill(tmp_inhomo_opes.begin(), tmp_inhomo_opes.end(), 0.0);
#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(int iorb = 0; iorb < norb; ++iorb){
      for(int jorb = 0; jorb < norb; ++jorb){
	//if(Bas.mval[iorb] == Bas.mval[jorb]){
	////////////////////////
	// omp parallelized
	////////////////////////
	for(int ik = llk; ik < ulk; ik++){
	  tmp_inhomo_opes[iorb * ksize + ik] += evr[iorb * norb + jorb] * opes_dt[jorb * ksize + ik];
	  //old2 tmp_inhomo_opes[iorb * ksize + ik] += inv_evr[iorb * norb + jorb] * opes_dt[jorb * ksize + ik];
	  //old3 tmp_inhomo_opes[iorb * ksize + ik] += inv_evr[jorb * norb + iorb] * opes_dt[jorb * ksize + ik];
	}
	//}
      }
      for(int ik = llk; ik < ulk; ik++){
	tmp_inhomo_opes[iorb * ksize + ik] *= dt * pade_iexp[iorb];
      }
    }
  }

  std::fill(inhomo_opes.begin(), inhomo_opes.end(), 0.0);
#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(int iorb = 0; iorb < norb; ++iorb){
      for(int jorb = 0; jorb < norb; ++jorb){
	//if(Bas.mval[iorb] == Bas.mval[jorb]){
	////////////////////////
	// omp parallelized
	////////////////////////
	for(int ik = llk; ik < ulk; ik++){
	  inhomo_opes[iorb * ksize + ik] += inv_evr[iorb * norb + jorb] * tmp_inhomo_opes[jorb * ksize + ik];
	  //old2 inhomo_opes[iorb * ksize + ik] += evr[iorb * norb + jorb] * tmp_inhomo_opes[jorb * ksize + ik];
	  //old3 inhomo_opes[iorb * ksize + ik] += evr[jorb * norb + iorb] * tmp_inhomo_opes[jorb * ksize + ik];
	}
	//}
      }
    }
  }



  ///////////////////////////////////
  // renormarize and sumation homo and inhomo parts
  ///////////////////////////////////

  if(renormarize == 1){
    static std::vector<dcomplex> old_den1((norb - Bas.ORMAS.ncore) * (norb - Bas.ORMAS.ncore), 0.0);
#pragma omp parallel default(shared)
    {
      int ithr = omp_get_thread_num();
      int llk, ulk;
      MPIP.omp_divide(ithr, 0, ksize, llk, ulk);
      int ncore = Bas.ORMAS.ncore;
      std::vector<dcomplex> old_rhok(norb * norb, 0.0);
      std::vector<dcomplex> new_rhok(norb * norb, 0.0);
      dcomplex old_norm;
      dcomplex new_norm;

      for(int ik = llk; ik < ulk; ik++){
	std::fill(old_rhok.begin(), old_rhok.end(), 0.0);
	std::fill(new_rhok.begin(), new_rhok.end(), 0.0);
	old_norm = 0.0;
	new_norm = 0.0;
	for(int iorb = 0; iorb < norb; ++iorb){
	  for(int jorb = 0; jorb < norb; ++jorb){
	    old_rhok[iorb * norb + jorb] = conj(opes[iorb * ksize + ik]) * opes[jorb * ksize + ik];
	    new_rhok[iorb * norb + jorb] = conj(homo_opes[iorb * ksize + ik]) * homo_opes[jorb * ksize + ik];
	  }
	}
	for(int iorb = 0; iorb < ncore; ++iorb){
	  old_norm += 2.0 * old_rhok[iorb * norb + iorb];
	  new_norm += 2.0 * new_rhok[iorb * norb + iorb];
	}
	for(int iorb = ncore; iorb < norb; ++iorb){
	  for(int jorb = ncore; jorb < norb; ++jorb){
	    old_norm += old_den1[(iorb - ncore) * (norb - ncore) + (jorb - ncore)] * old_rhok[iorb * norb + jorb];
	    new_norm +=     den1[(iorb - ncore) * (norb - ncore) + (jorb - ncore)] * new_rhok[iorb * norb + jorb];
	  }
	}
	if(std::abs(std::abs(new_norm) - std::abs(old_norm)) > 10E-15){
	  //if(std::abs(new_norm) - std::abs(old_norm) > 10E-15){
	  for(int iorb = 0; iorb < norb; ++iorb){
	    homo_opes[iorb * ksize + ik] *= std::sqrt(std::abs(old_norm / new_norm));
	  }
	}
      }
    }
    old_den1 = den1;
  }

#pragma omp parallel default(shared)
  {
    int ithr = omp_get_thread_num();
    int llk, ulk;
    MPIP.omp_divide(ithr, 0, ksize, llk, ulk);

    for(int iorb = 0; iorb < norb; ++iorb){
      ////////////////////////
      // omp parallelized
      ////////////////////////
      for(int ik = llk; ik < ulk; ik++){
	opes[iorb * ksize + ik] = homo_opes[iorb * ksize + ik] + inhomo_opes[iorb * ksize + ik];
	//test9 10  opes[iorb * ksize + ik] += inhomo_opes[iorb * ksize + ik];
      }
    }
  }


  return;
}
//////////////////////////////////////////////////////////////////////////
void Integ_EI::calc_pade1(std::vector<dcomplex> &eval, std::vector<dcomplex> &pade)
{
  // calc exp(x)
  
  int norder = 2;
  std::vector<double> coef_nume(norder + 1);
  coef_nume[0] = -60;
  coef_nume[1] = -24;
  coef_nume[2] = -3;

  int dorder = 3;
  std::vector<double> coef_deno(dorder + 1);
  coef_deno[0] = -60;
  coef_deno[1] = 36;
  coef_deno[2] = -9;
  coef_deno[3] = 1;

  dcomplex tmp_nume = 0.0;
  dcomplex tmp_deno = 0.0;
  for(int iorb = 0; iorb < norb; ++iorb){
    tmp_nume = 0.0;
    tmp_deno = 0.0;
    for(int j = 0; j < norder + 1; ++j){
      tmp_nume += coef_nume[j] * std::pow(eval[iorb], j);
    }
    for(int j = 0; j < dorder + 1; ++j){
      tmp_deno += coef_deno[j] * std::pow(eval[iorb], j);
    }
    pade[iorb] = tmp_nume / tmp_deno;
  }
  
  return;
}
//////////////////////////////////////////////////////////////////////////
void Integ_EI::calc_pade2(std::vector<dcomplex> &eval, std::vector<dcomplex> &pade)
{
  // calc (exp(x) - 1) / x
  
  int norder = 2;
  std::vector<double> coef_nume(norder + 1);
  coef_nume[0] = -60;
  coef_nume[1] = 6;
  coef_nume[2] = -1;

  int dorder = 3;
  std::vector<double> coef_deno(dorder + 1);
  coef_deno[0] = -60;
  coef_deno[1] = 36;
  coef_deno[2] = -9;
  coef_deno[3] = 1;

  dcomplex tmp_nume = 0.0;
  dcomplex tmp_deno = 0.0;
  for(int iorb = 0; iorb < norb; ++iorb){
    tmp_nume = 0.0;
    tmp_deno = 0.0;
    for(int j = 0; j < norder + 1; ++j){
      tmp_nume += coef_nume[j] * std::pow(eval[iorb], j);
    }
    for(int j = 0; j < dorder + 1; ++j){
      tmp_deno += coef_deno[j] * std::pow(eval[iorb], j);
    }
    pade[iorb] = tmp_nume / tmp_deno;
  }
  

  return;
}
//////////////////////////////////////////////////////////////////////////
void Integ_EI::inv_mat(std::vector<dcomplex> &inv_mat)
{
  //static int int dim = norb;
  //static int int lda= norb;
  //  static std::vector<int long> ipiv(norb);
  //static int int lwork = norb;
  //static std::vector<dcomplex>  work(norb*norb);
  //static int int info;
  
// Sato_TSURFF
  int norbl = norb;
  futil_gmatinv_(&norbl, &ZERO, &inv_mat[0], &inv_mat[0]);
//  zgetrf(&dim, &dim, &inv_mat[0], &lda, &ipiv[0], &info);
//  if(info != 0){
//    std::cout << "zgetrf_ failed at Integ_EI::prop" << std::endl;
//    abort();
//  }
//
//  zgetri(&dim, &inv_mat[0], &lda, &ipiv[0], &work[0], &lwork, &info);
//  if(info != 0){
//    std::cout << "zgetri_ failed at Integ_EI::prop" << std::endl;
//    abort();
//  }
// Sato_TSURFF
  
  return;
}		       
//////////////////////////////////////////////////////////////////////////
void Integ_EI::eigen_decompose(double dt, const std::vector<dcomplex> &v2xmat, std::vector<dcomplex> &eval,  
			       std::vector<dcomplex> &evl,  std::vector<dcomplex> &evr)
{
  static std::vector<dcomplex> dt_v2xmat(norb*norb);
  for(int iorb = 0; iorb < norb; ++iorb){
    for(int jorb = 0; jorb < norb; ++jorb){
      dt_v2xmat[iorb * norb + jorb] = v2xmat[iorb * norb + jorb] * dt;
    }
  }

  //static int int dim = norb;
  //static int int lda= norb;

  //static int int lwork = 2 * norb * norb;
  //static std::vector<dcomplex>  work(lwork);
  //static std::vector<double> rwork(2 * norb);
  //static int int info;

// Sato_TSURFF
//  std::cout << "surff::integ_EI: CHECK here!" << std::endl; abort();
//  zgeev("V", "V", &dim, &dt_v2xmat[0], &lda, &eval[0], &evl[0], &dim, &evr[0], &dim, &work[0], &lwork, &rwork[0], &info);
  int norbl = norb;
  bool dsc = false;
  futil_gdiag_comp_(&dsc, &norbl, &dt_v2xmat[0], &evl[0], &evr[0]);
  for (int iorb = 0; iorb < norb; ++iorb) eval[iorb] = dt_v2xmat[iorb*norb+iorb];
// Sato_TSURFF

  return;
}
//////////////////////////////////////////////////////////////////////////

//   // DEBUG previously in prop

//   for(int i = 0; i < 5; ++i){
//     for(int j = 0; j < 5; ++j){
//       dt_v2xmat[i * norb + j] = v2xmat[i * norb + j].real();
//     }
//   }


//   norb = 3;
//   dt_v2xmat.resize(norb*norb);

//   for(int i = 0; i < norb; ++i){
//     for(int j = 0; j < norb; ++j){
//       printf("%+20.10E %+20.10Ej  ", dt_v2xmat[i * norb + j].real(), dt_v2xmat[i * norb + j].imag());
//     }
//     printf("\n");
//   }
//   printf("\n");


//   static int int dim = norb;
//   static int int lda= norb;
//   static std::vector<dcomplex> eval(norb);
//   static std::vector<dcomplex> evl(norb*norb);
//   static std::vector<dcomplex> evr(norb*norb);

//   static int int lwork = norb * norb;
//   static std::vector<dcomplex>  work(lwork);
//   static std::vector<double> rwork(2 * norb);

//   static int int info;
//   zgeev("V", "V", &dim, &dt_v2xmat[0], &lda, &eval[0], &evl[0], &dim, &evr[0], &dim, &work[0], &lwork, &rwork[0], &info);
//   if(info != 0){
//     std::cout << "zzeev_ failed at Integ_EI::prop" << std::endl;
//     abort();
//   }


//   for(int i = 0; i < norb; ++i){
//     for(int j = 0; j < norb; ++j){
//       printf("%+20.10E  ", v2xmat[i * norb + j].real());
//     }
//     printf("\n");
//   }
//   printf("\n");


//   for(int i = 0; i < norb; ++i){
//     printf("%+20.10E  %+20.10Ej", eval[i].real(), eval[i].imag());
//     printf("\n");
//   }
//   printf("\n");

//   for(int i = 0; i < norb; ++i){
//     for(int j = 0; j < norb; ++j){
//       printf("%+20.10E  %+20.10Ej", evl[i * norb + j].real(),  evl[i * norb + j].imag());
//     }
//     printf("\n");
//   }
//   printf("\n");

//   for(int i = 0; i < norb; ++i){
//     for(int j = 0; j < norb; ++j){
//       printf("%+20.10E  %+20.10Ej", evr[i * norb + j].real(), evr[i * norb + j].imag());
//     }
//     printf("\n");
//   }
//   printf("\n");
//   abort();

//////////////////////////////////////////////////////////////////////////


