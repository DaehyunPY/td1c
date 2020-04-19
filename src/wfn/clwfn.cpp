////////////////////////////////////////////////////////////////////////
// Wavefunction
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
// Orimo_ECS
#include "surff.hpp"
// Orimo_ECS
////////////////////////////////////////////////////////////////////////
clwfn::clwfn()
{
// Sato_tSURFF
  dosurff = false;
// // Orimo_ECS
//   surff = NULL;
// // Orimo_ECS
// Sato_tSURFF
}
////////////////////////////////////////////////////////////////////////
// Sato_tSURFF
clwfn::clwfn(bool dosurff_)
{
  dosurff = dosurff_ && clcontrol::tsurff;
}
// Sato_tSURFF
////////////////////////////////////////////////////////////////////////
clwfn::clwfn(const clbas& Bas, const clwfn& tWfn)
{
  copy(Bas, tWfn);
}
////////////////////////////////////////////////////////////////////////
clwfn::clwfn(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
// Sato_tSURFF
  dosurff = false;
// Sato_tSURFF
  gen(MPIP, IO, Bas);
  if (IO.iprint > 0) print();
}
////////////////////////////////////////////////////////////////////////
clwfn::clwfn(const clmpi& MPIP, const clio& IO, const clbas& Bas, bool dosurff_)
{
// Sato_tSURFF
  dosurff = dosurff_ && clcontrol::tsurff;
// Sato_tSURFF
  gen(MPIP, IO, Bas);
  if (IO.iprint > 0) print();
}
////////////////////////////////////////////////////////////////////////
clwfn::~clwfn()
{
// Sato_tSURFF
  if (dosurff) delete surff;
// // Orimo_ECS
//   if(surff != NULL){
//     //std::cout << "surff1 : "  << surff << NULL << std::endl;
//     delete surff;
//   }
//   surff = NULL;
// // Orimo_ECS
// Sato_tSURFF
}
////////////////////////////////////////////////////////////////////////
//void clwfn::init()
//{
//  printf("clwfn::init() is nyi.");
//  abort();
//}
//////////////////////////////////////////////////////////////////////////
//void clwfn::final()
//{
//  printf("clwfn::init() is nyi.");
//  abort();
//}
////////////////////////////////////////////////////////////////////////
void clwfn::clear(const clbas& Bas)
{
  zclear_omp_(&size1, &wfn[0]);
  zclear_omp_(&sizeg, &wfng[0]);
  zclear_omp_(&size2, &wfn[size1]);
}
////////////////////////////////////////////////////////////////////////
void clwfn::clearo(const clbas& Bas)
{
  zclear_omp_(&size1, &wfn[0]);
}
////////////////////////////////////////////////////////////////////////
void clwfn::clearc(const clbas& Bas)
{
  zclear_omp_(&size2, &wfn[size1]);
}
////////////////////////////////////////////////////////////////////////
void clwfn::copy(const clbas& Bas, const clwfn& tWfn)
{
  adjust(Bas, tWfn);
  wfn = tWfn.wfn;
  wfng = tWfn.wfng;
}
////////////////////////////////////////////////////////////////////////
void clwfn::adjust(const clbas& Bas, const clwfn& tWfn)
{
  size1 = tWfn.size1;
  size2 = tWfn.size2;
  size = tWfn.size;
  sizeg = tWfn.sizeg;
  alloc(Bas);
// Sato_tSURFF
  if (dosurff != tWfn.dosurff) {
    std::cout << "clwfn: error in adjust." << std::endl;
    abort();
  }
// Sato_tSURFF
}
////////////////////////////////////////////////////////////////////////
void clwfn::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  read_info(IO);
  size1 = Bas.ORMAS.nfun * Bas.nbas;
  size2 = Bas.ORMAS.lcic;
  sizeg = Bas.ORMAS.nfun * Bas.ngrid;
  size = size1 + size2;

// Sato_tSURFF
// Orimo_ECS
  //surff
  if (dosurff) surff = new Surff(MPIP, IO, Bas);
  //surff
// Orimo_ECS
// Sato_tSURFF

  alloc(Bas);
}
////////////////////////////////////////////////////////////////////////
void clwfn::get_nradfc(const clmpi& MPIP, const clio& IO, const clbas& Bas)
{
  // determine Bas.nradfc
  // it is better to remove
  // const from clbas& Bas...
  wfn_get_nradfc_(&wfn[0]);
}
////////////////////////////////////////////////////////////////////////
void clwfn::alloc(const clbas& Bas)
{
  wfn.resize(size);
  wfng.resize(sizeg);
}
////////////////////////////////////////////////////////////////////////
void clwfn::print() const
{
  printf("# clwfn: size1 = %10ld\n", size1);
  printf("# clwfn: size2 = %10ld\n", size2);
  printf("# clwfn: size  = %10ld\n", size);
  printf("# clwfn: sizeg = %10ld\n", sizeg);
}
////////////////////////////////////////////////////////////////////////
void clwfn::read(const clmpi& Proc, const clio& IO, const clbas& Bas)
{
  // full initialization
  zclear_omp_(&size, &wfn[0]);
  zclear_omp_(&sizeg, &wfng[0]);

  FILE *fpo, *fpc;
  double wfnr, wfni;

  // read in orbitals
  long nrad0, lmax0, mmax0, nfun0;
  fpo = fopen(IO.orb.c_str(), "r");
  fscanf(fpo, "%ld\n", &nrad0);
  fscanf(fpo, "%ld\n", &lmax0);
  fscanf(fpo, "%ld\n", &mmax0);
  fscanf(fpo, "%ld\n", &nfun0);

  long ind;
  for (long ifun = 0; ifun < nfun0; ifun ++) {
    for (long l = 0; l <= lmax0; l ++) {
      for (long irad = 1; irad < nrad0; irad ++) {
//	fscanf(fpo, "%lf%lf\n", &wfnr, &wfni); 
	fscanf(fpo, "%le%le\n", &wfnr, &wfni); 
	if (ifun < Bas.ORMAS.nfun && l <= Bas.GAng.lmax1 && irad < Bas.GRad.nrad) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
          	                 + l * (Bas.GRad.nrad - 1) + irad - 1;
	  wfn[ind] = wfnr * RUNIT + wfni * IUNIT;
	}
      }
    }
  }
  fclose(fpo);

  //NEW
  bool symlm;
  IO.read_info("symlm", false, symlm);
  if (symlm) hprod_symlm_(&wfn[0]);
  //NEW

  // read in ci coefficients
  long lcic0;
  fpc = fopen(IO.cic.c_str(), "r");
  fscanf(fpc, "%ld\n", &lcic0);
  if (Bas.ORMAS.lcic != lcic0) {
    ormas_cic0_(&wfn[size1]);
    printf("clwfn::read-cic: cic was not read\n");
  } else {
    for (long idet = 0; idet < Bas.ORMAS.lcic; idet ++) {
      fscanf(fpc, "%le%le\n", &wfnr, &wfni);
      wfn[size1+idet] = wfnr * RUNIT + wfni * IUNIT;
    }
  }
  fclose(fpc);
}
////////////////////////////////////////////////////////////////////////
void clwfn::read_orb(const clmpi& Proc, const clio& IO, const clbas& Bas)
{
  FILE *fpo;
  double wfnr, wfni;

  // read in orbitals
  long nrad0, lmax0, mmax0, nfun0;
  fpo = fopen(IO.orb.c_str(), "r");
  fscanf(fpo, "%ld", &nrad0);
  fscanf(fpo, "%ld", &lmax0);
  fscanf(fpo, "%ld", &mmax0);
  fscanf(fpo, "%ld", &nfun0);

  long ind;
  for (long ifun = 0; ifun < nfun0; ifun ++) {
    for (long l = 0; l <= lmax0; l ++) {
      for (long irad = 1; irad < nrad0; irad ++) {
//	fscanf(fpo, "%lf%lf", &wfnr, &wfni); 
	fscanf(fpo, "%le%le", &wfnr, &wfni); 
//	if (l <= Bas.GAng.lmax1 && irad < Bas.GRad.nrad) {
	if (ifun < Bas.ORMAS.nfun && l <= Bas.GAng.lmax1 && irad < Bas.GRad.nrad) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
          	                 + l * (Bas.GRad.nrad - 1) + irad - 1;
	  wfn[ind] = wfnr * RUNIT + wfni * IUNIT;
	}
      }
    }
  }
  fclose(fpo);

  //NEW
  bool symlm;
  IO.read_info("symlm", false, symlm);
  if (symlm) hprod_symlm_(&wfn[0]);
  //NEW
}
////////////////////////////////////////////////////////////////////////
void clwfn::read_cic(const clmpi& Proc, const clio& IO, const clbas& Bas)
{
  FILE *fpc;
  double wfnr, wfni;

  // read in ci coefficients
  long lcic0;
  fpc = fopen(IO.cic.c_str(), "r");
  fscanf(fpc, "%ld", &lcic0);

  if (Bas.ORMAS.lcic != lcic0) {
    ormas_cic0_(&wfn[size1]);
    printf("clwfn::read-cic: cic was not read\n");
  } else {
    for (long idet = 0; idet < Bas.ORMAS.lcic; idet ++) {
      fscanf(fpc, "%le%le\n", &wfnr, &wfni);
      wfn[size1+idet] = wfnr * RUNIT + wfni * IUNIT;
    }
  }
  fclose(fpc);
}
////////////////////////////////////////////////////////////////////////
void clwfn::write(const clmpi& Proc, const clio& IO, const clbas& Bas) const
{
  FILE *fpo, *fpc;
  double wfnr, wfni;
  // write out orbitals
  fpo = fopen(IO.orb.c_str(), "w");
  fprintf(fpo, "%10ld\n", Bas.GRad.nrad);
  fprintf(fpo, "%10ld\n", Bas.GAng.lmax1);
  fprintf(fpo, "%10ld\n", Bas.GAng.mmax1);
  fprintf(fpo, "%10ld\n", Bas.ORMAS.nfun);

  long ind;
  for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
    for (long l = 0; l <= Bas.GAng.lmax1; l ++) {
      for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
	ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
          	               + l * (Bas.GRad.nrad - 1) + irad - 1;
	wfnr = real(wfn[ind]);
	wfni = imag(wfn[ind]);
//	fprintf(fpo, "%25.15lf%25.15lf\n", wfnr, wfni);
	fprintf(fpo, "%25.15le%25.15le\n", wfnr, wfni);
      }
    }
  }
  fclose(fpo);

  fpc = fopen(IO.cic.c_str(), "w");
  fprintf(fpc, "%10ld\n", Bas.ORMAS.lcic);
  for (long idet = 0; idet < Bas.ORMAS.lcic; idet ++) {
      wfnr = real(wfn[size1+idet]);
      wfni = imag(wfn[size1+idet]);
//    fprintf(fpc, "%25.15lf%25.15lf\n", wfnr, wfni); 
      fprintf(fpc, "%25.15le%25.15le\n", wfnr, wfni); 
  }
  fclose(fpc);
}
////////////////////////////////////////////////////////////////////////
void clwfn::read_info(const clio& IO)
{
  IO.read_info("orth_type", "schmidt", orth_type);
  IO.read_info("ci_normalize", false, ci_normalize);
}
////////////////////////////////////////////////////////////////////////
void clwfn::orth(const clbas& Bas)
{
  wfn_orth_cic_(&wfn[size1]);
  if (orth_type.compare("schmidt") == 0) {
    wfn_orth_orb_(&wfn[0]);
  } else if (orth_type.compare("symmetric") == 0) {
    wfn_orth_orb_symm_(&wfn[0]);
  }
}
////////////////////////////////////////////////////////////////////////
void clwfn::ortho(const clbas& Bas)
{
  if (orth_type.compare("schmidt") == 0) {
    wfn_orth_orb_(&wfn[0]);
  } else if (orth_type.compare("symmetric") == 0) {
    wfn_orth_orb_symm_(&wfn[0]);
  }
}
////////////////////////////////////////////////////////////////////////
void clwfn::orthc(const clbas& Bas)
{
  wfn_orth_cic_(&wfn[size1]);
}
////////////////////////////////////////////////////////////////////////
void clwfn::mask(const clbas& Bas)
{
  if (ci_normalize) {
    wfn_orth_cic_(&wfn[size1]);
  }
  if (Bas.GRad.ecs_flag == 1) {
    return;
  } else {
    wfn_mask_(&wfn[0]);
  }
}
////////////////////////////////////////////////////////////////////////
void clwfn::proj(const clbas& Bas, const clwfn& WfnP)
{
  wfn_proj_(&WfnP.wfn[0], &wfn[0]);
}
////////////////////////////////////////////////////////////////////////
void clwfn::projg(const clmpi& Proc, const clbas& Bas, const clwfn& WfnP)
{
  wfn_projg_(&WfnP.wfng[0], &wfng[0]);
}
////////////////////////////////////////////////////////////////////////
void clwfn::ladapt(const clmpi& Proc, const clbas& Bas)
{
  if (Bas.lconst) {
    wfn_ladapt_(&wfn[0], &wfn[size1]);
  } else if (Bas.lconst_core) {
    wfn_ladapt_core_(&wfn[0]);
  }
}
////////////////////////////////////////////////////////////////////////
void clwfn::madapt(const clmpi& Proc, const clbas& Bas)
{
  //  std::cout << "wfn:madapt: DISABLE madapt!" << std::endl;
  wfn_madapt_(&wfn[size1]);
}
////////////////////////////////////////////////////////////////////////
void clwfn::print_cic(std::string FNAME, const clbas& Bas) const
{
  ormas_cic_print_(&wfn[size1], FNAME.c_str(), FNAME.length());
}
////////////////////////////////////////////////////////////////////////
void clwfn::print_orb(long npt, std::string FNAME, const clbas& Bas) const
{
  FILE *fpo;
  fpo = fopen(FNAME.c_str(), "w");

  long ind, irad, num_dvr;
  double x_ll, x_ul, dpt;
  double point, orbr, orbi;
  dcomplex wfn_val;
  std::vector<dcomplex> bas_val(Bas.GRad.nmax + 1);

  for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
    fprintf(fpo, "# ifun = %10ld\n", ifun);
    for (long ife = 0; ife < Bas.GRad.nfe; ife ++) {
      x_ll = Bas.GRad.get_x0(ife);
      x_ul = Bas.GRad.get_x1(ife);
      dpt = (x_ul - x_ll) / npt;
      num_dvr = Bas.GRad.get_ndvr(ife);

      point = x_ll;
      while (point + dpt * HALF < x_ul) {

	for(long idvr = 0; idvr <= num_dvr; idvr ++) {
	  irad = Bas.GRad.mapf[ife] + idvr;
	  //bas_val[idvr] = Bas.GRad.get_val(irad, point);
	  bas_val[idvr] = Bas.GRad.get_val(ife, idvr, point);
	  //DEBUG
	  // if (ifun == 0) {
	  //	printf("%10ld%10ld%10ld%20.10f%20.10f%20.10f\n", ife, idvr, irad, point, bas_val[idvr].real(), bas_val[idvr].imag());
	  // }
	  //DEBUG
	}

	fprintf(fpo, "%20.10f", point);
	for (long l = 0; l <= Bas.GAng.lmax1; l ++) {
	  wfn_val = CZERO;
	  for(long idvr = 0; idvr <= num_dvr; idvr ++) {
	    //	  for(long idvr = 0; idvr < num_dvr; idvr ++) {
	    irad = Bas.GRad.mapf[ife] + idvr;
	    if (clcontrol::fedvr_normalized) {
	      if (irad > 0 && irad < Bas.GRad.nrad) {
		ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
             	                        + l * (Bas.GRad.nrad - 1) + irad - 1;
		wfn_val += wfn[ind] * bas_val[idvr] / sqrt(Bas.GRad.wrad[irad]);
	      }
	    } else {
	      if (irad > 0 && irad < Bas.GRad.nrad) {
		ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
             	                        + l * (Bas.GRad.nrad - 1) + irad - 1;
		wfn_val += wfn[ind] * bas_val[idvr];
	      }
	    }
	  }
//	  if (clcontrol::fedvr_normalized) {
//	    orbr = real(wfn_val / sqrt(Bas.GRad.wrad[irad]));
//	    orbi = imag(wfn_val / sqrt(Bas.GRad.wrad[irad]));
//	  } else {
	    orbr = real(wfn_val);
	    orbi = imag(wfn_val);
//	  }
	  fprintf(fpo, "%15.5e%15.5e", orbr, orbi);
	}
	fprintf(fpo, "\n");
	point += dpt;
      }
    }
    fprintf(fpo, "\n");
    fprintf(fpo, "\n");
  }

  fclose(fpo);
}
////////////////////////////////////////////////////////////////////////
void clwfn::print_wang(const clbas& Bas) const
{
  long ind;
  dcomplex tmp;

  long size = Bas.ORMAS.nfun * Bas.GAng.nsph1;
  std::vector<double> wlm(size);

  for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
    //    printf("%5ld", ifun);
    for (long l = get_abs(Bas.mval[ifun]); l <= Bas.GAng.lmax1; l ++) {
      tmp = CZERO;
      if (clcontrol::fedvr_normalized) {
	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1)
       	                        +  l * (Bas.GRad.nrad - 1) + irad - 1;
	  tmp += conj(wfn[ind]) * wfn[ind];
	}
      } else {
	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1)
       	                        +  l * (Bas.GRad.nrad - 1) + irad - 1;
	  tmp += conj(wfn[ind]) * wfn[ind] * Bas.GRad.wrad[irad];
	}
      }
      wlm[ifun * Bas.GAng.nsph1 + l] = sqrt(real(tmp));
    }
  }

//  printf("    l");
//  printf("    m");
//  for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
//    printf("%15ld", ifun);
//  }
//  printf("\n");
  for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
    printf("%10ld%10ld", ifun, Bas.mval[ifun]);
    for (long l = 0; l <= Bas.GAng.lmax1; l ++) {
      printf("%20.10e", wlm[ifun * Bas.GAng.nsph1 + l]);
    }
    printf("\n");
  }
}
////////////////////////////////////////////////////////////////////////
