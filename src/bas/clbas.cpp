////////////////////////////////////////////////////////////////////////
// Basis
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
// static members
////////////////////////////////////////////////////////////////////////
int clbas::num_bas = 0;
////////////////////////////////////////////////////////////////////////
clbas::clbas()
{
  num_bas ++;
  std::cout << "# clbas: num_bas = " << num_bas << std::endl;
}
////////////////////////////////////////////////////////////////////////
clbas::clbas(const clmpi& MPIP, const clio& IO)
{
  num_bas ++;
  std::cout << "# clbas: num_bas = " << num_bas << std::endl;
  read_info(IO);
  gen(MPIP, IO);
  if (IO.iprint > 0) print(IO);
}
////////////////////////////////////////////////////////////////////////
clbas::~clbas()
{
  std::cout << "~clbas" << std::endl;
  if (clcontrol::psp) bas_pp_final_();
  if (clcontrol::exact3j) sph_3j_final_();
}
////////////////////////////////////////////////////////////////////////
void clbas::read_info(const clio& IO)
{
  IO.read_info("znuc", znuc);
  IO.read_info("smul", smul);
  IO.read_info("ltot", ltot);
  IO.read_info("mtot", mtot);
  IO.read_info("lconst", false, lconst);
  IO.read_info("lconst_core", false, lconst_core);
  IO.read_info("psp_label", znuc, psp_label);
}
////////////////////////////////////////////////////////////////////////
void clbas::print(const clio& IO) const
{
  printf("# nbas  = %10d\n", nbas);
  printf("# ngrid = %10d\n", ngrid);
  if (IO.iprint > 1) GRad.print();
  if (IO.iprint > 1) GAng.print(IO);
}
////////////////////////////////////////////////////////////////////////
void clbas::gen(const clmpi& MPIP, const clio& IO)
{
  GRad.gen(MPIP, IO);
  GAng.gen(MPIP, IO);
  // c++/fortran bindings
  if (num_bas == 1) {
// tdcis-teramura
// Orimo_ECS
    bas_rad_bind_(&GRad.nfe, &GRad.nmax, &GRad.nrad, &GRad.nradfc, &GRad.mapf[0], &GRad.mapb[0],
		  &GRad.xrad[0], &GRad.wrad[0], &GRad.radp[0], &GRad.radk[0], &GRad.radk0[0],
		  &GRad.rmask, &GRad.mask[0],
		  &GRad.ecs_flag, &GRad.theta, &GRad.recs, &GRad.cwrad[0], &GRad.radkI_ecs, &GRad.irad_ecs,
		  &GRad.cxrad[0], &GRad.bra_wrad[0], &GRad.rdr[0], &GRad.wdw[0], 
		  &GRad.type_mkint1_sph, &GRad.type_mkint2_sph, &GRad.type_mkv2mf, &GRad.type_mkxmat_aa,
		  &GRad.inf_range, &GRad.irad_inf, &GRad.exp_factor, &GRad.switchoff, &GRad.irad_sw, &GRad.pot_type, &GRad.trunc_irad,
		  &GRad.nradgs);
//old    bas_rad_bind_(&GRad.nfe, &GRad.nmax, &GRad.nrad, &GRad.nradfc, &GRad.mapf[0], &GRad.mapb[0],
//old		  &GRad.xrad[0], &GRad.wrad[0], &GRad.radp[0], &GRad.radk[0], &GRad.radk0[0],
//old		  &GRad.rmask, &GRad.mask[0],
//old		  &GRad.ecs_flag, &GRad.theta, &GRad.recs, &GRad.cwrad[0], &GRad.radkI_ecs, &GRad.irad_ecs,
//old		  &GRad.cxrad[0], &GRad.bra_wrad[0]);
// Orimo_ECS
// tdcis-teramura
// tdcis-teramura
    bas_sph_bind_(&GAng.lmax1, &GAng.lmax2, &GAng.mmax1, &GAng.mmax2, &GAng.nsph1, 
		  &GAng.nlat, &GAng.nphi, &GAng.nang, &GAng.wgt_lat[0], &GAng.wgt_phi, 
		  &GAng.wgt_ang[0], &GAng.cost[0], &GAng.sint[0], 
		  &GAng.legf1[0], &GAng.legb1[0], &GAng.legf2[0], &GAng.legb2[0],
		  &GAng.lmax1gs, &GAng.lmax2gs);
// tdcis-teramura
  }
  if (clcontrol::exact3j) sph_3j_init_();

  ORMAS.gen(MPIP, IO);

  nbas = (GRad.nrad - 1) * GAng.nsph1;
  nbas2 = (GRad.nrad - 1) * GAng.nsph2;
  ngrid = (GRad.nrad - 1) * GAng.nang;

  // grid
  grid.resize(4 * ngrid);
  wgt.resize(ngrid);
  nval.resize(ORMAS.nfun);
  lval.resize(ORMAS.nfun);
  mval.resize(ORMAS.nfun);
  gen_grid();
  gen_wgt();
  read_lmval(IO);

  // lm coefficient
  alph_lm.resize((GAng.lmax1 + 1) * (2 * GAng.mmax1 + 1));

  // operator matrices
  int size0 = (2 * GRad.nmax + 1) * (GRad.nrad - 1);
  int size1 = (2 * GRad.nmax + 1) * (GRad.nrad - 1) * (GAng.lmax1 + 1);
  int size2 = (GRad.nmax + 1) * (GRad.nrad - 1) * (GAng.lmax2 + 1);
  pmat.resize(size0); // radial first derivative, transposed
  kmat.resize(size1); // kinetic
  tmat.resize(size1); // kinetic + nuclear
  d2ll.resize(size2); // lu-decomposed laplacian
  bas_zfac.resize  ((GRad.nrad - 1)*(GAng.lmax1 + 1)*(2*GAng.mmax1 + 1)); // for zprod
  bas_pzfac1.resize(                (GAng.lmax1 + 1)*(2*GAng.mmax1 + 1)); // for pzprod
  bas_pzfac2.resize((GRad.nrad - 1)*(GAng.lmax1 + 1)*(2*GAng.mmax1 + 1)); // for pzprod
  bas_azfac.resize ((GRad.nrad - 1)*(GAng.lmax1 + 1)*(2*GAng.mmax1 + 1)); // for aprod
  bas_d2fac1.resize(                 GAng.lmax2 + 1);  // for Poisson
  bas_d2fac2.resize((GRad.nrad - 1)*(GAng.lmax2 + 1)); // for Poisson
  bas_d2invr.resize((GRad.nrad - 1)*(GAng.lmax2 + 1)); // for Poisson
  bas_d2rpl0.resize((GRad.nrad - 1)*(GAng.lmax2 + 1)); // for Poisson
  bas_d2rpl1.resize((GRad.nrad - 1)*(GAng.lmax2 + 1)); // for Poisson
// Orimo_ECS
  bas_d2crpl1.resize((GRad.nrad - 1)*(GAng.lmax2 + 1)); // for Poisson
  int size3 = (3 * GRad.nmax + 1) * (GRad.nrad - 1) * (GAng.lmax2 + 1);
  int size4 = (GRad.nmax + 1) * (GRad.irad_ecs - 1) * (GAng.lmax2 + 1);
  d2ll_ecs.resize(size3);
  confd2ll.resize(size4); // lu-decomposed laplacian confined in irad_ecs
  ipiv_ecs.resize((GRad.nrad - 1)*(GAng.lmax2 + 1));
  d1mat.resize(size0); // ecs contour redial first  derivative
  d2mat.resize(size0); // ecs contour redial second derivative
// Orimo_ECS

  // c++/fortran bindings
  if (num_bas == 1) {
// Orimo_ECS
    bas_bas_bind_(&znuc, &smul, &ltot, &mtot, &nbas, &nbas2, &ngrid, &ORMAS.nfun, 
		  &nval[0], &lval[0], &mval[0], &grid[0], &wgt[0], &alph_lm[0], 
		  &pmat[0], &kmat[0], &tmat[0], &d2ll[0], &bas_zfac[0], &bas_pzfac1[0], 
		  &bas_pzfac2[0], &bas_azfac[0], &bas_d2fac1[0], &bas_d2fac2[0],
		  &bas_d2invr[0], &bas_d2rpl0[0], &bas_d2rpl1[0], 
		  &d2ll_ecs[0], &ipiv_ecs[0], &d1mat[0], &d2mat[0], &confd2ll[0], &bas_d2crpl1[0],
		  &psp_label);
//old    bas_bas_bind_(&znuc, &smul, &ltot, &mtot, &nbas, &nbas2, &ngrid, &ORMAS.nfun, 
//old		  &nval[0], &lval[0], &mval[0], &grid[0], &wgt[0], &alph_lm[0], 
//old		  &pmat[0], &kmat[0], &tmat[0], &d2ll[0], &bas_zfac[0], &bas_pzfac1[0], 
//old		  &bas_pzfac2[0], &bas_azfac[0], &bas_d2fac1[0], &bas_d2fac2[0],
//old		  &bas_d2invr[0], &bas_d2rpl0[0], &bas_d2rpl1[0], &psp_label);
// Orimo_ECS
  }

  if (clcontrol::psp) bas_pp_init_();
  bas_gen_alph_();
  bas_gen_pmat_(&pmat[0]);
  bas_gen_kmat_(&kmat[0]);
  bas_gen_tmat_(&tmat[0]);
  bas_gen_d2ll_(&d2ll[0]);
  bas_gen_zfac_();
  bas_gen_pzfac_();
  bas_gen_azfac_();
  bas_gen_d2fac_();
// Orimo_ECS
  //bas_gen_confd2ll_(&confd2ll[0]);
  //bas_gen_d2ll_ecs_(&d2ll_ecs[0]);
  //bas_gen_d2mat_(&d2mat[0]);
  bas_gen_d1mat_(&d1mat[0]);
// Orimo_ECS

  //DEBUG
  //  std::cout << "abort in clbas::gen for debug." << std::endl;
  //  abort();
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void clbas::read_lmval(const clio& IO)
//
// Read m value of each orbital
//
{
  std::string key = "orbital:";
  std::ifstream ifs(IO.inp.c_str());
  ifs.seekg(0, std::ios::beg);

  std::string line;
  std::stringstream ioss;
  int ttmp;

  while ( getline(ifs, line) && line.find(key,0) == std::string::npos ) {}
  if(! ifs.eof()) {
    for (int ifun = 0; ifun < ORMAS.nfun; ifun ++) {
      getline(ifs, line);
      ioss.str("");
      ioss.clear(std::stringstream::goodbit);
      ioss << line.c_str();
      ioss >> nval[ifun]
	   >> lval[ifun]
	   >> mval[ifun]
	   >> ttmp;
    }
    std::cout << "# lmval:" << std::endl;
    for (int ifun = 0; ifun < ORMAS.nfun; ifun ++) {
      printf("#%5d%5d%5d%5d\n", ifun, nval[ifun], lval[ifun], mval[ifun]);
    }
    //    abort();
  } else {
    std::cout<< "No lmval information." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clbas::gen_grid()
{
  double cosp, sinp, cost, sint, phi, rad;

  int igrid = 0;
  for (int iphi = 0; iphi < GAng.nphi; iphi ++) {
    phi = TWO * PI / GAng.nphi * iphi;
    cosp = cos(phi);
    sinp = sin(phi);
    for (int ilat = 0; ilat < GAng.nlat; ilat ++) {
      for (int irad = 1; irad < GRad.nrad; irad ++) {
	rad = GRad.xrad[irad];
	grid[            igrid] = rad * GAng.sint[ilat] * cosp; // x
	grid[    ngrid + igrid] = rad * GAng.sint[ilat] * sinp; // y
	grid[2 * ngrid + igrid] = rad * GAng.cost[ilat];        // z
	grid[3 * ngrid + igrid] = rad;                          // sqrt(x**2+y**2+z**2)
	igrid ++;
      }
    }
  }
  if (igrid != ngrid) {
    printf("Error in clbas_gen_grid: %10d != %10d\n", igrid, ngrid);
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clbas::gen_wgt()
{
  double wgt_ang;

  int igrid = 0;
  for (int iang = 0; iang < GAng.nang; iang ++) {
    for (int irad = 1; irad < GRad.nrad; irad ++) {
      wgt[igrid] = GAng.wgt_ang[iang];
      igrid ++;
    }
  }
  if (igrid != ngrid) {
    printf("Error in clgrid_gen_wgt: %10d != %10d\n", igrid, ngrid);
    abort();
  }
}
//////////////////////////////////////////////////////////////////////////
void clbas::ang2sph1(const clmpi& MPIP, 
		     const std::vector<dcomplex>& fang,
		     std::vector<dcomplex>& fsph) const
{
  bas_ang2sph1_(&fang[0], &fsph[0]);
}
//////////////////////////////////////////////////////////////////////////
void clbas::sph2ang1(const clmpi& MPIP,
		     const std::vector<dcomplex>& fsph,
		     std::vector<dcomplex>& fang) const
{
  bas_sph2ang1_(&fsph[0], &fang[0]);
}
//////////////////////////////////////////////////////////////////////////
void clbas::ang2sph2(const clmpi& MPIP,
		     const std::vector<dcomplex>& fang,
		     std::vector<dcomplex>& fsph) const
{
  bas_ang2sph2_(&fang[0], &fsph[0]);
}
//////////////////////////////////////////////////////////////////////////
void clbas::sph2ang2(const clmpi& MPIP,
		     const std::vector<dcomplex>& fsph,
		     std::vector<dcomplex>& fang) const
{
  bas_sph2ang2_(&fsph[0], &fang[0]);
}
//////////////////////////////////////////////////////////////////////////
void clbas::proj(const clbas& Bas2, const std::vector<dcomplex>& orb1, std::vector<dcomplex>& orb2) const
{
  int ife1, irad1, ind1, ind2;
  double rval, swgt, mat_proj;
  int nfun_tmp = std::max(ORMAS.nfun, Bas2.ORMAS.nfun);
  int lmax_tmp = std::max(GAng.lmax1, Bas2.GAng.lmax1);

  int size2 = Bas2.ORMAS.nfun * Bas2.nbas;
  zclear_omp_(&size2, &orb2[0]);

  for (int irad2 = 1; irad2 < Bas2.GRad.nrad; irad2 ++) {
    rval = Bas2.GRad.xrad[irad2];
    swgt = sqrt(Bas2.GRad.wrad[irad2]);
    ife1 = GRad.get_ife(rval);
    //    printf("%10d%20.10f%10d\n", irad2, Bas2.GRad.xrad[irad2], ife1);

    if (ife1 >= 0) {
      for (int m = std::max(LONE - ife1, LZERO); m < GRad.get_ndvr(ife1); m ++) {

	irad1 = GRad.mapf[ife1] + m;
	mat_proj = GRad.get_val(irad1, rval) * swgt;

	//	printf("%10d%10d%20.10f%20.10f\n", irad2, irad1, rval, GRad.xrad[irad1]);

	for (int ifun = 0; ifun < nfun_tmp; ifun ++) {
	  for (int l = 0; l <= lmax_tmp; l ++) {
	    ind1 = ifun * GAng.nsph1 * (GRad.nrad - 1) 
     	                        + l * (GRad.nrad - 1) + irad1 - 1;
	    ind2 = ifun * Bas2.GAng.nsph1 * (Bas2.GRad.nrad - 1) 
          	                     + l * (Bas2.GRad.nrad - 1) + irad2 - 1;
	    orb2[ind2] += orb1[ind1] * mat_proj;
	  }
	}
      }
    }

    if (ife1 < GRad.nfe - 1) {

      irad1 = GRad.mapf[ife1 + 1];
      mat_proj = GRad.get_val0(ife1, rval) * swgt;

      for (int ifun = 0; ifun < nfun_tmp; ifun ++) {
	for (int l = 0; l <= lmax_tmp; l ++) {
	  ind1 = ifun * GAng.nsph1 * (GRad.nrad - 1) 
     	                      + l * (GRad.nrad - 1) + irad1 - 1;
	  ind2 = ifun * Bas2.GAng.nsph1 * (Bas2.GRad.nrad - 1) 
          	                   + l * (Bas2.GRad.nrad - 1) + irad2 - 1;
	  orb2[ind2] += orb1[ind1] * mat_proj;
	}
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////
void clbas::ppgenkb(const std::vector<dcomplex>& orb) const
{
  bas_pp_genkb_(&orb[0]);
}
//////////////////////////////////////////////////////////////////////////
