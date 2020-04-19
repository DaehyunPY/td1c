////////////////////////////////////////////////////////////////////////
// Sph: basic
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
#include <fftw3.h>
//#include <shtns.h>
////////////////////////////////////////////////////////////////////////
clsph::clsph()
{
}
////////////////////////////////////////////////////////////////////////
clsph::clsph(const clmpi& MPIP, const clio& IO)
{
  gen(MPIP, IO);
}
////////////////////////////////////////////////////////////////////////
clsph::~clsph()
{
  std::cout << "~clsph" << std::endl;
}
////////////////////////////////////////////////////////////////////////
void clsph::print(const clio& IO) const
{
  //  printf("# sph.nlat = %10ld sph.nphi = %10ld sph.nang = %10ld\n", nlat, nphi, nang);

  printf("# SPH: Nodes and weights for %ld latitudes:\n", nlat);
  for (long ilat = 0; ilat < nlat; ilat ++) {
    double theta = acos(cost[ilat]) / PI * 180.0;
    printf("%10ld%20.10f%20.10f%20.10f%20.10f\n", ilat, theta, cost[ilat], sint[ilat], wgt_lat[ilat]);
  }

  long indf, indb;
  printf("# SPH: forward Legendre transform: legf1\n");
  for (long m = -mmax1; m <= mmax1; m ++) {
    printf("m = %10ld\n", m);
    for (long l = get_abs(m); l <= lmax1; l ++) {
  	printf("l: %10ld", l);
  	for (long ilat = 0; ilat < nlat; ilat ++) {
	  //  	  indf = (m * nlat + ilat) * (lmax1 + 1) + l;
	  indf = (mmax1 + m) * nlat * (lmax1 + 1)  // c++: indf[-mmax:mmax][1:nlat][0:lmax]
                             + ilat * (lmax1 + 1)  // f90: indf(0:lmax, 1:nlat, -mmax:mmax)
  	                                     + l;
  	  printf("%20.10f", legf1[indf]);
  	}
  	printf("\n");
    }
  }
  printf("# SPH: backward Legendre transform: legb1\n");
  for (long m = -mmax1; m <= mmax1; m ++) {
    printf("m = %10ld\n", m);
    for (long ilat = 0; ilat < nlat; ilat ++) {
      printf("ilat: %10ld", ilat);
      for (long l = get_abs(m); l <= lmax1; l ++) {
	//	indb = (m * (lmax1 + 1) + l) * nlat + ilat;
	indb = (mmax1 + m) * (lmax1 + 1) * nlat  // c++: indb[-mmax:mmax][0:lmax][1:nlat]
                                    + l  * nlat  // f90: indb(1:nlat, 0:lmax, -mmax:mmax)
	                                 + ilat;
	printf("%20.10f", legb1[indb]);
      }
      printf("\n");
    }
  }

  printf("# SPH: forward Legendre transform: legf2\n");
  for (long m = -mmax2; m <= mmax2; m ++) {
    printf("m = %10ld\n", m);
    for (long l = get_abs(m); l <= lmax2; l ++) {
  	printf("l: %10ld", l);
  	for (long ilat = 0; ilat < nlat; ilat ++) {
	  //  	  indf = (m * nlat + ilat) * (lmax2 + 1) + l;
	  indf = (mmax2 + m) * nlat * (lmax2 + 1)  // c++: indf[-mmax:mmax][1:nlat][0:lmax]
                             + ilat * (lmax2 + 1)  // f90: indf(0:lmax, 1:nlat, -mmax:mmax)
  	                                     + l;
  	  printf("%20.10f", legf2[indf]);
  	}
  	printf("\n");
    }
  }
  printf("# SPH: backward Legendre transform: legb2\n");
  for (long m = -mmax2; m <= mmax2; m ++) {
    printf("m = %10ld\n", m);
    for (long ilat = 0; ilat < nlat; ilat ++) {
      printf("ilat: %10ld", ilat);
      for (long l = get_abs(m); l <= lmax2; l ++) {
	//	indb = (m * (lmax2 + 1) + l) * nlat + ilat;
	indb = (mmax2 + m) * (lmax2 + 1) * nlat  // c++: indb[-mmax:mmax][0:lmax][1:nlat]
                                    + l  * nlat  // f90: indb(1:nlat, 0:lmax, -mmax:mmax)
	                                 + ilat;
	printf("%20.10f", legb2[indb]);
      }
      printf("\n");
    }
  }
}
////////////////////////////////////////////////////////////////////////
void clsph::gen(const clmpi& MPIP, const clio& IO)
{
  long m1max = read_mmax1(IO); 
  IO.read_info("lmax1", lmax1);
  IO.read_info("lmax2", lmax2);
  IO.read_info("mmax1", m1max,        mmax1);
  IO.read_info("mmax2", m1max * LTWO, mmax2);

  long nlat_def = lmax2 + 1;
  long nphi_def = 1;
  IO.read_info("nlat", nlat_def, nlat);
  IO.read_info("nphi", nphi_def, nphi);
  if (nphi != 1) {std::cout << "nphi should be 1!" << std::endl; abort();}
  nang = nlat * nphi;
  nsph1 = lmax1 + 1;
  nsph2 = lmax2 + 1;

  wgt_lat.resize(nlat);
  wgt_ang.resize(nang);
  cost.resize(nlat);
  sint.resize(nlat);

  long size1 = nlat * (lmax1 + 1) * (2 * mmax1 + 1);
  long size2 = nlat * (lmax2 + 1) * (2 * mmax2 + 1);
  legf1.resize(size1);
  legb1.resize(size1);
  legf2.resize(size2);
  legb2.resize(size2);

  //DEBUG printf("# SPH: SPH_THETA_CONTIGUOUS = %d\n", SHT_THETA_CONTIGUOUS);
  //DEBUG printf("# SPH: SPH_PHI_CONTIGUOUS   = %d\n", SHT_PHI_CONTIGUOUS);
  //DEBUG printf("# SPH: NSPAT_ALLOC = %d\n", NSPAT_ALLOC(shtns));

  sph_gen_(&lmax1, &lmax2, &mmax1, &mmax2, &nlat, &nphi,
	   &wgt_lat[0], &wgt_phi, &wgt_ang[0], &cost[0], &sint[0],
	   &legf1[0], &legb1[0], &legf2[0], &legb2[0]);
  if (IO.iprint > 4) print(IO);
}
////////////////////////////////////////////////////////////////////////
long clsph::read_mmax1(const clio& IO)
//
// Read m value of each orbital
//
{
  long norb;
  IO.read_info("nfun", norb);

  std::string key = "orbital:";
  std::ifstream ifs(IO.inp.c_str());
  ifs.seekg(0, std::ios::beg);

  std::string line;
  std::stringstream ioss;
  long ntmp, ltmp, mtmp, ttmp;

  long m1max = 0;

  while ( getline(ifs, line) && line.find(key,0) == std::string::npos ) {}
  if(! ifs.eof()) {
    for (long ifun = 0; ifun < norb; ifun ++) {
      getline(ifs, line);
      ioss.str("");
      ioss.clear(std::stringstream::goodbit);
      ioss << line.c_str();
      ioss >> ntmp
	   >> ltmp
	   >> mtmp
	   >> ttmp;
      m1max = std::max(m1max, get_abs(mtmp));
    }
  } else {
    std::cout<< "clbas::read_mmax1: No mval information." << std::endl;
    abort();
  }

  return m1max;
}
////////////////////////////////////////////////////////////////////////
