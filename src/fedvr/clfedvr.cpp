//2016/10/22 Yuki Orimo Changed
////////////////////////////////////////////////////////////////////////
// FEDVR: basic
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
clfedvr::clfedvr()
{
}
////////////////////////////////////////////////////////////////////////
clfedvr::clfedvr(const clmpi& MPIP, const clio& IO)
{
  gen(MPIP, IO);
}
////////////////////////////////////////////////////////////////////////
clfedvr::~clfedvr()
{
  std::cout << "~clfedvr" << std::endl;
//  x0.resize(0);
//  x1.resize(0);
//  dfe.resize(0);
//  ndvr.resize(0);
//  DVR.resize(0);
//  mapf.resize(0);
//  mapb.resize(0);
//  mask.resize(0);
//
//  // 2014/05/10
//  xrad.resize(0);
//  wrad.resize(0);
//  cxrad.resize(0);
//  cwrad.resize(0);
//  bra_wrad.resize(0);
//
//  radk.resize(0);
//  radp.resize(0);
//  //old  kinetic.resize(0, 0, 0, 0);
//  //old  nabla.resize(0, 0, 0, 0);
//  //old  nabla_new.resize(0);
//  //old  kloc.resize(0);
}
//////////////////////////////////////////////////////////////////////////
//void clfedvr::init()
//{
//  //  printf("# WARNING: clfedvr::init() does nothing.\n");
//  //  abort();
//}
//////////////////////////////////////////////////////////////////////////
//void clfedvr::final()
//{
//  //  printf("# WARNING: clfedvr::final() does nothing.\n");
//  //  abort();
//}
//////////////////////////////////////////////////////////////////////////
void clfedvr::clear()
{
//  x0.clear();
//  x1.clear();
//  ndvr.clear();
//  DVR.clear();
//  mapf.clear();
//  mapb.clear();
//  mask.clear();
//  xrad.clear();
//  wrad.clear();
//  radk.clear();
//  radp.clear();
}
////////////////////////////////////////////////////////////////////////
void clfedvr::print() const
{
  long i, m;
  long iB0, iB;
  double DVRX, DVRW, RADW, DMSK;
  std::string flag;

//DEBUG  std::cout << "# FE-DVR: local --> global:" << std::endl;
//DEBUG  for (i = 0; i < nfe; i ++) {
//DEBUG    iB0 = mapf[i];
//DEBUG    for(m = 0; m <= ndvr[i]; m ++) {
//DEBUG      if (i == 0 && m == 0) {
//DEBUG	flag = "xL";
//DEBUG      } else if (i == nfe - 1 && m == ndvr[i]) {
//DEBUG	flag = "xR";
//DEBUG      } else if (m == 0) {
//DEBUG	flag = " b";
//DEBUG      } else if (m == ndvr[i]) {
//DEBUG	flag = "!b";
//DEBUG      } else {
//DEBUG	flag = " o";
//DEBUG      }
//DEBUG
//DEBUG      iB = iB0 + m;
//DEBUG      DVRX = DVR[i].xpt[m];
//DEBUG      DVRW = DVR[i].wpt[m];
//DEBUG      RADW = wrad[iB];
//DEBUG      DMSK = mask[iB];
//DEBUG      printf( "%8ld %8ld (%8ld) %s %20.10f %20.10f %20.10f %20.10f\n", 
//DEBUG	      i, m, iB, flag.c_str(), DVRX, DVRW, RADW, DMSK );
//DEBUG    }
//DEBUG  }
//DEBUG
//DEBUG  std::cout << "# FE-DVR: global --> local:" << std::endl;
//DEBUG  for (iB = 0; iB <= nrad; iB ++) {
//DEBUG    i = mapb[iB];
//DEBUG    m = iB - mapf[i];
//DEBUG
//DEBUG    if (i == 0 && m == 0) {
//DEBUG      flag = "xL";
//DEBUG    } else if (i == nfe - 1 && m == ndvr[i]) {
//DEBUG      flag = "xR";
//DEBUG    } else if ( m == 0 ) {
//DEBUG      flag = " b";
//DEBUG    } else if ( m == ndvr[ i ] ) {
//DEBUG      flag = "!b";
//DEBUG    } else {
//DEBUG      flag = " o";
//DEBUG    }
//DEBUG
//DEBUG    DVRX = DVR[ i ].xpt[ m ];
//DEBUG    DVRW = DVR[ i ].wpt[ m ];
//DEBUG    RADW = wrad[ iB ];
//DEBUG    DMSK = mask[ iB ];
//DEBUG    printf( "%8ld (%8ld %8ld) %s %20.10f %20.10f %20.10f %20.10f\n", 
//DEBUG	    iB, i, m, flag.c_str(), DVRX, DVRW, RADW, DMSK );
//DEBUG  }

  //old  std::cout << "# FE-DVR: kinetic energy matrix:" << std::endl;
  //old  for (long ife = 0; ife < nfe; ife ++) {
  //old    long iB0 = mapf[ife];
  //old    printf("# element %5ld: n = %5ld, %10.5f <= x <= %10.5f\n", 
  //old	   ife, ndvr[ife], x0[ife], x1[ife]);
  //old    for (long m = 0; m <= ndvr[ife]; m ++) {
  //old      long mB = iB0 + m;
  //old      printf("%10ld", m);
  //old      for (long k = 0; k <= ndvr[ife]; k ++) {
  //old	long kB = iB0 + k;
  //old	printf("%13.5E", kinetic(kB, mB));
  //old      }
  //old      printf("\n");
  //old    }
  //old  }
  //old
  //old  std::cout << "# FE-DVR: nabla operator:" << std::endl;
  //old  for (long ife = 0; ife < nfe; ife ++) {
  //old    long iB0 = mapf[ife];
  //old    printf("# element %5ld: n = %5ld, %10.5f <= x <= %10.5f\n", 
  //old	   ife, ndvr[ife], x0[ife], x1[ife]);
  //old    for (long m = 0; m <= ndvr[ife]; m ++) {
  //old      long mB = iB0 + m;
  //old      printf("%10ld", m);
  //old      for (long k = 0; k <= ndvr[ife]; k ++) {
  //old	long kB = iB0 + k;
  //old	printf("%13.5E", nabla(kB, mB));
  //old      }
  //old      printf("\n");
  //old    }
  //old  }

//DEBUG  std::cout << "# FE-DVR: radk operator:" << std::endl;
//DEBUG  for (long ife = 0; ife < nfe; ife ++) {
//DEBUG    long iB0 = mapf[ife];
//DEBUG    printf("# element %5ld: n = %5ld, %10.5f <= x <= %10.5f\n", 
//DEBUG	   ife, ndvr[ife], x0[ife], x1[ife]);
//DEBUG    for (long m = 0; m <= ndvr[ife]; m ++) {
//DEBUG      long mB = iB0 + m;
//DEBUG      printf("%10ld", m);
//DEBUG      for (long k = 0; k <= ndvr[ife]; k ++) {
//DEBUG	long kB = iB0 + k;
//DEBUG	long kB_B = nmax + kB - mB;
//DEBUG	long kmB = (2 * nmax + 1) * mB + kB_B;
//DEBUG	printf("%13.5E", radk[kmB]);
//DEBUG      }
//DEBUG      printf("\n");
//DEBUG    }
//DEBUG  }
//DEBUG
//DEBUG  std::cout << "# FE-DVR: radp operator:" << std::endl;
//DEBUG  for (long ife = 0; ife < nfe; ife ++) {
//DEBUG    long iB0 = mapf[ife];
//DEBUG    printf("# element %5ld: n = %5ld, %10.5f <= x <= %10.5f\n", 
//DEBUG	   ife, ndvr[ife], x0[ife], x1[ife]);
//DEBUG    for (long m = 0; m <= ndvr[ife]; m ++) {
//DEBUG      long mB = iB0 + m;
//DEBUG      printf("%10ld", m);
//DEBUG      for (long k = 0; k <= ndvr[ife]; k ++) {
//DEBUG	long kB = iB0 + k;
//DEBUG	long kB_B = nmax + kB - mB;
//DEBUG	long kmB = (2 * nmax + 1) * mB + kB_B;
//DEBUG	printf("%13.5E", radp[kmB]);
//DEBUG      }
//DEBUG      printf("\n");
//DEBUG    }
//DEBUG  }
}
////////////////////////////////////////////////////////////////////////
void clfedvr::gen(const clmpi& MPIP, const clio& IO)
{
  //2016.03.08
  clear();

  IO.read_info("nfe", nfe);
  x0.resize(nfe);
  x1.resize(nfe);
  dfe.resize(nfe);
  ndvr.resize(nfe);
  DVR.resize(nfe);

  ///////////
  read_x(IO);
  ///////////

  mapf.resize(nfe);
  mapb.resize(nrad + 1);
  mask.resize(nrad + 1);
  xrad.resize(nrad + 1);
  wrad.resize(nrad + 1);
  cxrad.resize(nrad + 1);
  cwrad.resize(nrad + 1);
  bra_wrad.resize(nrad + 1);

  //old  kinetic.resize(nrad + 1, nrad + 1, nmax, nmax);
  //old  nabla.resize(nrad + 1, nrad + 1, nmax, nmax);
  //old  nabla_new.resize((2 * nmax + 1) * (nrad + 1));
  radk.resize((2 * nmax + 1) * (nrad + 1));
  radk0.resize((2 * nmax + 1) * (nrad + 1));
  radp.resize((2 * nmax + 1) * (nrad + 1));

  // debug
  printf("clfedvr nrad  = %10ld\n", nrad);   // 0 <= irad <= nrad
  printf("clfedvr nmax  = %10ld\n", nmax); 
  // debug

// Orimo_ECS
  gen_abc(IO);
  rdr.resize(nrad + 1);
  wdw.resize(nrad + 1);
// Orimo_ECS

  gen_map();
  gen_grid();
  gen_mask();
  gen_radk();
  gen_radp();

// Orimo_ECS
  set_potential(IO);
// Orimo_ECS

  if (IO.iprint > 0) print();
}
////////////////////////////////////////////////////////////////////////
void clfedvr::read_x(const clio& IO)
{
  std::string key = "grid:";
  std::ifstream ifs(IO.inp.c_str());
  ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;

  nmax = 0;
  while ( getline(ifs, line) && line.find(key,0) == std::string::npos ) {}
  if(! ifs.eof()) {

    nrad = 0;

    long ife = 0;
    long ndvr_tmp;
    long nfe_tmp;
    double x0_tmp;
    double x1_tmp;
    double dfe_tmp;
    while (ife < nfe) {
      getline(ifs, line);
      ioss.str("");
      ioss.clear( std::stringstream::goodbit );
      ioss << line.c_str();
      ioss >> x0_tmp
	   >> x1_tmp
	   >> ndvr_tmp
	   >> nfe_tmp;
      dfe_tmp = (x1_tmp - x0_tmp) / nfe_tmp;
      for (long ife_tmp = 0; ife_tmp < nfe_tmp; ife_tmp ++) {
	x0[ife + ife_tmp] = x0_tmp + dfe_tmp * ife_tmp;
	x1[ife + ife_tmp] = x0[ife + ife_tmp] + dfe_tmp;
	ndvr[ife + ife_tmp] = ndvr_tmp;
	dfe [ife + ife_tmp] = dfe_tmp;
      }
      nmax = std::max(nmax, ndvr_tmp);
      nrad += ndvr_tmp * nfe_tmp;
      ife += nfe_tmp;
    }
    //BUG  nrad -= 1;

    std::cout << "# grid:" << std::endl;
    std::cout << "# " << "nfe  = " << nfe << std::endl;
    for ( long ife = 0; ife < nfe; ife ++ ) {
      printf( "#%10ld %20.10f%20.10f%20.10f%20ld\n",
	      ife, x0[ ife ], x1[ ife ], dfe[ ife ], ndvr[ ife ] );
    }

  } else {

    std::cout<< "No grid information." << std::endl;
    abort();
  }

  rll = x0[ 0 ];
  rul = x1[ nfe - 1 ];
}
////////////////////////////////////////////////////////////////////////
void clfedvr::gen_map()
//
// Conversion map between local and global bases.
//
{
  long NB = 0;

  for (long ife = 0; ife < nfe; ife ++) {
    mapf[ife] = NB;
    //    for (long m = m0[ife]; m < ndvr[ife]; m ++) {
    for (long m = 0; m < ndvr[ife]; m ++) {
      mapb[NB] = ife;
      NB ++;
    }
  }
  mapb[nrad] = nfe - 1;

  // initialize nrad_fc by nrad
  nradfc = nrad - 1;

  if ( NB != nrad ) {
  //BUG  if (NB != nrad - 1) {
    std::cout << " nrad = " << nrad << " NB = " << NB << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clfedvr::gen_x()
//
// Generate left and right boundaries for each element
//
{
  double xp = rll;

  for ( long i = 0; i < nfe; i ++ ) {
    x0[ i ] = xp;
    x1[ i ] = x0[ i ] + dfe[ i ];
    xp = x1[ i ];
  }

  // debug
  std::cout << "# x0 and x1 for each element:" << std::endl;
  for ( long i = 0; i < nfe; i ++ ) {
    printf( "%10ld %20.10f %20.10f\n", i, x0[ i ], x1[ i ] );      
  }
  // debug
}
////////////////////////////////////////////////////////////////////////
void clfedvr::gen_mask()
//
// mask function of cos type.
//
{
  long i, m;
  const double pih = PI * HALF;
  const double oo4 = ONE / FOUR;
  const double oo8 = oo4 * HALF;;
  double absx, tmp, order;

  mask[0] = ZERO;
  mask[nrad] = ZERO;
  if (mask_cos_order > 0) {
    order = ONE / mask_cos_order;
  } else {
    order = ONE * get_abs(mask_cos_order);
  }
  //DEBUG
  //  printf("# clfedvr::gen_mask. order = %20.10f\n", order);
  //DEBUG

  for (long iB = 1; iB < nrad; iB ++) {
    i = mapb[iB];
    m = iB - mapf[i];
    absx = sqrt(DVR[i].xpt[m] * DVR[i].xpt[m]);
    if (absx <= rmask) {
      mask[iB] = ONE;
    } else {
      tmp = pih * (absx - rmask) / (rul - rmask);
      tmp = cos(tmp);
      mask[iB] = pow(tmp, order);
    }
    //DEBUG
    //    printf( " %5ld absx = %20.10f, mask = %20.10f\n", iB, absx, mask[ iB ] );
    //DEBUG
  }
}
////////////////////////////////////////////////////////////////////////
void clfedvr::gen_grid()
// Generate DVR nodes and weights within each element. 
// Local DVR of i'th element consists of (n + 1) quadrature nodes, 
// in which the 0th and nth points coincide the left and right edges 
// of the element.
// They also support local (n + 1) Lagrange interpolation polynomials.
// The left and right boundary points, as well as redundant bridges, 
// are included here.
{
  // Orimo_ECS
  // for (long i = 0; i < nfe; i ++) {
  //   DVR[i].gen(ndvr[i], 0, ndvr[i], x0[i], x1[i]);
  // }
  if (!inf_range) {
    for (long i = 0; i < nfe; i ++) {
      DVR[i].gen(ndvr[i], 0, ndvr[i], x0[i], x1[i]);
    }
  } else if(inf_range){
    for (long i = 0; i < nfe - 1; i ++) {
      DVR[i].gen(ndvr[i], 0, ndvr[i], x0[i], x1[i]);
    }
    int lfe = nfe -1;
    DVR[lfe].gen_radau(ndvr[lfe], 0, ndvr[lfe], x0[lfe], x1[lfe], 2.0*exp_factor);
  }
  // Orimo_ECS

  // Global arrays
  long i, m;
  //BUG  for (long iB = 0; iB < nrad; iB ++) {
  for (long iB = 0; iB <= nrad; iB ++) {
    i = mapb[iB];
    m = iB - mapf[i];
    xrad[iB] = DVR[i].xpt[m];
    wrad[iB] = DVR[i].wpt[m];
    if(xrad[iB] >= recs){
      cxrad[iB] = recs + exp(IUNIT * theta) * (xrad[iB] -recs);
      cwrad[iB] = DVR[i].wpt[m] * exp(IUNIT * theta);
    } else{
      cxrad[iB] = xrad[iB];
      cwrad[iB] = DVR[i].wpt[m];
    }
  }


  // Adjust weights for bridge functions, namely,
  // w^i_m = w^i_m, if (i|m) is not a bridge ( 1 <= m <= n-1 ), 
  // w^i_0 = w^(i-1)_n + w^i_0 = w^(i-1)_n
  //         ^^^^^^^^^^^^^^^^^
  for (long i = 1; i < nfe; i ++) {
    long iB = mapf[i];
    wrad[iB] += DVR[i - 1].wpt[ndvr[i - 1]];
    if(xrad[iB] > recs) 
      cwrad[iB] += DVR[i - 1].wpt[ndvr[i - 1]] * exp(IUNIT * theta);
    else 
      cwrad[iB] += DVR[i - 1].wpt[ndvr[i - 1]];
  }

// Orimo_ECS
//old  for(long irad = 0; irad < nrad + 1; irad++){
//old    bra_wrad[irad] = wrad[irad] / real(sqrt(std::conj(cwrad[irad]) * cwrad[irad]));
//old  }
//old
//old  if (ecs_flag == 1){
//old    irad_ecs = -1;
//old    for(int irad = 0; irad < nrad + 1; irad++){
//old      if(xrad[irad] == recs) irad_ecs = irad;
//old    }
//old    if(irad_ecs == -1){
//old      std::cout << "bad recs condition" << std::endl;
//old      abort();
//old    } 
//old  }
  rdr[0] = 1;
  wdw[0] = 1;
  bra_wrad[0] = 1;
  for(long irad = 1; irad < nrad + 1; irad++){
    bra_wrad[irad] = wrad[irad] / real(sqrt(std::conj(cwrad[irad]) * cwrad[irad]));
    rdr[irad] = cxrad[irad] / std::conj(cxrad[irad]);
    wdw[irad] = std::sqrt(cwrad[irad] / std::conj(cwrad[irad]));
  }

  if (ecs_flag == 1){
    double THV = pow(10., -10.);
    irad_ecs = -1;
    int test_boundary = -1;
    for(int irad = 0; irad < nrad + 1; irad++){
      if(fabs(xrad[irad] -recs) <= THV){
	std::cout << "error = " << abs(xrad[irad] -recs) << std::endl;
	irad_ecs = irad;
	irad_sw = irad;
	break;
      }
    }
    std::cout<< "THV = " << THV << std::endl;
    if(irad_ecs != -1){
      for(int ife = 0; ife < nfe; ife++){
	if(fabs(xrad[irad_ecs] - x0[ife]) <= THV){
	  test_boundary = 0;
	  break;
	}
      }
    }
    if(irad_ecs == -1 || test_boundary == -1){
      std::cout << "bad recs or rab condition" << std::endl;
      abort();
    } else {
      std::cout<< "irad_ecs = " << irad_ecs << std::endl;
      std::cout<< "irad_sw = " << irad_sw << std::endl;
    }
  }
// Orimo_ECS
}
////////////////////////////////////////////////////////////////////////
//oldvoid clfedvr::gen_kinetic()
//old{
//old  kinetic.clear();
//old
//old  double tmp;
//old  for (long i = 0; i < nfe; i ++) {
//old    long iB0 = mapf[i];
//old    for (long m = 0; m <= ndvr[i]; m ++) {
//old      long mB = iB0 + m;
//old      for (long k = 0; k <= ndvr[i]; k ++) {
//old	long kB = iB0 + k;
//old	long kmi = (ndvr[0] + 1) * (ndvr[0] + 1) * i
//old                                 + (ndvr[0] + 1) * m
//old                                                 + k;
//old	if (i == 0 && (m == 0 || k == 0)) {
//old	  kinetic(kB, mB) = ZERO;
//old	  kloc[kmi] = ZERO;
//old	} else if (i == nfe - 1 && (m == ndvr[i] || k == ndvr[i])) {
//old	  kinetic(kB, mB) = ZERO;
//old	  kloc[kmi] = ZERO;
//old	} else {
//old	  tmp = ZERO;
//old	  for (long l = 0; l <= ndvr[i]; l ++) {
//old	    tmp += DVR[i].wpt[ l ] * DVR[i].dshape(l, m) 
//old                                   * DVR[i].dshape(l, k);
//old	  }
//old	  tmp = HALF * tmp / sqrt(wrad[mB] * wrad[kB]);
//old	  //	  kin[i].darray[ k ][ m ] = tmp;
//old	  //	  kinetic(kB, mB) = tmp;
//old	  kinetic(kB, mB) += tmp;
//old	  kloc[kmi] = tmp;
//old	}
//old      }
//old    }
//old  }
//old}
//old////////////////////////////////////////////////////////////////////////
//oldvoid clfedvr::gen_nabla()
//old{
//old  nabla.clear();
//old
//old  double tmp;
//old  for (long i = 0; i < nfe; i ++) {
//old    long iB0 = mapf[i];
//old    for (long m = 0; m <= ndvr[i]; m ++) {
//old      long mB = iB0 + m;
//old      for (long k = 0; k <= ndvr[i]; k ++) {
//old	long kB = iB0 + k;
//old	if (i == 0 && (m == 0 || k == 0)) {
//old	  nabla(kB, mB) = ZERO;
//old	} else if (i == nfe - 1 && (m == ndvr[i] || k == ndvr[i])) {
//old	  nabla(kB, mB) = ZERO;
//old	} else {
//old	  tmp = DVR[i].wpt[m] * DVR[i].dshape(m, k);
//old	  nabla(kB, mB) += tmp / sqrt(wrad[mB] * wrad[kB]);
//old	}
//old      }
//old    }
//old  }
//old}
//old////////////////////////////////////////////////////////////////////////
//oldvoid clfedvr::gen_nabla_new()
//old{
//old  nabla_new.clear();
//old  //  fedvr_gen_nabla_new();
//old  double tmp;
//old  long kB_B, kmB;
//old  for (long i = 0; i < nfe; i ++) {
//old    long iB0 = mapf[i];
//old    for (long m = 0; m <= ndvr[i]; m ++) {
//old      long mB = iB0 + m;
//old      for (long k = 0; k <= ndvr[i]; k ++) {
//old	long kB = iB0 + k;
//old	kB_B = nmax + kB - mB;
//old//bug	kmB = (2 * nmax + 1) * kB + kB_B;
//old	kmB = (2 * nmax + 1) * mB + kB_B;
//old
//old	if (i == 0 && (m == 0 || k == 0)) {
//old	  nabla_new[kmB] = ZERO;
//old	} else if (i == nfe - 1 && (m == ndvr[i] || k == ndvr[i])) {
//old	  nabla_new[kmB] = ZERO;
//old	} else {
//old	  tmp = DVR[i].wpt[m] * DVR[i].dshape(m, k);
//old	  nabla_new[kmB] += tmp / sqrt(wrad[mB] * wrad[kB]);
//old	}
//old      }
//old    }
//old  }
//old}
////////////////////////////////////////////////////////////////////////
long clfedvr::get_irad(double rad) const
{
  long irad;
  double small = 1.E-10;

  irad = 1;
  while (rad + small > xrad[irad]) {
    irad ++;
    if (irad == nrad - 1) break;
  }

  //  printf("get_irad: irad = %ld", irad);
  return irad;
}
////////////////////////////////////////////////////////////////////////
long clfedvr::get_irad(long ife, long m) const
{
  return mapf[ife] + m;
}
//////////////////////////////////////////////////////////////////////////
long clfedvr::get_ndvr(long ife) const
{
  return ndvr[ife];
}
//////////////////////////////////////////////////////////////////////////
long clfedvr::get_ife(double rval) const
{
  long ife = -1;
  double small = 1.E-10;
  for (long i = 0; i < nfe; i ++) {
    if (rval + small < x1[i]) {
      ife = i;
      break;
    }
  }
  return ife;
}
////////////////////////////////////////////////////////////////////////
double clfedvr::get_x0(long ife) const
{
  return xrad[mapf[ife]];
}
////////////////////////////////////////////////////////////////////////
double clfedvr::get_x1(long ife) const
{
  return xrad[mapf[ife]+ndvr[ife]];
}
//////////////////////////////////////////////////////////////////////////
double clfedvr::get_val(long irad, double rval) const
{
  long ife = mapb[irad];
  long m = irad - mapf[ife];
//  printf("%10ld%20.10f%10ld%10ld%20.10f%20.10f%20.10f\n", irad, rval, ife, m, 
//	 DVR[ife].get_val(m, rval), sqrt(wrad[irad]),
//	 DVR[ife].get_val(m, rval) / sqrt(wrad[irad]));
//return DVR[ife].get_val(m, rval) / sqrt(wrad[irad]);
  return DVR[ife].get_val(m, rval);
}
//////////////////////////////////////////////////////////////////////////
double clfedvr::get_val(long ife, long idvr, double rval) const
{
  long irad = mapf[ife] + idvr;
//return DVR[ife].get_val(idvr, rval) / sqrt(wrad[irad]);
  return DVR[ife].get_val(idvr, rval);
}
//////////////////////////////////////////////////////////////////////////
double clfedvr::get_val0(long ife, double rval) const
{
  long irad = mapf[ife + 1];
//return DVR[ife].get_val(ndvr[ife], rval) / sqrt(wrad[irad]);
  return DVR[ife].get_val(ndvr[ife], rval);
}
////////////////////////////////////////////////////////////////////////
void clfedvr::gen_radk()
{
  //bug  radk.clear();
  radk.assign(radk.size(), ZERO);
  radk0.assign(radk0.size(), ZERO);

  double tmp, tmp0;
  long lm, lk, iB0, mB, kB, kB_B, kmB;

  for (long i = 0; i < nfe; i ++) {
    iB0 = mapf[i];
    for (long m = 0; m <= ndvr[i]; m ++) {
      mB = iB0 + m;
      for (long k = 0; k <= ndvr[i]; k ++) {
	kB = iB0 + k;
	kB_B = nmax + kB - mB;
	kmB = (2 * nmax + 1) * mB + kB_B;

	if (i == 0 && (m == 0 || k == 0)) {
	  radk[kmB] = ZERO;
	  radk0[kmB] = ZERO;
	} else if (i == nfe - 1 && (m == ndvr[i] || k == ndvr[i])) {
	  radk[kmB] = ZERO;
	  radk0[kmB] = ZERO;
	} else {
	  tmp = ZERO;
	  tmp0 = ZERO;
	  if (i == 0) {
	    for (long l = 1; l <= ndvr[i]; l ++) {
	      lm = (ndvr[i] + 1) * m + l;
	      lk = (ndvr[i] + 1) * k + l;
	      tmp0 += DVR[i].wpt[ l ] * DVR[i].dshape[lm] 
                                      * DVR[i].dshape[lk];
	    }
	  } else {
	    for (long l = 0; l <= ndvr[i]; l ++) {
	      lm = (ndvr[i] + 1) * m + l;
	      lk = (ndvr[i] + 1) * k + l;
	      tmp0 += DVR[i].wpt[ l ] * DVR[i].dshape[lm] 
                                      * DVR[i].dshape[lk];
	    }
	  }
	  for (long l = 0; l <= ndvr[i]; l ++) {
	    lm = (ndvr[i] + 1) * m + l;
	    lk = (ndvr[i] + 1) * k + l;
	    tmp += DVR[i].wpt[ l ] * DVR[i].dshape[lm] 
                                   * DVR[i].dshape[lk];
	  }

	  // for radkI_ecs, when r <= R0 = recs 
	  if(xrad[mB] == recs && xrad[kB] == recs && xrad[iB0] < recs){
	    radkI_ecs = 0; // initialization
	    for (long l = 0; l <= ndvr[i]; l ++) {
	      lm = (ndvr[i] + 1) * m + l;
	      lk = (ndvr[i] + 1) * k + l;
	      radkI_ecs += DVR[i].wpt[ l ] * DVR[i].dshape[lm] 
	                                   * DVR[i].dshape[lk]; 
	    }
	  }
	  // for radkI_ecs, when r >= Ro = recs
	  if(xrad[mB] == recs && xrad[kB] == recs && xrad[iB0] >= recs){
	    for (long l = 0; l <= ndvr[i]; l ++) {
	      lm = (ndvr[i] + 1) * m + l;
	      lk = (ndvr[i] + 1) * k + l;
	      radkI_ecs += DVR[i].wpt[ l ] * exp(IUNIT * theta) 
		* exp( (-2.0) * IUNIT * theta ) * DVR[i].dshape[lm] * DVR[i].dshape[lk]; 
	    }
	     radkI_ecs = HALF * radkI_ecs / cwrad[mB]; // normalization like below \/:: cwrad[mB] = sqrt(cwrad[mB] * cwrad[kB])  
	  }

	  if (clcontrol::fedvr_normalized) {
	    radk[kmB] += HALF * tmp / sqrt(wrad[mB] * wrad[kB]);
	    radk0[kmB] += HALF * tmp0 / sqrt(wrad[mB] * wrad[kB]);
	  } else {
	    radk[kmB] += HALF * tmp;
	    radk0[kmB] += HALF * tmp0;
	  }
	}
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////
void clfedvr::gen_radp()
{
  //bug  radp.clear();
  radp.assign(radp.size(), ZERO);

  double tmp;
  long iB0, mB, kB, kB_B, kmB, mk;

  for (long i = 0; i < nfe; i ++) {
    iB0 = mapf[i];
    for (long m = 0; m <= ndvr[i]; m ++) {
      mB = iB0 + m;
      for (long k = 0; k <= ndvr[i]; k ++) {
	kB = iB0 + k;
	kB_B = nmax + kB - mB;
	kmB = (2 * nmax + 1) * mB + kB_B;
	mk = (ndvr[i] + 1) * k + m;

	if (i == 0 && (m == 0 || k == 0)) {
	  radp[kmB] = ZERO;
	} else if (i == nfe - 1 && (m == ndvr[i] || k == ndvr[i])) {
	  radp[kmB] = ZERO;
	} else {
	  tmp = DVR[i].wpt[m] * DVR[i].dshape[mk];
	  if (clcontrol::fedvr_normalized) {
	    radp[kmB] += tmp / sqrt(wrad[mB] * wrad[kB]);
	  } else {
	    radp[kmB] += tmp;
	  }
	}
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////
// Orimo_ECS
void clfedvr::gen_abc(const clio& IO)
{
  IO.read_info("abc_type", "mask", abc_type);

  double recs_def = rul * 0.8;
  double rmask_def = rul * 0.8;

  IO.read_info("recs", recs_def, recs); 
  double theta_def = 30.0; // (degree) = 1/6*PI (rad) 
  IO.read_info("theta", theta_def, theta);
  theta = theta / 180.0 * PI;

  IO.read_info("rmask", rmask_def, rmask);
  IO.read_info("mask_type", "cos", mask_type);
  IO.read_info("mask_cos_order", LFOUR, mask_cos_order);

  double rab_def = -1;
  IO.read_info("rab", rab_def, rab);
  if(rab != -1){
    recs = rab;
    rmask = rab;
  }

  if(abc_type == "ECS" || abc_type == "ecs"){
    abc_type = "ecs";
    ecs_flag = 1;
    IO.read_info("retreive_flag", 0, retreive_flag);
    rmask = rul + 1; // reason for [+1] is to avoid 0 division
  } else if(abc_type == "mask"){
    ecs_flag = 0;
    recs = rul;
    theta = 0;
  } else if(abc_type == "ecs_mask" ||abc_type == "ECS_mask" || abc_type == "mask_ecs" || abc_type == "mask_ECS"){
      abc_type = "ecs_mask";
      ecs_flag = 1;
  }

  if(!(abc_type == "mask" || abc_type == "ecs" ||abc_type == "ecs_mask")){
    std::cout << "improper Absornig Boudary Condition type!" << std::endl;
    abort();
  }  
  
  
  if(ecs_flag == 1){
    IO.read_info("type_mkint1_sph", 4, type_mkint1_sph);
    IO.read_info("type_mkint2_sph", 4, type_mkint2_sph);
    IO.read_info("type_mkv2mf", 5, type_mkv2mf);
    IO.read_info("type_mkxmat_aa", 4, type_mkxmat_aa);
  }else if(ecs_flag == 0){
    IO.read_info("type_mkint1_sph", -1, type_mkint1_sph);
    IO.read_info("type_mkint2_sph", -1, type_mkint2_sph);
    IO.read_info("type_mkv2mf", -1, type_mkv2mf);
    IO.read_info("type_mkxmat_aa", -1, type_mkxmat_aa);
    IO.read_info("irad_ecs", nrad, irad_ecs);
  }
  
  double def_exp_factor = HALF;
  IO.read_info("inf_range", false, inf_range);
  IO.read_info("exp_factor", def_exp_factor, exp_factor);
  irad_inf = 0;
  for(int ife = 0; ife < nfe-1; ife++)  irad_inf += ndvr[ife];
  if(inf_range) std::cout<< "# irad_inf = " << irad_inf <<std::endl;
  else irad_inf += ndvr[nfe-1];
  
  std::string sw_info;
  if(ecs_flag == 0){
    IO.read_info("switchoff", "no", sw_info);
    IO.read_info("irad_sw", nrad-1, irad_sw);
  }else if(ecs_flag == 1){
    IO.read_info("switchoff", "mkrho2", sw_info);

  }
  if(sw_info == "no") switchoff = 0;
  else if(sw_info == "mkrho2") switchoff = 1;
  else if(sw_info == "mkv2mf") switchoff = 2;

//OLD   IO.read_info("abc_type", "mask", abc_type);
//OLD   if(!(abc_type == "mask" || abc_type == "ECS" || abc_type == "ecs")){
//OLD     std::cout << "improper Absornig Boudary Condition type!" << std::endl;
//OLD     abort();
//OLD   }
//OLD
//OLD   double recs_def = rul * 0.8;
//OLD   double rmask_def = rul * 0.8;
//OLD
//OLD   IO.read_info("recs", recs_def, recs); 
//OLD   double theta_def = 30.0; // (degree) = 1/6*PI (rad) 
//OLD   IO.read_info("theta", theta_def, theta);
//OLD   theta = theta / 180.0 * PI;
//OLD
//OLD   IO.read_info("rmask", rmask_def, rmask);
//OLD   IO.read_info("mask_type", "cos", mask_type);
//OLD   IO.read_info("mask_cos_order", LFOUR, mask_cos_order);
//OLD
//OLD   double rab_def = -1;
//OLD   IO.read_info("rab", rab_def, rab);
//OLD   if(rab != -1){
//OLD     recs = rab;
//OLD     rmask = rab;
//OLD   }
//OLD
//OLD   if(abc_type == "ECS" || abc_type == "ecs"){
//OLD     abc_type = "ecs";
//OLD     ecs_flag = 1;
//OLD     rmask = rul + 1; // reason of [+1] is to avoid 0 division
//OLD   } else if(abc_type == "mask"){
//OLD     ecs_flag = 0;
//OLD     recs = rul;
//OLD     theta = 0;
//OLD   }

  return;
}
// Orimo_ECS
////////////////////////////////////////////////////////////////////////
// Orimo_ECS
void clfedvr::set_potential(const clio& IO)
{
  trunc_irad = nrad;
  std::string pot_str;
  IO.read_info("pot_type", "atom", pot_str);
  if(pot_str == "atom") {
    pot_type = 1;
  } else{
    std::cout << "pot_type option is not implemented in this version" << std::endl;
    abort();
  }
//   if(pot_str == "atom") {
//     pot_type = 1;
//   }
//   else if(pot_str == "harmonic"){
//      pot_type = 2;
//   }
//   else if(pot_str == "yukawa"){
//     pot_type = 3;
//   }
//   else if(pot_str == "all_truncation"){
//     pot_type = 4;
//     double trunc_rrad = -1.0;
//     trunc_irad = -1;
//     IO.read_info("truncation_radius", trunc_rrad);
//     trunc_irad = get_irad(trunc_rrad) - 1;
//     if(trunc_irad < 0){
//       std::cout << "bad truncation radius." << std::endl;
//       abort();
//     }
//     std::cout << "# trunc_irad: " << trunc_irad << std::endl;
//     std::cout << "# trunc_rrad: " << trunc_rrad << std::endl;
//   }
//   else {
//     std::cout << "bad pot_type" << std::endl;
//     abort();
//   }
  
  return;
}
// Orimo_ECS
////////////////////////////////////////////////////////////////////////
