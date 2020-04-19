//2015/10/22 Yuki Orimo Changed
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
void guess(const clmpi& Proc, const clio& IO, const clbas& Bas, clwfn& Wfn)
{
  clock_t time0 = clock();

  std::string guess_type;
  IO.read_info("guess_type", guess_type);

  // initialization
  Wfn.clear(Bas);
  //old Wfn.wfn[Wfn.size1] = RUNIT;
  Bas.ORMAS.cic0(&Wfn.wfn[Wfn.size1]);
  //  ormas_cic0_(&Wfn.wfn[Wfn.size1]);
  //DEBUG
  //  std::string fname = IO.name + ".ci0";;
  //  Wfn.print_cic(fname, Bas);
  //  abort();
  //DEBUG

  if (guess_type.compare("card") == 0) {
    guess_card(Proc, IO, Bas, Wfn);
    Wfn.orth(Bas);
  } else if (guess_type.compare("read") == 0) {
    Wfn.read(Proc, IO, Bas);
  } else if (guess_type.compare("read_card") == 0) {
    guess_card(Proc, IO, Bas, Wfn);
    Wfn.read_orb(Proc, IO, Bas);
    Wfn.read_cic(Proc, IO, Bas);
    Wfn.orth(Bas);
  } else if (guess_type.compare("read_orb") == 0) {
    Wfn.read_orb(Proc, IO, Bas);
    Wfn.orth(Bas);
  } else if (guess_type.compare("read_proj") == 0) {
    std::string inp_proj;
    inp_proj = IO.name; 
    inp_proj += "_proj.inp";
    clio IO_Proj(inp_proj);
    clbas Bas_Proj(Proc, IO_Proj);
    clwfn Wfn_Proj(Proc, IO_Proj, Bas_Proj);
    Wfn_Proj.read(Proc, IO_Proj, Bas_Proj);

    Bas_Proj.proj(Bas, Wfn_Proj.wfn, Wfn.wfn);
    zcopy_omp_(&Wfn.size2, &Wfn_Proj.wfn[Wfn_Proj.size1], &Wfn.wfn[Wfn.size1]);
    Wfn.orth(Bas);
//nyi  } else if (guess_type.compare("g09") == 0) {
//nyi
//nyi    guess_g09(Proc, IO, Bas, Wfn);
//nyi    Wfn.orth(Bas);
  } else {
    std::cout << "Bad guess_type." << std::endl;
    abort();
  }

  // smaller box for frozen-core orbitals
  Wfn.get_nradfc(Proc, IO, Bas);

  // save
  Wfn.write(Proc, IO, Bas);

  printf("lm-components of guess orbitals:\n");
  Wfn.print_wang(Bas);

  bool print_orb;
  long print_orb_nstep;
  std::string fname = IO.name + "_guess.orbital";
  IO.read_info("print_orb", false, print_orb);
  if (print_orb) {
    IO.read_info("print_orb_nstep", LONE*10, print_orb_nstep);
    Wfn.print_orb(print_orb_nstep, fname, Bas);
  }

  clock_t time1 = clock();
  std::cout << "# guess:" << (double)(time1 - time0) / CLOCKS_PER_SEC << std::endl;
}
////////////////////////////////////////////////////////////////////////
double guess_rfun_hlike(long n, long l, const long& znuc, const double& xrad, const double& wrad)
{
  double rval = xrad * znuc;
  double polr, expr;

  ///// gaussians /////
  //  expr = exp(-rval * rval); //1s
  //  polr = rval;              //1s
  ///// scaled slater /////
  //  expr = exp(-rval * 1.2); //1s
  //  polr = rval;             //1s

  ///// exact hydrogenic orbitals /////
  if (n == 1) { // 1s, or [1p], [1d], ...
    expr = exp(-rval);
    polr = rval;      
  } else if (n == 2 && l == 0) { // 2s
    expr = exp(-rval / 2.0);
    polr = rval * (2.0 - rval);
  } else if (n == 2 && l >= 1) { // 2p, or [2d], [2f], ...
    expr = exp(-rval / 2.0);
    polr = rval * rval;
  } else if (n == 3 && l == 0) { // 3s
    expr = exp(-rval / 3.0);
    polr = rval * (27.0 + rval * (-18.0 + 2.0 * rval));
  } else if (n == 3 && l == 1) { // 3p
    expr = exp(-rval / 3.0);
    polr = rval * rval * (6.0 - rval);
  } else if (n == 3 && l >= 2) { // 3d, or [3f], [3g], ...
    expr = exp(-rval / 3.0);
    polr = rval * rval * rval;
  } else if (n == 4 && l == 0) { // 4s
    expr = exp(-rval / 4.0);
    polr = rval * (192.0 + rval * (-144.0 + rval * (24.0 - rval)));
  } else if (n == 4 && l == 1) { // 4p
    expr = exp(-rval / 4.0);
    polr = rval * rval * (80.0 + rval * (- 20.0 + rval));
  } else if (n == 4 && l == 2) { // 4d
    expr = exp(-rval / 4.0);
    polr = rval * rval * rval * (12.0 - rval);
  } else if (n == 4 && l >= 3) { // 4f, or [4g], [4h], ...
    expr = exp(-rval / 4.0);
    polr = rval * rval * rval * rval;
  } else if (n >= 5) {           // 5g, 6h, 7i, etc...
    double dn = static_cast<double>(n);
    expr = exp(-rval / dn);
    polr = pow(rval, dn);
  } else {
    printf("guess_rfun: n > MAXN.\n");
    abort();
  }

  double rfun;
  if (clcontrol::fedvr_normalized) {
    rfun = polr * expr * sqrt(wrad);
  } else {
    rfun = polr * expr;
  }
  return rfun;
}
////////////////////////////////////////////////////////////////////////
dcomplex guess_rfun_hlike_ecs(long n, long l, const long& znuc, const double& xrad,
			      const dcomplex& cwrad, const clbas& Bas)
{
  dcomplex rval;
  if(xrad < Bas.GRad.recs)
    rval = xrad * znuc;
  if(xrad >= Bas.GRad.recs)
    rval = ( exp(IUNIT * Bas.GRad.theta) * (xrad - Bas.GRad.recs) + Bas.GRad.recs ) * (double)znuc;

  dcomplex polr, expr;

  ///// gaussians /////
  //  expr = exp(-rval * rval); //1s
  //  polr = rval;              //1s
  ///// scaled slater /////
  //  expr = exp(-rval * 1.2); //1s
  //  polr = rval;             //1s

  ///// exact hydrogenic orbitals /////
  if (n == 1) { // 1s, or [1p], [1d], ...
    expr = exp(-rval);
    polr = rval;      
  } else if (n == 2 && l == 0) { // 2s
    expr = exp(-rval / 2.0);
    polr = rval * (2.0 - rval);
  } else if (n == 2 && l >= 1) { // 2p, or [2d], [2f], ...
    expr = exp(-rval / 2.0);
    polr = rval * rval;
  } else if (n == 3 && l == 0) { // 3s
    expr = exp(-rval / 3.0);
    polr = rval * (27.0 + rval * (-18.0 + 2.0 * rval));
  } else if (n == 3 && l == 1) { // 3p
    expr = exp(-rval / 3.0);
    polr = rval * rval * (6.0 - rval);
  } else if (n == 3 && l >= 2) { // 3d, or [3f], [3g], ...
    expr = exp(-rval / 3.0);
    polr = rval * rval * rval;
  } else if (n == 4 && l == 0) { // 4s
    expr = exp(-rval / 4.0);
    polr = rval * (192.0 + rval * (-144.0 + rval * (24.0 - rval)));
  } else if (n == 4 && l == 1) { // 4p
    expr = exp(-rval / 4.0);
    polr = rval * rval * (80.0 + rval * (- 20.0 + rval));
  } else if (n == 4 && l == 2) { // 4d
    expr = exp(-rval / 4.0);
    polr = rval * rval * rval * (12.0 - rval);
  } else if (n == 4 && l >= 3) { // 4f, or [4g], [4h], ...
    expr = exp(-rval / 4.0);
    polr = rval * rval * rval * rval;
  } else if (n >= 5) {           // 5g, 6h, 7i, etc...
    double dn = static_cast<double>(n);
    expr = exp(-rval / dn);
    polr = pow(rval, dn);
  } else {
    printf("guess_rfun: n > MAXN.\n");
    abort();
  }

  dcomplex rfun;
  if (clcontrol::fedvr_normalized) {
    rfun = polr * expr * sqrt(cwrad);
  } else {
    rfun = polr * expr;
  }
  return rfun;
}
////////////////////////////////////////////////////////////////////////
double guess_rfun_slater(double neff, double xi, double xrad, double wrad)
{
  double rval = xi * xrad;
  double polr, expr;

  expr = exp(-rval);
  polr = pow(xrad, neff);

  double rfun;
  if (clcontrol::fedvr_normalized) {
    rfun = polr * expr * sqrt(wrad);
  } else {
    rfun = polr * expr;
  }
  return rfun;
}
////////////////////////////////////////////////////////////////////////
dcomplex guess_rfun_slater_ecs(double neff, double xi, double xrad, dcomplex cwrad,const clbas& Bas)
{

  dcomplex rval = xi * xrad;
  if(xrad > Bas.GRad.recs)
    rval = xi * ( exp( IUNIT * Bas.GRad.theta) * (xrad - Bas.GRad.recs) + Bas.GRad.recs);

  dcomplex polr, expr;

  expr = exp(-rval);
  polr = pow(xrad, neff);

  dcomplex rfun;
  if (clcontrol::fedvr_normalized) {
    rfun = polr * expr * sqrt(cwrad);
  } else {
    rfun = polr * expr;
  }
  return rfun;;
}
////////////////////////////////////////////////////////////////////////
void guess_card(const clmpi& Proc, const clio& IO, const clbas& Bas, clwfn& Wfn)
{
  double neff, xi, faccs;
  long n0, l0, m0, ind;
  long nel = Bas.ORMAS.neltot[2];

  std::string guess_rad_type;
  IO.read_info("guess_rad_type", "slater", guess_rad_type);

  for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {

    guess_card_getnlm(IO, ifun, n0, l0, m0);
    if (m0 != Bas.mval[ifun]) {
      std::cout << "inconsistent mval." << std::endl;
      abort();
    }

    // Condon Shortley factor
    faccs = ONE;
    //    if (clcontrol::docs1 && Bas.mval[ifun] > 0 && Bas.mval[ifun] % 2 == 1) faccs = -ONE;

    if (l0 >= 0) {
      printf("guess_card: ifun = %3ld n = %3ld l = %3ld m = %3ld\n", ifun, n0, l0, m0);
      if (guess_rad_type.compare("slater") == 0) {
	guess_aufbau_getneff(n0, l0, m0, nel, Bas.znuc, neff, xi);
	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
         	                + l0 * (Bas.GRad.nrad - 1) + irad - 1;
	  Wfn.wfn[ind] = guess_rfun_slater(neff, xi, Bas.GRad.xrad[irad], Bas.GRad.wrad[irad]);
	  // Condon Shortley factor
	  Wfn.wfn[ind] = Wfn.wfn[ind] * faccs;
	}
      } else if (guess_rad_type.compare("hlike") == 0) {
	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
         	               +  l0 * (Bas.GRad.nrad - 1) + irad - 1;
	  Wfn.wfn[ind] = guess_rfun_hlike(n0, l0, Bas.znuc, Bas.GRad.xrad[irad], Bas.GRad.wrad[irad]);
	  // Condon Shortley factor
	  Wfn.wfn[ind] = Wfn.wfn[ind] * faccs;
	}
      } else if (guess_rad_type.compare("hlike_ecs") == 0) {
	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
	    +  l0 * (Bas.GRad.nrad - 1) + irad - 1;
	  Wfn.wfn[ind] = guess_rfun_hlike_ecs(n0, l0, Bas.znuc, Bas.GRad.xrad[irad], Bas.GRad.cwrad[irad], Bas);
	  // Condon Shortley factor
	  Wfn.wfn[ind] = Wfn.wfn[ind] * faccs;
	}
      } else if (guess_rad_type.compare("slater_ecs") == 0) {
	guess_aufbau_getneff(n0, l0, m0, nel, Bas.znuc, neff, xi);
	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
	    + l0 * (Bas.GRad.nrad - 1) + irad - 1;
	  Wfn.wfn[ind] = guess_rfun_slater_ecs(neff, xi, Bas.GRad.xrad[irad], Bas.GRad.cwrad[irad], Bas);
	  // Condon Shortley factor
	  Wfn.wfn[ind] = Wfn.wfn[ind] * faccs;
	}
      }  
    } else {
      std::cout << "guess_card: bad l0." << std::endl;
      abort();
    }
  }
}
////////////////////////////////////////////////////////////////////////
void guess_card_getnlm(const clio& IO, long ifun, long& n, long& l, long& m)
// Read guess information
{
  std::string key = "orbital:";
  std::ifstream ifs(IO.inp.c_str());
  ifs.seekg(0, std::ios::beg);

  std::string line;
  std::stringstream ioss;

  while (getline(ifs, line) && line.find(key,0) == std::string::npos) {}
  if(! ifs.eof()) {
    long jfun = 0;
    while (jfun < ifun) {
      getline(ifs, line);
      jfun ++;
    }

    getline(ifs, line);
    ioss.str("");
    ioss.clear(std::stringstream::goodbit);
    ioss << line.c_str();
    ioss >> n
	 >> l
	 >> m;

  } else {

    std::cout<< "No guess information." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void guess_aufbau_getneff(long n, long l, long m, long nel, double znuc, double& neff, double& xi)
{
  double scon;

  if (n == 1) {
    if (nel == 1) {
      scon = 0.0;
    } else {
      scon = 0.3;
    }
  } else if (n == 2) {
    if (nel == 1) {
      scon = 0.0;
    } else if (nel == 2) {
      scon = 0.85 * 2;
    } else {
      scon = 0.85 * 2 + 0.35 * (std::min(nel - 2, (long) 8) - 1);
    }
  } else if (n == 3) {
    if (l == 0 || l == 1) {
      if (nel == 1) {
	scon = 0.0;
      } else if (nel <= 2) {
	scon = nel * 1.0;
      } else if (nel <= 10) {
	scon = 2.0 + (nel - 2) * 0.85;
      } else {
	scon = 2.0 + 8 * 0.85 + 0.35 * (std::min(nel - 10, (long) 8) - 1);
      }
    } else if (l == 2) {
      if (nel == 1) {
	scon = 0.0;
      }	else if (nel <= 18) {
	scon = nel * 1.0;
      } else if (nel <= 20) {
	scon = 18.0;
      } else {
	scon = 18.0 + 0.35 * (std::min(nel - 20, (long) 10) - 1);
      }
    }
  } else if (n == 4) {
    if (l == 0 || l == 1) {
      if (nel == 1) {
	scon = 0.0;
      } else if (nel <= 10) {
	scon = nel * 1.0;
      } else if (nel <= 18) {
	scon = 10.0 + (nel - 10) * 0.85;
      } else if (nel <= 20) {
	scon = 18.0 + (nel - 18) * 0.35;
      } else if (nel <= 30) {
	scon = 18.0 + 2 * 0.35 + (nel - 20) * 0.85;
      } else {
	scon = 18.0 + 2 * 0.35 + 10 * 0.85 + 0.35 * (std::min(nel - 30, (long) 6) - 1);
      }
    } else if (l == 2) {
      if (nel == 1) {
	scon = 0.0;
      } else if (nel <= 36) {
	scon = nel * 1.0;
      } else if (nel <= 38) {
	scon = 36.0;
      } else {
	scon = 36.0 + 0.35 * (std::min(nel - 38, (long) 10) - 1);
      }
    } else {
      printf("guess_aufbau_getneff: l > 2 nyi.\n");
      abort();      
    }
  } else if (n == 5) {
    if (l == 0 || l == 1) {
      if (nel == 1) {
	scon = 0.0;
      } else if (nel <= 18) {
	scon = nel * 1.0;
      } else if (nel <= 20) {
	scon = 10.0 + (nel - 20) * 0.85;
      } else if (nel <= 30) {
	scon = 10.0 + 2 * 0.85 + (nel - 20) * ONE;
      } else if (nel <= 36) {
	scon = 10.0 + 2 * 0.85 + 10.0 + (nel - 30) * 0.85;
      } else if (nel <= 38) {
	scon = 10.0 + 2 * 0.85 + 10.0 + 6 * 0.85 + (nel - 36) * 0.35;
      } else if (nel <= 48) {
	scon = 10.0 + 2 * 0.85 + 10.0 + 6 * 0.85 + 2 * 0.35 + (nel - 38) * 0.85;
      }	else if (nel <= 54) {
	scon = 10.0 + 2 * 0.85 + 10.0 + 6 * 0.85 + 2 * 0.35 + 10 * 0.85 + 0.35 * (std::min(nel - 48, (long) 6) - 1);
      }
    } else {
      printf("guess_aufbau_getneff: nl > 5sp nyi.\n");
      abort();      
    }
  } else {
    printf("guess_aufbau_getneff: n > 5sp nyi.\n");
    abort();
  }

  neff = double(n);
  xi = (znuc - scon) / neff;

  printf("guess_aufbau_getneff: nlm = %5ld%5ld%5ld. neff = %10.5f xi = %10.5f\n", n, l, m, neff, xi);
}
////////////////////////////////////////////////////////////////////////
//void guess_g09(const clmpi& Proc, const clio& IO, const clbas& Bas, clwfn& Wfn)
//{
//  clgbas Gauss(IO);
//  if (Bas.ORMAS.nfun > Gauss.ngfun) {
//    printf("no enough functions in fck.\n");
//    abort();
//  }
//  if (Bas.GAng.lmax < Gauss.glmax) {
//    printf("too large lmax in fck.\n");
//    abort();
//  }
//
//  long l, l0, lm, lmp, lmm, mu, icmo, ind, pind, mind;
//  double cp, ap, rpt, wpt, rsph;
//  dcomplex cfacp, cfacm;
//
//  std::vector<double> rfun(Bas.GRad.nrad - 1);
//  long bdone = 0;
//  long pdone = 0;
//
//  for (long ishell = 0; ishell < Gauss.nshell; ishell ++) {
//    // radial wavefunction
//    l = get_abs(Gauss.type[ishell]);
//    for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
//      rfun[irad - 1] = ZERO;
//    }
//    for (long iprm = pdone; iprm < pdone + Gauss.nprm[ishell]; iprm ++) {
//      cp = Gauss.cont[iprm];
//      ap = Gauss.alph[iprm];
//      for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
//	rpt = Bas.GRad.xrad[irad];
//	rfun[irad - 1] += cp * exp(-ap * rpt * rpt);
//      }
//    }
//    pdone += Gauss.nprm[ishell];
//
//    // fedvr coefficients of radial distribution function
//    for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
//      rpt = Bas.GRad.xrad[irad];
//      wpt = Bas.GRad.wrad[irad];
//      rfun[irad - 1] *= rpt * sqrt(wpt) * pow(rpt, l);
//    }
//
//    // lcao expansion
//    l0 = l * (l + 1);
//    if (l != 1) {
//      mu = 0;
//      // m = 0
//      for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
//	icmo = ifun * Gauss.ngbas + bdone + mu;
//	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
//	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//         	                + l0 * (Bas.GRad.nrad - 1) + irad - 1;
//	  Wfn.wfn[ind] += rfun[irad - 1] * Gauss.cmo[icmo];
//	}      
//      }
//      mu ++;
//
//      // m = +k, -k
//      for (long m = 1; m <= l; m ++) {
//	lmp = l0 + m;
//	lmm = l0 - m;
//	for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
//	  cfacp = ONE / sqrt(TWO) * pow(-ONE, 1);
//	  cfacm = ONE / sqrt(TWO);
//	  icmo = ifun * Gauss.ngbas + bdone + mu;
//	  for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
//	    rsph = rfun[irad - 1] * Gauss.cmo[icmo];
//	    pind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//         	                  + lmp * (Bas.GRad.nrad - 1) + irad - 1;
//	    mind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//         	                  + lmm * (Bas.GRad.nrad - 1) + irad - 1;
//	    Wfn.wfn[pind] += rsph * cfacp;
//	    Wfn.wfn[mind] += rsph * cfacm;
//	  }      
//
//	  cfacp = ONE / sqrt(TWO) * pow(-ONE, 1) * IUNIT;
//	  cfacm = ONE / sqrt(TWO) * (-IUNIT);
//	  icmo = ifun * Gauss.ngbas + bdone + mu + 1;
//	  for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
//	    rsph = rfun[irad - 1] * Gauss.cmo[icmo];
//	    pind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//         	                  + lmp * (Bas.GRad.nrad - 1) + irad - 1;
//	    mind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//         	                  + lmm * (Bas.GRad.nrad - 1) + irad - 1;
//	    Wfn.wfn[pind] += rsph * cfacp;
//	    Wfn.wfn[mind] += rsph * cfacm;
//	  }
//	}
//	mu += 2;
//      }
//    } else {
//      mu = 0;
//      // p_x
//      for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
//	icmo = ifun * Gauss.ngbas + bdone + mu;
//	cfacp = ONE / sqrt(TWO) * pow(-ONE, 1);
//	cfacm = ONE / sqrt(TWO);
//	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
//	  rsph = rfun[irad - 1] * Gauss.cmo[icmo];
//	  pind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//	                   + (l0 + 1) * (Bas.GRad.nrad - 1) + irad - 1;
//	  mind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//       	                   + (l0 - 1) * (Bas.GRad.nrad - 1) + irad - 1;
//	  Wfn.wfn[pind] += rsph * cfacp;
//	  Wfn.wfn[mind] += rsph * cfacm;
//	}      
//      }
//      mu ++;
//      // p_y
//      for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
//	icmo = ifun * Gauss.ngbas + bdone + mu;
//	cfacp = ONE / sqrt(TWO) * pow(-ONE, 1) * IUNIT;
//	cfacm = ONE / sqrt(TWO) * (-IUNIT);
//	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
//	  rsph = rfun[irad - 1] * Gauss.cmo[icmo];
//	  pind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//	                   + (l0 + 1) * (Bas.GRad.nrad - 1) + irad - 1;
//	  mind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//       	                   + (l0 - 1) * (Bas.GRad.nrad - 1) + irad - 1;
//	  Wfn.wfn[pind] += rsph * cfacp;
//	  Wfn.wfn[mind] += rsph * cfacm;
//	}      
//      }
//      mu ++;
//      // p_z
//      for (long ifun = 0; ifun < Bas.ORMAS.nfun; ifun ++) {
//	icmo = ifun * Gauss.ngbas + bdone + mu;
//	for (long irad = 1; irad < Bas.GRad.nrad; irad ++) {
//	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
//         	                + l0 * (Bas.GRad.nrad - 1) + irad - 1;
//	  Wfn.wfn[ind] += rfun[irad - 1] * Gauss.cmo[icmo];
//	}      
//      }
//      mu ++;
//    }
//    bdone += 2 * get_abs(Gauss.type[ishell]) + 1;
//  }
//  // debug
//  printf("guess_g09: pdone = %5ld\n", pdone);
//  printf("guess_g09: bdone = %5ld\n", bdone);
//  // debug
//}
////////////////////////////////////////////////////////////////////////
