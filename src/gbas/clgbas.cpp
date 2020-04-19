////////////////////////////////////////////////////////////////////////
// Gaussian basis function
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
clgbas::clgbas(const clio& IO)
{
  fck = IO.name;
  fck += ".fck";

  read_fck_info("Number of basis functions", ngbas);
  read_fck_info("Number of independent functions", ngfun);
  read_fck_info("Highest angular momentum", glmax);
  read_fck_info("Number of contracted shells", nshell);
  read_fck_info("Number of primitive shells", nshelp);

  type.resize(nshell);
  nprm.resize(nshell);
  alph.resize(nshelp);
  cont.resize(nshelp);
  cmo.resize(ngbas * ngbas);

  read_fck_array("Shell types", nshell, type);
  read_fck_array("Number of primitives per shell", nshell, nprm);
  read_fck_array("Primitive exponents", nshelp, alph);
  read_fck_array("Contraction coefficients", nshelp, cont);
  read_fck_array("Alpha MO coefficients", ngbas * ngbas, cmo);

  overlap();
  normalize1();
  overlap();
  normalize2();
  overlap();
}
////////////////////////////////////////////////////////////////////////
clgbas::~clgbas()
{
  type.resize(0);
  nprm.resize(0);
  alph.resize(0);
  cont.resize(0);
  cmo.resize(0);
}
////////////////////////////////////////////////////////////////////////
void clgbas::read_fck_info(const std::string key, std::string& val) const
{
  /* Read a bool val labeled by key from fck. */

  std::ifstream _ifs( fck.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;
  long ind_key;

  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "I" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> val;

    std::cout << "# " << key.c_str() << " = " << val << std::endl;

  } else {

    std::cout << "# " << "Input error: " << key << " not found." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clgbas::read_fck_info(const std::string key, bool& val) const
{
  /* Read a bool val labeled by key from fck. */

  std::ifstream _ifs( fck.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;
  long ind_key;

  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "I" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> val;

    std::cout << "# " << key.c_str() << " = " << val << std::endl;

  } else {

    std::cout << "# " << "Input error: " << key << " not found." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clgbas::read_fck_info(const std::string key, long& val) const
{
  /* Read an integer val labeled by key from fck. */

  std::ifstream _ifs( fck.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;
  long ind_key;

  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "I" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> val;

    std::cout << "# " << key.c_str() << " = " << val << std::endl;

  } else {

    std::cout << "# " << "Input error: " << key << " not found." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clgbas::read_fck_info(const std::string key, double& val) const
{
  /* Read an double val labeled by key from fck. */

  std::ifstream _ifs( fck.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;
  long ind_key;

  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "I" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> val;

    std::cout << "# " << key.c_str() << " = " << val << std::endl;

  } else {

    std::cout << "Input error: " << key << " not found." << std::endl;
    abort();

  }

}
////////////////////////////////////////////////////////////////////////
void clgbas::read_fck_array(const std::string key, long nval, std::vector<long>& val) const
{
  /* Read integer array val labeled by key from fck. */

  std::ifstream _ifs( fck.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;
  long ind_key, nval1;
  const long ncol = 6;

  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "=" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> nval1;

    if (nval == nval1) {

      long ndo;
      long done = 0;
      while ( done < nval ) {
	getline(_ifs, line);
	ioss.str("");
	ioss.clear(std::stringstream::goodbit);
	ioss << line.c_str();

	ndo = std::min(nval - done, ncol);
	for (long ido = 0; ido < ndo; ido ++) {
	  ioss >> val[done + ido];
	}
	done += ndo;
      }

      std::cout << "# " << key.c_str() << ": N = " << nval << std::endl;
      done = 0;
      while ( done < nval ) {
	ndo = std::min(nval - done, ncol);
	for (long ido = 0; ido < ndo; ido ++) {
	  //	  std::cout << val[done + ido];
	  printf("%12ld", val[done + ido]);
	}
	std::cout << std::endl;
	done += ndo;
      }

    } else {

      std::cout << "# " << "Input error: " << key << " not matched." << std::endl;
      abort();
    }

  } else {

    std::cout << "# " << "Input error: " << key << " not found." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clgbas::read_fck_array(const std::string key, long nval, std::vector<double>& val) const
{
  /* Read integer array val labeled by key from fck. */

  std::ifstream _ifs( fck.c_str() );
  _ifs.seekg( 0, std::ios::beg );

  std::string line;
  std::stringstream ioss;
  long ind_key, nval1;
  const long ncol = 5;

  while ( getline(_ifs, line) && line.find(key,0) == std::string::npos ) {}

  if(! _ifs.eof()) {

    ind_key = line.find( "=" ) + 1;

    ioss.str("");
    ioss << line.substr( ind_key );
    ioss >> nval1;

    if (nval == nval1) {

      long ndo;
      long done = 0;
      while ( done < nval ) {
	getline(_ifs, line);
	ioss.str("");
	ioss.clear(std::stringstream::goodbit);
	ioss << line.c_str();

	ndo = std::min(nval - done, ncol);
	for (long ido = 0; ido < ndo; ido ++) {
	  ioss >> val[done + ido];
	}
	done += ndo;
      }

      long expd;
      double data;      
      std::cout << "# " << key.c_str() << ": N = " << nval << std::endl;
      done = 0;
      while ( done < nval ) {
	ndo = std::min(nval - done, ncol);
	for (long ido = 0; ido < ndo; ido ++) {
	  //	  data = val[done + ido];
	  //	  expd = (long) log10(get_abs(data));
	  //	  printf("%.3fe%ld", data / pow(10,expd), expd);
	  printf("%16.8E", val[done + ido]);
	}
	std::cout << std::endl;
	done += ndo;
      }

    } else {

      std::cout << "# " << "Input error: " << key << " not matched." << std::endl;
      abort();
    }

  } else {

    std::cout << "# " << "Input error: " << key << " not found." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clgbas::normalize1()
{
  long l, lfac;
  double cpi, api, cpj, apj, norm, normp, cont2, alph2, fac1, fac2;

  long pdone = 0;
  for (long ishell = 0; ishell < nshell; ishell ++) {
    l = get_abs(type[ishell]);
    lfac = 1;
    for (long k = 0; k <= l; k ++) {
      lfac *= 2 * k + 1;
    }

    norm = ZERO;
    for (long iprm = pdone; iprm < pdone + nprm[ishell]; iprm ++) {
      cpi = cont[iprm];
      api = alph[iprm];
      alph2 = api * TWO;
      fac1 = HALF * sqrt(PI / alph2);
      alph2 *= TWO;
      fac2 = lfac / pow(alph2, l + 1);
      normp = fac1 * fac2;
      printf("iprm = %5ld normp = %20.10f\n", iprm, normp);
      cont[iprm] /= sqrt(normp);
    }
    // 4PI ///////////
    //    norm *= FOUR * PI;
    //    norm *= TWO * PI;
    //////////////////
//    norm = ONE / sqrt(norm);
//    for (long iprm = pdone; iprm < pdone + nprm[ishell]; iprm ++) {
//      cont[iprm] *= norm;
//      //DEBUG if (ishell == 0) {
//      //DEBUG 	cont[iprm] *= 18.2369;
//      //DEBUG }
//      //DEBUG
//    }
    pdone += nprm[ishell];
  }
}
////////////////////////////////////////////////////////////////////////
void clgbas::normalize2()
{
  long l, lfac;
  double cpi, api, cpj, apj, norm, normp, cont2, alph2, fac1, fac2;

  long pdone = 0;
  for (long ishell = 0; ishell < nshell; ishell ++) {
    l = get_abs(type[ishell]);
    lfac = 1;
    for (long k = 0; k <= l; k ++) {
      lfac *= 2 * k + 1;
    }

    norm = ZERO;
    for (long iprm = pdone; iprm < pdone + nprm[ishell]; iprm ++) {
      cpi = cont[iprm];
      api = alph[iprm];
      for (long jprm = pdone; jprm < pdone + nprm[ishell]; jprm ++) {
	cpj = cont[jprm];
	apj = alph[jprm];
	cont2 = cpi * cpj;
	alph2 = api + apj;
	fac1 = HALF * sqrt(PI / alph2);
	alph2 *= TWO;
	fac2 = lfac / pow(alph2, l + 1);
	norm += cont2 * fac1 * fac2;
      }
    }
    // 4PI ///////////
    //    norm *= FOUR * PI;
    //    norm *= TWO * PI;
    //////////////////
    norm = ONE / sqrt(norm);
    for (long iprm = pdone; iprm < pdone + nprm[ishell]; iprm ++) {
      cont[iprm] *= norm;
      //DEBUG if (ishell == 0) {
      //DEBUG 	cont[iprm] *= 18.2369;
      //DEBUG }
      //DEBUG
    }
    pdone += nprm[ishell];
  }
}
////////////////////////////////////////////////////////////////////////
void clgbas::overlap()
{
  long l, lfac;
  double cpi, api, cpj, apj, ovlp, cont2, alph2, fac1, fac2;

  long bdone1 = 0;
  long pdone1 = 0;
  for (long ishell = 0; ishell < nshell; ishell ++) {
    l = get_abs(type[ishell]);
    lfac = 1;
    for (long k = 0; k <= l; k ++) {
      lfac *= 2 * k + 1;
    }

    long bdone2 = 0;
    long pdone2 = 0;
    for (long jshell = 0; jshell < nshell; jshell ++) {
      ovlp = ZERO;
      if (get_abs(type[jshell]) == l) {
	for (long iprm = pdone1; iprm < pdone1 + nprm[ishell]; iprm ++) {
	  cpi = cont[iprm];
	  api = alph[iprm];
	  for (long jprm = pdone2; jprm < pdone2 + nprm[jshell]; jprm ++) {
	    cpj = cont[jprm];
	    apj = alph[jprm];
	    cont2 = cpi * cpj;
	    alph2 = api + apj;
	    fac1 = HALF * sqrt(PI / alph2);
	    alph2 *= TWO;
	    fac2 = lfac / pow(alph2, l + 1);
	    ovlp += cont2 * fac1 * fac2;
	  }
	}
	// 4PI ///////////
	//	ovlp *= FOUR * PI;
	//	ovlp *= TWO * PI;
	//////////////////
      }
      // debug
      printf("clgbas::overlap: %5ld%5ld%20.10f\n", ishell, jshell, ovlp);
      // debug
      pdone2 += nprm[jshell];
    }
    pdone1 += nprm[ishell];
  }
}
////////////////////////////////////////////////////////////////////////
