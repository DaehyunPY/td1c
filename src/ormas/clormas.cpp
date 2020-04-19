////////////////////////////////////////////////////////////////////////
// ORMAS parameters
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
clormas::clormas()
{
}
////////////////////////////////////////////////////////////////////////
clormas::clormas(const clmpi& MPIP, const clio& IO)
{
  gen(MPIP, IO);
}
////////////////////////////////////////////////////////////////////////
clormas::~clormas()
{
  std::cout << "~clormas" << std::endl;
  ormas_final_();
}
////////////////////////////////////////////////////////////////////////
void clormas::gen(const clmpi& MPIP, const clio& IO)
{
  long nela, nelb;

  IO.read_info("donly", false, donly);
  IO.read_info("den2_abonly", false, den2_abonly);
  IO.read_info("tdcc", false, tdcc);
  IO.read_info("cc_type", 0, cc_type);
  IO.read_info("cc_code", "manual", cc_code);

  IO.read_info("nela", nela);
  IO.read_info("nelb", nelb);
  IO.read_info("nfun", nfun);
  IO.read_info("nblock", nblock);
  type_block.resize(nblock);
  nfun_block.resize(nblock);
  IO.read_info("type_block", type_block);
  IO.read_info("nfun_block", nfun_block);

  IO.read_info("nsub", LONE, nsub);
  norb_sub.resize(nsub);
  lorb_sub.resize(nsub*LTWO);
  min_sub.resize(nsub);
  max_sub.resize(nsub);
  if (nsub > 1) {
    IO.read_info("norb_sub", norb_sub);
    IO.read_info("min_sub", min_sub);
    IO.read_info("max_sub", max_sub);
  }
  IO.read_info("thradfc", -ONE, thradfc);
  IO.read_info("fab_den2", false, fab_den2);

  neltot[0] = nela;
  neltot[1] = nelb;
  neltot[2] = neltot[0] + neltot[1];

  mval.resize(nfun);
  froz.resize(nfun);
  read_mval(IO);
//printf("mval @ clormas::gen: %10ld\n", mval[0]);
//printf("mval @ clormas::gen: %10ld\n", mval[1]);
  ormas_bind_(&nblock, &type_block[0], &nfun_block[0], &nfcore2, &nfcore1, &nfcore, &ndcore, &ncore, &nact,
	      &nocc, &nvir, &nfun, &mval[0], &froz[0], nelcore, nelact, neltot, &nsub, &norb_sub[0], &lorb_sub[0],
	      &min_sub[0], &max_sub[0], &nstr_alph, &nstr_beta, &ndet, &lcic, &ndetx, &thradfc, &fab_den2, &donly,
	      &den2_abonly, &tdcc, &cc_type, cc_code.c_str(), cc_code.length());

  long mtot, iprint;
  IO.read_info("mtot", mtot);
  IO.read_info("iprint", LZERO, iprint);
  ormas_init_(&iprint, &mtot);
  check_orb(IO);


  //DEBUG
  //DEBUG  std::cout << "abort in clormas::gen for debug." << std::endl;
  //DEBUG  abort();
  //DEBUG
}
////////////////////////////////////////////////////////////////////////
void clormas::cic0(dcomplex* cic) const
{
  ormas_cic0_(&cic[0]);
}
////////////////////////////////////////////////////////////////////////
void clormas::check_orb(const clio& IO) const
{
  std::string key = "orbital:";
  std::ifstream ifs(IO.inp.c_str());
  ifs.seekg(0, std::ios::beg);

  std::string line;
  std::stringstream ioss;
  long ntmp, ltmp, mtmp, orb_type;
  long nfcore2x = 0;
  long nfcore1x = 0;
  long ndcorex = 0;
  long nactx = 0;
  long nvirx = 0;

  while ( getline(ifs, line) && line.find(key,0) == std::string::npos ) {}
  if(! ifs.eof()) {
    for (long ifun = 0; ifun < nfun; ifun ++) {
      getline(ifs, line);
      ioss.str("");
      ioss.clear(std::stringstream::goodbit);
      ioss << line.c_str();
      ioss >> ntmp
	   >> ltmp
	   >> mtmp
	   >> orb_type;
      if (orb_type == -2) {
	nfcore2x ++;
      } else if (orb_type == -1) {
	nfcore1x ++;
      } else if (orb_type == 0) {
	ndcorex ++;
      } else if (orb_type == 1) {
	nactx ++;
      } else if (orb_type == 2) {
	nvirx ++;
      }
    }
    if (nfcore2x != nfcore2 ||
	nfcore1x != nfcore1 ||
	ndcorex != ndcore ||
	nactx != nact ||
	nvirx != nvir) {
      std::cout<< "Bad orbital information." << std::endl;
      abort();
    }
    printf(" # norb : %10ld%10ld%10ld%10ld%10ld\n", nfcore2, nfcore1, ndcore, nact, nvir);

  } else {
    std::cout<< "No orbital information." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
void clormas::read_mval(const clio& IO)
//
// Read m value of each orbital
//
{
  std::string key = "orbital:";
  std::ifstream ifs(IO.inp.c_str());
  ifs.seekg(0, std::ios::beg);

  std::string line;
  std::stringstream ioss;
  long ntmp, ltmp, ttmp, afroz;

  while ( getline(ifs, line) && line.find(key,0) == std::string::npos ) {}
  if(! ifs.eof()) {
    for (long ifun = 0; ifun < nfun; ifun ++) {
      getline(ifs, line);
      ioss.str("");
      ioss.clear(std::stringstream::goodbit);
      ioss << line.c_str();
//      ioss >> ntmp
//	   >> ltmp
//	   >> mval[ifun]
//	   >> ttmp;
//      printf("# clormas::read_mval: f%10ld%10ld%10ld%10ld\n", ntmp,ltmp,mval[ifun],ttmp);

      mval[ifun] = 0;
      froz[ifun] = 0;
      ioss >> ntmp
	   >> ltmp
	   >> mval[ifun]
	   >> ttmp
	   >> froz[ifun];
      printf("# clormas::read_mval: %10ld%10ld%10ld%10ld%10ld\n", ntmp,ltmp,mval[ifun],ttmp,froz[ifun]);
    }
  } else {
    std::cout<< "No lmval information." << std::endl;
    abort();
  }
}
////////////////////////////////////////////////////////////////////////
