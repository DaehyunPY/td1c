////////////////////////////////////////////////////////////////////////
// Job control parameters
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
// static members
////////////////////////////////////////////////////////////////////////
int clcontrol::num_control = 0;
bool clcontrol::fedvr_normalized = true;
int clcontrol::icomp;
int clcontrol::igauge;
int clcontrol::oorot_type;
int clcontrol::split_type;
int clcontrol::iprojfc;
int clcontrol::type_dcx;
bool clcontrol::jfc_implicit;
bool clcontrol::xfc_implicit;
int clcontrol::h1rat_maxcyc;
double clcontrol::h1rat_thresh;
int clcontrol::maxipx;
int clcontrol::maxipd;
double clcontrol::radipx;
bool clcontrol::docs1;
bool clcontrol::docs2;
bool clcontrol::sae;
bool clcontrol::psp;
int clcontrol::psp_type;
int clcontrol::dft_type;
// teramura-tdcis
bool clcontrol::istdcis;
bool clcontrol::tdcis_rvg;
bool clcontrol::print_xrad;
// teramura-tdcis
int clcontrol::reg_type;
double clcontrol::throcc1;
double clcontrol::throcc2;
double clcontrol::throcc3;
int clcontrol::xact2_type;
int clcontrol::xact2_maxitr;
double clcontrol::xact2_thresh;
int clcontrol::ncut_occ3;
bool clcontrol::exact3j;
bool clcontrol::cionly;
int clcontrol::rrad_type;
double clcontrol::rrad_rion;
double clcontrol::rrad_r2in;
double clcontrol::rrad_r2out;
bool clcontrol::tsurff;
////////////////////////////////////////////////////////////////////////
clcontrol::clcontrol()
{
}
////////////////////////////////////////////////////////////////////////
clcontrol::~clcontrol()
{
  if (dft_type != 0) {
    control_xc_final_();
  }
}
////////////////////////////////////////////////////////////////////////
clcontrol::clcontrol(const clio& IO)
{
  gen(IO);
}
////////////////////////////////////////////////////////////////////////
void clcontrol::gen(const clio& IO)
{
  num_control ++;
  std::cout << "clcontrol::gen: num_control = " << num_control << std::endl;
  if (num_control != 1) abort();

  // from clhprod...
  if (IO.job_type.compare("td") == 0) {
    icomp = 1;
  } else {
    icomp = 0;
  }

  std::string gauge;
  IO.read_info("gauge", "length1", gauge);
  if (gauge.compare("length1") == 0 || gauge.compare("length2") == 0) {
    igauge = 0;
  } else if (gauge.compare("velocity1") == 0 || gauge.compare("velocity2") == 0) {
    igauge = 1;
  } else {
    std::cout << "clcontrol::gen: bad gauge." << std::endl;
    abort();
  }

  IO.read_info("oorot_type", LZERO, oorot_type);
  if (oorot_type < 0 || oorot_type > 2) {
    std::cout << "clcontrol::gen: bad oorot_type." << std::endl;
    std::cout << "clcontrol::gen: oorot_type = " << oorot_type << std::endl;
    abort();
  } 

  IO.read_info("split_type", LONE, split_type);
  if (split_type < 0 || split_type > 2) {
    std::cout << "clcontrol::gen: bad split_type." << std::endl;
    std::cout << "clcontrol::gen: split_type = " << split_type << std::endl;
    abort();
  } 

  if (oorot_type != 0 && split_type == 0) {
    std::cout << "clcontrol::gen: bad oorot/split_type." << std::endl;
    abort();
  } 

  IO.read_info("iprojfc", 2, iprojfc);
  if (iprojfc != 1 && iprojfc != 2) {
    std::cout << "clcontrol::gen: bad projfc." << std::endl;
    std::cout << "clcontrol::gen: projfc = " << iprojfc << std::endl;
    abort();
  } 

  IO.read_info("type_dcx", LZERO, type_dcx);
  IO.read_info("maxipx", LONE, maxipx);
  IO.read_info("maxipd", LONE, maxipd);
  IO.read_info("radipx", 20.0, radipx);
  IO.read_info("jfc_implicit", true, jfc_implicit);
  IO.read_info("xfc_implicit", true, xfc_implicit);
  IO.read_info("h1rat_maxcyc", 20, h1rat_maxcyc);
  IO.read_info("h1rat_thresh", (double) 1.0E-10, h1rat_thresh);
  IO.read_info("docs1", false, docs1);
  IO.read_info("docs2", false, docs2);
  IO.read_info("sae", false, sae);
  IO.read_info("psp", false, psp);
  IO.read_info("psp_type", LONE, psp_type);
  IO.read_info("dft_type", LZERO, dft_type);
  IO.read_info("reg_type", 0, reg_type);

// teramura-tdcis
  IO.read_info("istdcis", false, istdcis);
  IO.read_info("tdcis_rvg", false, tdcis_rvg);
  if(igauge == 0 && tdcis_rvg == true) {
    std::cout << "clcontrol::gen: bad rvg and igauge" << std::endl;
    std::cout << "clcontrol::gen: rvg = "<< tdcis_rvg
       << ", igauge = " << igauge << std::endl;
    abort();
  }
  IO.read_info("print_xrad", false, print_xrad);
// teramura-tdcis

//  throcc1 is for D^{-1} in Q terms
//  throcc2 is for (2-D)^{-1} in core-virtual terms
//  throcc3 is for A matrix in active-active terms
  IO.read_info("throcc1", 1.0E-15, throcc1);
  IO.read_info("throcc2", 1.0E-15, throcc2);
  IO.read_info("throcc3", 1.0E-15, throcc3);
//  IO.read_info("throcc1", throcc1);
//  IO.read_info("throcc2", throcc2);
//  IO.read_info("throcc3", throcc3);

//  Algorithm for active-active rotation solver
  IO.read_info("xact2_type", 0, xact2_type);
  IO.read_info("xact2_maxitr", 100, xact2_maxitr);
  IO.read_info("xact2_thresh", 1.0E-15, xact2_thresh);

  IO.read_info("ncut_occ3", 0, ncut_occ3);
  IO.read_info("exact3j", false, exact3j);
  IO.read_info("cionly", false, cionly);
  IO.read_info("rrad_type", -1, rrad_type);
  IO.read_info("rrad_rion", 20E0, rrad_rion);
  IO.read_info("rrad_r2in", 0E0, rrad_r2in);
  IO.read_info("rrad_r2out", 1E+10, rrad_r2out);
// Sato_tSURFF
  IO.read_info("tsurff", false, tsurff);
// Sato_tSURFF

  if (! jfc_implicit && xfc_implicit) {
    std::cout << "clcontrol::gen: jfc_implicit=F && xfc_implicit=T." << std::endl;
    abort();
  }

  if (dft_type != 0) {
//    int x_type, c_type;
//    IO.read_info("x_type", 1, x_type); // LDAx
//    IO.read_info("c_type", 9, c_type); // LDAc (Perdew-Zunger 81)
    control_xc_init_(&dft_type);
  }

  // c++/fortran bindings
// tdcis-teramura
  control_bind_(&icomp, &igauge, &oorot_type, &split_type, &iprojfc, &type_dcx, 
		&jfc_implicit, &xfc_implicit, &h1rat_maxcyc, &h1rat_thresh, 
		&docs1, &docs2, &sae, &psp, &psp_type, &dft_type, &reg_type, &throcc1, &throcc2,
		&throcc3, &xact2_type, &xact2_maxitr, &xact2_thresh, &ncut_occ3, &exact3j, 
		&cionly, &istdcis, &tdcis_rvg);
// tdcis-teramura
}
//////////////////////////////////////////////////////////////////////////
