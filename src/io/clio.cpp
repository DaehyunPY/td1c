////////////////////////////////////////////////////////////////////////
// IO: basics
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
clio::clio()
{
  name = "test";
  gen();
}
////////////////////////////////////////////////////////////////////////
clio::clio(std::string inp_)
{
  gen(inp_);
}
////////////////////////////////////////////////////////////////////////
clio::clio(int argc, char** argv)
{
  gen(argc, argv);
}
////////////////////////////////////////////////////////////////////////
clio::~clio()
{
  if (nprint_ene > 0) fclose(fp_ene);
  if (nprint_dip > 0) fclose(fp_op1);
  if (nprint_opx > 0) fclose(fp_op0);
  if (nprint_opx > 0) fclose(fp_opd);
  if (nprint_opx > 0) fclose(fp_opx);
  if (nprint_opn > 0) {
    fclose(fp_dipn);
    fclose(fp_veln);
    fclose(fp_accn);
  }
  if (nprint_ion > 0) fclose(fp_ipx);
  if (nprint_ipd > 0) fclose(fp_ipd);
  if (nprint_rrad > 0) fclose(fp_rrad);
  if (nprint_full > 0) {
    fclose(fp_torb);
    fclose(fp_tcic);
  }
}
////////////////////////////////////////////////////////////////////////
void clio::print() const
{
  printf("# name      = %s\n", name.c_str());
  printf("# inp  file = %s\n", inp.c_str());
  printf("# orbp file = %s\n", orbp.c_str());
  printf("# orb  file = %s\n", orb.c_str());
  printf("# cic  file = %s\n", cic.c_str());
  printf("# rhok file = %s\n", rhok.c_str());
  printf("# rhokz file = %s\n", rhokz.c_str());
  if (nprint_ene > 0) printf("# ene file = %s\n", ene.c_str());
  if (nprint_dip > 0) printf("# op1 file = %s\n", op1.c_str());
  if (nprint_opx > 0) printf("# op0 file = %s\n", op0.c_str());
  if (nprint_opx > 0) printf("# opd file = %s\n", opd.c_str());
  if (nprint_opx > 0) printf("# opx file = %s\n", opx.c_str());
  if (nprint_opn > 0) {
    printf("# dipn file = %s\n", dipn.c_str());
    printf("# veln file = %s\n", veln.c_str());
    printf("# accn file = %s\n", accn.c_str());
  }
  if (nprint_ion > 0) printf("# ipx file = %s\n", ipx.c_str());
  if (nprint_ipd > 0) printf("# ipd file = %s\n", ipd.c_str());
  if (nprint_rrad > 0) printf("# rrad file = %s\n", rrad.c_str());
  if (nprint_full > 0) {
    printf("# torb file = %s\n", torb.c_str());
    printf("# tcic file = %s\n", tcic.c_str());
  }
}
////////////////////////////////////////////////////////////////////////
void clio::gen(std::string inp_)
{
  name = inp_.substr(0, inp_.find_last_of("."));
  gen();
}
////////////////////////////////////////////////////////////////////////
void clio::gen(int argc, char** argv)
{
  if ( argc < 1 ) {
    printf( "Specify input file name." );
    exit( 1 );
  }
  inp = argv[1];
  name = inp.substr(0, inp.find_last_of("."));
  gen();
}
////////////////////////////////////////////////////////////////////////
void clio::gen()
{
  inp = name; inp += ".inp";
  orb = name; orb += ".orb";
  cic = name; cic += ".cic";
  ene = name; ene += ".ene";
  op1 = name; op1 += ".op1";
  op0 = name; op0 += ".op0";
  opd = name; opd += ".opd";
  opx = name; opx += ".opx";
  dipn = name; dipn += ".dipn";
  veln = name; veln += ".veln";
  accn = name; accn += ".accn";
  ipx = name; ipx += ".ipx";
  ipd = name; ipd += ".ipd";
  orbp = name; orbp += ".orbp";
  rhok = name; rhok += ".rhok";
  rhokz = name; rhokz += ".rhokz";
  rrad = name; rrad += ".rrad";
  torb = name; torb += ".torb";
  tcic = name; tcic += ".tcic";

  read_info("job_type", "init", job_type);
  read_info("iprint", LONE, iprint);
  if (job_type.compare("init") == 0) {
    read_info("nprint_ene", -LONE, nprint_ene);
    read_info("nprint_dip", -LONE, nprint_dip);
    read_info("nprint_opx", -LONE, nprint_opx);
    read_info("nprint_opn", -LONE, nprint_opn);
    read_info("nprint_ion", -LONE, nprint_ion);
    read_info("nprint_ipd", -LONE, nprint_ipd);
    read_info("nprint_norm",-LONE, nprint_norm);
    read_info("nprint_rrad", -LONE, nprint_rrad);
    read_info("nprint_full", -LONE, nprint_full);
  } else {
    read_info("nprint_dip", LONE, nprint_dip);
    read_info("nprint_ene", LONE, nprint_ene);
    read_info("nprint_ion", LONE, nprint_ion);
    read_info("nprint_opx", -LONE, nprint_opx);
    read_info("nprint_opn", -LONE, nprint_opn);
    read_info("nprint_ipd", -LONE, nprint_ipd);
    read_info("nprint_norm",-LONE, nprint_norm);
    read_info("nprint_rrad", -LONE, nprint_rrad);
    read_info("nprint_full", -LONE, nprint_full);
  }

  if (nprint_ene > 0) fp_ene = fopen(ene.c_str(), "w");
  if (nprint_dip > 0) fp_op1 = fopen(op1.c_str(), "w");
  if (nprint_opx > 0) fp_op0 = fopen(op0.c_str(), "w");
  if (nprint_opx > 0) fp_opd = fopen(opd.c_str(), "w");
  if (nprint_opx > 0) fp_opx = fopen(opx.c_str(), "w");
  if (nprint_opn > 0) {
    fp_dipn = fopen(dipn.c_str(), "w");
    fp_veln = fopen(veln.c_str(), "w");
    fp_accn = fopen(accn.c_str(), "w");
  }
  if (nprint_ion > 0) fp_ipx = fopen(ipx.c_str(), "w");
  if (nprint_ipd > 0) fp_ipd = fopen(ipd.c_str(), "w");
  if (nprint_rrad > 0) fp_rrad = fopen(rrad.c_str(), "w");
  if (nprint_full > 0) {
    fp_torb = fopen(torb.c_str(), "w");
    fp_tcic = fopen(tcic.c_str(), "w");
  }
  print();
  io_bind_(name.c_str(), name.length());
}
////////////////////////////////////////////////////////////////////////
