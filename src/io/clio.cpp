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
  if (nprint_op1 > 0) fclose(fp_op1);
  if (nprint_opx > 0) fclose(fp_op0);
  if (nprint_opx > 0) fclose(fp_opd);
  if (nprint_opx > 0) fclose(fp_opx);
  if (nprint_op1tr > 0) {
    fclose(fp_atrP);
    fclose(fp_atrQ);
  }
  if (nprint_opn > 0) {
    fclose(fp_dipn);
    fclose(fp_veln);
    fclose(fp_accn);
  }

  if (nprint_ipx > 0) fclose(fp_ipx);
  if (nprint_ipd > 0) fclose(fp_ipd);
  if (nprint_opipx > 0) {
    fclose(fp_dipipx);
    fclose(fp_velipx);
    fclose(fp_accipx);
  }
  if (nprint_opipd > 0) {
    fclose(fp_dipipd);
    fclose(fp_velipd);
    fclose(fp_accipd);
  }
  if (nprint_rrad > 0) fclose(fp_rrad);

  if (nprint_full > 0 || nread_full > 0) {
    fclose(fp_torb);
    fclose(fp_tcic);
  }
// tdcis-teramura
  if (nprint_ci0 > 0) fclose(fp_ci0);
// tdcis-teramura
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
  if (nprint_op1 > 0) printf("# op1 file = %s\n", op1.c_str());
  if (nprint_opx > 0) printf("# op0 file = %s\n", op0.c_str());
  if (nprint_opx > 0) printf("# opd file = %s\n", opd.c_str());
  if (nprint_opx > 0) printf("# opx file = %s\n", opx.c_str());
  if (nprint_op1tr > 0) {
    printf("# atrP file = %s\n", atrP.c_str());
    printf("# atrQ file = %s\n", atrQ.c_str());
  }
  if (nprint_opn > 0) {
    printf("# dipn file = %s\n", dipn.c_str());
    printf("# veln file = %s\n", veln.c_str());
    printf("# accn file = %s\n", accn.c_str());
  }
  if (nprint_ipx > 0) printf("# ipx file = %s\n", ipx.c_str());
  if (nprint_ipd > 0) printf("# ipd file = %s\n", ipd.c_str());
  if (nprint_opipx > 0) {
    printf("# dipipx file = %s\n", dipipx.c_str());
    printf("# velipx file = %s\n", velipx.c_str());
    printf("# accipx file = %s\n", accipx.c_str());
  }
  if (nprint_opipd > 0) {
    printf("# dipipd file = %s\n", dipipd.c_str());
    printf("# velipd file = %s\n", velipd.c_str());
    printf("# accipd file = %s\n", accipd.c_str());
  }
  if (nprint_rrad > 0) printf("# rrad file = %s\n", rrad.c_str());
  if (nprint_full > 0 || nread_full > 0) {
    printf("# torb file = %s\n", torb.c_str());
    printf("# tcic file = %s\n", tcic.c_str());
  }
// tdcis-teramura
  if (nprint_ci0 > 0) printf("# ci0  file = %s\n", ci0.c_str());
// tdcis-teramura
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
  atrP = name; atrP += ".atrP";
  atrQ = name; atrQ += ".atrQ";
  dipn = name; dipn += ".dipn";
  veln = name; veln += ".veln";
  accn = name; accn += ".accn";
  ipx = name; ipx += ".ipx";
  ipd = name; ipd += ".ipd";
  dipipx = name; dipipx += ".dipipx";
  velipx = name; velipx += ".velipx";
  accipx = name; accipx += ".accipx";
  dipipd = name; dipipd += ".dipipd";
  velipd = name; velipd += ".velipd";
  accipd = name; accipd += ".accipd";
  orbp = name; orbp += ".orbp";
  rhok = name; rhok += ".rhok";
  rhokz = name; rhokz += ".rhokz";
  rrad = name; rrad += ".rrad";
  torb = name; torb += ".torb";
  tcic = name; tcic += ".tcic";
// tdcis-teramura
  ci0 = name; ci0 += ".ci0";
// tdcis-teramura

  read_info("job_type", "init", job_type);
  read_info("iprint", LONE, iprint);
  if (job_type.compare("init") == 0) {
    read_info("nprint_ene", -LONE, nprint_ene);
    read_info("nprint_op1", -LONE, nprint_op1);
    read_info("nprint_opx", -LONE, nprint_opx);
    read_info("nprint_opn", -LONE, nprint_opn);
    read_info("nprint_op1tr", -LONE, nprint_op1tr);
    read_info("nprint_ipx", -LONE, nprint_ipx);
    read_info("nprint_ipd", -LONE, nprint_ipd);
    read_info("nprint_opipx", -LONE, nprint_opipx);
    read_info("nprint_opipd", -LONE, nprint_opipd);
    read_info("nprint_norm",-LONE, nprint_norm);
    read_info("nprint_rrad", -LONE, nprint_rrad);
    read_info("nprint_full", -LONE, nprint_full);
    read_info("nread_full", -LONE, nread_full);
  } else {
    read_info("nprint_op1", LONE, nprint_op1);
    read_info("nprint_ene", LONE, nprint_ene);
    read_info("nprint_ipx", LONE, nprint_ipx);
    read_info("nprint_opipx", -LONE, nprint_opipx);
    read_info("nprint_opipd", -LONE, nprint_opipd);
    read_info("nprint_opx", -LONE, nprint_opx);
    read_info("nprint_opn", -LONE, nprint_opn);
    read_info("nprint_op1tr", -LONE, nprint_op1tr);
    read_info("nprint_ipd", -LONE, nprint_ipd);
    read_info("nprint_norm",-LONE, nprint_norm);
    read_info("nprint_rrad", -LONE, nprint_rrad);
    read_info("nprint_full", -LONE, nprint_full);
    read_info("nread_full", -LONE, nread_full);
// tdcis-teramra
    read_info("nprint_every", -LONE, nprint_every);
    read_info("nprint_ci0", -LONE, nprint_ci0);
// tdcis-teramra
  }

  if (nprint_ene > 0) fp_ene = fopen(ene.c_str(), "w");
  if (nprint_op1 > 0) fp_op1 = fopen(op1.c_str(), "w");
  if (nprint_opx > 0) fp_op0 = fopen(op0.c_str(), "w");
  if (nprint_opx > 0) fp_opd = fopen(opd.c_str(), "w");
  if (nprint_opx > 0) fp_opx = fopen(opx.c_str(), "w");
  if (nprint_op1tr > 0) {
    fp_atrP = fopen(atrP.c_str(), "w");
    fp_atrQ = fopen(atrQ.c_str(), "w");
  }
  if (nprint_opn > 0) {
    fp_dipn = fopen(dipn.c_str(), "w");
    fp_veln = fopen(veln.c_str(), "w");
    fp_accn = fopen(accn.c_str(), "w");
  }
  if (nprint_ipx > 0) fp_ipx = fopen(ipx.c_str(), "w");
  if (nprint_ipd > 0) fp_ipd = fopen(ipd.c_str(), "w");
  if (nprint_opipx > 0) {
    fp_dipipx = fopen(dipipx.c_str(), "w");
    fp_velipx = fopen(velipx.c_str(), "w");
    fp_accipx = fopen(accipx.c_str(), "w");
  }
  if (nprint_opipd > 0) {
    fp_dipipd = fopen(dipipd.c_str(), "w");
    fp_velipd = fopen(velipd.c_str(), "w");
    fp_accipd = fopen(accipd.c_str(), "w");
  }
  if (nprint_rrad > 0) fp_rrad = fopen(rrad.c_str(), "w");
  if (nprint_full > 0) {
    fp_torb = fopen(torb.c_str(), "w");
    fp_tcic = fopen(tcic.c_str(), "w");
  } else if (nread_full > 0) {
    fp_torb = fopen(torb.c_str(), "r");
    fp_tcic = fopen(tcic.c_str(), "r");
  }
// tdcis-teramura
  if (nprint_ci0 > 0) fp_ci0 = fopen(ci0.c_str(), "w");
// tdcis-teramura
  print();
  io_bind_(name.c_str(), name.length());
}
////////////////////////////////////////////////////////////////////////
