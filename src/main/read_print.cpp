////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
void read_print_read(const clmpi& Proc, const clio& IO, const clbas& Bas, 
		     const clfield& Field, clhprod& HPW, clwfn& Wfn)
{
  double wfnr, wfni;
  // read in orbitals
  fprintf(IO.fp_torb, "###\n");
  fprintf(IO.fp_torb, "%10d\n", Bas.GRad.nrad);
  fprintf(IO.fp_torb, "%10d\n", Bas.GAng.lmax1);
  fprintf(IO.fp_torb, "%10d\n", Bas.GAng.mmax1);
  fprintf(IO.fp_torb, "%10d\n", Bas.ORMAS.nfun);

  char line[256];
  int nrad0, lmax0, mmax0, nfun0;
  fscanf(IO.fp_torb, "%s\n", line);
  fscanf(IO.fp_torb, "%d\n", &nrad0);
  fscanf(IO.fp_torb, "%d\n", &lmax0);
  fscanf(IO.fp_torb, "%d\n", &mmax0);
  fscanf(IO.fp_torb, "%d\n", &nfun0);
  int ind;
  for (int ifun = 0; ifun < nfun0; ifun ++) {
    for (int l = 0; l <= lmax0; l ++) {
      for (int irad = 1; irad < nrad0; irad ++) {
	fscanf(IO.fp_torb, "%le%le\n", &wfnr, &wfni); 
	if (ifun < Bas.ORMAS.nfun && l <= Bas.GAng.lmax1 && irad < Bas.GRad.nrad) {
	  ind = ifun * Bas.GAng.nsph1 * (Bas.GRad.nrad - 1) 
          	                 + l * (Bas.GRad.nrad - 1) + irad - 1;
	  Wfn.wfn[ind] = wfnr * RUNIT + wfni * IUNIT;
	}
      }
    }
  }
  //  fscanf(IO.fp_torb, "\n", line);
  //  fscanf(IO.fp_torb, "\n", line);
  fscanf(IO.fp_torb, "\n");
  fscanf(IO.fp_torb, "\n");

  // read in ci coefficients
  int lcic0;
  fscanf(IO.fp_tcic, "%s\n", line);
  fscanf(IO.fp_tcic, "%d\n", &lcic0);
  if (Bas.ORMAS.lcic != lcic0) {
    printf("clwfn::read_print_read: lcic != lcic0\n");
    abort();
  } else {
    for (int idet = 0; idet < Bas.ORMAS.lcic; idet ++) {
      fscanf(IO.fp_tcic, "%le%le\n", &wfnr, &wfni);
      Wfn.wfn[Wfn.size1+idet] = wfnr * RUNIT + wfni * IUNIT;
    }
  }
  fscanf(IO.fp_tcic, "\n");
  fscanf(IO.fp_tcic, "\n");
}
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  printf("# read_print: start    \n"); 
  system("date");

  clmpi Proc;
  clio IO(argc, argv);
  clcontrol CTRL(IO);
  clbas Bas(Proc, IO);
  clwfn Wfn(Proc, IO, Bas, true);
  guess(Proc, IO, Bas, Wfn);

  clfield Field(Proc, IO);
  Field.init();

  clhprod HPW(Proc, IO, Bas, Field);
  HPW.set_wfn0(Proc, IO, Bas, Wfn);
  if (IO.nprint_op1tr > 0) {
    HPW.op1tr_init(IO, Bas, Wfn);
    HPW.op1tr_printp_tag(IO, Bas);
    HPW.op1tr_printq_tag(IO, Bas);
  }

  //debug
  //std::cout << "for debug 1." << std::endl;
  //exit(1);
  //debug

  while (! Field.finished()) {
    if (IO.nread_full > 0 && Field.step % IO.nread_full == 0) {
      read_print_read(Proc,IO,Bas,Field,HPW,Wfn);
      tdse_print(Proc, IO, Bas, Field, HPW, Wfn);
    }
    //debug
    //std::cout << "for debug 2." << std::endl;
    //exit(1);
    //debug
    Field.step ++;
    Field.time += Field.dtime;
  }

  if (IO.nprint_op1tr > 0) HPW.op1tr_final();

  system("date");
  printf("# read_print: end      \n"); 
}
////////////////////////////////////////////////////////////////////////
