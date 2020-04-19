////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
void calprop_print_rrad(const clmpi& Proc, const clio& IO, const clbas& Bas, 
			const clfield& Field, clhprod& HPW, const clwfn& Wfn)
{
  double wgt, lfield[9];
  Field.get_value(lfield);
  if (IO.nprint_rrad > 0 && Field.step % IO.nprint_rrad == 0) {
    double drrad, rradp, rradm;
    std::vector<dcomplex> rrad(Bas.GRad.nrad + 1);
    std::vector<dcomplex> rradpw((Bas.GRad.nrad + 1)*(Bas.GAng.lmax1 + 1));
    HPW.mkrradx(Wfn, rrad, rradpw);
    //    HPW.mkrrad(Wfn, rrad, rradpw);
    //    HPW.mkrrad1(Wfn, rrad, rradpw);
    fprintf(IO.fp_rrad, "###\n");
    for (int irad = 1; irad < Bas.GRad.nrad; irad ++) {
      wgt = Bas.GRad.wrad[irad];
      //      drrad = rrad[irad].real();
      drrad = rrad[irad].real() * wgt;
      fprintf(IO.fp_rrad, " %20.10e %20.10e", Bas.GRad.xrad[irad], drrad);
      for (int l = 0; l <= Bas.GAng.lmax1; l++) {
	//	drrad = rradpw[(Bas.GRad.nrad+1)*l+irad].real();
	drrad = rradpw[(Bas.GRad.nrad+1)*l+irad].real() * wgt;
	fprintf(IO.fp_rrad, "%20.10e", drrad);
      }
      fprintf(IO.fp_rrad, "\n");
    }
    fprintf(IO.fp_rrad, "\n");
    fprintf(IO.fp_rrad, "\n");
  }
}
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  printf("# td1c: start    "); 
  system("date");

  clmpi Proc;
  clio IO(argc, argv);
  clcontrol CTRL(IO);
  clbas Bas(Proc, IO);
  clwfn Wfn(Proc, IO, Bas);

  guess(Proc, IO, Bas, Wfn);

  clfield Field(Proc, IO);
  clhprod HPW(Proc, IO, Bas, Field);
  clpes PES1(Proc, IO, Bas, HPW);

  bool print_orb;
  int print_orb_nstep;
  IO.read_info("print_orb", false, print_orb);

  Field.step = ZERO;
  Field.time = ZERO;
  //  tdse_print(Proc, IO, Bas, Field, HPW, Wfn);
  calprop_print_rrad(Proc, IO, Bas, Field, HPW, Wfn);

  system("date");
  printf("# td1c: end      "); 
}
////////////////////////////////////////////////////////////////////////
