////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  printf("# pes: start    "); 
  system("date");

  clmpi Proc;
  clio IO(argc, argv);
  clcontrol CTRL(IO);
  clbas Bas(Proc, IO);
  clwfn Wfn(Proc, IO, Bas);

  guess(Proc, IO, Bas, Wfn);

  clfield Field(Proc, IO);
  clhprod HPW(Proc, IO, Bas, Field);

  FILE *fppes;
  double kval;

  clpes PES1(Proc, IO, Bas, HPW);

  //  PES1.spec1_kz(Proc, IO, Bas, HPW, Wfn);
  //  fppes = fopen(IO.rhokz.c_str(), "w");
  PES1.spec1_k(Proc, IO, Bas, HPW, Wfn);
  fppes = fopen(IO.rhok.c_str(), "w");

  for (int ik = 0; ik <= PES1.pes_numk; ik ++) {
    kval = PES1.pes_k_min + ik * PES1.pes_k_step;
    fprintf(fppes, "%10d %20.10f %20.10f\n", ik, kval, PES1.pes_rhok[ik].real());
  }
  fclose(fppes);

  system("date");
  printf("# pes: end      "); 
}
////////////////////////////////////////////////////////////////////////
