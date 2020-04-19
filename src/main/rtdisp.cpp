////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
void rtdisp(const clmpi&, const clio&, const clbas&, clwfn&);
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  printf("# rtdisp: start    "); 
  system("date");

  clmpi Proc;
  clio IO(argc, argv);
  clcontrol CTRL(IO);
  clbas Bas(Proc, IO);
  clwfn Wfn(Proc, IO, Bas);

  guess(Proc, IO, Bas, Wfn);
  init(Proc, IO, Bas, Wfn);
  rtdisp(Proc, IO, Bas, Wfn);
  Wfn.write(Proc, IO, Bas);

  system("date");
  printf("# rtdisp: end      "); 
}
////////////////////////////////////////////////////////////////////////
void rtdisp(const clmpi& Proc, const clio& IO, const clbas& Bas, clwfn& Wfn)
{
  if (IO.job_type.compare("td") != 0) return;
  clock_t time0 = clock();

  clfield Field(Proc, IO);
  clhprod HPW(Proc, IO, Bas, Field);
  HPW.set_wfn0(Proc, IO, Bas, Wfn);
  cletdrb H12P(Proc, IO, Bas, Field, HPW);

  typhys Phys;
  int nfreq;
  double lfield[9], resp, xfreq, ufreq, test, thresh, dipw, fieldw;

  IO.read_info("rtdisp_thresh", 1.E-15, thresh);
  IO.read_info("rtdisp_nfreq", 10, nfreq);
  double pol1[nfreq];

  Field.step = 0;
  Field.time = ZERO;
  test = thresh*10;
  for (int ifreq = nfreq; ifreq >= 1; ifreq --) {
    pol1[ifreq-1] += ZERO;
  }

  while (! Field.finished()) {
    Field.get_value(lfield);
    HPW.dipole(Proc, Bas, lfield, Wfn, Phys);

    double maxresp = ZERO;
    for (int ifreq = nfreq; ifreq >= 1; ifreq --) {
      xfreq = cos((2*ifreq-1)*PI/(4*nfreq));
      ufreq = xfreq / sqrt(ONE - xfreq*xfreq);

      //resp = -exp(-ufreq*(Field.step-1)*Field.dtime) * Phys.dip[2] / Field.Fields[0].famp;
      //resp = -exp(-ufreq*Field.time) * Phys.dip[2] / Field.Fields[0].famp * Field.dtime;

      dipw = -exp(-ufreq*Field.time) * Phys.dip[2];
      //fieldw = Field.Fields[0].famp * (1.D+0 + exp(-ufreq*Field.dtime)) * HALF;
      fieldw = Field.Fields[0].famp * (1.D+0 + exp(-ufreq*Field.dtime));
      //resp = dipw / fieldw;
      resp = dipw / fieldw * Field.dtime;

      pol1[ifreq-1] += resp;
      if (ifreq == nfreq) test = resp;
    }

    if (IO.nprint_op1 > 0 && Field.step % IO.nprint_op1 == 0) {
      fprintf(IO.fp_op1, "%10d %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n", 
	      Field.step, Field.time, Field.ncyc(), lfield[5], Phys.dip[2], test, pol1[nfreq-1]);
    }

    H12P.prop(Proc, Bas, Field, HPW, Wfn);
    Wfn.mask(Bas);
    Field.step ++;
    Field.time += Field.dtime;
  }

  FILE *fp_alph;
  std::string alph;
  alph = IO.name; alph += ".alph";
  fp_alph = fopen(alph.c_str(), "w");
  for (int ifreq = nfreq; ifreq >= 1; ifreq --) {
    xfreq = cos((2*ifreq-1)*PI/(4*nfreq));
    ufreq = xfreq / sqrt(ONE - xfreq*xfreq);
    fprintf(fp_alph, "%20.10e %20.10e\n", ufreq, pol1[ifreq-1]);
  }
  fclose(fp_alph);
}
////////////////////////////////////////////////////////////////////////
