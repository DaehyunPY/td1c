////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  clmpi Proc;
  clio IO(argc, argv);
  clcontrol CTRL(IO);
  clfield Field(Proc, IO);
  double lfield[9], dfield[9];
  Field.step = ZERO;
  Field.time = ZERO;
  while (! Field.finished()) {
    Field.get_value(lfield);
    Field.get_der(dfield);
    printf("%10d %20.10e %20.10e %20.10e %20.10e %20.10e\n",
	   Field.step, Field.time, Field.ncyc(), lfield[5], lfield[8], dfield[2]);

    Field.step ++;
    Field.time += Field.dtime;
  }
}
////////////////////////////////////////////////////////////////////////
