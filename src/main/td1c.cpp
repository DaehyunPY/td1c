////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  printf("# td1c: start    \n"); 
  system("date");

  clmpi Proc;
  clio IO(argc, argv);
  clcontrol CTRL(IO);
  clbas Bas(Proc, IO);
  clwfn Wfn(Proc, IO, Bas, true);

  guess(Proc, IO, Bas, Wfn);
  init(Proc, IO, Bas, Wfn);
  tdse(Proc, IO, Bas, Wfn);

  system("date");
  printf("# td1c: end      \n"); 
}
////////////////////////////////////////////////////////////////////////
