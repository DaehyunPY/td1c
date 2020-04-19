////////////////////////////////////////////////////////////////////////
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <iostream>
extern "C"
{
  void test_zaxpy_omp_();
}
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  printf("# td1c: start    \n"); 
  system("date");

  test_zaxpy_omp_();

//  clmpi Proc;
//  clio IO(argc, argv);
//  clcontrol CTRL(IO);
//  clbas Bas(Proc, IO);
//// Sato_tSURFF
//  //clwfn Wfn(Proc, IO, Bas);
//  clwfn Wfn(Proc, IO, Bas, true);
//// Sato_tSURFF
//
//  guess(Proc, IO, Bas, Wfn);
//  init(Proc, IO, Bas, Wfn);
//  tdse(Proc, IO, Bas, Wfn);

  system("date");
  printf("# td1c: end      \n"); 
}
////////////////////////////////////////////////////////////////////////
