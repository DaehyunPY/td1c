////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  clmpi Proc;
  printf("# main: start    "); print_date();
  if ( argc < 2 ) {
    printf( "syntax: proj INP1 INP2" );
    exit( 1 );
  }

  std::string inp1 = argv[1];
  std::string inp2 = argv[2];

  clio IO1(inp1);
  clbas Bas1(Proc, IO1);
  clwfn Wfn1(Proc, IO1, Bas1);
  guess(Proc, IO1, Bas1, Wfn1);

  clio IO2(inp2);
  clbas Bas2(Proc, IO2);
  clwfn Wfn2(Proc, IO2, Bas2);
  Bas1.proj(Bas2, Wfn1.wfn, Wfn2.wfn);
  Wfn2.print(Bas2);

  printf("# main: end      "); print_date();
}
////////////////////////////////////////////////////////////////////////
