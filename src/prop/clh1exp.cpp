////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
clh1exp::clh1exp()
{
  std::cout << "clh1exp" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clh1exp::~clh1exp()
{
  std::cout << "~clh1exp" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clh1exp::clh1exp(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		 const clfield& Field, const clhprod& HPW)	     
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
void clh1exp::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		  const clfield& Field, const clhprod& HPW)
{
  std::cout << "clh1itr::gen" << std::endl;
  clh1prop::gen_basic(MPIP, IO, Bas, Field, HPW);

  if (projfc) {
    std::cout << "clh1exp::gen: projfc is nyi" << std::endl;
    abort();
  }

  IO.read_info("h1exp_maxcyc", h1exp_maxcyc);
  IO.read_info("h1exp_thresh", h1exp_thresh);
}
////////////////////////////////////////////////////////////////////////
void clh1exp::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		  const clfield& Field, double dtime, const clhprod& HPW)
{
  std::cout << "clh1exp::gen_2 nyi." << std::endl;
  abort();
}
////////////////////////////////////////////////////////////////////////
void clh1exp::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field, 
		   long STEP, clhprod& HPW, clwfn& Wfn)
{
  double time;
  double lfield[9];

  if (STEP == 1) {
    time = Field.time;
  } else {
    time = Field.time + Field.dtime * HALF;
  }

  double dt2 = HALF * Field.dtime;
  Field.get_value(time + dt2 * HALF, lfield);  
  h1exp_prop_(&Field.td_type, &clcontrol::split_type, &Field.igauge, lfield, 
	      &dt2, &h1exp_maxcyc, &h1exp_thresh, &Wfn.wfn[0]);
//  h1exp_prop1_(&Field.td_type, &HPW.split_type, &Field.igauge, lfield, 
//	       &dt2, &h1exp_maxcyc, &Wfn.wfn[0]);
}
////////////////////////////////////////////////////////////////////////
void clh1exp::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  std::cout << "clh1exp::prop_2 nyi." << std::endl;
  abort();
}
//////////////////////////////////////////////////////////////////////////
