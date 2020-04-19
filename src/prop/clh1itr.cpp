////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
clh1itr::clh1itr()
{
  std::cout << "clh1itr" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clh1itr::~clh1itr()
{
  std::cout << "~clh1itr" << std::endl;
}
////////////////////////////////////////////////////////////////////////
clh1itr::clh1itr(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		 const clfield& Field, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, HPW);
}
////////////////////////////////////////////////////////////////////////
clh1itr::clh1itr(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		 const clfield& Field, double dtime, const clhprod& HPW)
{
  gen(MPIP, IO, Bas, Field, dtime, HPW);
}
////////////////////////////////////////////////////////////////////////
void clh1itr::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		  const clfield& Field, const clhprod& HPW)
{
  std::cout << "clh1itr::gen" << std::endl;
  clh1prop::gen_basic(MPIP, IO, Bas, Field, HPW);

  IO.read_info("h1itr_maxcyc", (long) 20, h1itr_maxcyc);
  IO.read_info("h1itr_thresh", 1.0E-10, h1itr_thresh);

  cnpiv.resize((Bas.GRad.nrad - 1) *  (Bas.GAng.lmax1 + 1));
  cninv.resize((3 * Bas.GRad.nmax + 1) *  (Bas.GRad.nrad - 1) *  (Bas.GAng.lmax1 + 1));

  double dt2 = HALF * Field.dtime;
  h1itr_gen_cninv_(&Field.td_type, &dt2, &cnpiv[0], &cninv[0]);
}
////////////////////////////////////////////////////////////////////////
void clh1itr::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		  const clfield& Field, double dtime, const clhprod& HPW)
{
  std::cout << "clh1itr::gen" << std::endl;
  clh1prop::gen_basic(MPIP, IO, Bas, Field, HPW);

  IO.read_info("h1itr_maxcyc", h1itr_maxcyc);
  IO.read_info("h1itr_thresh", h1itr_thresh);

  cnpiv.resize((Bas.GRad.nrad - 1) *  (Bas.GAng.lmax1 + 1));
  cninv.resize((3 * Bas.GRad.nmax + 1) *  (Bas.GRad.nrad - 1) *  (Bas.GAng.lmax1 + 1));

  double dt2 = HALF * dtime;
  h1itr_gen_cninv_(&Field.td_type, &dt2, &cnpiv[0], &cninv[0]);
}
////////////////////////////////////////////////////////////////////////
void clh1itr::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field, 
		   long STEP, clhprod& HPW, clwfn& Wfn)
{
  if (clcontrol::split_type == 0) return;

  double time;
  double lfield[9];

  if (STEP == 1) {
    time = Field.time;
  } else {
    time = Field.time + Field.dtime * HALF;
  }

  long iprojfc = 0;
  if (projfc) iprojfc = 1;

//  if (projfc && Field.vgauge) {
//    std::cout << "clh1itr::prop: projfc is nyi for velocity gauge" << std::endl;
//    abort();
//  }

  double dt2 = HALF * Field.dtime;
  Field.get_value(time + dt2 * HALF, lfield);  
  h1itr_prop_(&Field.td_type, &clcontrol::split_type, &Field.igauge, &iprojfc, 
	      lfield, &dt2, &h1itr_maxcyc, &h1itr_thresh, &cnpiv[0], 
	      &cninv[0], &Wfn.wfn[0]);
//DONOT WORK  double lfield1[9];
//DONOT WORK  double lfield2[9];
//DONOT WORK  Field.get_value(time, lfield1);  
//DONOT WORK  Field.get_value(time + dt2, lfield2);
//DONOT WORK  h1itr_prop2_(&Field.td_type, &clcontrol::split_type, &Field.igauge, lfield1, lfield2,
//DONOT WORK	       &dt2, &h1itr_maxcyc, &h1itr_thresh, &cnpiv[0], &cninv[0], &Wfn.wfn[0]);
}
////////////////////////////////////////////////////////////////////////
void clh1itr::prop(const clmpi& Proc, const clbas& Bas, const clfield& Field,
		   double time, double dtime, clhprod& HPW, clwfn& Wfn)
{
  if (clcontrol::split_type == 0) return;

  double lfield[9];

  long iprojfc = 0;
  if (projfc) iprojfc = 1;

//  if (projfc && Field.vgauge) {
//    std::cout << "clh1itr::prop: projfc is nyi for velocity gauge" << std::endl;
//    abort();
//  }

  double dt2 = HALF * dtime;
  Field.get_value(time + dt2 * HALF, lfield);  
  h1itr_prop_(&Field.td_type, &clcontrol::split_type, &Field.igauge, &iprojfc, 
	      lfield, &dt2, &h1itr_maxcyc, &h1itr_thresh, &cnpiv[0], 
	      &cninv[0], &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
