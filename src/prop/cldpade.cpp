////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
////////////////////////////////////////////////////////////////////////
cldpade::cldpade()
{
}
////////////////////////////////////////////////////////////////////////
cldpade::~cldpade()
{
}
////////////////////////////////////////////////////////////////////////
cldpade::cldpade(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		 double dtime, dcomplex alpha, int td_type, int split_type)
	     
{
  gen(MPIP, IO, Bas, dtime, alpha, td_type, split_type);
}
////////////////////////////////////////////////////////////////////////
void cldpade::gen(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
		  double dtime, dcomplex alpha, int td_type, int split_type)
{
  dpade_alpha = alpha;
  dpade_dtime = dtime;

  IO.read_info("h1itr_maxcyc", 20, dpade_maxcyc);
  IO.read_info("h1itr_thresh", (double) 1.0E-10, dpade_thresh);
  IO.read_info("h1itr_midpt", true, dpade_midpt);
  tpiv.resize((Bas.GRad.nrad - 1) * (Bas.GAng.lmax1 + 1));
  tinv.resize((3 * Bas.GRad.nmax + 1) * (Bas.GRad.nrad - 1) * (Bas.GAng.lmax1 + 1));
  dpade_gen_(&td_type, &dtime, &alpha, &tinv[0], &tpiv[0]);
}
////////////////////////////////////////////////////////////////////////
//void cldpade::gen2(const clmpi& MPIP, const clio& IO, const clbas& Bas, 
//		   dcomplex d0, dcomplex d1)
//{
//  dpade_coeff0 = d0;
//  dpade_coeff1 = d1;
//
//  IO.read_info("h1itr_maxcyc", 20, dpade_maxcyc);
//  IO.read_info("h1itr_thresh", (double) 1.0E-10, dpade_thresh);
//  IO.read_info("h1itr_midpt", true, dpade_midpt);
//
//  tpiv.resize((Bas.GRad.nrad - 1) * (Bas.GAng.lmax1 + 1));
//  tinv.resize((3 * Bas.GRad.nmax + 1) * (Bas.GRad.nrad - 1) * (Bas.GAng.lmax1 + 1));
//  dpade_gen2_(&dpade_coeff0, &dpade_coeff1, &tinv[0], &tpiv[0]);
//}
//////////////////////////////////////////////////////////////////////////
void cldpade::prod(const clmpi& Proc, const clbas& Bas, double time,
		   const clfield& Field, clhprod& HPW, clwfn& Wfn) const
{
  if (clcontrol::iprojfc !=2 && Field.vgauge) {
    std::cout << "cldpade::prod: projfc for vgauge is nyi." << std::endl;
    abort();
  }

  double lfield[9];
  Field.get_value(time, lfield);

  dpade_prod_(&clcontrol::split_type, &Field.igauge, &Field.td_type, &LZERO, lfield, &dpade_dtime,
	      &dpade_alpha, &dpade_maxcyc, &dpade_thresh, &tpiv[0], &tinv[0], &Wfn.wfn[0]);
	      
}
//////////////////////////////////////////////////////////////////////////
//void cldpade::prod2(const clmpi& Proc, const clbas& Bas, double time,
//		    const clfield& Field, clhprod& HPW, clwfn& Wfn) const
//{
//  if (clcontrol::iprojfc != 2) {
//    std::cout << "cldpade::prod2: projfc is nyi." << std::endl;
//    abort();
//  }
//
//  double lfield[9];
//  Field.get_value(time, lfield);
//
//  dpade_prod2_(&clcontrol::split_type, &Field.igauge, lfield, &dpade_coeff1,
//	       &dpade_maxcyc, &dpade_thresh, &tpiv[0], &tinv[0], &Wfn.wfn[0]);
//}
//////////////////////////////////////////////////////////////////////////
void cldpade::cnic(const clmpi& Proc, const clbas& Bas, double time,
		   const clfield& Field, clhprod& HPW, clwfn& Wfn) const
{
  if (clcontrol::iprojfc != 2 && Field.vgauge) {
    std::cout << "cldpade::prod: projfc for vgauge is nyi." << std::endl;
    abort();
  }

  double lfield[9];
  Field.get_value(time, lfield);

  dpade_cnic_(&clcontrol::split_type, &Field.igauge, &Field.td_type, &LZERO, lfield, &dpade_dtime,
	      &dpade_alpha, &dpade_maxcyc, &dpade_thresh, &tpiv[0], &tinv[0], &Wfn.wfn[0]);	      
}
//////////////////////////////////////////////////////////////////////////
