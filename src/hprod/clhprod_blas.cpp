////////////////////////////////////////////////////////////////////////
// Hamiltonian
////////////////////////////////////////////////////////////////////////
#include "td1c.hpp"
#include "wrapper.hpp"
//////////////////////////////////////////////////////////////////////////
void clhprod::clear(const clmpi& MPIP, const clbas& Bas, clwfn& Wfn)
{
  zclear_omp_(&Wfn.size, &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::clearo(const clmpi& MPIP, const clbas& Bas, clwfn& Wfn)
{
  zclear_omp_(&Wfn.size1, &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::clearc(const clmpi& MPIP, const clbas& Bas, clwfn& Wfn)
{
  zclear_omp_(&Wfn.size2, &Wfn.wfn[Wfn.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::clearog(const clmpi& MPIP, const clbas& Bas, clwfn& Wfn)
{
  zclear_omp_(&Wfn.sizeg, &Wfn.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void clhprod::copy(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  zcopy_omp_(&Wfn1.size, &Wfn1.wfn[0], &Wfn2.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::copyo(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  zcopy_omp_(&Wfn1.size1, &Wfn1.wfn[0], &Wfn2.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::copyc(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  zcopy_omp_(&Wfn1.size2, &Wfn1.wfn[Wfn1.size1], &Wfn2.wfn[Wfn2.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::copyog(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  zcopy_omp_(&Wfn1.sizeg, &Wfn1.wfng[0], &Wfn2.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void clhprod::scal(const clmpi& MPIP, const clbas& Bas, dcomplex fac, clwfn& Wfn)
{
  zscal_omp_(&Wfn.size, &fac, &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::scalo(const clmpi& MPIP, const clbas& Bas, dcomplex fac, clwfn& Wfn)
{
  zscal_omp_(&Wfn.size1, &fac, &Wfn.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::scalc(const clmpi& MPIP, const clbas& Bas, dcomplex fac, clwfn& Wfn)
{
  zscal_omp_(&Wfn.size2, &fac, &Wfn.wfn[Wfn.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::scalog(const clmpi& MPIP, const clbas& Bas, dcomplex fac, clwfn& Wfn)
{
  zscal_omp_(&Wfn.sizeg, &fac, &Wfn.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void clhprod::xpy(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  zxpy_omp_(&Wfn1.size, &Wfn1.wfn[0], &Wfn2.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xpyo(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  zxpy_omp_(&Wfn1.size1, &Wfn1.wfn[0], &Wfn2.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xpyog(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  zxpy_omp_(&Wfn1.sizeg, &Wfn1.wfng[0], &Wfn2.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xpyc(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  zxpy_omp_(&Wfn1.size2, &Wfn1.wfn[Wfn1.size1], &Wfn2.wfn[Wfn2.size1]);
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void clhprod::xpyz(const clmpi& MPIP, const clbas& Bas, 
		   const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  zxpyz_omp_(&Wfn1.size, &Wfn1.wfn[0], &Wfn2.wfn[0], &Wfn3.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xpyzo(const clmpi& MPIP, const clbas& Bas, 
		    const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  zxpyz_omp_(&Wfn1.size1, &Wfn1.wfn[0], &Wfn2.wfn[0], &Wfn3.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xpyzc(const clmpi& MPIP, const clbas& Bas,
		    const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  zxpyz_omp_(&Wfn1.size2, &Wfn1.wfn[Wfn1.size1], &Wfn2.wfn[Wfn2.size1], &Wfn3.wfn[Wfn3.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xpyzog(const clmpi& MPIP, const clbas& Bas, 
		     const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  zxpyz_omp_(&Wfn1.sizeg, &Wfn1.wfng[0], &Wfn2.wfng[0], &Wfn3.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void clhprod::xmy(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  dcomplex fac = - RUNIT;
  zaxpy_omp_(&Wfn1.size, &fac, &Wfn1.wfn[0], &Wfn2.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xmyo(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  dcomplex fac = - RUNIT;
  zaxpy_omp_(&Wfn1.size1, &fac, &Wfn1.wfn[0], &Wfn2.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xmyc(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  dcomplex fac = - RUNIT;
  zaxpy_omp_(&Wfn1.size2, &fac, &Wfn1.wfn[Wfn1.size1], &Wfn2.wfn[Wfn2.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xmyog(const clmpi& MPIP, const clbas& Bas, const clwfn& Wfn1, clwfn& Wfn2)
{
  dcomplex fac = - RUNIT;
  zaxpy_omp_(&Wfn1.sizeg, &fac, &Wfn1.wfng[0], &Wfn2.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void clhprod::xmyz(const clmpi& MPIP, const clbas& Bas, 
		   const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  dcomplex fac = - RUNIT;
  //  zaxpyz_omp_(&Wfn1.size, &fac, &Wfn1.wfn[0], &Wfn2.wfn[0], &Wfn3.wfn[0]);
  zaxpyz_omp_(&Wfn1.size, &fac, &Wfn2.wfn[0], &Wfn1.wfn[0], &Wfn3.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xmyzo(const clmpi& MPIP, const clbas& Bas, 
		    const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  dcomplex fac = - RUNIT;
  //  zaxpyz_omp_(&Wfn1.size1, &fac, &Wfn1.wfn[0], &Wfn2.wfn[0], &Wfn3.wfn[0]);
  zaxpyz_omp_(&Wfn1.size1, &fac, &Wfn2.wfn[0], &Wfn1.wfn[0], &Wfn3.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xmyzc(const clmpi& MPIP, const clbas& Bas, 
		    const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  dcomplex fac = - RUNIT;
  //  zaxpyz_omp_(&Wfn1.size2, &fac, &Wfn1.wfn[Wfn1.size1], &Wfn2.wfn[Wfn2.size1], &Wfn3.wfn[Wfn3.size1]);
  zaxpyz_omp_(&Wfn1.size2, &fac, &Wfn2.wfn[Wfn2.size1], &Wfn1.wfn[Wfn1.size1], &Wfn3.wfn[Wfn3.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::xmyzog(const clmpi& MPIP, const clbas& Bas, 
		     const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  dcomplex fac = - RUNIT;
  //  zaxpyz_omp_(&Wfn1.sizeg, &fac, &Wfn1.wfng[0], &Wfn2.wfng[0], &Wfn3.wfng[0]);
  zaxpyz_omp_(&Wfn1.sizeg, &fac, &Wfn2.wfng[0], &Wfn1.wfng[0], &Wfn3.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void clhprod::axpy(const clmpi& MPIP, const clbas& Bas, dcomplex fac, 
		   const clwfn& Wfn1, clwfn& Wfn2)
{
  zaxpy_omp_(&Wfn1.size, &fac, &Wfn1.wfn[0], &Wfn2.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::axpyo(const clmpi& MPIP, const clbas& Bas, dcomplex fac, 
		    const clwfn& Wfn1, clwfn& Wfn2)
{
  zaxpy_omp_(&Wfn1.size1, &fac, &Wfn1.wfn[0], &Wfn2.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::axpyc(const clmpi& MPIP, const clbas& Bas, dcomplex fac, 
		    const clwfn& Wfn1, clwfn& Wfn2)
{
  zaxpy_omp_(&Wfn1.size2, &fac, &Wfn1.wfn[Wfn1.size1], &Wfn2.wfn[Wfn2.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::axpyog(const clmpi& MPIP, const clbas& Bas, dcomplex fac, 
		     const clwfn& Wfn1, clwfn& Wfn2)
{
  zaxpy_omp_(&Wfn1.sizeg, &fac, &Wfn1.wfng[0], &Wfn2.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void clhprod::axpyz(const clmpi& MPIP, const clbas& Bas, dcomplex fac, 
		    const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  zaxpyz_omp_(&Wfn1.size, &fac, &Wfn1.wfn[0], &Wfn2.wfn[0], &Wfn3.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::axpyzo(const clmpi& MPIP, const clbas& Bas, dcomplex fac, 
		     const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  zaxpyz_omp_(&Wfn1.size1, &fac, &Wfn1.wfn[0], &Wfn2.wfn[0], &Wfn3.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::axpyzc(const clmpi& MPIP, const clbas& Bas, dcomplex fac, 
		    const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  zaxpyz_omp_(&Wfn1.size2, &fac, &Wfn1.wfn[Wfn1.size1], &Wfn2.wfn[Wfn2.size1], &Wfn3.wfn[Wfn3.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::axpyzog(const clmpi& MPIP, const clbas& Bas, dcomplex fac, 
		      const clwfn& Wfn1, const clwfn& Wfn2, clwfn& Wfn3)
{
  zaxpyz_omp_(&Wfn1.sizeg, &fac, &Wfn1.wfng[0], &Wfn2.wfng[0], &Wfn3.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
void clhprod::axpbyz(const clmpi& MPIP, const clbas& Bas, 
		     dcomplex fac1, const clwfn& Wfn1, 
		     dcomplex fac2, const clwfn& Wfn2, 
		     clwfn& Wfn3)
{
  zaxpbyz_omp_(&Wfn1.size, &fac1, &Wfn1.wfn[0], &fac2, &Wfn2.wfn[0], &Wfn3.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::axpbyzo(const clmpi& MPIP, const clbas& Bas, 
		      dcomplex fac1, const clwfn& Wfn1, 
		      dcomplex fac2, const clwfn& Wfn2, 
		      clwfn& Wfn3)
{
  zaxpbyz_omp_(&Wfn1.size1, &fac1, &Wfn1.wfn[0], &fac2, &Wfn2.wfn[0], &Wfn3.wfn[0]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::axpbyzc(const clmpi& MPIP, const clbas& Bas, 
		      dcomplex fac1, const clwfn& Wfn1, 
		      dcomplex fac2, const clwfn& Wfn2, 
		      clwfn& Wfn3)
{
  zaxpbyz_omp_(&Wfn1.size2, &fac1, &Wfn1.wfn[Wfn1.size1], &fac2, &Wfn2.wfn[Wfn2.size1], &Wfn3.wfn[Wfn3.size1]);
}
//////////////////////////////////////////////////////////////////////////
void clhprod::axpbyzog(const clmpi& MPIP, const clbas& Bas, 
		       dcomplex fac1, const clwfn& Wfn1, 
		       dcomplex fac2, const clwfn& Wfn2, 
		       clwfn& Wfn3)
{
  zaxpbyz_omp_(&Wfn1.sizeg, &fac1, &Wfn1.wfng[0], &fac2, &Wfn2.wfng[0], &Wfn3.wfng[0]);
}
//////////////////////////////////////////////////////////////////////////
