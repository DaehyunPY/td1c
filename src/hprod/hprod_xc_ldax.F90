!     subroutine hprod_xc_ldax (ideriv,npt,rhoa1,sigmaaa1,zk,vrhoa,vsigmaaa,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
      subroutine hprod_xc_ldax (ideriv,npt,rhoa1,zk,vrhoa)

      use, intrinsic :: iso_c_binding
!
!     P.A.M. Dirac
!     Proceedings of the Cambridge Philosophical Society, 26 (1930) 376
!
!     CITATION:
!
!     Functionals were obtained from the Density Functional Repository 
!     as developed and distributed by the Quantum Chemistry Group, 
!     CCLRC Daresbury Laboratory, Daresbury, Cheshire, WA4 4AD 
!     United Kingdom. Contact Huub van Dam (h.j.j.vandam@dl.ac.uk) or 
!     Paul Sherwood for further information.
!
!     COPYRIGHT:
!
!     Users may incorporate the source code into software packages and
!     redistribute the source code provided the source code is not
!     changed in anyway and is properly cited in any documentation or
!     publication related to its use.
!
!     ACKNOWLEDGEMENT:
!
!     The source code was generated using Maple 8 through a modified
!     version of the dfauto script published in:
!
!        R. Strange, F.R. Manby, P.J. Knowles
!        Automatic code generation in density functional theory
!        Comp. Phys. Comm. 136 (2001) 310-318.
!
      implicit real*8 (a-h,o-z)
      integer(c_int) ideriv,npt
      real*8 rhoa1(npt)
!SATO real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt) !SATO,vsigmaaa(npt)
!SATO real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = rho**(1.D0/3.D0)
      zk(i) = -0.7385587663820224D0*t1*rho
      vrhoa(i) = -0.9847450218426965D0*t1
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
!      do i=1,npt
!      rho = dmax1(0.D0,rhoa1(i))
!      if(rho.gt.tol) then
!      t1 = rho**(1.D0/3.D0)
!      zk(i) = -0.7385587663820224D0*t1*rho
!      vrhoa(i) = -0.9847450218426965D0*t1
!      t5 = t1**2
!      v2rhoa2(i) = -0.6564966812284644D0/t5
!      else ! rho
!      zk(i) = 0.0d0
!      vrhoa(i) = 0.0d0
!      v2rhoa2(i) = 0.0d0
!      endif ! rho
!      enddo
      
      endif ! ideriv
      return
      end subroutine hprod_xc_ldax
