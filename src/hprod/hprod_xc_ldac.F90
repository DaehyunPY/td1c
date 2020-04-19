      real*8 function piecewise_ldac(o,a,b)
!
!     Fix for a Maple problem. The Maple Fortran generator is not
!     able to transform its piecewise command to Fortran. Therefore
!     this function was written to take care of this.
!
      implicit none
      logical :: o
      real*8 :: a, b
      if (o) then
         piecewise_ldac = a
      else
         piecewise_ldac = b
      endif
      return
      end function piecewise_ldac
!
!    subroutine rks_c_pz81(ideriv,npt,rhoa1,sigmaaa1,zk,vrhoa,vsigmaaa,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
subroutine hprod_xc_ldac(ideriv,npt,rhoa1,zk,vrhoa)
!
      use, intrinsic :: iso_c_binding
!
!     J.P. Perdew, and A. Zunger
!     Self-interaction correction to density-functional approximations 
!     for many-electron systems
!     Phys. Rev. B23 (1981) 5048-5079
!
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
!SATO      real*8 sigmaaa1(npt)
      real*8 zk(npt),vrhoa(npt) !SATO,vsigmaaa(npt)
!SATO      real*8 v2rhoa2(npt),v2rhoasigmaaa(npt),v2sigmaaa2(npt)
      parameter(tol=1.0d-20)
      logical t4
      
      if(ideriv.eq.0) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t5 = t1**(1.D0/6.D0)
      t11 = dlog(t3)
      t17 = piecewise_ldac(1.D0 .le. t3,-0.1423D0/(1.D0 &
    & +0.8292885914166397D0*t5+0.20682485366586D0*t2),0.311D-1*t11 &
    & -0.48D-1+0.12407009817988D-2*t2*t11-0.719606569443304D-2*t2)
      zk(i) = rho*t17
      else ! rho
      zk(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.1) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.8292885914166397D0*t5+0.20682485366586D0*t2
      t11 = dlog(t3)
      t17 = piecewise_ldac(t4,-0.1423D0/t8,0.311D-1*t11-0.48D-1 &
    & +0.12407009817988D-2*t2*t11-0.719606569443304D-2*t2)
      zk(i) = rho*t17
      t18 = piecewise_ldac(t4,0.D0,0.D0)
      t20 = t8**2
      t22 = t5**2
      t23 = t22**2
      t26 = rho**2
      t27 = 1/t26
      t30 = t2**2
      t31 = 1/t30
      t32 = t31*t27
      t45 = piecewise_ldac(t4,0.1423D0/t20*(-0.1382147652361066D0/t23/t5 &
    & *t27-0.6894161788861999D-1*t32),-0.1036666666666667D-1*t1 &
    & -0.4135669939329333D-3*t31*t11*t27-0.4135669939329333D-3*t2*t1 &
    & +0.2398688564811013D-2*t32)
      vrhoa(i) = rho*t18+t17+rho*t45
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
      endif ! rho
      enddo
      
      else if(ideriv.eq.2) then
      
      do i=1,npt
      rho = dmax1(0.D0,rhoa1(i))
      if(rho.gt.tol) then
      t1 = 1/rho
      t2 = t1**(1.D0/3.D0)
      t3 = 0.6203504908994D0*t2
      t4 = 1.D0 .le. t3
      t5 = t1**(1.D0/6.D0)
      t8 = 1.D0+0.8292885914166397D0*t5+0.20682485366586D0*t2
      t10 = 0.1423D0/t8
      t11 = dlog(t3)
      t13 = t2*t11
      t17 = piecewise_ldac(t4,-t10,0.311D-1*t11-0.48D-1+0.12407009817988D-2 &
    & *t13-0.719606569443304D-2*t2)
      zk(i) = rho*t17
      t18 = piecewise_ldac(t4,0.D0,0.D0)
      t19 = rho*t18
      t20 = t8**2
      t21 = 1/t20
      t22 = t5**2
      t23 = t22**2
      t24 = t23*t5
      t25 = 1/t24
      t26 = rho**2
      t27 = 1/t26
      t30 = t2**2
      t31 = 1/t30
      t32 = t31*t27
      t34 = -0.1382147652361066D0*t25*t27-0.6894161788861999D-1*t32
      t38 = t31*t11
      t45 = piecewise_ldac(t4,0.1423D0*t21*t34,-0.1036666666666667D-1*t1 &
    & -0.4135669939329333D-3*t38*t27-0.4135669939329333D-3*t2*t1 &
    & +0.2398688564811013D-2*t32)
      vrhoa(i) = t19+t17+rho*t45
!SATO      t53 = (-0.843D-1/(1.D0+0.1101176160755631D1*t5
!SATO     #+0.1619735131738333D0*t2)+t10)*t27
!SATO      t59 = (-0.1555D-1*t11+0.211D-1-0.80645563816922D-3*t13
!SATO     #+0.421838333811592D-2*t2)*t27
!SATO      t61 = piecewise_ldac(t4,0.1709920934161366D1*t53,0.1709920934161366D1
!SATO     #*t59)
!SATO      t68 = t34**2
!SATO      t73 = t26**2
!SATO      t74 = 1/t73
!SATO      t78 = 1/t26/rho
!SATO      t82 = 1/t30/t1
!SATO      t83 = t82*t74
!SATO      t85 = t31*t78
!SATO      t102 = piecewise_ldac(t4,-0.2846D0/t20/t8*t68+0.1423D0*t21*(
!SATO     #-0.1151789710300888D0/t24/t1*t74+0.2764295304722132D0*t25*t78
!SATO     #-0.4596107859241333D-1*t83+0.13788323577724D0*t85
!SATO     #),0.1036666666666667D-1*t27-0.2757113292886222D-3*t82*t11*t74
!SATO     #-0.4521665800333405D-2*t85+0.8271339878658667D-3*t38*t78
!SATO     #+0.4135669939329333D-3*t2*t27+0.1599125709874009D-2*t83)
!SATO      t107 = piecewise_ldac(t4,-0.1709920934161366D1*t53,
!SATO     #-0.1709920934161366D1*t59)
!SATO      v2rhoa2(i) = rho*t61+4.D0*t18+4.D0*t19+4.D0*t45+2.D0*rho*t102+rho*t107
      else ! rho
      zk(i) = 0.0d0
      vrhoa(i) = 0.0d0
!SATO      v2rhoa2(i) = 0.0d0
      endif ! rho
      enddo
      
      endif ! ideriv
      return
end subroutine hprod_xc_ldac
