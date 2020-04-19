!######################################################################
! test FC use nfcore_tdcis instead of nfcore
!######################################################################
subroutine hprod_mfprod2_tdcis(phi, chi, v2, v20, g1wfn)

  ! dynamical-core-eXchange
  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nfcore, ndcore, ncore, nfun
  use mod_const, only : zero, two, ctwo, pi
  use mod_hprod, only : tdcis_eig
  use mod_ormas, only : nfcore_tdcis

  implicit none
  complex(c_double_complex), intent(in) :: phi(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: chi(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v20(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfn(1:ngrid, 1:nfun)
  integer(c_int) :: ifun, jfun, igrid, llg, ulg

  !$omp parallel default(shared) private(llg, ulg, ifun, jfun)
  call util_omp_disp(1, ngrid, llg, ulg)
  ! do ifun = nfcore + 1, nfun
  do ifun = nfcore_tdcis + 1, nfun
     do igrid = llg, ulg
        g1wfn(igrid, ifun) = g1wfn(igrid, ifun) - tdcis_eig(ifun)*chi(igrid,ifun)
     end do
     ! acting on frozen core
     ! do jfun = 1, nfcore
     do jfun = 1, nfcore_tdcis
        do igrid = llg, ulg
           g1wfn(igrid, ifun) = g1wfn(igrid, ifun) &
                + v20(igrid, jfun, jfun)*chi(igrid, ifun)*2d0 &
                - v2(igrid, jfun, ifun)*phi(igrid, jfun) 
        end do
     end do
     ! acting on active
     ! do jfun = nfcore + 1, nfun
     do jfun = nfcore_tdcis + 1, nfun
        do igrid = llg, ulg
           g1wfn(igrid, ifun) = g1wfn(igrid, ifun) &
                 + v20(igrid, jfun, jfun)*chi(igrid, ifun)*2d0 &
                 - v20(igrid, jfun, ifun)*chi(igrid, jfun) &
                 + v2(igrid, jfun, jfun)*phi(igrid, ifun)*2d0 &
                 - v2(igrid, jfun, ifun)*phi(igrid, jfun) 
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod2_tdcis
!######################################################################
 
