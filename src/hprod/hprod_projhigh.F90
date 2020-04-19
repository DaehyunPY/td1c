!######################################################################
subroutine hprod_projhigh(wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : projhigh,projhigh_orbs,projhigh_ncut

  implicit none
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex) :: tmp
  integer(c_long) :: lll,ull,ifun,l,jcut,irad

  if (.not. projhigh) return

  !$omp parallel default(shared) private(lll,ull,tmp)
  !###########################
  call util_omp_disp(0, lmax1, lll, ull)
  do ifun = nfcore + 1, nfun
     do l = lll, ull
        if (l >= abs(mval(ifun))) then
           do jcut = 1, projhigh_ncut(l)
              tmp= 0d0
              do irad = 1, nrad - 1
                 tmp = tmp + projhigh_orbs(irad, jcut, l) * wfn(irad, l, ifun)
              end do
              do irad = 1, nrad - 1
                 wfn(irad, l, ifun) = wfn(irad, l, ifun) - projhigh_orbs(irad, jcut, l)*tmp
              end do
           end do
        end if
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_projhigh
!######################################################################
