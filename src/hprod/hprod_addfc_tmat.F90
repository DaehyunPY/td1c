!######################################################################
subroutine hprod_addfc_tmat()
!
  use, intrinsic :: iso_c_binding
  use mod_bas, only : tmat
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, ndvr
  use mod_const, only : pi
  use mod_hprod, only : v2jfc
  use mod_ormas, only : nfcore
  use mod_control, only : jfc_implicit, exact3j

  implicit none
  integer(c_int) :: irad, l

  if (nfcore == 0) return

  if (jfc_implicit) then
     if (.not. exact3j) then
        do l = 0, lmax1
           do irad = 1, nrad - 1
              tmat(1+ndvr, irad, l) = tmat(1+ndvr, irad, l) + v2jfc(irad)
           end do
        end do
     else
        do l = 0, lmax1
           do irad = 1, nrad - 1
              tmat(1+ndvr, irad, l) = tmat(1+ndvr, irad, l) + v2jfc(irad) / (2d0*pi)
           end do
        end do
     end if
  end if

end subroutine hprod_addfc_tmat
!######################################################################
