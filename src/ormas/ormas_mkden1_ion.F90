!######################################################################
subroutine ormas_mkden1_ion(smat, cic, den1)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit
  use mod_ormas, only : nact
  use mod_ormas, only : nstr_alph
  use mod_ormas, only : nstr_beta

  implicit none
  complex(c_double_complex), intent(inout) :: smat(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(out) :: den1(1:nact, 1:nact)

  integer(c_long) :: iact, jact, kact, lact
  complex(c_double_complex), allocatable :: den2(:,:,:,:)
  allocate(den2(1:nact, 1:nact, 1:nact, 1:nact))

  if (nact == 0) return
  call ormas_mkden1(cic, den1)
  call ormas_mkden2(cic, den1, den2)

  smat(1:nact, 1:nact) = -smat(1:nact, 1:nact)
  do iact = 1, nact
     smat(iact, iact) = smat(iact, iact) + runit
  end do

  do iact = 1, nact
     do jact = 1, nact
        den1(iact, jact) = czero
        do kact = 1, nact
           do lact = 1, nact
              den1(iact, jact) = den1(iact, jact) + den2(iact, jact, kact, lact) * smat(lact, kact)
           end do
        end do
     end do
  end do

  deallocate(den2)

end subroutine ormas_mkden1_ion
!######################################################################
subroutine ormas_mkden1_smat(smat, cic, den1)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit
  use mod_ormas, only : nact
  use mod_ormas, only : nstr_alph
  use mod_ormas, only : nstr_beta

  implicit none
  complex(c_double_complex), intent(in) :: smat(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(out) :: den1(1:nact, 1:nact)

  integer(c_long) :: iact, jact, kact, lact
  complex(c_double_complex), allocatable :: den2(:,:,:,:)
  allocate(den2(1:nact, 1:nact, 1:nact, 1:nact))

  if (nact == 0) return
  call ormas_mkden1(cic, den1)
  call ormas_mkden2(cic, den1, den2)

  do iact = 1, nact
     do jact = 1, nact
        den1(iact, jact) = czero
        do kact = 1, nact
           do lact = 1, nact
              den1(iact, jact) = den1(iact, jact) + den2(iact, jact, kact, lact) * smat(lact, kact)
           end do
        end do
     end do
  end do

  deallocate(den2)

end subroutine ormas_mkden1_smat
!######################################################################
