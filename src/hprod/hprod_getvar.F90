!#######################################################################
subroutine hprod_getden1(output)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_hprod, only : den1

  implicit none
  complex(c_double_complex), intent(out) :: output(1:nact,1:nact)
  output(:,:) = den1(:,:)

end subroutine hprod_getden1
!#######################################################################
subroutine hprod_getint1(output)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_hprod, only : int1e

  implicit none
  complex(c_double_complex), intent(out) :: output(1:nact,1:nact)
  output(:,:) = int1e(:,:)

end subroutine hprod_getint1
!#######################################################################
subroutine hprod_getint2(output)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_hprod, only : int2e

  implicit none
  complex(c_double_complex), intent(out) :: output(1:nact,1:nact,1:nact,1:nact)
  output(:,:,:,:) = int2e(:,:,:,:)

end subroutine hprod_getint2
!#######################################################################
