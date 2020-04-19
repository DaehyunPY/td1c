!################################################################################
module mod_h1rat

  use, intrinsic :: iso_c_binding

  implicit none

  real(c_double) :: h1rat_eref
  complex(c_double_complex), allocatable :: h1rat_int1e(:,:)
  complex(c_double_complex), allocatable :: h1rat_int2e(:,:,:,:)

end module mod_h1rat
!################################################################################
