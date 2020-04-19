!######################################################################
subroutine ccnew_mktcc1(cc1, tcc1)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_ccnew, only : norb1

  implicit none
  complex(kind(0d0)), intent(in) :: cc1((norb1+1):nact, 1:norb1, 1:*)
  complex(kind(0d0)), intent(out) :: tcc1((norb1+1):nact, 1:norb1, 1:*)

  integer(c_int) :: h1, p2

  do h1 = 1, norb1
  do p2 = norb1 + 1, nact
     tcc1(p2, h1, 1) = + cc1(p2, h1, 1)
     tcc1(p2, h1, 2) = + cc1(p2, h1, 1)
  end do
  end do

end subroutine ccnew_mktcc1
!######################################################################
