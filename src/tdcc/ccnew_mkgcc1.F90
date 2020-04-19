!######################################################################
subroutine ccnew_mkgcc1(cc1, gcc1)

  use, intrinsic :: iso_c_binding
  use omp_mod
  use wfn_mod, only : nact, norb1

  implicit none
  complex(kind(0d0)), intent(in) :: cc1(1:norb1, (norb1+1):nact, 1:*)
  complex(kind(0d0)), intent(out) :: gcc1(1:norb1, (norb1+1):nact, 1:*)

  integer(c_int) :: h1, p2

  do h1 = 1, norb1
  do p2 = norb1 + 1, nact
     gcc1(h1, p2, 1) = + cc1(h1, p2, 1)
     gcc1(h1, p2, 2) = + cc1(h1, p2, 1)
  end do
  end do

end subroutine ccnew_mkgcc1
!######################################################################
