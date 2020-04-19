!######################################################################
subroutine ccnew_mkgcc2(cc2, gcc2)

  use, intrinsic :: iso_c_binding
  use omp_mod
  use wfn_mod, only : nact, norb1

  implicit none
  complex(kind(0d0)), intent(in) :: cc2(1:norb1, 1:norb1, (norb1+1):nact, (norb1+1):nact, 1:*)
  complex(kind(0d0)), intent(out) :: gcc2(1:norb1, 1:norb1, (norb1+1):nact, (norb1+1):nact, 1:*)

  integer(c_int) :: h1, h2, p3, p4

  gcc2(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:6)=0d0
  do h1 = 1, norb1
  do h2 = 1, norb1
  do p3 = norb1 + 1, nact
  do p4 = norb1 + 1, nact
     if (h1.ne.h2 .and. p3.ne.p4) then
        gcc2(h1, h2, p3, p4, 1) = + cc2(h1, h2, p3, p4, 1)  ! 1 (h+h+|p+p+)
        gcc2(h1, h2, p3, p4, 2) = + cc2(h1, h2, p3, p4, 1)  ! 2 (h-h-|p-p-)
     end if
     gcc2(h1, h2, p3, p4, 3) = + cc2(h1, h2, p3, p4, 2)  ! 3 (h+h-|p+p-)
     gcc2(h1, h2, p3, p4, 4) = + cc2(h2, h1, p4, p3, 2)  ! 4 (h-h+|p-p+)
     gcc2(h1, h2, p3, p4, 5) = - cc2(h1, h2, p4, p3, 2)  ! 5 (h+h-|p-p+)
     gcc2(h1, h2, p3, p4, 6) = - cc2(h2, h1, p3, p4, 2)  ! 6 (h-h+|p+p-)
  end do
  end do
  end do
  end do

end subroutine ccnew_mkgcc2
!######################################################################
