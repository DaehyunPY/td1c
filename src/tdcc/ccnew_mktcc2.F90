!######################################################################
subroutine ccnew_mktcc2(cc2, tcc2)

  use, intrinsic :: iso_c_binding
  use omp_mod
  use wfn_mod, only : nact, norb1

  implicit none
  complex(kind(0d0)), intent(in) :: cc2((norb1+1):nact, (norb1+1):nact, 1:norb1, 1:norb1, 1:*)
  complex(kind(0d0)), intent(out) :: tcc2((norb1+1):nact, (norb1+1):nact, 1:norb1, 1:norb1, 1:*)

  integer(c_int) :: h1, h2, p3, p4

  tcc2((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:6)=0d0
  do h1 = 1, norb1
  do h2 = 1, norb1
  do p3 = norb1 + 1, nact
  do p4 = norb1 + 1, nact
     if (h1.ne.h2 .and. p3.ne.p4) then
        tcc2(p3, p4, h1, h2, 1) = + cc2(p3, p4, h1, h2, 1)   ! 1 (p+p+|h+h+)
        tcc2(p3, p4, h1, h2, 2) = + cc2(p3, p4, h1, h2, 1)   ! 2 (p-p-|h-h-)
     end if
     tcc2(p3, p4, h1, h2, 3) = + cc2(p3, p4, h1, h2, 2)   ! 3 (p+p-|h+h-)
     tcc2(p3, p4, h1, h2, 4) = + cc2(p4, p3, h2, h1, 2)   ! 4 (p-p+|h-h+)
     tcc2(p3, p4, h1, h2, 5) = - cc2(p3, p4, h2, h1, 2)   ! 5 (p+p-|h-h+)
     tcc2(p3, p4, h1, h2, 6) = - cc2(p4, p3, h1, h2, 2)   ! 6 (p-p+|h+h-)
  end do
  end do
  end do
  end do

end subroutine ccnew_mktcc2
!######################################################################
