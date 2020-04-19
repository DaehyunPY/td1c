!######################################################################
subroutine ccnew_mktcc3(cc3, tcc3)

  !ASSUMING RHF reference

  use, intrinsic :: iso_c_binding
  use omp_mod
  use wfn_mod, only : nact, norb1

  implicit none
  complex(kind(0d0)),intent(in) ::   cc3((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:*)
  complex(kind(0d0)),intent(out) :: tcc3((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:*)

  integer(c_int) :: h1,h2,h3,p4,p5,p6

!  call ccnew_tcc3_clean(cc3)

  tcc3((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:20)=0d0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
  do p4 = norb1 + 1,nact
  do p5 = norb1 + 1,nact
  do p6 = norb1 + 1,nact
     if (h1.ne.h2 .and. h2.ne.h3 .and. h3.ne.h1 .and. &
         p4.ne.p5 .and. p5.ne.p6 .and. p6.ne.p4) then
        tcc3(p4,p5,p6,h1,h2,h3, 1) = + cc3(p4,p5,p6,h1,h2,h3,1) ! (+++|+++)
        tcc3(p4,p5,p6,h1,h2,h3, 2) = + cc3(p4,p5,p6,h1,h2,h3,1) ! (---|---)
     end if

     if (p4.ne.p5) then
        if (h1.ne.h2) tcc3(p4,p5,p6,h1,h2,h3, 3) = + cc3(p4,p5,p6,h1,h2,h3,2) ! (++-|++-)
        if (h3.ne.h1) tcc3(p4,p5,p6,h1,h2,h3, 4) = - cc3(p4,p5,p6,h1,h3,h2,2) ! (++-|+-+)
        if (h2.ne.h3) tcc3(p4,p5,p6,h1,h2,h3, 5) = - cc3(p4,p5,p6,h3,h2,h1,2) ! (++-|-++)

        if (h1.ne.h2) tcc3(p4,p5,p6,h1,h2,h3,12) = + cc3(p4,p5,p6,h1,h2,h3,2) ! (--+|--+)
        if (h3.ne.h1) tcc3(p4,p5,p6,h1,h2,h3,13) = - cc3(p4,p5,p6,h1,h3,h2,2) ! (--+|-+-)
        if (h2.ne.h3) tcc3(p4,p5,p6,h1,h2,h3,14) = - cc3(p4,p5,p6,h3,h2,h1,2) ! (--+|+--)
     end if

     if (p6.ne.p4) then
        if (h1.ne.h2) tcc3(p4,p5,p6,h1,h2,h3, 6) = - cc3(p4,p6,p5,h1,h2,h3,2) ! (+-+|++-)
        if (h3.ne.h1) tcc3(p4,p5,p6,h1,h2,h3, 7) = + cc3(p4,p6,p5,h1,h3,h2,2) ! (+-+|+-+)
        if (h2.ne.h3) tcc3(p4,p5,p6,h1,h2,h3, 8) = + cc3(p4,p6,p5,h3,h2,h1,2) ! (+-+|-++)

        if (h1.ne.h2) tcc3(p4,p5,p6,h1,h2,h3,15) = - cc3(p4,p6,p5,h1,h2,h3,2) ! (-+-|--+)
        if (h3.ne.h1) tcc3(p4,p5,p6,h1,h2,h3,16) = + cc3(p4,p6,p5,h1,h3,h2,2) ! (-+-|-+-)
        if (h2.ne.h3) tcc3(p4,p5,p6,h1,h2,h3,17) = + cc3(p4,p6,p5,h3,h2,h1,2) ! (-+-|+--)
     end if

     if (p5.ne.p6) then
        if (h1.ne.h2) tcc3(p4,p5,p6,h1,h2,h3, 9) = - cc3(p6,p5,p4,h1,h2,h3,2) ! (-++|++-)
        if (h3.ne.h1) tcc3(p4,p5,p6,h1,h2,h3,10) = + cc3(p6,p5,p4,h1,h3,h2,2) ! (-++|+-+)
        if (h2.ne.h3) tcc3(p4,p5,p6,h1,h2,h3,11) = + cc3(p6,p5,p4,h3,h2,h1,2) ! (-++|-++)

        if (h1.ne.h2) tcc3(p4,p5,p6,h1,h2,h3,18) = - cc3(p6,p5,p4,h1,h2,h3,2) ! (+--|--+)
        if (h3.ne.h1) tcc3(p4,p5,p6,h1,h2,h3,19) = + cc3(p6,p5,p4,h1,h3,h2,2) ! (+--|-+-)
        if (h2.ne.h3) tcc3(p4,p5,p6,h1,h2,h3,20) = + cc3(p6,p5,p4,h3,h2,h1,2) ! (+--|+--)
     end if

!old     tcc3(p4,p5,p6,h1,h2,h3, 1) = + cc3(p4,p5,p6,h1,h2,h3,1) ! (+++|+++)
!old     tcc3(p4,p5,p6,h1,h2,h3, 2) = + cc3(p4,p5,p6,h1,h2,h3,1) ! (---|---)
!old
!old     tcc3(p4,p5,p6,h1,h2,h3, 3) = + cc3(p4,p5,p6,h1,h2,h3,2) ! (++-|++-)
!old     tcc3(p4,p5,p6,h1,h2,h3, 4) = - cc3(p4,p5,p6,h1,h3,h2,2) ! (++-|+-+)
!old     tcc3(p4,p5,p6,h1,h2,h3, 5) = + cc3(p4,p5,p6,h3,h2,h1,2) ! (++-|-++)
!old     tcc3(p4,p5,p6,h1,h2,h3, 6) = - cc3(p4,p6,p5,h1,h2,h3,2) ! (+-+|++-)
!old     tcc3(p4,p5,p6,h1,h2,h3, 7) = + cc3(p4,p6,p5,h1,h3,h2,2) ! (+-+|+-+)
!old     tcc3(p4,p5,p6,h1,h2,h3, 8) = - cc3(p4,p6,p5,h3,h2,h1,2) ! (+-+|-++)
!old     tcc3(p4,p5,p6,h1,h2,h3, 9) = + cc3(p6,p5,p4,h1,h2,h3,2) ! (-++|++-)
!old     tcc3(p4,p5,p6,h1,h2,h3,10) = - cc3(p6,p5,p4,h1,h3,h2,2) ! (-++|+-+)
!old     tcc3(p4,p5,p6,h1,h2,h3,11) = + cc3(p6,p5,p4,h3,h2,h1,2) ! (-++|-++)
!old
!old     tcc3(p4,p5,p6,h1,h2,h3,12) = + cc3(p4,p5,p6,h1,h2,h3,2) ! (--+|--+)
!old     tcc3(p4,p5,p6,h1,h2,h3,13) = - cc3(p4,p5,p6,h1,h3,h2,2) ! (--+|-+-)
!old     tcc3(p4,p5,p6,h1,h2,h3,14) = + cc3(p4,p5,p6,h3,h2,h1,2) ! (--+|+--)
!old     tcc3(p4,p5,p6,h1,h2,h3,15) = - cc3(p4,p6,p5,h1,h2,h3,2) ! (-+-|--+)
!old     tcc3(p4,p5,p6,h1,h2,h3,16) = + cc3(p4,p6,p5,h1,h3,h2,2) ! (-+-|-+-)
!old     tcc3(p4,p5,p6,h1,h2,h3,17) = - cc3(p4,p6,p5,h3,h2,h1,2) ! (-+-|+--)
!old     tcc3(p4,p5,p6,h1,h2,h3,18) = + cc3(p6,p5,p4,h1,h2,h3,2) ! (+--|--+)
!old     tcc3(p4,p5,p6,h1,h2,h3,19) = - cc3(p6,p5,p4,h1,h3,h2,2) ! (+--|-+-)
!old     tcc3(p4,p5,p6,h1,h2,h3,20) = + cc3(p6,p5,p4,h3,h2,h1,2) ! (+--|+--)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccnew_mktcc3
!######################################################################
