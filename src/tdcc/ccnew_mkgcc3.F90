!######################################################################
subroutine ccnew_mkgcc3(cc3,gcc3)

  !ASSUMING RHF reference

  use, intrinsic :: iso_c_binding
  use omp_mod
  use wfn_mod,only : nact,norb1

  implicit none
  complex(kind(0d0)),intent(in) ::   cc3(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:*)
  complex(kind(0d0)),intent(out) :: gcc3(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:*)

  integer(c_int) :: h1,h2,h3,p1,p2,p3

!  call ccnew_gcc3_clean(cc3)

  gcc3(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:20)=0d0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
  do p1 = norb1 + 1,nact
  do p2 = norb1 + 1,nact
  do p3 = norb1 + 1,nact
     if (h1.ne.h2 .and. h2.ne.h3 .and. h3.ne.h1 .and. &
         p1.ne.p2 .and. p2.ne.p3 .and. p3.ne.p1) then
        gcc3(h1,h2,h3,p1,p2,p3, 1) = + cc3(h1,h2,h3,p1,p2,p3,1)  ! 1 (+++|+++)
        gcc3(h1,h2,h3,p1,p2,p3, 2) = + cc3(h1,h2,h3,p1,p2,p3,1)  ! 2 (---|---)
     end if

     if (p1.ne.p2) then
        if (h1.ne.h2) gcc3(h1,h2,h3,p1,p2,p3, 3) = + cc3(h1,h2,h3,p1,p2,p3,2)  ! 3 (++-|++-)
        if (h3.ne.h1) gcc3(h1,h2,h3,p1,p2,p3, 4) = - cc3(h1,h3,h2,p1,p2,p3,2)  ! 3 (+-+|++-)
        if (h2.ne.h3) gcc3(h1,h2,h3,p1,p2,p3, 5) = - cc3(h3,h2,h1,p1,p2,p3,2)  ! 3 (-++|++-)

        if (h1.ne.h2) gcc3(h1,h2,h3,p1,p2,p3,12) = + cc3(h1,h2,h3,p1,p2,p3,2)  ! 3 (++-|++-)
        if (h3.ne.h1) gcc3(h1,h2,h3,p1,p2,p3,13) = - cc3(h1,h3,h2,p1,p2,p3,2)  ! 3 (+-+|++-)
        if (h2.ne.h3) gcc3(h1,h2,h3,p1,p2,p3,14) = - cc3(h3,h2,h1,p1,p2,p3,2)  ! 3 (-++|++-)
     end if

     if (p3.ne.p1) then
        if (h1.ne.h2) gcc3(h1,h2,h3,p1,p2,p3, 6) = - cc3(h1,h2,h3,p1,p3,p2,2)  ! 3 (++-|+-+)
        if (h3.ne.h1) gcc3(h1,h2,h3,p1,p2,p3, 7) = + cc3(h1,h3,h2,p1,p3,p2,2)  ! 3 (+-+|+-+)
        if (h2.ne.h3) gcc3(h1,h2,h3,p1,p2,p3, 8) = + cc3(h3,h2,h1,p1,p3,p2,2)  ! 3 (-++|+-+)

        if (h1.ne.h2) gcc3(h1,h2,h3,p1,p2,p3,15) = - cc3(h1,h2,h3,p1,p3,p2,2)  ! 3 (++-|+-+)
        if (h3.ne.h1) gcc3(h1,h2,h3,p1,p2,p3,16) = + cc3(h1,h3,h2,p1,p3,p2,2)  ! 3 (+-+|+-+)
        if (h2.ne.h3) gcc3(h1,h2,h3,p1,p2,p3,17) = + cc3(h3,h2,h1,p1,p3,p2,2)  ! 3 (-++|+-+)
     end if

     if (p2.ne.p3) then
        if (h1.ne.h2) gcc3(h1,h2,h3,p1,p2,p3, 9) = - cc3(h1,h2,h3,p3,p2,p1,2)  ! 3 (++-|-++)
        if (h3.ne.h1) gcc3(h1,h2,h3,p1,p2,p3,10) = + cc3(h1,h3,h2,p3,p2,p1,2)  ! 3 (+-+|-++)
        if (h2.ne.h3) gcc3(h1,h2,h3,p1,p2,p3,11) = + cc3(h3,h2,h1,p3,p2,p1,2)  ! 3 (-++|-++)

        if (h1.ne.h2) gcc3(h1,h2,h3,p1,p2,p3,18) = - cc3(h1,h2,h3,p3,p2,p1,2)  ! 3 (++-|-++)
        if (h3.ne.h1) gcc3(h1,h2,h3,p1,p2,p3,19) = + cc3(h1,h3,h2,p3,p2,p1,2)  ! 3 (+-+|-++)
        if (h2.ne.h3) gcc3(h1,h2,h3,p1,p2,p3,20) = + cc3(h3,h2,h1,p3,p2,p1,2)  ! 3 (-++|-++)
     end if
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccnew_mkgcc3
!######################################################################
subroutine ccnew_mkgcc3_bug(cc3,gcc3)

  use omp_mod
  use wfn_mod,only : nact,norb1

  implicit none
  complex(kind(0d0)),intent(in) ::   cc3(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:*)
  complex(kind(0d0)),intent(out) :: gcc3(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:*)

  integer(c_int) :: h1,h2,h3,p1,p2,p3

  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
  do p1 = norb1 + 1,nact
  do p2 = norb1 + 1,nact
  do p3 = norb1 + 1,nact
     gcc3(h1,h2,h3,p1,p2,p3, 1) = + cc3(h1,h2,h3,p1,p2,p3,1)  ! 1 (+++|+++)
     gcc3(h1,h2,h3,p1,p2,p3, 2) = + cc3(h1,h2,h3,p1,p2,p3,1)  ! 2 (---|---)

     gcc3(h1,h2,h3,p1,p2,p3, 3) = + cc3(h1,h2,h3,p1,p2,p3,2)  ! 3 (++-|++-)
     gcc3(h1,h2,h3,p1,p2,p3, 4) = - cc3(h1,h2,h3,p1,p3,p2,2)  ! 3 (++-|+-+)
     gcc3(h1,h2,h3,p1,p2,p3, 5) = - cc3(h1,h2,h3,p3,p2,p1,2)  ! 3 (++-|-++)

     gcc3(h1,h2,h3,p1,p2,p3, 6) = - cc3(h1,h3,h2,p1,p2,p3,2)  ! 3 (+-+|++-)
     gcc3(h1,h2,h3,p1,p2,p3, 7) = + cc3(h1,h3,h2,p1,p3,p2,2)  ! 3 (+-+|+-+)
     gcc3(h1,h2,h3,p1,p2,p3, 8) = + cc3(h1,h3,h2,p3,p2,p1,2)  ! 3 (+-+|-++)

     gcc3(h1,h2,h3,p1,p2,p3, 9) = - cc3(h3,h2,h1,p1,p2,p3,2)  ! 3 (-++|++-)
     gcc3(h1,h2,h3,p1,p2,p3,10) = + cc3(h3,h2,h1,p1,p3,p2,2)  ! 3 (-++|+-+)
     gcc3(h1,h2,h3,p1,p2,p3,11) = + cc3(h3,h2,h1,p3,p2,p1,2)  ! 3 (-++|-++)

     gcc3(h1,h2,h3,p1,p2,p3,12) = + cc3(h1,h2,h3,p1,p2,p3,2)  ! 3 (++-|++-)
     gcc3(h1,h2,h3,p1,p2,p3,13) = - cc3(h1,h2,h3,p1,p3,p2,2)  ! 3 (++-|+-+)
     gcc3(h1,h2,h3,p1,p2,p3,14) = - cc3(h1,h2,h3,p3,p2,p1,2)  ! 3 (++-|-++)

     gcc3(h1,h2,h3,p1,p2,p3,15) = - cc3(h1,h3,h2,p1,p2,p3,2)  ! 3 (+-+|++-)
     gcc3(h1,h2,h3,p1,p2,p3,16) = + cc3(h1,h3,h2,p1,p3,p2,2)  ! 3 (+-+|+-+)
     gcc3(h1,h2,h3,p1,p2,p3,17) = + cc3(h1,h3,h2,p3,p2,p1,2)  ! 3 (+-+|-++)

     gcc3(h1,h2,h3,p1,p2,p3,18) = - cc3(h3,h2,h1,p1,p2,p3,2)  ! 3 (-++|++-)
     gcc3(h1,h2,h3,p1,p2,p3,19) = + cc3(h3,h2,h1,p1,p3,p2,2)  ! 3 (-++|+-+)
     gcc3(h1,h2,h3,p1,p2,p3,20) = + cc3(h3,h2,h1,p3,p2,p1,2)  ! 3 (-++|-++)
  end do
  end do
  end do
  end do
  end do
  end do

end subroutine ccnew_mkgcc3_bug
!######################################################################
