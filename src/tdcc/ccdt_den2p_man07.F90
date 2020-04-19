!##########################################################
subroutine ccdt_den2p_man07(i0ab,i0ba,work1,work2,work3)

!  i0 ( a b i j )_ytt + = +1/2 * Sum ( k l ) * i1 ( k l i j )_yt * t ( a b k l )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1(1:norb1,1:norb1,1:norb1,1:norb1),&
       work2(1:norb1,1:norb1,1:norb1,1:norb1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den2p_man07_1(work1,work2)
  !call tdcc_fillooooaa(work1) ! not needed

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,tmp)
  !$omp do
  do icc = 1, noovvab
     i = h1_oovvab(icc)
     j = h2_oovvab(icc)
     a = p3_oovvab(icc)
     b = p4_oovvab(icc)
     tmp = 0d0
     do k = 1,norb1
        do l = 1,norb1
           tmp = tmp + work2(k,l,i,j)*t2inp(a,b,k,l,spin_t2ab)
        end do
     end do
     i0ab(a,b,i,j) = i0ab(a,b,i,j) + tmp
     i0ba(a,b,j,i) = i0ba(a,b,j,i) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man07
!##########################################################
subroutine ccdt_den2p_man07_1(i1aa,i1ab)

! i1 ( i j k l )_yt + = +1/2 * Sum ( a b ) * y ( i j a b )_y * t ( a b k l )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,1:norb1,1:norb1), &
       i1ab(1:norb1,1:norb1,1:norb1,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  ! not needed
  i1aa = 0d0

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, nooooab
     i = h1_ooooab(icc)
     j = h2_ooooab(icc)
     k = h3_ooooab(icc)
     l = h4_ooooab(icc)
     do a = norb1+1,nact
        do b = norb1+1,nact
           i1ab(i,j,k,l) = i1ab(i,j,k,l) + g2inp(i,j,a,b,spin_g2ab)*t2inp(a,b,k,l,spin_t2ab)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man07_1
!##########################################################
