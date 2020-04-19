!##########################################################
subroutine ccdt_den2p_man14(i0ab,i0ba,work1,work2,work3)

! i0 ( a b c i )_ytt + = -1/4 * Sum ( j k ) 
!  * i1 ( j k c i )_yt * t ( a b j k )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1(1:norb1,1:norb1,(norb1+1):nact,1:norb1), &
       work2(1:norb1,1:norb1,(norb1+1):nact,1:norb1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den2p_man14_1(work1,work2)
  !call tdcc_filloovoaa(work1) ! i1aa not needed

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, novvvab
     a = p3_ovvvab(icc)
     b = p4_ovvvab(icc)
     c = p2_ovvvab(icc)
     i = h1_ovvvab(icc)
     tmp = 0d0
     do j = 1,norb1
     do k = 1,norb1
        tmp = tmp - work2(j,k,c,i)*t2inp(a,b,j,k,spin_t2ab)
     end do
     end do
     i0ab(a,b,c,i) = i0ab(a,b,c,i) + tmp
     i0ba(b,a,c,i) = i0ba(b,a,c,i) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man14
!##########################################################
subroutine ccdt_den2p_man14_1(i1aa,i1ab)

! i1 ( i j a k )_yt + = +1 * Sum ( l b c ) 
!  * y ( i j l a b c )_y * t ( b c k l )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,(norb1+1):nact,1:norb1), &
       i1ab(1:norb1,1:norb1,(norb1+1):nact,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  ! i1aa not needed
  i1aa = 0d0

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooovab
     i = h1_ooovab(icc)
     j = h2_ooovab(icc)
     a = p4_ooovab(icc)
     k = h3_ooovab(icc)
     do l = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                + g3inp(j,l,i,b,c,a,spin_g3aab)*t2inp(b,c,k,l,spin_t2aa)
        end do

        do c = norb1+1,nact
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                + g3inp(i,l,j,a,c,b,spin_g3aab)*t2inp(b,c,k,l,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man14_1
!##########################################################
