!##########################################################
subroutine ccdt_l2p_man11(i0,work1,work2,work3)

!11: i0 ( i j a b )_ytv + = -1/4 * Sum ( k l ) * i1 ( i j k l )_yt * v ( k l a b )_v 2

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1:norb1,1:norb1,1:norb1,1:norb1)
  complex(kind(0d0)),intent(inout) :: work2(1:norb1,1:norb1,1:norb1,1:norb1)
  complex(kind(0d0)),intent(inout) :: work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
!  call ccdt_l2p_11_1(1,1,1,1,work1)
!  call ccdt_l2p_11_2(1,1,1,1,work1)
!  call ccdt_l2p_11_1(1,2,1,2,work2)
!  call ccdt_l2p_11_2(1,2,1,2,work2)
!  work1 = -0.5d0*work1
!  work2 = -0.5d0*work2
  call ccdt_l2p_man11_1(work1,work2)
  call tdcc_fillooooaa(work1)

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     do k = 1,norb1
     do l = 1,k-1
        i0(i,j,a,b,1) = i0(i,j,a,b,1) + work1(i,j,k,l)*int2x(k,l,a,b,spin_int2aa)
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     do k = 1,norb1
     do l = 1,norb1
        i0(i,j,a,b,2) = i0(i,j,a,b,2) + work2(i,j,k,l)*int2x(k,l,a,b,spin_int2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man11
!##########################################################
subroutine ccdt_l2p_man11_1(i1aa,i1ab)

! 11-1
! i1 ( i j k l )_yt + = 1/2 * Sum ( a b ) 
!  * t ( a b k l )_t * y ( i j a b )_y 0
! 11-2
! i1 ( i j k l )_yt + = 1/6 * Sum ( m a b c ) 
!  * t ( a b c m k l )_t * y ( i j m a b c )_y 0
!
! Demanding: ORDER-8
!
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,cc_rank
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1aa(1:norb1,1:norb1,1:norb1,1:norb1)
  complex(kind(0d0)),intent(inout) :: i1ab(1:norb1,1:norb1,1:norb1,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooooaa
     i = h1_ooooaa(icc)
     j = h2_ooooaa(icc)
     k = h3_ooooaa(icc)
     l = h4_ooooaa(icc)
     ! 11-1
     do a = norb1+1,nact
     do b = norb1+1,a-1
        i1aa(i,j,k,l) = i1aa(i,j,k,l) + t2inp(a,b,k,l,spin_t2aa)*g2inp(i,j,a,b,spin_g2aa)
     end do
     end do

     ! 11-2
     if (cc_rank < 3) cycle
     do m = 1,norb1
     do a = norb1+1,nact
     do b = norb1+1,a-1
        do c = norb1+1,b-1
           i1aa(i,j,k,l) = i1aa(i,j,k,l) &
                + t3inp(a,b,c,k,l,m,spin_t3aaa)*g3inp(i,j,m,a,b,c,spin_g3aaa)
        end do
        do c = norb1+1,nact
           i1aa(i,j,k,l) = i1aa(i,j,k,l) &
                + t3inp(a,b,c,k,l,m,spin_t3aab)*g3inp(i,j,m,a,b,c,spin_g3aab)
        end do
     end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nooooab
     i = h1_ooooab(icc)
     j = h2_ooooab(icc)
     k = h3_ooooab(icc)
     l = h4_ooooab(icc)
     ! 11-1
     do a = norb1+1,nact
     do b = norb1+1,nact
        i1ab(i,j,k,l) = i1ab(i,j,k,l) + t2inp(a,b,k,l,spin_t2ab)*g2inp(i,j,a,b,spin_g2ab)
     end do
     end do

     ! 11-2
     if (cc_rank < 3) cycle
     do m = 1,norb1
     do a = norb1+1,nact
     do b = norb1+1,a-1
        do c = norb1+1,nact
           i1ab(i,j,k,l) = i1ab(i,j,k,l) &
                + t3inp(a,b,c,k,m,l,spin_t3aab)*g3inp(i,m,j,a,b,c,spin_g3aab) &
                + t3inp(a,b,c,l,m,k,spin_t3aab)*g3inp(j,m,i,a,b,c,spin_g3aab)
        end do
     end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man11_1
!##########################################################
