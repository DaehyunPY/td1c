!##########################################################
subroutine ccdt_l2p_man10(i0,work1,work2,work3)

!9:  i0 ( i j a b )_ytv + = 1/2 * P( 2 ) * Sum ( c ) * i1 ( c a )_yt * v ( i j b c )_v 2
!10: i0 ( i j a b )_ytv + = 1/2 * P( 2 ) * Sum ( k ) * i1 ( i k )_yt * v ( j k a b )_v 2

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1((norb1+1):nact,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work2(1:norb1,1:norb1)
  complex(kind(0d0)),intent(inout) :: work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  call ccdt_l2p_man10_1(work1) ! automatic 9-1,2
  work2 = 0d0
  call ccdt_l2p_man10_2(work2) ! automatic 10-1,2

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     ! diagram 9
     do c = norb1+1,nact
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             - work1(c,a)*int2x(i,j,c,b,spin_int2aa) &
             - work1(c,b)*int2x(i,j,a,c,spin_int2aa)
     end do

     ! diagram 10
     do k = 1,norb1
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             - work2(i,k)*int2x(k,j,a,b,spin_int2aa) &
             - work2(j,k)*int2x(i,k,a,b,spin_int2aa)
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     ! diagram 9
     do c = norb1+1,nact
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             - work1(c,a)*int2x(i,j,c,b,spin_int2ab) &
             - work1(c,b)*int2x(i,j,a,c,spin_int2ab)
     end do

     ! diagram 10
     do k = 1,norb1
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             - work2(i,k)*int2x(k,j,a,b,spin_int2ab) &
             - work2(j,k)*int2x(i,k,a,b,spin_int2ab)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man10
!##########################################################
subroutine ccdt_l2p_man10_1(i1)

! 9-1
! i1 ( a b )_yt + = 1/2 * Sum ( i j c ) 
!  * t ( a c i j )_t * y ( i j b c )_y 0
! 9-2
! i1 ( a b )_yt + = 1/12 * Sum ( i j k c d ) 
!  * t ( c d a i j k )_t * y ( i j k b c d )_y 0
!
! Demanding: ORDER-7
!

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,cc_rank
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nvv
     a = p1_vv(icc)
     b = p2_vv(icc)
     ! 9-1
     do c = norb1+1,nact
     do i = 1,norb1
        do j = 1,i-1
           i1(a,b) = i1(a,b) + t2inp(a,c,i,j,spin_t2aa)*g2inp(i,j,b,c,spin_g2aa)
        end do
        do j = 1,norb1
           i1(a,b) = i1(a,b) + t2inp(a,c,i,j,spin_t2ab)*g2inp(i,j,b,c,spin_g2ab)
        end do
     end do
     end do

     if (cc_rank < 3) cycle
     ! 9-2
     do c = norb1+1,nact
     do d = norb1+1,c-1
     do i = 1,norb1
     do j = 1,i-1
     do k = 1,j-1
        i1(a,b) = i1(a,b) + t3inp(a,c,d,i,j,k,spin_t3aaa)*g3inp(i,j,k,b,c,d,spin_g3aaa)
     end do
     end do
     end do
     end do
     end do

     do c = norb1+1,nact
     do d = norb1+1,c-1
     do i = 1,norb1
     do j = 1,norb1
     do k = 1,j-1
        i1(a,b) = i1(a,b) + t3inp(c,d,a,j,k,i,spin_t3aab)*g3inp(j,k,i,c,d,b,spin_g3aab)
     end do
     end do
     end do
     end do
     end do

     do c = norb1+1,nact
     do d = norb1+1,nact
     do i = 1,norb1
     do j = 1,i-1
     do k = 1,norb1
        i1(a,b) = i1(a,b) + t3inp(a,c,d,i,j,k,spin_t3aab)*g3inp(i,j,k,b,c,d,spin_g3aab)
     end do
     end do
     end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man10_1
!##########################################################
subroutine ccdt_l2p_man10_2(i1)

!10-1
! i1 ( i j )_yt + = 1/2 * Sum ( k a b ) 
!  * t ( a b j k )_t * y ( i k a b )_y 0
!10-2
! i1 ( i j )_yt + = 1/12 * Sum ( k l a b c ) 
!  * t ( a b c j k l )_t * y ( i k l a b c )_y 0
!
! Demanding: ORDER-7
!
  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,cc_rank
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, noo
     i = h1_oo(icc)
     j = h2_oo(icc)
     !10-1
     do k = 1,norb1
     do a = norb1+1,nact
     do b = norb1+1,a-1
        i1(i,j) = i1(i,j) + t2inp(a,b,j,k,spin_t2aa)*g2inp(i,k,a,b,spin_g2aa)
     end do
     end do
     end do
     do k = 1,norb1
     do a = norb1+1,nact
     do b = norb1+1,nact
        i1(i,j) = i1(i,j) + t2inp(a,b,j,k,spin_t2ab)*g2inp(i,k,a,b,spin_g2ab)
     end do
     end do
     end do

     if (cc_rank < 3) cycle
     !10-2
     do k = 1,norb1
     do l = 1,k-1
     do a = norb1+1,nact
     do b = norb1+1,a-1
     do c = norb1+1,b-1
        i1(i,j) = i1(i,j) + t3inp(a,b,c,j,k,l,spin_t3aaa)*g3inp(i,k,l,a,b,c,spin_g3aaa)
     end do
     end do
     end do
     end do
     end do

     do k = 1,norb1
     do l = 1,k-1
     do a = norb1+1,nact
     do b = norb1+1,nact
     do c = norb1+1,b-1
        i1(i,j) = i1(i,j) + t3inp(b,c,a,k,l,j,spin_t3aab)*g3inp(k,l,i,b,c,a,spin_g3aab)
     end do
     end do
     end do
     end do
     end do

     do k = 1,norb1
     do l = 1,norb1
     do a = norb1+1,nact
     do b = norb1+1,a-1
     do c = norb1+1,nact
        i1(i,j) = i1(i,j) + t3inp(a,b,c,j,k,l,spin_t3aab)*g3inp(i,k,l,a,b,c,spin_g3aab)
     end do
     end do
     end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man10_2
!##########################################################
