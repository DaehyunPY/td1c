!##########################################################
subroutine ccdt_l2p_man03(i0,work1,work2,work3)

!2:  i0 ( i j a b )_yf + = -1 * P( 2 ) * Sum ( k ) * y ( i k a b )_y * i1 ( j k )_f 2
!3:  i0 ( i j a b )_yf + = 1 * P( 2 ) * Sum ( c ) * y ( i j a c )_y * i1 ( c b )_f 2

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1:norb1,1:norb1)
  complex(kind(0d0)),intent(inout) :: work2((norb1+1):nact,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
!  call ccdt_l2p_2_1(1,1,work1)
!  call ccdt_l2p_2_2(1,1,work1)
  call ccdt_l2p_man03_1(work1) ! automatic 2
  work2 = 0d0
!  call ccdt_l2p_3_1(1,1,work2)
!  call ccdt_l2p_3_2(1,1,work2)
  call ccdt_l2p_man03_2(work2) ! automatic 3

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     ! diagram 2
     do k = 1,norb1
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             - g2inp(i,k,a,b,spin_g2aa)*work1(j,k) &
             + g2inp(j,k,a,b,spin_g2aa)*work1(i,k)
     end do

     ! diagram 3
     do c = norb1+1,nact
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             + g2inp(i,j,a,c,spin_g2aa)*work2(c,b) &
             - g2inp(i,j,b,c,spin_g2aa)*work2(c,a)
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     ! diagram 2
     do k = 1,norb1
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             - g2inp(i,k,a,b,spin_g2ab)*work1(j,k) &
             - g2inp(k,j,a,b,spin_g2ab)*work1(i,k)
     end do

     ! diagram 3
     do c = norb1+1,nact
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             + g2inp(i,j,a,c,spin_g2ab)*work2(c,b) &
             + g2inp(i,j,c,b,spin_g2ab)*work2(c,a)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man03
!##########################################################
subroutine ccdt_l2p_man03_1(i1)

!2-1 i1 ( i j )_f + = 1 * f ( i j )_f 0
!2-2 i1 ( i j )_vt + = -1/2 * Sum ( k a b ) * t ( a b k j )_t * v ( i k a b )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,fock,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,e)
  !$omp do
  do icc = 1, noo
     i = h1_oo(icc)
     j = h2_oo(icc)
     i1(i,j) = i1(i,j) + fock(i,j,spin_focka)
     do k = 1,norb1
     do a = norb1+1,nact
        do b = norb1+1,a-1
           i1(i,j) = i1(i,j) - t2inp(a,b,k,j,spin_t2aa)*int2x(i,k,a,b,spin_int2aa)
        end do
        do b = norb1+1,nact
           i1(i,j) = i1(i,j) + t2inp(a,b,j,k,spin_t2ab)*int2x(i,k,a,b,spin_int2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man03_1
!##########################################################
subroutine ccdt_l2p_man03_2(i1)

!3-1 i1 ( a b )_f + = 1 * f ( a b )_f 0
!3-2 i1 ( a b )_vt + = -1/2 * Sum ( i j c ) * t ( a c i j )_t * v ( i j b c )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,fock,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,e)
  !$omp do
  do icc = 1, nvv
     a = p1_vv(icc)
     b = p2_vv(icc)
     i1(a,b) = i1(a,b) + fock(a,b,spin_focka)
     do c = norb1+1,nact
     do i = 1,norb1
        do j = 1,i-1
           i1(a,b) = i1(a,b) - t2inp(a,c,i,j,spin_t2aa)*int2x(i,j,b,c,spin_int2aa)
        end do
        do j = 1,norb1
           i1(a,b) = i1(a,b) - t2inp(a,c,i,j,spin_t2ab)*int2x(i,j,b,c,spin_int2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man03_2
!##########################################################
