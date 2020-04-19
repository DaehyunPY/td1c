!##########################################################
subroutine ccdt_t2p_man01(i0,work1,work2,work3)

!2: i0 ( a b i j )_tf + = -1 * P( 2 ) * Sum ( k ) * t ( a b i k )_t * i1 ( k j )_f 2
!3: i0 ( a b i j )_tf + =  1 * P( 2 ) * Sum ( c ) * t ( a c i j )_t * i1 ( b c )_f 2

  use mod_ormas,only : nact
  use mod_cc,only : int2x,norb1,ncc2aa,ncc2ab,t2inp,t3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,1:norb1),work2((norb1+1):nact,(norb1+1):nact),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_t2p_man01_1(work1) ! automatic 2_1,2
  call ccdt_t2p_man01_2(work2) ! automatic 3_1,2

  !$omp parallel default(shared) private(a,b,i,j,tmp)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     tmp = 0d0
     ! diagram 2
     do k = 1,norb1
        tmp = tmp - t2inp(a,b,i,k,spin_t2aa)*work1(k,j) &
                  + t2inp(a,b,j,k,spin_t2aa)*work1(k,i)
     end do
     ! diagram 3
     do c = norb1+1,nact
        tmp = tmp + t2inp(a,c,i,j,spin_t2aa)*work2(b,c) &
                  - t2inp(b,c,i,j,spin_t2aa)*work2(a,c)
     end do
     i0(a,b,i,j,1) = i0(a,b,i,j,1) + tmp
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     tmp = 0d0
     ! diagram 2
     do k = 1,norb1
        tmp = tmp - t2inp(a,b,i,k,spin_t2ab)*work1(k,j) &
                  - t2inp(a,b,k,j,spin_t2ab)*work1(k,i)
     end do
     ! diagram 3
     do c = norb1+1,nact
        tmp = tmp + t2inp(a,c,i,j,spin_t2ab)*work2(b,c) &
                  + t2inp(c,b,i,j,spin_t2ab)*work2(a,c)
     end do
     i0(a,b,i,j,2) = i0(a,b,i,j,2) + tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t2p_man01
!##########################################################
subroutine ccdt_t2p_man01_1(i1)

! i1 ( i j )_f + = 1 * f ( i j )_f 0
! i1 ( i j )_vt + = 1/2 * Sum ( k a b ) * t ( a b j k )_t * v ( i k a b )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,fock,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, noo
     i = h1_oo(icc)
     j = h2_oo(icc)
     tmp = fock(i,j,spin_focka)
     do k = 1,norb1
     do a = norb1+1,nact
        do b = norb1+1,a-1
           tmp = tmp + t2inp(a,b,j,k,spin_t2aa)*int2x(i,k,a,b,spin_int2aa)
        end do
        do b = norb1+1,nact
           tmp = tmp + t2inp(a,b,j,k,spin_t2ab)*int2x(i,k,a,b,spin_int2ab)
        end do
     end do
     end do
     i1(i,j) = i1(i,j) + tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t2p_man01_1
!##########################################################
subroutine ccdt_t2p_man01_2(i1)

! i1 ( a b )_f + = 1 * f ( a b )_f 0
! i1 ( a b )_vt + = -1/2 * Sum ( i j c ) * t ( a c i j )_t * v ( i j b c )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,fock,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, nvv
     a = p1_vv(icc)
     b = p2_vv(icc)
     tmp = fock(a,b,spin_focka)
     do c = norb1+1,nact
     do i = 1,norb1
        do j = 1,i-1
           tmp = tmp - t2inp(a,c,i,j,spin_t2aa)*int2x(i,j,b,c,spin_int2aa)
        end do
        do j = 1,norb1
           tmp = tmp - t2inp(a,c,i,j,spin_t2ab)*int2x(i,j,b,c,spin_int2ab)
        end do
     end do
     end do
     i1(a,b) = i1(a,b) + tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t2p_man01_2
!##########################################################
