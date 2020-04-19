!##########################################################
subroutine ccdt_l2p_man08(i0,work1,work2,work3)

!8:  i0 ( i j a b )_yv + = -1/2 * P( 2 ) * Sum ( k c d ) * y ( i j k a c d )_y * i1 ( c d k b )_v 4

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work2((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
!  call ccdt_l2p_new_8_1(1,1,1,1,work1)
!  call ccdt_l2p_new_8_2(1,1,1,1,work1)
!  call ccdt_l2p_new_8_1(1,2,1,2,work2)
!  call ccdt_l2p_new_8_2(1,2,1,2,work2)
!  call ccdt_l2p_new_8_3(work1,work2)
!  call ccdt_l2p_new_8_4(work1,work2)
  call ccdt_l2p_man08_1(work1,work2)
  call tdcc_fillvvovaa(work1)

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     do k = 1,norb1
     do c = norb1+1,nact
        do d = norb1+1,c-1
           i0(i,j,a,b,1) = i0(i,j,a,b,1) &
                - g3inp(i,j,k,a,c,d,spin_g3aaa)*work1(c,d,k,b) &
                + g3inp(i,j,k,b,c,d,spin_g3aaa)*work1(c,d,k,a)
        end do
        do d = norb1+1,nact
           i0(i,j,a,b,1) = i0(i,j,a,b,1) &
                + g3inp(i,j,k,a,c,d,spin_g3aab)*work2(d,c,k,b) &
                - g3inp(i,j,k,b,c,d,spin_g3aab)*work2(d,c,k,a)
        end do
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
     do c = norb1+1,nact
        do d = norb1+1,c-1
           i0(i,j,a,b,2) = i0(i,j,a,b,2) &
                - g3inp(j,k,i,c,d,a,spin_g3aab)*work1(c,d,k,b) &
                + g3inp(i,k,j,d,c,b,spin_g3aab)*work1(c,d,k,a)
        end do
        do d = norb1+1,nact
           i0(i,j,a,b,2) = i0(i,j,a,b,2) &
                + g3inp(i,k,j,a,c,d,spin_g3aab)*work2(c,d,k,b) &
                - g3inp(k,j,i,b,c,d,spin_g3aab)*work2(c,d,k,a)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man08
!##########################################################
subroutine ccdt_l2p_man08_1(i1aa,i1ab)

! 8-1
! i1 ( a b i c )_v + = 1 * v ( a b i c )_v 0
! 8-2
! i1 ( a b i c )_vt + = 1/2 * Sum ( j k ) * t ( a b j k )_t * v ( j k i c )_v 0
! 8-3
! original: not antisymmetric
! i1 ( a b i c )_vt + = 2 * Sum ( j d ) * t ( a d j i )_t * v ( j b c d )_v 0
! revised: force antisymmetricity
! i1 ( a b i c )_vt + = Sum ( j d ) * t ( a d j i )_t * v ( j b c d )_v 0
!                         - Sum ( j d ) * t ( b d j i )_t * v ( j a c d )_v 0

! 8-4
! i1 ( a b i c )_vt + = 1/2 * Sum ( j k d ) 
!  * t ( a d b j k i )_t * v ( j k c d )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,t3inp,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       i1ab((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novvvaa
     a = p3_ovvvaa(icc)
     b = p4_ovvvaa(icc)
     i = h1_ovvvaa(icc)
     c = p2_ovvvaa(icc)
     ! 8-1
     i1aa(a,b,i,c) = i1aa(a,b,i,c) + int2x(a,b,i,c,spin_int2aa)

     ! 8-2
     do j = 1,norb1
        do k = 1,j-1
           i1aa(a,b,i,c) = i1aa(a,b,i,c) + t2inp(a,b,j,k,spin_t2aa)*int2x(j,k,i,c,spin_int2aa)
        end do
     end do

     ! 8-3
     do j = 1,norb1
     do d = norb1+1,nact
        i1aa(a,b,i,c) = i1aa(a,b,i,c) &
             + t2inp(a,d,j,i,spin_t2aa)*int2x(j,b,c,d,spin_int2aa) &
             + t2inp(a,d,i,j,spin_t2ab)*int2x(b,j,c,d,spin_int2ab) &
             - t2inp(b,d,j,i,spin_t2aa)*int2x(j,a,c,d,spin_int2aa) &
             - t2inp(b,d,i,j,spin_t2ab)*int2x(a,j,c,d,spin_int2ab)
     end do
     end do

     ! 8-4
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1aa(a,b,i,c) = i1aa(a,b,i,c) + t3inp(a,d,b,j,k,i,spin_t3aaa)*int2x(j,k,c,d,spin_int2aa)
        end do
        do k = 1,norb1
           i1aa(a,b,i,c) = i1aa(a,b,i,c) + t3inp(a,b,d,j,i,k,spin_t3aab)*int2x(j,k,c,d,spin_int2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, novvvab
     a = p3_ovvvab(icc)
     b = p4_ovvvab(icc)
     i = h1_ovvvab(icc)
     c = p2_ovvvab(icc)
     ! 8-1
     i1ab(a,b,i,c) = i1ab(a,b,i,c) + int2x(a,b,i,c,spin_int2ab)

     ! 8-2
     do j = 1,norb1
        do k = 1,norb1
           i1ab(a,b,i,c) = i1ab(a,b,i,c) + t2inp(a,b,j,k,spin_t2ab)*int2x(j,k,i,c,spin_int2ab)
        end do
     end do

     ! 8-3
     do j = 1,norb1
     do d = norb1+1,nact
        i1ab(a,b,i,c) = i1ab(a,b,i,c) &
             - t2inp(a,d,j,i,spin_t2aa)*int2x(j,b,d,c,spin_int2ab) &
             - t2inp(a,d,i,j,spin_t2ab)*int2x(j,b,c,d,spin_int2aa) &
             - t2inp(b,d,j,i,spin_t2ab)*int2x(j,a,c,d,spin_int2ab)
     end do
     end do

     ! 8-4
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1ab(a,b,i,c) = i1ab(a,b,i,c) + t3inp(d,b,a,j,k,i,spin_t3aab)*int2x(j,k,c,d,spin_int2aa)
        end do
        do k = 1,norb1
           i1ab(a,b,i,c) = i1ab(a,b,i,c) + t3inp(a,d,b,k,i,j,spin_t3aab)*int2x(j,k,c,d,spin_int2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man08_1
!##########################################################
