!##########################################################
subroutine ccdt_l3p_man07(i0,work1,work2,work3)

! i0 ( i j k a b c )_yv + = -1 * P( 9 ) * Sum ( l d ) 
!  * y ( i j l a b d )_y * i1 ( k d l c )_v 2

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: &
       work1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       work2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       work3(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  work3 = 0d0
!  call ccdt_l3p_7_1(1,1,1,1,work1)
!  call ccdt_l3p_7_2(1,1,1,1,work1)
!  call ccdt_l3p_7_1(1,2,1,2,work2)
!  call ccdt_l3p_7_2(1,2,1,2,work2)
!  call ccdt_l3p_7_1(1,2,2,1,work3)
!  call ccdt_l3p_7_2(1,2,2,1,work3)
  call ccdt_l3p_man07_1(work1,work2,work3)

  !##################################################
  !$omp parallel default(shared) private(a,b,c,i,j,k)
  !$omp do
  do icc = 1,ncc3aaa
     a = p1_cc3aaa(icc)
     b = p2_cc3aaa(icc)
     c = p3_cc3aaa(icc)
     i = h1_cc3aaa(icc)
     j = h2_cc3aaa(icc)
     k = h3_cc3aaa(icc)
     do l = 1,norb1
     do d = norb1+1,nact
        i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
             - g3inp(i,j,l,a,b,d,spin_g3aaa)*work1(k,d,l,c) &
             - g3inp(i,j,l,a,b,d,spin_g3aab)*work3(k,d,l,c) &
             + g3inp(k,j,l,a,b,d,spin_g3aaa)*work1(i,d,l,c) &
             + g3inp(k,j,l,a,b,d,spin_g3aab)*work3(i,d,l,c) &
             + g3inp(i,k,l,a,b,d,spin_g3aaa)*work1(j,d,l,c) &
             + g3inp(i,k,l,a,b,d,spin_g3aab)*work3(j,d,l,c) &
             + g3inp(i,j,l,c,b,d,spin_g3aaa)*work1(k,d,l,a) &
             + g3inp(i,j,l,c,b,d,spin_g3aab)*work3(k,d,l,a) &
             + g3inp(i,j,l,a,c,d,spin_g3aaa)*work1(k,d,l,b) &
             + g3inp(i,j,l,a,c,d,spin_g3aab)*work3(k,d,l,b) &
             - g3inp(k,j,l,c,b,d,spin_g3aaa)*work1(i,d,l,a) &
             - g3inp(k,j,l,c,b,d,spin_g3aab)*work3(i,d,l,a) &
             - g3inp(i,k,l,c,b,d,spin_g3aaa)*work1(j,d,l,a) &
             - g3inp(i,k,l,c,b,d,spin_g3aab)*work3(j,d,l,a) &
             - g3inp(k,j,l,a,c,d,spin_g3aaa)*work1(i,d,l,b) &
             - g3inp(k,j,l,a,c,d,spin_g3aab)*work3(i,d,l,b) &
             - g3inp(i,k,l,a,c,d,spin_g3aaa)*work1(j,d,l,b) &
             - g3inp(i,k,l,a,c,d,spin_g3aab)*work3(j,d,l,b)
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc3aab
     a = p1_cc3aab(icc)
     b = p2_cc3aab(icc)
     c = p3_cc3aab(icc)
     i = h1_cc3aab(icc)
     j = h2_cc3aab(icc)
     k = h3_cc3aab(icc)
     do l = 1,norb1
     do d = norb1+1,nact
        i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
             - g3inp(i,j,l,a,b,d,spin_g3aaa)*work3(k,d,l,c) &
             - g3inp(i,j,l,a,b,d,spin_g3aab)*work1(k,d,l,c) &
             - g3inp(l,j,k,a,b,d,spin_g3aab)*work2(i,d,l,c) &
             - g3inp(i,l,k,a,b,d,spin_g3aab)*work2(j,d,l,c) &
             - g3inp(i,j,l,d,b,c,spin_g3aab)*work2(k,d,l,a) &
             - g3inp(i,j,l,a,d,c,spin_g3aab)*work2(k,d,l,b) &
             - g3inp(j,l,k,b,d,c,spin_g3aab)*work1(i,d,l,a) &
             - g3inp(k,l,j,c,d,b,spin_g3aab)*work3(i,d,l,a) &
             - g3inp(i,l,k,d,b,c,spin_g3aab)*work1(j,d,l,a) &
             - g3inp(l,k,i,c,d,b,spin_g3aab)*work3(j,d,l,a) &
             - g3inp(l,j,k,a,d,c,spin_g3aab)*work1(i,d,l,b) &
             - g3inp(k,l,j,d,c,a,spin_g3aab)*work3(i,d,l,b) &
             - g3inp(i,l,k,a,d,c,spin_g3aab)*work1(j,d,l,b) &
             - g3inp(k,l,i,c,d,a,spin_g3aab)*work3(j,d,l,b)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man07
!##########################################################
subroutine ccdt_l3p_man07_1(i1aa,i1ab,i1ba)

! 7-1
!     i1 ( i a j b )_v + = 1 * v ( i a j b )_v 0
! 7-2
!     i1 ( i a j b )_vt + = 1 * Sum ( k c ) * t ( a c k j )_t * v ( i k b c )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,int2x,t2inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       i1ab(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       i1ba(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novovab
     i = h1_ovovab(icc)
     a = p2_ovovab(icc)
     j = h3_ovovab(icc)
     b = p4_ovovab(icc)
     ! 7-1
     i1aa(i,a,j,b) = i1aa(i,a,j,b) + int2x(i,a,j,b,spin_int2aa)
     i1ab(i,a,j,b) = i1ab(i,a,j,b) + int2x(i,a,j,b,spin_int2ab)
     i1ba(i,a,j,b) = i1ba(i,a,j,b) - int2x(i,a,b,j,spin_int2ab)

     ! 7-2
     do k = 1,norb1
     do c = norb1+1,nact
        i1aa(i,a,j,b) = i1aa(i,a,j,b) &
          + t2inp(a,c,k,j,spin_t2aa)*int2x(i,k,b,c,spin_int2aa) &
          - t2inp(a,c,j,k,spin_t2ab)*int2x(i,k,b,c,spin_int2ab)
        i1ab(i,a,j,b) = i1ab(i,a,j,b) &
          - t2inp(a,c,k,j,spin_t2ab)*int2x(i,k,c,b,spin_int2ab)
        i1ba(i,a,j,b) = i1ba(i,a,j,b) &
          + t2inp(a,c,k,j,spin_t2aa)*int2x(i,k,b,c,spin_int2ab) &
          - t2inp(a,c,j,k,spin_t2ab)*int2x(i,k,b,c,spin_int2aa)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_man07_1
!##########################################################
!##########################################################
