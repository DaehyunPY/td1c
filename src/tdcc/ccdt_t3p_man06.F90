!##########################################################
subroutine ccdt_t3p_man06(i0,work1,work2,work3)

! i0 ( a b c i j k )_vt + = -1 * P( 9 ) * Sum ( l d ) 
!  * t ( a b d i j l )_t * i1 ( l c k d )_v 2

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,act1_ll,act1_ul
  use mod_cc, only : t3inp,ncc3aaa,ncc3aab,norb1
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa 
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab 
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       work2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       work3(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  work3 = 0d0
  call ccdt_t3p_man06_1(work1,work2,work3)

  !$omp parallel default(shared) private(a,b,c,i,j,k)
  !$omp do
  do icc = 1, ncc3aaa
     !aaaaaa
     a = p1_cc3aaa(icc)
     b = p2_cc3aaa(icc)
     c = p3_cc3aaa(icc)
     i = h1_cc3aaa(icc)
     j = h2_cc3aaa(icc)
     k = h3_cc3aaa(icc)
     do l = act1_ll,norb1
     do d = norb1+1,act1_ul
        i0(a,b,c,i,j,k,1) = i0(a,b,c,i,j,k,1) &
             + t3inp(a,b,d,i,j,l,spin_t3aaa) * work1(l,c,k,d) + t3inp(a,b,d,i,j,l,spin_t3aab) * work3(l,c,k,d) &
             + t3inp(a,b,d,l,j,k,spin_t3aaa) * work1(l,c,i,d) - t3inp(a,b,d,k,j,l,spin_t3aab) * work3(l,c,i,d)&
             + t3inp(a,b,d,i,l,k,spin_t3aaa) * work1(l,c,j,d) - t3inp(a,b,d,i,k,l,spin_t3aab) * work3(l,c,j,d)&
             + t3inp(d,b,c,i,j,l,spin_t3aaa) * work1(l,a,k,d) - t3inp(c,b,d,i,j,l,spin_t3aab) * work3(l,a,k,d)&
             + t3inp(a,d,c,i,j,l,spin_t3aaa) * work1(l,b,k,d) - t3inp(a,c,d,i,j,l,spin_t3aab) * work3(l,b,k,d)&
             + t3inp(d,b,c,l,j,k,spin_t3aaa) * work1(l,a,i,d) + t3inp(c,b,d,k,j,l,spin_t3aab) * work3(l,a,i,d) &
             + t3inp(a,d,c,l,j,k,spin_t3aaa) * work1(l,b,i,d) + t3inp(a,c,d,k,j,l,spin_t3aab) * work3(l,b,i,d) &
             + t3inp(d,b,c,i,l,k,spin_t3aaa) * work1(l,a,j,d) + t3inp(c,b,d,i,k,l,spin_t3aab) * work3(l,a,j,d) &
             + t3inp(a,d,c,i,l,k,spin_t3aaa) * work1(l,b,j,d) + t3inp(a,c,d,i,k,l,spin_t3aab) * work3(l,b,j,d)
     end do
     end do
  end do
  !$omp end do

  !$omp do
  do icc = 1, ncc3aab
     ! aabaab
     a = p1_cc3aab(icc)
     b = p2_cc3aab(icc)
     c = p3_cc3aab(icc)
     i = h1_cc3aab(icc)
     j = h2_cc3aab(icc)
     k = h3_cc3aab(icc)
     do l = act1_ll,norb1
     do d = norb1+1,act1_ul
        i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
             + t3inp(a,b,d,i,j,l,spin_t3aaa) * work3(l,c,k,d) + t3inp(a,b,d,i,j,l,spin_t3aab) * work1(l,c,k,d) &
             + t3inp(a,b,d,l,j,k,spin_t3aab) * work2(l,c,i,d) &
             + t3inp(a,b,d,i,l,k,spin_t3aab) * work2(l,c,j,d) &
             + t3inp(d,b,c,i,j,l,spin_t3aab) * work2(l,a,k,d) &
             + t3inp(a,d,c,i,j,l,spin_t3aab) * work2(l,b,k,d) &
             + t3inp(d,b,c,l,j,k,spin_t3aab) * work1(l,a,i,d) + t3inp(c,d,b,k,l,j,spin_t3aab) * work3(l,a,i,d) &
             + t3inp(a,d,c,l,j,k,spin_t3aab) * work1(l,b,i,d) + t3inp(d,c,a,k,l,j,spin_t3aab) * work3(l,b,i,d) &
             + t3inp(d,b,c,i,l,k,spin_t3aab) * work1(l,a,j,d) + t3inp(d,c,b,k,l,i,spin_t3aab) * work3(l,a,j,d) &
             + t3inp(a,d,c,i,l,k,spin_t3aab) * work1(l,b,j,d) + t3inp(c,d,a,k,l,i,spin_t3aab) * work3(l,b,j,d)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man06
!##########################################################
subroutine ccdt_t3p_man06_1(i1aa,i1ab,i1ba)

! i1 ( i a j b )_v + = -1 * v ( i a j b )_v 0
! i1 ( i a j b )_vt + = Sum ( k c ) 
!  * t ( a c j k )_t * v ( i k b c )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       i1ab(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       i1ba(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novovab
     i = h1_ovovab(icc)
     a = p2_ovovab(icc)
     j = h3_ovovab(icc)
     b = p4_ovovab(icc)
     i1aa(i,a,j,b) = i1aa(i,a,j,b) - int2x(i,a,j,b,spin_int2aa)
     i1ab(i,a,j,b) = i1ab(i,a,j,b) - int2x(i,a,j,b,spin_int2ab)
     i1ba(i,a,j,b) = i1ba(i,a,j,b) + int2x(i,a,b,j,spin_int2ab)
     do k = 1,norb1
     do c = norb1+1,nact
        i1aa(i,a,j,b) = i1aa(i,a,j,b) &
             + t2inp(a,c,j,k,spin_t2aa)*int2x(i,k,b,c,spin_int2aa) &
             + t2inp(a,c,j,k,spin_t2ab)*int2x(i,k,b,c,spin_int2ab)
        i1ab(i,a,j,b) = i1ab(i,a,j,b) &
             + t2inp(c,a,j,k,spin_t2ab)*int2x(i,k,c,b,spin_int2ab)
        i1ba(i,a,j,b) = i1ba(i,a,j,b) &
             + t2inp(a,c,j,k,spin_t2ab)*int2x(i,k,b,c,spin_int2aa) &
             + t2inp(a,c,j,k,spin_t2aa)*int2x(i,k,b,c,spin_int2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man06_1
!##########################################################
