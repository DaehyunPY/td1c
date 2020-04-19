!##########################################################
subroutine ccdt_t3p_man04(i0,work1,work2,work3)

! i0 ( a b c i j k )_tf + = 1 * P( 3 ) * Sum ( d ) 
!  * t ( a b d i j k )_t * i1 ( c d )_f 2

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact,act1_ul
  use mod_cc, only : t3inp,ncc3aaa,ncc3aab,norb1
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa 
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab 
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: work1((norb1+1):nact,(norb1+1):nact),work2(1),work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  call ccdt_t3p_man04_1(work1)

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
     do d = norb1+1,act1_ul
        i0(a,b,c,i,j,k,1) = i0(a,b,c,i,j,k,1) &
             + t3inp(a,b,d,i,j,k,spin_t3aaa) * work1(c,d) &
             + t3inp(d,b,c,i,j,k,spin_t3aaa) * work1(a,d) &
             + t3inp(a,d,c,i,j,k,spin_t3aaa) * work1(b,d)
     end do
  end do
  !$omp end do

  !$omp do
  do icc = 1, ncc3aab
     !aabaab
     a = p1_cc3aab(icc)
     b = p2_cc3aab(icc)
     c = p3_cc3aab(icc)
     i = h1_cc3aab(icc)
     j = h2_cc3aab(icc)
     k = h3_cc3aab(icc)
     do d = norb1+1,act1_ul
        i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
             + t3inp(a,b,d,i,j,k,spin_t3aab) * work1(c,d) &
             + t3inp(d,b,c,i,j,k,spin_t3aab) * work1(a,d) &
             + t3inp(a,d,c,i,j,k,spin_t3aab) * work1(b,d)
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man04
!##########################################################
subroutine ccdt_t3p_man04_1(i1)

! i1 ( a b )_f + = 1 * f ( a b )_f 0
! i1 ( a b )_vt + = -1/2 * Sum ( i j c ) 
!  * t ( a c i j )_t * v ( i j b c )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,fock,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
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

end subroutine ccdt_t3p_man04_1
!##########################################################
