!##########################################################
subroutine ccdt_t3p_man03(i0,work1,work2,work3)

! i0 ( a b c i j k )_tf + = -1 * P( 3 ) * Sum ( l ) 
!  * t ( a b c i j l )_t * i1 ( l k )_f 2

  use mod_ormas, only : nact
  use mod_cc, only : t3inp,ncc3aaa,ncc3aab,norb1
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa 
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab 
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1:norb1,1:norb1),work2(1),work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  call ccdt_t3p_man03_1(work1)

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
     do l = 1,norb1
        i0(a,b,c,i,j,k,1) = i0(a,b,c,i,j,k,1) &
             - t3inp(a,b,c,i,j,l,spin_t3aaa) * work1(l,k) &
             - t3inp(a,b,c,l,j,k,spin_t3aaa) * work1(l,i) &
             - t3inp(a,b,c,i,l,k,spin_t3aaa) * work1(l,j)
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
     do l = 1,norb1
        i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
             - t3inp(a,b,c,i,j,l,spin_t3aab) * work1(l,k) &
             - t3inp(a,b,c,l,j,k,spin_t3aab) * work1(l,i) &
             - t3inp(a,b,c,i,l,k,spin_t3aab) * work1(l,j)
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man03
!##########################################################
subroutine ccdt_t3p_man03_1(i1)

! i1 ( i j )_f + = 1 * f ( i j )_f 0
! i1 ( i j )_vt + = 1/2 * Sum ( k a b ) * t ( a b j k )_t * v ( i k a b )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,fock,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, noo
     i = h1_oo(icc)
     j = h2_oo(icc)
     i1(i,j) = i1(i,j) + fock(i,j,spin_focka)
     do k = 1,norb1
     do a = norb1+1,nact
        do b = norb1+1,a-1
           i1(i,j) = i1(i,j) + t2inp(a,b,j,k,spin_t2aa)*int2x(i,k,a,b,spin_int2aa)
        end do
        do b = norb1+1,nact
           i1(i,j) = i1(i,j) + t2inp(a,b,j,k,spin_t2ab)*int2x(i,k,a,b,spin_int2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man03_1
!##########################################################
