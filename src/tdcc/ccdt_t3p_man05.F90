!##########################################################
subroutine ccdt_t3p_man05(i0,work1,work2,work3)

! i0 ( a b c i j k )_vt + = 1/2 * P( 3 ) * Sum ( l m )
!  * t ( a b c i l m )_t * i1 ( l m j k )_v 2

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
       work1(1:norb1,1:norb1,1:norb1,1:norb1), &
       work2(1:norb1,1:norb1,1:norb1,1:norb1),work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_t3p_man05_1(work1,work2)
  call tdcc_fillooooaa(work1)

  !$omp parallel default(shared) private(a,b,c,i,j,k)
  !$omp do
  do icc = 1,ncc3aaa
     !aaaaaa
     a = p1_cc3aaa(icc)
     b = p2_cc3aaa(icc)
     c = p3_cc3aaa(icc)
     i = h1_cc3aaa(icc)
     j = h2_cc3aaa(icc)
     k = h3_cc3aaa(icc)
     do l = act1_ll,norb1
        do m = act1_ll,l-1
           i0(a,b,c,i,j,k,1) = i0(a,b,c,i,j,k,1) &
                + t3inp(a,b,c,i,l,m,spin_t3aaa) * work1(l,m,j,k) &
                + t3inp(a,b,c,l,j,m,spin_t3aaa) * work1(l,m,i,k) &
                + t3inp(a,b,c,l,m,k,spin_t3aaa) * work1(l,m,i,j)
        end do
     end do
  end do
  !$omp end do

  !$omp do
  do icc = 1,ncc3aab
     ! aabaab
     a = p1_cc3aab(icc)
     b = p2_cc3aab(icc)
     c = p3_cc3aab(icc)
     i = h1_cc3aab(icc)
     j = h2_cc3aab(icc)
     k = h3_cc3aab(icc)
     do l = act1_ll,norb1
        do m = act1_ll,norb1
           i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
                + t3inp(a,b,c,i,l,m,spin_t3aab) * work2(l,m,j,k) &
                + t3inp(a,b,c,l,j,m,spin_t3aab) * work2(l,m,i,k)
        end do
        do m = act1_ll,l-1
           i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
                + t3inp(a,b,c,l,m,k,spin_t3aab) * work1(l,m,i,j)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man05
!##########################################################
subroutine ccdt_t3p_man05_1(i1aa,i1ab)

! i1 ( i j k l )_v + = 1 * v ( i j k l )_v 0
! i1 ( i j k l )_vt + = 1/2 * Sum ( a b ) 
!  * t ( a b k l )_t * v ( i j a b )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,1:norb1,1:norb1), &
       i1ab(1:norb1,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooooaa
     i = h1_ooooaa(icc)
     j = h2_ooooaa(icc)
     k = h3_ooooaa(icc)
     l = h4_ooooaa(icc)
     i1aa(i,j,k,l) = i1aa(i,j,k,l) + int2x(i,j,k,l,spin_int2aa)
     do a = norb1+1,nact
        do b = norb1+1,a-1
           i1aa(i,j,k,l) = i1aa(i,j,k,l) + t2inp(a,b,k,l,spin_t2aa)*int2x(i,j,a,b,spin_int2aa)
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
     i1ab(i,j,k,l) = i1ab(i,j,k,l) + int2x(i,j,k,l,spin_int2ab)
     do a = norb1+1,nact
        do b = norb1+1,nact
           i1ab(i,j,k,l) = i1ab(i,j,k,l) + t2inp(a,b,k,l,spin_t2ab)*int2x(i,j,a,b,spin_int2ab)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man05_1
!##########################################################
