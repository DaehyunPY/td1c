!##########################################################
subroutine ccdt_l3p_man06(i0,work1,work2,work3)

! i0 ( i j k a b c )_yv + = -1/2 * P( 3 ) * Sum ( l m ) 
!  * y ( i m l a b c )_y * i1 ( j k m l )_v 2

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: &
       work1(1:norb1,1:norb1,1:norb1,1:norb1), &
       work2(1:norb1,1:norb1,1:norb1,1:norb1),work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  !call ccdt_l3p_6_1(1,1,1,1,work1)
  !call ccdt_l3p_6_2(1,1,1,1,work1)
  !call ccdt_l3p_6_1(1,2,1,2,work2)
  !call ccdt_l3p_6_2(1,2,1,2,work2)
  !work1 = -work1
  !work2 = -work2
  call ccdt_l3p_man06_1(work1,work2)
  call tdcc_fillooooaa(work1)

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
        do m = 1,l-1
           i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
                + g3inp(i,l,m,a,b,c,spin_g3aaa)*work1(j,k,l,m) &
                - g3inp(j,l,m,a,b,c,spin_g3aaa)*work1(i,k,l,m) &
                - g3inp(k,l,m,a,b,c,spin_g3aaa)*work1(j,i,l,m)
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
        do m = 1,l-1
           i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
                - g3inp(l,m,k,a,b,c,spin_g3aab)*work1(j,i,l,m)
        end do
     end do
     do l = 1,norb1
        do m = 1,norb1
           i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
                + g3inp(i,l,m,a,b,c,spin_g3aab)*work2(j,k,l,m) &
                - g3inp(j,l,m,a,b,c,spin_g3aab)*work2(i,k,l,m)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man06
!##########################################################
subroutine ccdt_l3p_man06_1(i1aa,i1ab)

! 6-1
!     i1 ( i j k l )_v + = v ( i j k l )_v 0
! 6-2
!     i1 ( i j k l )_vt + = 1/2 * Sum ( a b ) * t ( a b k l )_t * v ( i j a b )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,int2x,t2inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,1:norb1,1:norb1), &
       i1ab(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooooaa
     i = h1_ooooaa(icc)
     j = h2_ooooaa(icc)
     k = h3_ooooaa(icc)
     l = h4_ooooaa(icc)
     ! 6-1
     i1aa(i,j,k,l) = i1aa(i,j,k,l) + int2x(i,j,k,l,spin_int2aa)

     ! 6-2
     do a = norb1+1,nact
        do b = norb1+1,a-1
           i1aa(i,j,k,l) = i1aa(i,j,k,l) &
                + t2inp(a,b,k,l,spin_t2aa)*int2x(i,j,a,b,spin_int2aa)
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
     ! 6-1
     i1ab(i,j,k,l) = i1ab(i,j,k,l) + int2x(i,j,k,l,spin_int2ab)

     ! 6-2
     do a = norb1+1,nact
        do b = norb1+1,nact
           i1ab(i,j,k,l) = i1ab(i,j,k,l) &
                + t2inp(a,b,k,l,spin_t2ab)*int2x(i,j,a,b,spin_int2ab)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_man06_1
!##########################################################
