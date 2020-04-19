!##########################################################
subroutine ccdt_l3p_man04(i0,work1,work2,work3)

! i0 ( i j k a b c )_yf + = -1 * P( 3 ) * Sum ( l ) 
!  * y ( i j l a b c )_y * i1 ( k l )_f 2

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: work1(1:norb1,1:norb1),work2(1),work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
!  call ccdt_l3p_4_1(1,1,work1)
!  call ccdt_l3p_4_2(1,1,work1)
  call ccdt_l3p_man04_1(work1)

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
     do l = act1_ll,norb1
        i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
             - g3inp(i,j,l,a,b,c,spin_g3aaa)*work1(k,l) &
             + g3inp(k,j,l,a,b,c,spin_g3aaa)*work1(i,l) &
             + g3inp(i,k,l,a,b,c,spin_g3aaa)*work1(j,l)
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
     do l = act1_ll,norb1
        i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
             - g3inp(i,j,l,a,b,c,spin_g3aab)*work1(k,l) &
             - g3inp(l,j,k,a,b,c,spin_g3aab)*work1(i,l) &
             - g3inp(i,l,k,a,b,c,spin_g3aab)*work1(j,l)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man04
!##########################################################
subroutine ccdt_l3p_man04_1(i1)

! 4-1
! i1 ( i j )_f + = 1 * f ( i j )_f 0
! 4-2
! i1 ( i j )_vt + = -1/2 * Sum ( k a b ) * t ( a b k j )_t * v ( i k a b )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : t2inp,fock,int2x,norb1
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, noo
     i = h1_oo(icc)
     j = h2_oo(icc)
     ! 4-1
     i1(i,j) = i1(i,j) + fock(i,j,spin_focka)

     ! 4-2
     do k = 1,norb1
     do a = norb1+1,nact
        do b = norb1+1,a-1
           i1(i,j) = i1(i,j) &
                - t2inp(a,b,k,j,spin_t2aa)*int2x(i,k,a,b,spin_int2aa)
        end do
        do b = norb1+1,nact
           i1(i,j) = i1(i,j) &
                + t2inp(a,b,j,k,spin_t2ab)*int2x(i,k,a,b,spin_int2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_man04_1
!##########################################################
