!##########################################################
subroutine ccdt_l3p_man02(i0,work1,work2,work3)

! i0 ( i j k a b c )_yv + = -1 * P( 9 ) * Sum ( l ) * y ( i l a b )_y * v ( j k l c )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: work1(1),work2(1),work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

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
     do l = 1, norb1
        i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
             - g2inp(i,l,a,b,spin_g2aa)*int2x(j,k,l,c,spin_int2aa) &
             + g2inp(j,l,a,b,spin_g2aa)*int2x(i,k,l,c,spin_int2aa) &
             + g2inp(k,l,a,b,spin_g2aa)*int2x(j,i,l,c,spin_int2aa) &
             + g2inp(i,l,c,b,spin_g2aa)*int2x(j,k,l,a,spin_int2aa) &
             + g2inp(i,l,a,c,spin_g2aa)*int2x(j,k,l,b,spin_int2aa) &
             - g2inp(j,l,c,b,spin_g2aa)*int2x(i,k,l,a,spin_int2aa) &
             - g2inp(j,l,a,c,spin_g2aa)*int2x(i,k,l,b,spin_int2aa) &
             - g2inp(k,l,c,b,spin_g2aa)*int2x(j,i,l,a,spin_int2aa) &
             - g2inp(k,l,a,c,spin_g2aa)*int2x(j,i,l,b,spin_int2aa)
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
     do l = 1, norb1
        i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
             - g2inp(i,l,a,b,spin_g2aa)*int2x(j,k,l,c,spin_int2ab) &
             + g2inp(j,l,a,b,spin_g2aa)*int2x(i,k,l,c,spin_int2ab) &
!ZERO        + g2inp(k,l,a,b,spin_g2aa)*int2x(j,i,l,c,spin_int2aa) &
             + g2inp(i,l,b,c,spin_g2ab)*int2x(j,k,a,l,spin_int2ab) &
             - g2inp(i,l,a,c,spin_g2ab)*int2x(j,k,b,l,spin_int2ab) &
             - g2inp(j,l,b,c,spin_g2ab)*int2x(i,k,a,l,spin_int2ab) &
             + g2inp(j,l,a,c,spin_g2ab)*int2x(i,k,b,l,spin_int2ab) &
             - g2inp(k,l,c,b,spin_g2ab)*int2x(j,i,l,a,spin_int2aa) &
             + g2inp(l,k,a,c,spin_g2ab)*int2x(j,i,l,b,spin_int2aa)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man02
!##########################################################
