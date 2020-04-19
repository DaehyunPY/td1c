!##########################################################
subroutine bccdt_l3p_man01(i0,work1,work2,work3)

!1: i0 ( i j k a b c )_yv + = 1 * P( 9 ) * y ( i a )_y * v ( j k b c )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g1inp,t3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: work1(1),work2(1),work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

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
     i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
          + g1inp(i,a,spin_g1a)*int2x(j,k,b,c,spin_int2aa) &
          - g1inp(j,a,spin_g1a)*int2x(i,k,b,c,spin_int2aa) &
          - g1inp(k,a,spin_g1a)*int2x(j,i,b,c,spin_int2aa) &
          - g1inp(i,b,spin_g1a)*int2x(j,k,a,c,spin_int2aa) &
          - g1inp(i,c,spin_g1a)*int2x(j,k,b,a,spin_int2aa) &
          + g1inp(j,b,spin_g1a)*int2x(i,k,a,c,spin_int2aa) &
          + g1inp(k,b,spin_g1a)*int2x(j,i,a,c,spin_int2aa) &
          + g1inp(j,c,spin_g1a)*int2x(i,k,b,a,spin_int2aa) &
          + g1inp(k,c,spin_g1a)*int2x(j,i,b,a,spin_int2aa)
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
     i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
          + g1inp(i,a,spin_g1a)*int2x(j,k,b,c,spin_int2ab) &
          - g1inp(j,a,spin_g1a)*int2x(i,k,b,c,spin_int2ab) &
          - g1inp(i,b,spin_g1a)*int2x(j,k,a,c,spin_int2ab) &
          + g1inp(j,b,spin_g1a)*int2x(i,k,a,c,spin_int2ab) &
          + g1inp(k,c,spin_g1a)*int2x(j,i,b,a,spin_int2aa)
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine bccdt_l3p_man01
!##########################################################
