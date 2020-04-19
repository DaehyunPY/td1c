!##########################################################
subroutine ccdt_l3p_1e_man01(i0)

! i0 ( i j k a b c )_yf + = 1 * P( 9 ) 
!  * y ( i j a b )_y * f ( k c )_f 0

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)

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
          + g2inp(i,j,a,b,spin_g2aa)*fock(k,c,spin_focka) &
          - g2inp(k,j,a,b,spin_g2aa)*fock(i,c,spin_focka) &
          - g2inp(i,k,a,b,spin_g2aa)*fock(j,c,spin_focka) &
          - g2inp(i,j,c,b,spin_g2aa)*fock(k,a,spin_focka) &
          - g2inp(i,j,a,c,spin_g2aa)*fock(k,b,spin_focka) &
          + g2inp(k,j,c,b,spin_g2aa)*fock(i,a,spin_focka) &
          + g2inp(i,k,c,b,spin_g2aa)*fock(j,a,spin_focka) &
          + g2inp(k,j,a,c,spin_g2aa)*fock(i,b,spin_focka) &
          + g2inp(i,k,a,c,spin_g2aa)*fock(j,b,spin_focka)
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
          + g2inp(i,j,a,b,spin_g2aa)*fock(k,c,spin_focka) &
          + g2inp(k,j,c,b,spin_g2ab)*fock(i,a,spin_focka) &
          - g2inp(i,k,b,c,spin_g2ab)*fock(j,a,spin_focka) &
          - g2inp(j,k,a,c,spin_g2ab)*fock(i,b,spin_focka) &
          + g2inp(i,k,a,c,spin_g2ab)*fock(j,b,spin_focka)
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_1e_man01
!##########################################################
