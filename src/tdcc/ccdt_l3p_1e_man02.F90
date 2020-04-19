!##########################################################
subroutine ccdt_l3p_1e_man02(i0)

! i0 ( i j k a b c )_yf + = -1 * P( 3 ) * Sum ( l ) 
!  * y ( i j l a b c )_y * f ( k l )_f 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)

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
     do l = 1,norb1
        i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
             - g3inp(i,j,l,a,b,c,spin_g3aaa)*fock(k,l,spin_focka) &
             + g3inp(k,j,l,a,b,c,spin_g3aaa)*fock(i,l,spin_focka) &
             + g3inp(i,k,l,a,b,c,spin_g3aaa)*fock(j,l,spin_focka)
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
        i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
             - g3inp(i,j,l,a,b,c,spin_g3aab)*fock(k,l,spin_focka) &
             - g3inp(l,j,k,a,b,c,spin_g3aab)*fock(i,l,spin_focka) &
             - g3inp(i,l,k,a,b,c,spin_g3aab)*fock(j,l,spin_focka)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_1e_man02
!##########################################################
