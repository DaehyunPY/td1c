!##########################################################
subroutine ccdt_l2p_1e_man02(i0)

! i0 ( i j a b )_yf + = 1 * P( 2 ) * Sum ( c ) 
!  * y ( i j a c )_y * f ( c b )_f 0

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     do c = norb1+1,nact
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             + g2inp(i,j,a,c,spin_g2aa)*fock(c,b,spin_focka) &
             - g2inp(i,j,b,c,spin_g2aa)*fock(c,a,spin_focka)
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     do c = norb1+1,nact
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             + g2inp(i,j,a,c,spin_g2ab)*fock(c,b,spin_focka) &
             + g2inp(i,j,c,b,spin_g2ab)*fock(c,a,spin_focka)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_1e_man02
!##########################################################
