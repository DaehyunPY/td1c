!##########################################################
subroutine bccdt_l2p_man01(i0,work1,work2,work3)

!2: i0 ( i j a b )_yf + = 1 * P( 4 ) * y ( i a )_y * f ( j b )_f 0
!3: i0 ( i j a b )_yv + = -1 * P( 2 ) * Sum ( k ) * y ( k a )_y * v ( i j k b )_v 0
!4: i0 ( i j a b )_yv + = -1 * P( 2 ) * Sum ( c ) * y ( i c )_y * v ( j c a b )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,g1inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     ! diagram 2
     i0(i,j,a,b,1) = i0(i,j,a,b,1) &
          + g1inp(i,a,spin_g1a)*fock(j,b,spin_focka) &
          - g1inp(j,a,spin_g1a)*fock(i,b,spin_focka) &
          - g1inp(i,b,spin_g1a)*fock(j,a,spin_focka) &
          + g1inp(j,b,spin_g1a)*fock(i,a,spin_focka)
     ! diagram 3
     do k = 1,norb1
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             - g1inp(k,a,spin_g1a)*int2x(i,j,k,b,spin_int2aa) &
             + g1inp(k,b,spin_g1a)*int2x(i,j,k,a,spin_int2aa)
     end do
     ! diagram 4
     do c = norb1+1,nact
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             - g1inp(i,c,spin_g1a)*int2x(j,c,a,b,spin_int2aa) &
             + g1inp(j,c,spin_g1a)*int2x(i,c,a,b,spin_int2aa)
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     ! diagram 2
     i0(i,j,a,b,2) = i0(i,j,a,b,2) &
          + g1inp(i,a,spin_g1a)*fock(j,b,spin_focka) &
          + g1inp(j,b,spin_g1a)*fock(i,a,spin_focka)
     ! diagram 3
     do k = 1,norb1
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             - g1inp(k,a,spin_g1a)*int2x(i,j,k,b,spin_int2ab) &
             - g1inp(k,b,spin_g1a)*int2x(i,j,a,k,spin_int2ab)
     end do
     ! diagram 4
     do c = norb1+1,nact
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             + g1inp(i,c,spin_g1a)*int2x(c,j,a,b,spin_int2ab) &
             + g1inp(j,c,spin_g1a)*int2x(i,c,a,b,spin_int2ab)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine bccdt_l2p_man01
!##########################################################
