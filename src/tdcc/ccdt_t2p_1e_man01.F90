!##########################################################
subroutine ccdt_t2p_1e_man01(i0)

!i0 ( a b i j )_tf + = -1 * P( 2 ) * Sum ( k ) 
!  * t ( a b i k )_t * f ( k j )_f 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:2)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     do k = 1,norb1
        i0(a,b,i,j,1) = i0(a,b,i,j,1) &
             - t2inp(a,b,i,k,spin_t2aa)*fock(k,j,spin_focka) &
             + t2inp(a,b,j,k,spin_t2aa)*fock(k,i,spin_focka)
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     do k = 1,norb1
        i0(a,b,i,j,2) = i0(a,b,i,j,2) &
             - t2inp(a,b,i,k,spin_t2ab)*fock(k,j,spin_focka) &
             - t2inp(a,b,k,j,spin_t2ab)*fock(k,i,spin_focka)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_t2p_1e_man01
!##########################################################
