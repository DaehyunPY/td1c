!##########################################################
subroutine ccdt_l2p_man06(i0,work1,work2,work3)

!1:  i0 ( i j a b )_v + = 1 * v ( i j a b )_v 0
!6:  i0 ( i j a b )_yv + = 1/2 * Sum ( c d ) * y ( i j c d )_y * v ( c d a b )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     ! diagram 1
     i0(i,j,a,b,1) = i0(i,j,a,b,1) + int2x(i,j,a,b,spin_int2aa)

     ! diagram 6
     do c = norb1+1,nact
     do d = norb1+1,c-1
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             + g2inp(i,j,c,d,spin_g2aa)*int2x(c,d,a,b,spin_int2aa)
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     ! diagram 1
     i0(i,j,a,b,2) = i0(i,j,a,b,2) + int2x(i,j,a,b,spin_int2ab)

     ! diagram 6
     do c = norb1+1,nact
     do d = norb1+1,nact
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             + g2inp(i,j,c,d,spin_g2ab)*int2x(c,d,a,b,spin_int2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man06
!##########################################################
