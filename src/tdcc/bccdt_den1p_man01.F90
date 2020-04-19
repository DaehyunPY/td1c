!##########################################################
subroutine bccdt_den1p_man01(i0,work1,work2)

!1:  i0 ( i a )_y + = +1 * y ( i a )_y 0
!2:  i0 ( a i )_yt + = +1 * Sum ( j b ) * y ( j b )_y * t ( b a j i )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g1inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, nov
     a = p2_ov(icc)
     i = h1_ov(icc)
     ! diagram 1
     i0(i,a) = i0(i,a) + g1inp(i,a,spin_g1a)

     ! diagram 2
     do j = 1,norb1
     do b = norb1+1,nact
        i0(a,i) = i0(a,i) &
             + g1inp(j,b,spin_g1a)*t2inp(b,a,j,i,spin_t2aa) &
             + g1inp(j,b,spin_g1a)*t2inp(b,a,j,i,spin_t2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine bccdt_den1p_man01
!##########################################################
