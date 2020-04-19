!##########################################################
subroutine bccdt_bmat2p_man01(i0,work1,work2)

!  i0 ( i a )_ydt + = +1 * Sum ( j b ) 
!  * y ( j b )_y * dt ( a b i j )_dt 0

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_ormas,only : nact,ncore
  use mod_cc,only : norb1,g1inp,dt2inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1)

  integer(c_int) :: icc,jcc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,jcc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, nov
     i = h1_ov(icc)
     a = p2_ov(icc)
     do jcc = 1, nov
        j = h1_ov(jcc)
        b = p2_ov(jcc)
        i0(i,a) = i0(i,a) + g1inp(j,b,spin_g1a)*( &
             dt2inp(a,b,i,j,spin_t2aa) &
           + dt2inp(a,b,i,j,spin_t2ab))
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine bccdt_bmat2p_man01
!##########################################################
