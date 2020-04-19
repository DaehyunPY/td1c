!##########################################################
subroutine ccdt_den2p_man08(i0ab,i0ba,work1,work2,work3)

!  i0 ( i a b j )_yt + = +1 * Sum ( k c ) * y ( k i c b )_y * t ( c a k j )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,tmp)
  !$omp do
  do icc = 1, novovab
     i = h1_ovovab(icc)
     a = p2_ovovab(icc)
     b = p4_ovovab(icc)
     j = h3_ovovab(icc)
     do k = 1,norb1
     do c = norb1+1,nact
        i0ab(i,a,b,j) = i0ab(i,a,b,j) &
             + g2inp(k,i,c,b,spin_g2aa)*t2inp(c,a,k,j,spin_t2ab) &
             + g2inp(k,i,c,b,spin_g2ab)*t2inp(c,a,k,j,spin_t2aa)
        i0ba(i,a,b,j) = i0ba(i,a,b,j) &
             + g2inp(i,k,c,b,spin_g2ab)*t2inp(c,a,j,k,spin_t2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man08
!##########################################################
