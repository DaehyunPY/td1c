!##########################################################
subroutine ccdt_den1p_man01(i0,work1,work2)

! i0 ( a b )_yt + = +1/2 * Sum ( i j c ) 
!  * y ( i j c b )_y * t ( c a i j )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, nvv
     a = p1_vv(icc)
     b = p2_vv(icc)
     do c = norb1+1,nact
        do i = 1,norb1
           do j = 1,i-1
              i0(a,b) = i0(a,b) + g2inp(i,j,b,c,spin_g2aa)*t2inp(a,c,i,j,spin_t2aa)
           end do
           do j = 1,norb1
              i0(a,b) = i0(a,b) + g2inp(i,j,b,c,spin_g2ab)*t2inp(a,c,i,j,spin_t2ab)
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den1p_man01
!##########################################################
