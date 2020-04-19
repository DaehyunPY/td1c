!##########################################################
subroutine ccdt_den1p_man05(i0,work1,work2)

!1: i0 ( a b )_yt + = +1/2 * Sum ( i j c ) 
!  * y ( i j c b )_y * t ( c a i j )_t 0
!5: i0 ( a b )_yt + = +1/12 * Sum ( i j k c d ) 
!  * y ( i j k c d b )_y * t ( c d a i j k )_t 0

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
     ! diagram 1
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

     ! diagram 5
     do i = 1,norb1
     do c = norb1+1,nact
        do d = norb1+1,c-1
        do j = 1,i-1
        do k = 1,j-1
           i0(a,b) = i0(a,b) + g3inp(i,j,k,c,d,b,spin_g3aaa)*t3inp(c,d,a,i,j,k,spin_t3aaa)
        end do
        end do
        end do

        do d = norb1+1,c-1
        do j = 1,norb1
        do k = 1,j-1
           i0(a,b) = i0(a,b) + g3inp(j,k,i,c,d,b,spin_g3aab)*t3inp(c,d,a,j,k,i,spin_t3aab)
        end do
        end do
        end do

        do d = norb1+1,nact
        do j = 1,norb1
        do k = 1,j-1
           i0(a,b) = i0(a,b) + g3inp(k,j,i,c,b,d,spin_g3aab)*t3inp(c,a,d,k,j,i,spin_t3aab)
        end do
        end do
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den1p_man05
!##########################################################
