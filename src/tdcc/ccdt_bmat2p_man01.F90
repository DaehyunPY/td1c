!##########################################################
subroutine ccdt_bmat2p_man01(i0,work1,work2)

!  i0 ( i a )_ydt + = +1/4 * Sum ( j k b c ) 
!  * y ( j k b c )_y * dt ( a b c i j k )_dt 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,dt3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, nov
     i = h1_ov(icc)
     a = p2_ov(icc)
     do j = 1,norb1
     do b = norb1+1,nact
        do k = 1,j-1
        do c = norb1+1,b-1
           i0(i,a) = i0(i,a) + g2inp(j,k,b,c,spin_g2aa)*( &
                dt3inp(a,b,c,i,j,k,spin_t3aaa) &
              + dt3inp(b,c,a,j,k,i,spin_t3aab))
        end do
        end do
        do k = 1,norb1
        do c = norb1+1,nact
           i0(i,a) = i0(i,a) &
                + g2inp(j,k,b,c,spin_g2ab)*dt3inp(a,b,c,i,j,k,spin_t3aab)
        end do
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_bmat2p_man01
!##########################################################
