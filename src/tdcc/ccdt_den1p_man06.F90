!##########################################################
subroutine ccdt_den1p_man06(i0,work1,work2)

!2: i0 ( i j )_yt + = -1/2 * Sum ( k a b ) 
!  * y ( k i a b )_y * t ( a b k j )_t 0
!6: i0 ( i j )_yt + = -1/12 * Sum ( k l a b c ) 
!  * y ( k l i a b c )_y * t ( a b c k l j )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, noo
     i = h1_oo(icc)
     j = h2_oo(icc)
     ! diagram 2
     do k = 1,norb1
        do a = norb1+1,nact
           do b = norb1+1,a-1
              i0(i,j) = i0(i,j) - g2inp(i,k,a,b,spin_g2aa)*t2inp(a,b,j,k,spin_t2aa)
           end do
           do b = norb1+1,nact
              i0(i,j) = i0(i,j) - g2inp(i,k,a,b,spin_g2ab)*t2inp(a,b,j,k,spin_t2ab)
           end do
        end do
     end do

     ! diagram 6
     do a = norb1+1,nact
     do k = 1,norb1
        do l = 1,k-1
        do b = norb1+1,a-1
        do c = norb1+1,b-1
           i0(i,j) = i0(i,j) - g3inp(k,l,i,a,b,c,spin_g3aaa)*t3inp(a,b,c,k,l,j,spin_t3aaa)
        end do
        end do
        end do

        do l = 1,k-1
        do b = norb1+1,nact
        do c = norb1+1,b-1
           i0(i,j) = i0(i,j) - g3inp(k,l,i,b,c,a,spin_g3aab)*t3inp(b,c,a,k,l,j,spin_t3aab)
        end do
        end do
        end do

        do l = 1,norb1
        do b = norb1+1,nact
        do c = norb1+1,b-1
           i0(i,j) = i0(i,j) - g3inp(k,i,l,c,b,a,spin_g3aab)*t3inp(c,b,a,k,j,l,spin_t3aab)
        end do
        end do
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den1p_man06
!##########################################################
