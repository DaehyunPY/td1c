!##########################################################
subroutine ccdt_amat_IIJB(i0,work1,work2,work3)

! i0 ( i2 i1 | b1 j1 )_yt + = 1/2 * Sum ( k c d )
!  * y ( i2 k c d )_y * t ( b1 c d i1 j1 k )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  !$omp parallel default(shared) private(i1,i2,j1,j2,a1,a2,b1,b2,i,j,k,l,m,a,b,c,d,e,tmp1,tmp2)

  !$omp do collapse(4)
  do j1 = act1_ll,norb1
  do b1 = norb1+1,act1_ul
  do i1 = act1_ll,norb1
  do i2 = 1,act1_ll-1
     tmp1 = 0d0
     tmp2 = 0d0
     do k = act1_ll,norb1
     do c = norb1+1,act1_ul
        do d = norb1+1,c-1
           tmp1 = tmp1 - g2inp(i2,k,c,d,spin_g2aa)*t3inp(b1,c,d,j1,i1,k,spin_t2aa)
           tmp2 = tmp2 - g2inp(i2,k,c,d,spin_g2aa)*t3inp(c,d,b1,i1,k,j1,spin_t2ab)
        end do
        do d = norb1+1,act1_ul
           tmp1 = tmp1 - g2inp(i2,k,c,d,spin_g2ab)*t3inp(b1,c,d,j1,i1,k,spin_t2ab)
           tmp2 = tmp2 - g2inp(i2,k,c,d,spin_g2ab)*t3inp(b1,d,c,j1,k,i1,spin_t2ab)
        end do
     end do
     end do
     i0(i2,i1,b1,j1,1) = i0(i2,i1,b1,j1,1) + tmp1
     i0(i2,i1,b1,j1,2) = i0(i2,i1,b1,j1,2) + tmp2
  end do
  end do
  end do
  end do
  !$omp end do

  !$omp end parallel
end subroutine ccdt_amat_IIJB
!##########################################################
