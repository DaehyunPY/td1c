!##########################################################
subroutine ccdt_amat_AABB(i0,work1,work2,work3)

! i0 ( a2 a1 | b1 a2 )_yt + = 1/2 * Sum ( j k c )
!  * y ( j k a1 c )_y * t ( b1 c j k )_t 0

! i0 ( a1 a2 | b2 a1 )_yt + = 1/2 * Sum ( j k c )
!  * y ( j k a2 c )_y * t ( b2 c j k )_t 0

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

  !$omp do collapse(2)
  do a1 = norb1+1,act1_ul
  do b1 = norb1+1,act1_ul
     tmp1 = 0d0
     do c = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           tmp1 = tmp1 + g2inp(j,k,a1,c,spin_g2aa)*t2inp(b1,c,j,k,spin_t2aa)
        end do
        do k = 1,norb1
           tmp1 = tmp1 + g2inp(j,k,a1,c,spin_g2ab)*t2inp(b1,c,j,k,spin_t2ab)
        end do
     end do
     end do
     do a2 = act1_ul+1,nact
        i0(a2,a1,b1,a2,1) = i0(a2,a1,b1,a2,1) + tmp1
     end do
  end do
  end do
  !$omp end do

  !$omp do collapse(2)
  do a2 = act1_ul+1,nact
  do b2 = act1_ul+1,nact
     tmp1 = 0d0
     do c = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           tmp1 = tmp1 + g2inp(j,k,a2,c,spin_g2aa)*t2inp(b2,c,j,k,spin_t2aa)
        end do
        do k = 1,norb1
           tmp1 = tmp1 + g2inp(j,k,a2,c,spin_g2ab)*t2inp(b2,c,j,k,spin_t2ab)
        end do
     end do
     end do
     do a1 = norb1+1,act1_ul
        i0(a1,a2,b2,a1,1) = i0(a1,a2,b2,a1,1) + tmp1
     end do
  end do
  end do
  !$omp end do

  !$omp end parallel
end subroutine ccdt_amat_AABB
!##########################################################
