!##########################################################
subroutine ccdt_amat_IIJJ(i0,work1,work2,work3)

! i0 ( i2 i1 | j1 i2 )_yt + = 1/2 * Sum ( k b ) 
!  * y ( i2 k a b )_y * t ( a b j2 k )_t 0

! i0 ( i1 i2 | j2 i1 )_yt + = 1/2 * Sum ( k b ) 
!  * y ( i1 k a b )_y * t ( a b j1 k )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(out) :: work1(1),work2(1),work3(1)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  !$omp parallel default(shared) private(i1,i2,j1,j2,a1,a2,b1,b2,i,j,k,l,m,a,b,c,d,e,tmp1,tmp2)
  !$omp do collapse(2)
  do i2 = 1,act1_ll-1
  do j2 = 1,act1_ll-1
     tmp1 = 0d0
     do k = 1,norb1
     do a = norb1+1,nact
        do b = norb1+1,a-1
           tmp1 = tmp1 + g2inp(i2,k,a,b,spin_g2aa)*t2inp(a,b,j2,k,spin_t2aa)
        end do
        do b = norb1+1,nact
           tmp1 = tmp1 + g2inp(i2,k,a,b,spin_g2ab)*t2inp(a,b,j2,k,spin_t2ab)
        end do
     end do
     end do
     do i1 = act1_ll,norb1
        i0(i2,i1,i1,j2,1) = i0(i2,i1,i1,j2,1) + tmp1
     end do
  end do
  end do
  !$omp end do

  !$omp do collapse(2)
  do i1 = act1_ll,norb1
  do j1 = act1_ll,norb1
     tmp1 = 0d0
     do k = 1,norb1
     do a = norb1+1,nact
        do b = norb1+1,a-1
           tmp1 = tmp1 + g2inp(i1,k,a,b,spin_g2aa)*t2inp(a,b,j1,k,spin_t2aa)
        end do
        do b = norb1+1,nact
           tmp1 = tmp1 + g2inp(i1,k,a,b,spin_g2ab)*t2inp(a,b,j1,k,spin_t2ab)
        end do
     end do
     end do
     do i2 = 1,act1_ll-1
        i0(i1,i2,i2,j1,1) = i0(i1,i2,i2,j1,1) + tmp1
     end do
  end do
  end do
  !$omp end do

  !$omp end parallel
end subroutine ccdt_amat_IIJJ
!##########################################################
