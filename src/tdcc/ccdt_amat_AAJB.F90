!##########################################################
subroutine ccdt_amat_AAJB(i0,work1,work2,work3)

! i0 ( a1 a2 | b1 j1 )_yt + = 1/2 * Sum ( c k l )
!  * y ( k l a2 c )_y * t ( b1 a1 c j1 k l )_t 0

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
  do a2 = act1_ul+1,nact
  do a1 = norb1+1,act1_ul
     tmp1 = 0d0
     tmp2 = 0d0
     do c = norb1+1,act1_ul
     do k = act1_ll,norb1
        do l = act1_ll,k-1
           tmp1 = tmp1 + g2inp(k,l,a2,c,spin_g2aa)*t3inp(b1,a1,c,j1,k,l,spin_t2aa)
           tmp2 = tmp2 + g2inp(k,l,a2,c,spin_g2aa)*t3inp(a1,c,b1,k,l,j1,spin_t2ab)
        end do
        do l = act1_ll,norb1
           tmp1 = tmp1 + g2inp(k,l,a2,c,spin_g2ab)*t3inp(b1,a1,c,j1,k,l,spin_t2ab)
           tmp2 = tmp2 + g2inp(k,l,a2,c,spin_g2ab)*t3inp(b1,c,a1,j1,l,k,spin_t2ab)
        end do
     end do
     end do
     i0(a1,a2,b1,j1,1) = i0(a1,a2,b1,j1,1) + tmp1
     i0(a1,a2,b1,j1,2) = i0(a1,a2,b1,j1,2) + tmp2
  end do
  end do
  end do
  end do
  !$omp end do

  !$omp end parallel
end subroutine ccdt_amat_AAJB
!##########################################################
