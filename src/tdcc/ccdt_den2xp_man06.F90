!##########################################################
subroutine ccdt_den2xp_man06(i0,work1,work2,work3)

! i0 ( a b i j )_ytt + = Sum ( k l ) 
!  * i1 ( k l i j )_yt * t ( a b k l )_t 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,1:norb1,1:norb1,1:norb1), &
       work2(1:norb1,1:norb1,1:norb1,1:norb1), &
       work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den2xp_man06_1(work1,work2)

  !$omp parallel default(shared) private(i,j,k,l,a,b,c,d)
  !$omp do collapse(4)
  do a = norb1+1,nact
  do b = norb1+1,nact
  do i = 1,norb1
  do j = 1,norb1
     do k = 1,norb1
        do l = 1,k-1
           i0(a,b,i,j,1) = i0(a,b,i,j,1) &
                + work1(k,l,i,j)*t2inp(a,b,k,l,spin_t2aa)
        end do
     end do
     do k = 1,norb1
        do l = 1,norb1
           i0(a,b,i,j,2) = i0(a,b,i,j,2) &
                + work2(k,l,i,j)*t2inp(a,b,k,l,spin_t2ab)
        end do
     end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2xp_man06
!##########################################################
subroutine ccdt_den2xp_man06_1(i1,i2)

! i1 ( k l i j )_yt + = Sum ( c d ) 
!  * y ( k l c d )_y * t ( c d i j )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) ::&
       i1(1:norb1,1:norb1,1:norb1,1:norb1), &
       i2(1:norb1,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared)
  !$omp do collapse(4)
  do k = 1,norb1
  do l = 1,norb1
  do i = 1,norb1
  do j = 1,norb1
     do c = norb1+1,nact
        do d = norb1+1,c-1
           i1(k,l,i,j) = i1(k,l,i,j) &
                + g2inp(k,l,c,d,spin_g2aa)*t2inp(c,d,i,j,spin_t2aa)
        end do
     end do
     do c = norb1+1,nact
        do d = norb1+1,nact
           i2(k,l,i,j) = i2(k,l,i,j) &
                + g2inp(k,l,c,d,spin_g2ab)*t2inp(c,d,i,j,spin_t2ab)
        end do
     end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2xp_man06_1
!##########################################################
