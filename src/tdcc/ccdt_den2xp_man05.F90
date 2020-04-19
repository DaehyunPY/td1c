!##########################################################
subroutine ccdt_den2xp_man05(i0,work1,work2,work3)

! i0 ( a b i j )_ytt + = -1 * Sum ( k c ) 
!  * i1 ( k b c i )_yt * t ( a c j k )_t 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1), &
       work2(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1), &
       work3(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  work3 = 0d0
  call ccdt_den2xp_man05_1(work1,work2,work3)

  !$omp parallel default(shared) private(i,j,k,l,a,b,c,d)
  !$omp do collapse(4)
  do a = norb1+1,nact
  do b = norb1+1,nact
  do i = 1,norb1
  do j = 1,norb1
     do c = norb1+1,nact
     do k = 1,norb1
        i0(a,b,i,j,1) = i0(a,b,i,j,1) &
             - work1(k,b,c,i)*t2inp(a,c,j,k,spin_t2aa) &
             - work2(k,b,c,i)*t2inp(a,c,j,k,spin_t2ab)
        i0(a,b,i,j,2) = i0(a,b,i,j,2) &
             + work3(k,b,c,i)*t2inp(a,c,k,j,spin_t2ab)
     end do
     end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2xp_man05
!##########################################################
subroutine ccdt_den2xp_man05_1(i1,i2,i3)

! i1 ( k b c i )_yt + = Sum ( l d ) 
!  * y ( k l c d )_y * t ( b d i l )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) ::&
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1), &
       i2(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1), &
       i3(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared)
  !$omp do collapse(4)
  do k = 1,norb1
  do b = norb1+1,nact
  do c = norb1+1,nact
  do i = 1,norb1
     do l = 1,norb1
        do d = norb1+1,nact
           i1(k,b,c,i) = i1(k,b,c,i) &
                + g2inp(k,l,c,d,spin_g2aa)*t2inp(b,d,i,l,spin_t2aa) &
                + g2inp(k,l,c,d,spin_g2ab)*t2inp(b,d,i,l,spin_t2ab)
           i2(k,b,c,i) = i2(k,b,c,i) &
                + g2inp(k,l,c,d,spin_g2ab)*t2inp(b,d,i,l,spin_t2aa) &
                + g2inp(k,l,c,d,spin_g2aa)*t2inp(b,d,i,l,spin_t2ab)
           i3(k,b,c,i) = i3(k,b,c,i) &
                + g2inp(k,l,d,c,spin_g2ab)*t2inp(b,d,l,i,spin_t2ab)
        end do
     end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2xp_man05_1
!##########################################################
