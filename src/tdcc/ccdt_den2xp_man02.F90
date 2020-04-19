!##########################################################
subroutine ccdt_den2xp_man02(i0,work1,work2)

! i0 ( a b i j )_ytt + = -1/2 * Sum ( k ) 
!  * i1 ( k j )_yt * t ( a b i k )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1),work2(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  call ccdt_den2xp_man02_1(work1)

  !$omp parallel default(shared) private(i,j,k,l,a,b,c,d)
  !$omp do collapse(4)
  do a = norb1+1,nact
  do b = norb1+1,nact
  do i = 1,norb1
  do j = 1,norb1
     do k = 1,norb1
        i0(a,b,i,j,1) = i0(a,b,i,j,1) &
             - work1(k,j)*t2inp(a,b,i,k,spin_t2aa)
        i0(a,b,i,j,2) = i0(a,b,i,j,2) &
             - work1(k,j)*t2inp(a,b,i,k,spin_t2ab)
     end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2xp_man02
!##########################################################
subroutine ccdt_den2xp_man02_1(i1)

! i1 ( i j )_yt + = +1 * Sum ( k a b ) 
!  * y ( i k a b )_y * t ( a b j k )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared)
  !$omp do collapse(2)
  do i = 1,norb1
  do j = 1,norb1
     do k = 1,norb1
     do a = norb1+1,nact
        do b = norb1+1,a-1
           i1(i,j) = i1(i,j) + &
                g2inp(i,k,a,b,spin_g2aa)*t2inp(a,b,j,k,spin_t2aa)
        end do
        do b = norb1+1,nact
           i1(i,j) = i1(i,j) + &
                g2inp(i,k,a,b,spin_g2ab)*t2inp(a,b,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2xp_man02_1
!##########################################################
