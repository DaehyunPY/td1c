!##########################################################
subroutine ccdt_den2xp_man04(i0,work1,work2)

! i0 ( a b i j )_ytt + = -1/4 * Sum ( k l c ) 
!  * i1 ( k l c i )_yt * t ( b c a j k l )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1,(norb1+1):nact,1:norb1),&
       work2(1:norb1,1:norb1,(norb1+1):nact,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den2xp_man04_1(work1,work2)

  !$omp parallel default(shared) private(i,j,k,l,a,b,c,d)
  !$omp do collapse(4)
  do a = norb1+1,nact
  do b = norb1+1,nact
  do i = 1,norb1
  do j = 1,norb1
     do c = norb1+1,nact
     do k = 1,norb1
        do l = 1,k-1
           i0(a,b,i,j,1) = i0(a,b,i,j,1) &
                - work1(k,l,c,i)*t3inp(b,c,a,j,k,l,spin_t3aaa)
           i0(a,b,i,j,2) = i0(a,b,i,j,2) &
                - work1(k,l,c,i)*t3inp(a,c,b,l,k,j,spin_t3aab)
        end do
        do l = 1,norb1
           i0(a,b,i,j,1) = i0(a,b,i,j,1) &
                - work2(k,l,c,i)*t3inp(b,a,c,j,l,k,spin_t3aab)
           i0(a,b,i,j,2) = i0(a,b,i,j,2) &
                - work2(k,l,c,i)*t3inp(b,c,a,j,k,l,spin_t3aab)
        end do
     end do
     end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2xp_man04
!##########################################################
subroutine ccdt_den2xp_man04_1(i1aa,i1ab)

! i1 ( i j a k )_yt + = +1 * Sum ( l b c ) 
!  * y ( i j l a b c )_y * t ( b c k l )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,(norb1+1):nact,1:norb1),&
       i1ab(1:norb1,1:norb1,(norb1+1):nact,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared)
  !$omp do collapse(4)
  do i = 1,norb1
  do j = 1,norb1
  do a = norb1+1,nact
  do k = 1,norb1
     do l = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i1aa(i,j,a,k) = i1aa(i,j,a,k) &
                + g3inp(i,j,l,a,b,c,spin_g3aaa)*t2inp(b,c,k,l,spin_t2aa)
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                + g3inp(j,l,i,b,c,a,spin_g3aab)*t2inp(b,c,k,l,spin_t2aa)
        end do
        do c = norb1+1,nact
           i1aa(i,j,a,k) = i1aa(i,j,a,k) &
                + g3inp(i,j,l,a,b,c,spin_g3aab)*t2inp(b,c,k,l,spin_t2ab)
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                + g3inp(i,l,j,a,c,b,spin_g3aab)*t2inp(b,c,k,l,spin_t2ab)
        end do
     end do
     end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2xp_man04_1
!##########################################################
