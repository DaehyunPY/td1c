!##########################################################
subroutine ccdt_l2p_man13(i0,work1,work2,work3)

!13: i0 ( i j a b )_ytv + = -1/2 * Sum ( k c ) * i1 ( k c a b )_yt * v ( i j k c )_v 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_l2p_man13_1(work1,work2)
  call tdcc_fillovvvaa(work1)

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     do k = 1,norb1
     do c = norb1+1,nact
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             + work1(k,c,a,b)*int2x(i,j,k,c,spin_int2aa)
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     do k = 1,norb1
     do c = norb1+1,nact
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             + work2(k,c,a,b)*int2x(i,j,k,c,spin_int2ab) &
             + work2(k,c,b,a)*int2x(i,j,c,k,spin_int2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man13
!##########################################################
subroutine ccdt_l2p_man13_1(i1aa,i1ab)

! i1 ( i a b c )_yt + = 1/2 * Sum ( j k d ) 
!  * t ( a d j k )_t * y ( i j k b c d )_y 0
!
! Demanding: ORDER-7
!
  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       i1ab(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novvvaa
     i = h1_ovvvaa(icc)
     a = p2_ovvvaa(icc)
     b = p3_ovvvaa(icc)
     c = p4_ovvvaa(icc)
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + t2inp(a,d,j,k,spin_t2aa)*g3inp(i,j,k,b,c,d,spin_g3aaa)
        end do
        do k = 1,norb1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + t2inp(a,d,j,k,spin_t2ab)*g3inp(i,j,k,b,c,d,spin_g3aab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, novvvab
     i = h1_ovvvab(icc)
     a = p2_ovvvab(icc)
     b = p3_ovvvab(icc)
     c = p4_ovvvab(icc)
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                + t2inp(a,d,j,k,spin_t2aa)*g3inp(j,k,i,c,d,b,spin_g3aab)
        end do
        do k = 1,norb1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                + t2inp(a,d,j,k,spin_t2ab)*g3inp(i,k,j,b,d,c,spin_g3aab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man13_1
!##########################################################
