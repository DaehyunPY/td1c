!##########################################################
subroutine ccdt_l2p_man18(i0,work1,work2,work3)

!18: i0 ( i j a b )_ytv + = 1/4 * P( 4 ) * Sum ( c k ) * i1 ( i c k a )_yt * v ( j k b c )_v 1

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work3(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  work3 = 0d0
  call ccdt_l2p_man18_1(work1,work2,work3)

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     do c = norb1+1,nact
     do k = 1,norb1
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             + work1(i,c,k,a)*int2x(j,k,b,c,spin_int2aa) &
             + work3(i,c,k,a)*int2x(j,k,b,c,spin_int2ab) &
             - work1(j,c,k,a)*int2x(i,k,b,c,spin_int2aa) &
             - work3(j,c,k,a)*int2x(i,k,b,c,spin_int2ab) &
             - work1(i,c,k,b)*int2x(j,k,a,c,spin_int2aa) &
             - work3(i,c,k,b)*int2x(j,k,a,c,spin_int2ab) &
             + work1(j,c,k,b)*int2x(i,k,a,c,spin_int2aa) &
             + work3(j,c,k,b)*int2x(i,k,a,c,spin_int2ab)
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
     do c = norb1+1,nact
     do k = 1,norb1
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             + work1(i,c,k,a)*int2x(j,k,b,c,spin_int2ab) &
             + work3(i,c,k,a)*int2x(j,k,b,c,spin_int2aa) &
             + work2(j,c,k,a)*int2x(i,k,c,b,spin_int2ab) &
             + work2(i,c,k,b)*int2x(k,j,a,c,spin_int2ab) &
             + work3(j,c,k,b)*int2x(i,k,a,c,spin_int2aa) &
             + work1(j,c,k,b)*int2x(i,k,a,c,spin_int2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man18
!##########################################################
subroutine ccdt_l2p_man18_1(i1aa,i1ab,i1ba)

! i1 ( i a j b )_yt + = 1 * Sum ( k l c d ) 
!  * t ( c d a k l j )_t * y ( i k l b c d )_y 0
!
! Demanding: ORDER-8
!
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       i1ab(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       i1ba(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novovab
     i = h1_ovovab(icc)
     a = p2_ovovab(icc)
     j = h3_ovovab(icc)
     b = p4_ovovab(icc)
     do k = 1,norb1
     do l = 1,k-1
     do c = norb1+1,nact
     do d = norb1+1,c-1
        i1aa(i,a,j,b) = i1aa(i,a,j,b) &
             + t3inp(c,d,a,k,l,j,spin_t3aaa)*g3inp(i,k,l,b,c,d,spin_g3aaa) &
             + t3inp(c,d,a,k,l,j,spin_t3aab)*g3inp(k,l,i,c,d,b,spin_g3aab)
        i1ba(i,a,j,b) = i1ba(i,a,j,b) &
             + t3inp(c,d,a,k,l,j,spin_t3aab)*g3inp(i,k,l,b,c,d,spin_g3aaa) &
             + t3inp(c,d,a,k,l,j,spin_t3aaa)*g3inp(k,l,i,c,d,b,spin_g3aab)
     end do
     end do
     end do
     end do

     do k = 1,norb1
     do l = 1,norb1
     do c = norb1+1,nact
     do d = norb1+1,nact
        i1aa(i,a,j,b) = i1aa(i,a,j,b) &
             + t3inp(a,c,d,j,k,l,spin_t3aab)*g3inp(i,k,l,b,c,d,spin_g3aab)
        i1ba(i,a,j,b) = i1ba(i,a,j,b) &
             + t3inp(d,a,c,l,j,k,spin_t3aab)*g3inp(i,k,l,b,c,d,spin_g3aab)
     end do
     end do
     end do
     end do

     do k = 1,norb1
     do l = 1,norb1
     do c = norb1+1,nact
     do d = norb1+1,c-1
        i1ab(i,a,j,b) = i1ab(i,a,j,b) &
             + t3inp(c,d,a,k,j,l,spin_t3aab)*g3inp(i,k,l,d,c,b,spin_g3aab)
     end do
     end do
     end do
     end do

     do k = 1,norb1
     do l = 1,k-1
     do c = norb1+1,nact
     do d = norb1+1,nact
        i1ab(i,a,j,b) = i1ab(i,a,j,b) &
             + t3inp(d,a,c,k,l,j,spin_t3aab)*g3inp(k,l,i,d,b,c,spin_g3aab)
     end do
     end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man18_1
!##########################################################
