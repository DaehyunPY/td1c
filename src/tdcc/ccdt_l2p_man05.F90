!##########################################################
subroutine ccdt_l2p_man05(i0,work1,work2,work3)

!5:  i0 ( i j a b )_yv + = -1 * P( 4 ) * Sum ( k c ) * y ( i k a c )_y * i1 ( j c k b )_v 2

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
!  call ccdt_l2p_5_1(1,1,1,1,work1)
!  call ccdt_l2p_5_2(1,1,1,1,work1)
!  call ccdt_l2p_5_1(1,2,1,2,work2)
!  call ccdt_l2p_5_2(1,2,1,2,work2)
!  call ccdt_l2p_5_1(2,1,1,2,work3)
!  call ccdt_l2p_5_2(2,1,1,2,work3)
  call ccdt_l2p_man05_1(work1,work2,work3)

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
             - g2inp(i,k,a,c,spin_g2aa)*work1(j,c,k,b) &
             - g2inp(i,k,a,c,spin_g2ab)*work3(j,c,k,b) &
             + g2inp(j,k,a,c,spin_g2aa)*work1(i,c,k,b) &
             + g2inp(j,k,a,c,spin_g2ab)*work3(i,c,k,b) &
             + g2inp(i,k,b,c,spin_g2aa)*work1(j,c,k,a) &
             + g2inp(i,k,b,c,spin_g2ab)*work3(j,c,k,a) &
             - g2inp(j,k,b,c,spin_g2aa)*work1(i,c,k,a) &
             - g2inp(j,k,b,c,spin_g2ab)*work3(i,c,k,a)
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
             - g2inp(i,k,a,c,spin_g2aa)*work3(j,c,k,b) &
             - g2inp(i,k,a,c,spin_g2ab)*work1(j,c,k,b) &
             - g2inp(k,j,a,c,spin_g2ab)*work2(i,c,k,b) &
             - g2inp(i,k,c,b,spin_g2ab)*work2(j,c,k,a) &
             - g2inp(j,k,b,c,spin_g2ab)*work1(i,c,k,a) &
             - g2inp(j,k,b,c,spin_g2aa)*work3(i,c,k,a)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man05
!##########################################################
subroutine ccdt_l2p_man05_1(i1aa,i1ab,i1ba)

!5-1 i1 ( i a j b )_v + = 1 * v ( i a j b )_v 0
!5-2 i1 ( i a j b )_vt + = 1 * Sum ( k c ) * t ( a c k j )_t * v ( i k b c )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,fock,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact),&
       i1ab(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact),&
       i1ba(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novovab
     i = h1_ovovab(icc)
     a = p2_ovovab(icc)
     j = h3_ovovab(icc)
     b = p4_ovovab(icc)
     i1aa(i,a,j,b) = i1aa(i,a,j,b) + int2x(i,a,j,b,spin_int2aa)
     i1ab(i,a,j,b) = i1ab(i,a,j,b) + int2x(i,a,j,b,spin_int2ab)
     i1ba(i,a,j,b) = i1ba(i,a,j,b) - int2x(i,a,b,j,spin_int2ab)
     do k = 1,norb1
     do c = norb1+1,nact
        i1aa(i,a,j,b) = i1aa(i,a,j,b) &
             + t2inp(a,c,k,j,spin_t2aa)*int2x(i,k,b,c,spin_int2aa) &
             - t2inp(a,c,j,k,spin_t2ab)*int2x(i,k,b,c,spin_int2ab)
        i1ab(i,a,j,b) = i1ab(i,a,j,b) &
             - t2inp(a,c,k,j,spin_t2ab)*int2x(i,k,c,b,spin_int2ab)
        i1ba(i,a,j,b) = i1ba(i,a,j,b) &
             - t2inp(c,a,k,j,spin_t2ab)*int2x(i,k,b,c,spin_int2aa) &
             + t2inp(a,c,k,j,spin_t2aa)*int2x(i,k,b,c,spin_int2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man05_1
!##########################################################
