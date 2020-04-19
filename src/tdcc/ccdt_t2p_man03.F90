!##########################################################
subroutine ccdt_t2p_man03(i0,work1,work2,work3)

!5: i0 ( a b i j )_vt + = -1 * P( 4 ) * Sum ( k c ) 
!  * t ( a c i k )_t * i1 ( k b j c )_v 2

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,t3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       work2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       work3(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  complex(kind(0d0)) :: tmp

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  work3 = 0d0
  call ccdt_t2p_man03_1(work1,work2,work3)

  !$omp parallel default(shared) private(a,b,i,j,tmp)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     tmp = 0d0
     ! diagram 5
     do k = 1,norb1
     do c = norb1+1,nact
        tmp = tmp - t2inp(a,c,i,k,spin_t2aa)*work1(k,b,j,c) &
                  - t2inp(a,c,i,k,spin_t2ab)*work3(k,b,j,c) &
                  + t2inp(a,c,j,k,spin_t2aa)*work1(k,b,i,c) &
                  + t2inp(a,c,j,k,spin_t2ab)*work3(k,b,i,c) &
                  + t2inp(b,c,i,k,spin_t2aa)*work1(k,a,j,c) &
                  + t2inp(b,c,i,k,spin_t2ab)*work3(k,a,j,c) &
                  - t2inp(b,c,j,k,spin_t2aa)*work1(k,a,i,c) &
                  - t2inp(b,c,j,k,spin_t2ab)*work3(k,a,i,c)
     end do
     end do
     i0(a,b,i,j,1) = i0(a,b,i,j,1) + tmp
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     tmp = 0d0
     ! diagram 5
     do k = 1,norb1
     do c = norb1+1,nact
        tmp = tmp - t2inp(a,c,i,k,spin_t2aa)*work3(k,b,j,c) &
                  - t2inp(a,c,i,k,spin_t2ab)*work1(k,b,j,c) &
                  - t2inp(a,c,k,j,spin_t2ab)*work2(k,b,i,c) &
                  - t2inp(c,b,i,k,spin_t2ab)*work2(k,a,j,c) &
                  - t2inp(b,c,j,k,spin_t2ab)*work1(k,a,i,c) &
                  - t2inp(b,c,j,k,spin_t2aa)*work3(k,a,i,c)
     end do
     end do
     i0(a,b,i,j,2) = i0(a,b,i,j,2) + tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t2p_man03
!##########################################################
subroutine ccdt_t2p_man03_1(i1aa,i1ab,i1ba)

! i1 ( i a j b )_vt + = -1/2 * Sum ( k c ) 
!  * t ( a c j k )_t * v ( i k b c )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact),&
       i1ab(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact),&
       i1ba(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  complex(kind(0d0)) :: tmp_aa,tmp_ab,tmp_ba

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp_aa,tmp_ab,tmp_ba)
  !$omp do
  do icc = 1, novovab
     i = h1_ovovab(icc)
     a = p2_ovovab(icc)
     j = h3_ovovab(icc)
     b = p4_ovovab(icc)
     tmp_aa = 0d0
     tmp_ab = 0d0
     tmp_ba = 0d0
     do k = 1,norb1
     do c = norb1+1,nact
        tmp_aa = tmp_aa &
             + t2inp(a,c,j,k,spin_t2aa)*int2x(i,k,b,c,spin_int2aa) &
             + t2inp(a,c,j,k,spin_t2ab)*int2x(i,k,b,c,spin_int2ab)
        tmp_ab = tmp_ab &
             + t2inp(c,a,j,k,spin_t2ab)*int2x(i,k,c,b,spin_int2ab)
        tmp_ba = tmp_ba &
             + t2inp(a,c,j,k,spin_t2ab)*int2x(i,k,b,c,spin_int2aa) &
             + t2inp(a,c,j,k,spin_t2aa)*int2x(i,k,b,c,spin_int2ab)
     end do
     end do
     i1aa(i,a,j,b) = i1aa(i,a,j,b) + int2x(i,a,j,b,spin_int2aa) - tmp_aa*0.5d0
     i1ab(i,a,j,b) = i1ab(i,a,j,b) + int2x(i,a,j,b,spin_int2ab) - tmp_ab*0.5d0
     i1ba(i,a,j,b) = i1ba(i,a,j,b) - int2x(i,a,b,j,spin_int2ab) - tmp_ba*0.5d0
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t2p_man03_1
!##########################################################
