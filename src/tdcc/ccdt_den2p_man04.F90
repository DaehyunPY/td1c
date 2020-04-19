!##########################################################
subroutine ccdt_den2p_man04(i0ab,i0ba,work1,work2,work3)

!  i0 ( a b i j )_ytt + = P( i j ) * P( a b ) * Sum ( k c ) * i1 ( k a c i )_yt * t ( b c j k )_t 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1), &
       work2(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1), &
       work3(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)
  complex(kind(0d0)) :: tmp

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  work3 = 0d0
  call ccdt_den2p_man04_1(work1,work2,work3)

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,tmp)
  !$omp do
  do icc = 1, noovvab
     i = h1_oovvab(icc)
     j = h2_oovvab(icc)
     a = p3_oovvab(icc)
     b = p4_oovvab(icc)
     tmp = 0d0
     do k = 1,norb1
     do c = norb1+1,nact
        tmp = tmp &
             + work1(k,a,c,i)*t2inp(b,c,j,k,spin_t2ab) &
             + work2(k,a,c,i)*t2inp(b,c,j,k,spin_t2aa) &
             + work3(k,a,c,j)*t2inp(c,b,i,k,spin_t2ab) &
             + work3(k,b,c,i)*t2inp(a,c,k,j,spin_t2ab) &
             + work2(k,b,c,j)*t2inp(a,c,i,k,spin_t2aa) &
             + work1(k,b,c,j)*t2inp(a,c,i,k,spin_t2ab)
     end do
     end do
     i0ab(a,b,i,j) = i0ab(a,b,i,j) + tmp
     i0ba(a,b,j,i) = i0ba(a,b,j,i) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man04
!##########################################################
subroutine ccdt_den2p_man04_1(i1aa,i1ab,i1ba)

! i1 ( i a b j )_yt + = +1/2 * Sum ( k c ) * y ( k i c b )_y * t ( c a k j )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1), &
       i1ab(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1), &
       i1ba(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)
  complex(kind(0d0)) :: tmp_aa,tmp_ab,tmp_ba
  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(i,a,b,j,k,c,tmp_aa,tmp_ab,tmp_ba)
  !$omp do
  do icc = 1, novovab
     i = h1_ovovab(icc)
     a = p2_ovovab(icc)
     b = p4_ovovab(icc)
     j = h3_ovovab(icc)
     tmp_aa = 0d0
     tmp_ab = 0d0
     tmp_ba = 0d0
     do k = 1,norb1
     do c = norb1+1,nact
        tmp_aa = tmp_aa &
             + g2inp(k,i,c,b,spin_g2aa)*t2inp(c,a,k,j,spin_t2aa) &
             + g2inp(k,i,c,b,spin_g2ab)*t2inp(c,a,k,j,spin_t2ab)
        tmp_ab = tmp_ab &
             + g2inp(k,i,c,b,spin_g2aa)*t2inp(c,a,k,j,spin_t2ab) &
             + g2inp(k,i,c,b,spin_g2ab)*t2inp(c,a,k,j,spin_t2aa)
        tmp_ba = tmp_ba &
             + g2inp(i,k,c,b,spin_g2ab)*t2inp(c,a,j,k,spin_t2ab)
     end do
     end do
     i1aa(i,a,b,j) = i1aa(i,a,b,j) + tmp_aa*0.5d0
     i1ab(i,a,b,j) = i1ab(i,a,b,j) + tmp_ab*0.5d0
     i1ba(i,a,b,j) = i1ba(i,a,b,j) + tmp_ba*0.5d0
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man04_1
!##########################################################
