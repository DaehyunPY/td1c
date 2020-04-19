!##########################################################
subroutine ccdt_l2p_man12(i0,work1,work2,work3)

!12: i0 ( i j a b )_ytf + = 1/2 * P( 2 ) * Sum ( k ) * i1 ( i j k a )_yt * f ( k b )_f 1

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  complex(kind(0d0)),intent(inout) :: work2(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  complex(kind(0d0)),intent(inout) :: work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_l2p_man12_1(work1,work2)
  call tdcc_filloovoaa(work1)

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     do k = 1, norb1
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             + work1(i,j,a,k)*fock(k,b,spin_focka) &
             - work1(i,j,b,k)*fock(k,a,spin_focka)
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
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             + work2(i,j,a,k)*fock(k,b,spin_focka) &
             + work2(j,i,b,k)*fock(k,a,spin_focka)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man12
!##########################################################
subroutine ccdt_l2p_man12_1(i1aa,i1ab)

! i1 ( i j a k )_yt + = -1/2 * Sum ( l b c ) 
!  * t ( b c k l )_t * y ( i j l a b c )_y 0
!
! Demanding: ORDER-7
!
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,(norb1+1):nact,1:norb1), &
       i1ab(1:norb1,1:norb1,(norb1+1):nact,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooovaa
     i = h1_ooovaa(icc)
     j = h2_ooovaa(icc)
     a = p4_ooovaa(icc)
     k = h3_ooovaa(icc)
     do l = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i1aa(i,j,a,k) = i1aa(i,j,a,k) &
                - t2inp(b,c,k,l,spin_t2aa)*g3inp(i,j,l,a,b,c,spin_g3aaa)
        end do
        do c = norb1+1,nact
           i1aa(i,j,a,k) = i1aa(i,j,a,k) &
                - t2inp(b,c,k,l,spin_t2ab)*g3inp(i,j,l,a,b,c,spin_g3aab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nooovab
     i = h1_ooovab(icc)
     j = h2_ooovab(icc)
     a = p4_ooovab(icc)
     k = h3_ooovab(icc)
     do l = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                - t2inp(b,c,k,l,spin_t2aa)*g3inp(j,l,i,b,c,a,spin_g3aab)
        end do
        do c = norb1+1,nact
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                - t2inp(b,c,k,l,spin_t2ab)*g3inp(i,l,j,a,c,b,spin_g3aab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man12_1
!##########################################################
