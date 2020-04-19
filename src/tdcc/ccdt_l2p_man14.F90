!##########################################################
subroutine ccdt_l2p_man14(i0,work1,work2,work3)

!14: i0 ( i j a b )_ytv + = -1/2 * P( 4 ) * Sum ( k l ) * i1 ( i k l a )_yt * v ( j l k b )_v 1

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
!  call ccdt_l2p_14_1(1,1,1,1,work1)
!  call ccdt_l2p_14_1(1,2,1,2,work2)
!  work1 = -0.5d0*work1
!  work2 = -0.5d0*work2
  call ccdt_l2p_man14_1(work1,work2)
  call tdcc_fillooovaa(work1)

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     do k = 1,norb1
     do l = 1,norb1
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             + work1(i,k,l,a)*int2x(j,l,k,b,spin_int2aa) &
             + work2(k,i,l,a)*int2x(j,l,b,k,spin_int2ab) &
             - work1(j,k,l,a)*int2x(i,l,k,b,spin_int2aa) &
             - work2(k,j,l,a)*int2x(i,l,b,k,spin_int2ab) &
             - work1(i,k,l,b)*int2x(j,l,k,a,spin_int2aa) &
             - work2(k,i,l,b)*int2x(j,l,a,k,spin_int2ab) &
             + work1(j,k,l,b)*int2x(i,l,k,a,spin_int2aa) &
             + work2(k,j,l,b)*int2x(i,l,a,k,spin_int2ab)
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
     do l = 1,norb1
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             + work2(k,i,l,a)*int2x(l,j,k,b,spin_int2aa) &
             + work1(k,i,l,a)*int2x(l,j,k,b,spin_int2ab) &
             - work2(j,k,l,a)*int2x(i,l,k,b,spin_int2ab) &
             - work2(i,k,l,b)*int2x(j,l,k,a,spin_int2ab) &
             + work2(k,j,l,b)*int2x(l,i,k,a,spin_int2aa) &
             + work1(k,j,l,b)*int2x(l,i,k,a,spin_int2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man14
!##########################################################
subroutine ccdt_l2p_man14_1(i1aa,i1ab)

! i1 ( i j k a )_yt + = 1/2 * Sum ( l b c ) 
!  * t ( b c l k )_t * y ( i l j a b c )_y 0
!
! Demanding: ORDER-7
!
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,1:norb1,(norb1+1):nact), &
       i1ab(1:norb1,1:norb1,1:norb1,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooovaa
     i = h1_ooovaa(icc)
     j = h2_ooovaa(icc)
     k = h3_ooovaa(icc)
     a = p4_ooovaa(icc)
     do l = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i1aa(i,j,k,a) = i1aa(i,j,k,a) &
                + t2inp(b,c,l,k,spin_t2aa)*g3inp(i,l,j,a,b,c,spin_g3aaa)
        end do
        do c = norb1+1,nact
           i1aa(i,j,k,a) = i1aa(i,j,k,a) &
                + t2inp(b,c,l,k,spin_t2ab)*g3inp(i,j,l,a,c,b,spin_g3aab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nooovab
     i = h1_ooovab(icc)
     j = h2_ooovab(icc)
     k = h3_ooovab(icc)
     a = p4_ooovab(icc)
     do l = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i1ab(i,j,k,a) = i1ab(i,j,k,a) &
                + t2inp(b,c,l,k,spin_t2aa)*g3inp(i,l,j,b,c,a,spin_g3aab)
        end do
        do c = norb1+1,nact
           i1ab(i,j,k,a) = i1ab(i,j,k,a) &
                + t2inp(b,c,l,k,spin_t2ab)*g3inp(l,j,i,a,b,c,spin_g3aab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man14_1
!##########################################################
