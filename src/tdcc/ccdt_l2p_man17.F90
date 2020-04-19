!##########################################################
subroutine ccdt_l2p_man17(i0,work1,work2,work3)

!17: i0 ( i j a b )_vty + = 1/12 * Sum ( c d ) * i1 ( c d a b )_yt * v ( i j c d )_v 1

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work2((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  complex(kind(0d0)),intent(inout) :: work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_l2p_man17_1(work1,work2)
  call tdcc_fillvvvvaa(work1)

  !##################################################
  !$omp parallel default(shared) private(a,b,i,j)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     do c = norb1+1,nact
     do d = norb1+1,c-1
        i0(i,j,a,b,1) = i0(i,j,a,b,1) &
             + work1(c,d,a,b)*int2x(i,j,c,d,spin_int2aa)
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
     do d = norb1+1,nact
        i0(i,j,a,b,2) = i0(i,j,a,b,2) &
             + work2(c,d,a,b)*int2x(i,j,c,d,spin_int2ab)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_man17
!##########################################################
subroutine ccdt_l2p_man17_1(i1aa,i1ab)

! i1 ( c d a b )_yt + = 1 * Sum ( e i j k ) 
!  * y ( i j k a b e )_y * t ( e c d i j k )_t 0
!
! Demanding: ORDER-8
!
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       i1ab((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nvvvvaa
     c = p1_vvvvaa(icc)
     d = p2_vvvvaa(icc)
     a = p3_vvvvaa(icc)
     b = p4_vvvvaa(icc)
     do e = norb1+1,nact
        do i = 1,norb1
           do j = 1,i-1
              do k = 1,j-1
                 i1aa(c,d,a,b) = i1aa(c,d,a,b) &
                      + g3inp(i,j,k,a,b,e,spin_g3aaa)*t3inp(e,c,d,i,j,k,spin_t3aaa)
              end do
           end do
        end do

        do k = 1,norb1
           do i = 1,norb1
              do j = 1,i-1
                 i1aa(c,d,a,b) = i1aa(c,d,a,b) &
                      + g3inp(i,j,k,a,b,e,spin_g3aab)*t3inp(c,d,e,i,j,k,spin_t3aab)
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nvvvvab
     c = p1_vvvvab(icc)
     d = p2_vvvvab(icc)
     a = p3_vvvvab(icc)
     b = p4_vvvvab(icc)
     do e = norb1+1,nact
        do k = 1,norb1
           do i = 1,norb1
              do j = 1,i-1
                 i1ab(c,d,a,b) = i1ab(c,d,a,b) &
                      - g3inp(i,j,k,a,e,b,spin_g3aab)*t3inp(e,c,d,i,j,k,spin_t3aab) &
                      - g3inp(i,j,k,b,e,a,spin_g3aab)*t3inp(e,d,c,i,j,k,spin_t3aab)
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_l2p_man17_1
!##########################################################
