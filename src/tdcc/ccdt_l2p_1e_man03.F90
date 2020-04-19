!##########################################################
subroutine ccdt_l2p_1e_man03(i0,work1,work2)

! i0 ( i j a b )_fty + = 1/2 * P( 2 ) * Sum ( c k l ) 
!  * y ( i k l a b c )_y * i1 ( j c k l )_ft 1

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,(norb1+1):nact,1:norb1,1:norb1),&
       work2(1:norb1,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_l2p_1e_man03_1(work1,work2)
  call tdcc_fillovooaa(work1)

!ok  work1 = 0d0
!ok  work2 = 0d0
!ok  call ccdt_l2_1e_3_1(1,1,1,1,work1)
!ok  call ccdt_l2_1e_3_1(1,2,1,2,work2)

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
        do l = 1,k-1
           i0(i,j,a,b,1) = i0(i,j,a,b,1) &
                + g3inp(i,k,l,a,b,c,spin_g3aaa)*work1(j,c,k,l) &
                - g3inp(j,k,l,a,b,c,spin_g3aaa)*work1(i,c,k,l)
        end do
        do l = 1,norb1
           i0(i,j,a,b,1) = i0(i,j,a,b,1) &
                + g3inp(i,k,l,a,b,c,spin_g3aab)*work2(j,c,k,l) &
                - g3inp(j,k,l,a,b,c,spin_g3aab)*work2(i,c,k,l)
        end do
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
        do l = 1,k-1
           i0(i,j,a,b,2) = i0(i,j,a,b,2) &
                + g3inp(k,l,i,b,c,a,spin_g3aab)*work1(j,c,k,l) &
                - g3inp(l,k,j,a,c,b,spin_g3aab)*work1(i,c,k,l)
        end do
        do l = 1,norb1
           i0(i,j,a,b,2) = i0(i,j,a,b,2) &
                + g3inp(i,l,k,a,c,b,spin_g3aab)*work2(j,c,k,l) &
                - g3inp(j,k,l,c,b,a,spin_g3aab)*work2(i,c,l,k)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l2p_1e_man03
!##########################################################
subroutine ccdt_l2p_1e_man03_1(i1aa,i1ab)

! i1 ( i a j k )_ft + = -1 * Sum ( b ) 
!  * t ( b a j k )_t * f ( i b )_f 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,fock
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,1:norb1,1:norb1), &
       i1ab(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooovaa
     i = h3_ooovaa(icc)
     a = p4_ooovaa(icc)
     j = h1_ooovaa(icc)
     k = h2_ooovaa(icc)
     do b = norb1+1,nact
        i1aa(i,a,j,k) = i1aa(i,a,j,k) &
             - t2inp(b,a,j,k,spin_t2aa)*fock(i,b,spin_focka)
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nooovab
     i = h3_ooovab(icc)
     a = p4_ooovab(icc)
     j = h1_ooovab(icc)
     k = h2_ooovab(icc)
     do b = norb1+1,nact
        i1ab(i,a,j,k) = i1ab(i,a,j,k) &
             - t2inp(b,a,j,k,spin_t2ab)*fock(i,b,spin_focka)
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l2p_1e_man03_1
!##########################################################
