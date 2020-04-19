!##########################################################
subroutine ccdt_t3p_1e_man03(i0,work1,work2)

! i0 ( a b c i j k )_ftt + = P( 9 ) * Sum ( l ) 
!  * t ( a b i l )_t * i1 ( l c j k )_ft 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2)
  complex(kind(0d0)) :: &
       work1(1:norb1,(norb1+1):nact,1:norb1,1:norb1),&
       work2(1:norb1,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_t3p_1e_man03_1(work1,work2)
  call tdcc_fillovooaa(work1)

  !##################################################
  !$omp parallel default(shared) private(a,b,c,i,j,k)
  !$omp do
  do icc = 1,ncc3aaa
     a = p1_cc3aaa(icc)
     b = p2_cc3aaa(icc)
     c = p3_cc3aaa(icc)
     i = h1_cc3aaa(icc)
     j = h2_cc3aaa(icc)
     k = h3_cc3aaa(icc)
     do l = 1,norb1
        i0(a,b,c,i,j,k,1) = i0(a,b,c,i,j,k,1) &
             + t2inp(a,b,i,l,spin_t2aa)*work1(l,c,j,k) &
             - t2inp(a,b,j,l,spin_t2aa)*work1(l,c,i,k) &
             - t2inp(a,b,k,l,spin_t2aa)*work1(l,c,j,i) &
             - t2inp(c,b,i,l,spin_t2aa)*work1(l,a,j,k) &
             - t2inp(a,c,i,l,spin_t2aa)*work1(l,b,j,k) &
             + t2inp(c,b,j,l,spin_t2aa)*work1(l,a,i,k) &
             + t2inp(c,b,k,l,spin_t2aa)*work1(l,a,j,i) &
             + t2inp(a,c,j,l,spin_t2aa)*work1(l,b,i,k) &
             + t2inp(a,c,k,l,spin_t2aa)*work1(l,b,j,i)
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc3aab
     a = p1_cc3aab(icc)
     b = p2_cc3aab(icc)
     c = p3_cc3aab(icc)
     i = h1_cc3aab(icc)
     j = h2_cc3aab(icc)
     k = h3_cc3aab(icc)
     do l = 1,norb1
        i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
             + t2inp(a,b,i,l,spin_t2aa)*work2(l,c,j,k) &
             - t2inp(a,b,j,l,spin_t2aa)*work2(l,c,i,k) &
             !- t2inp(a,b,k,l,spin_t2aa)*work1(l,c,j,i) &
             - t2inp(b,c,i,l,spin_t2ab)*work2(l,a,k,j) &
             + t2inp(a,c,i,l,spin_t2ab)*work2(l,b,k,j) &
             + t2inp(b,c,j,l,spin_t2ab)*work2(l,a,k,i) &
             + t2inp(c,b,k,l,spin_t2ab)*work1(l,a,j,i) &
             - t2inp(a,c,j,l,spin_t2ab)*work2(l,b,k,i) &
             - t2inp(a,c,l,k,spin_t2ab)*work1(l,b,j,i)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_t3p_1e_man03
!##########################################################
subroutine ccdt_t3p_1e_man03_1(i1aa,i1ab)

! i1 ( i a j k )_ft + = Sum ( b ) 
!  * t ( a b j k )_t * f ( i b )_f 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g3inp,fock
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,1:norb1,1:norb1), &
       i1ab(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooovaa
     i = h3_ooovaa(icc)
     a = p4_ooovaa(icc)
     j = h1_ooovaa(icc)
     k = h2_ooovaa(icc)
     do b = norb1+1,nact
        i1aa(i,a,j,k) = i1aa(i,a,j,k) + t2inp(a,b,j,k,spin_t2aa)*fock(i,b,spin_focka)
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
        i1ab(i,a,j,k) = i1ab(i,a,j,k) - t2inp(b,a,j,k,spin_t2ab)*fock(i,b,spin_focka)
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_t3p_1e_man03_1
!##########################################################
