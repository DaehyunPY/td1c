!##########################################################
subroutine ccdt_t3p_1e_man01(i0,work1,work2)

! i0 ( a b c i j k )_tf + = -1 * P( 3 ) * Sum ( l ) 
!  * t ( a b c i j l )_t * f ( l k )_f 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2)
  complex(kind(0d0)) :: work1(1),work2(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

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
             - t3inp(a,b,c,i,j,l,spin_t3aaa)*fock(l,k,spin_focka) &
             + t3inp(a,b,c,k,j,l,spin_t3aaa)*fock(l,i,spin_focka) &
             + t3inp(a,b,c,i,k,l,spin_t3aaa)*fock(l,j,spin_focka)
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
             - t3inp(a,b,c,i,j,l,spin_t3aab)*fock(l,k,spin_focka) &
             - t3inp(a,b,c,l,j,k,spin_t3aab)*fock(l,i,spin_focka) &
             - t3inp(a,b,c,i,l,k,spin_t3aab)*fock(l,j,spin_focka)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_t3p_1e_man01
!##########################################################
