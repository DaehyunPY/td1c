!##########################################################
subroutine ccdt_t3p_1e_man02(i0,work1,work2)

! i0 ( a b c i j k )_tf + = 1 * P( 3 ) * Sum ( d ) 
!  * t ( a b d i j k )_t * f ( c d )_f 0

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2)
  complex(kind(0d0)) :: work1(1),work2(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

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
     do d = norb1+1,nact
        i0(a,b,c,i,j,k,1) = i0(a,b,c,i,j,k,1) &
             + t3inp(a,b,d,i,j,k,spin_t3aaa)*fock(c,d,spin_focka) &
             - t3inp(c,b,d,i,j,k,spin_t3aaa)*fock(a,d,spin_focka) &
             - t3inp(a,c,d,i,j,k,spin_t3aaa)*fock(b,d,spin_focka)
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
     do d = norb1+1,nact
!     do l = 1,norb1
        i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
             + t3inp(a,b,d,i,j,k,spin_t3aab)*fock(c,d,spin_focka) &
             + t3inp(d,b,c,i,j,k,spin_t3aab)*fock(a,d,spin_focka) &
             + t3inp(a,d,c,i,j,k,spin_t3aab)*fock(b,d,spin_focka)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_t3p_1e_man02
!##########################################################
