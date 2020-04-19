!##########################################################
subroutine ccdt_l3p_man03(i0,work1,work2,work3)

! i0 ( i j k a b c )_yv + = -1 * P( 9 ) * Sum ( d ) 
!  * y ( i j a d )_y * v ( k d b c )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: work1(1),work2(1),work3(1)

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
        i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
             - g2inp(i,j,a,d,spin_g2aa)*int2x(k,d,b,c,spin_int2aa) &
             + g2inp(k,j,a,d,spin_g2aa)*int2x(i,d,b,c,spin_int2aa) &
             + g2inp(i,k,a,d,spin_g2aa)*int2x(j,d,b,c,spin_int2aa) &
             + g2inp(i,j,b,d,spin_g2aa)*int2x(k,d,a,c,spin_int2aa) &
             + g2inp(i,j,c,d,spin_g2aa)*int2x(k,d,b,a,spin_int2aa) &
             - g2inp(k,j,b,d,spin_g2aa)*int2x(i,d,a,c,spin_int2aa) &
             - g2inp(i,k,b,d,spin_g2aa)*int2x(j,d,a,c,spin_int2aa) &
             - g2inp(k,j,c,d,spin_g2aa)*int2x(i,d,b,a,spin_int2aa) &
             - g2inp(i,k,c,d,spin_g2aa)*int2x(j,d,b,a,spin_int2aa)
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
        i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
             + g2inp(i,j,a,d,spin_g2aa)*int2x(d,k,b,c,spin_int2ab) &
             - g2inp(j,k,a,d,spin_g2ab)*int2x(i,d,b,c,spin_int2ab) &
             + g2inp(i,k,a,d,spin_g2ab)*int2x(j,d,b,c,spin_int2ab) &
             - g2inp(i,j,b,d,spin_g2aa)*int2x(d,k,a,c,spin_int2ab) &
             + g2inp(j,k,b,d,spin_g2ab)*int2x(i,d,a,c,spin_int2ab) &
             - g2inp(i,k,b,d,spin_g2ab)*int2x(j,d,a,c,spin_int2ab) &
             - g2inp(j,k,d,c,spin_g2ab)*int2x(i,d,b,a,spin_int2aa) &
             + g2inp(i,k,d,c,spin_g2ab)*int2x(j,d,b,a,spin_int2aa)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man03
!##########################################################
