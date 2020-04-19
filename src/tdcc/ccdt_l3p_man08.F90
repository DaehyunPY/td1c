!##########################################################
subroutine ccdt_l3p_man08(i0,work1,work2,work3)

! 1
! i0 ( i j k a b c )_yf + = 1 * P( 9 ) * y ( i j a b )_y * f ( k c )_f 0
! 2
! i0 ( i j k a b c )_yv + = -1 * P( 9 ) * Sum ( l ) 
!  * y ( i l a b )_y * v ( j k l c )_v 0
! 3
! i0 ( i j k a b c )_yv + = -1 * P( 9 ) * Sum ( d ) 
!  * y ( i j a d )_y * v ( k d b c )_v 0
! 8
! i0 ( i j k a b c )_yv + = 1/2 * P( 3 ) * Sum ( d e ) 
!  * y ( i j k a d e )_y * v ( d e b c )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: work1(1),work2(1),work3(1)

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
     ! diagram 1
     i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
          + g2inp(i,j,a,b,spin_g2aa)*fock(k,c,spin_focka) &
          - g2inp(k,j,a,b,spin_g2aa)*fock(i,c,spin_focka) &
          - g2inp(i,k,a,b,spin_g2aa)*fock(j,c,spin_focka) &
          - g2inp(i,j,c,b,spin_g2aa)*fock(k,a,spin_focka) &
          - g2inp(i,j,a,c,spin_g2aa)*fock(k,b,spin_focka) &
          + g2inp(k,j,c,b,spin_g2aa)*fock(i,a,spin_focka) &
          + g2inp(k,j,a,c,spin_g2aa)*fock(i,b,spin_focka) &
          + g2inp(i,k,c,b,spin_g2aa)*fock(j,a,spin_focka) &
          + g2inp(i,k,a,c,spin_g2aa)*fock(j,b,spin_focka)
     ! diagram 2
     do l = 1, norb1
        i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
             - g2inp(i,l,a,b,spin_g2aa)*int2x(j,k,l,c,spin_int2aa) &
             + g2inp(j,l,a,b,spin_g2aa)*int2x(i,k,l,c,spin_int2aa) &
             + g2inp(k,l,a,b,spin_g2aa)*int2x(j,i,l,c,spin_int2aa) &
             + g2inp(i,l,c,b,spin_g2aa)*int2x(j,k,l,a,spin_int2aa) &
             + g2inp(i,l,a,c,spin_g2aa)*int2x(j,k,l,b,spin_int2aa) &
             - g2inp(j,l,c,b,spin_g2aa)*int2x(i,k,l,a,spin_int2aa) &
             - g2inp(j,l,a,c,spin_g2aa)*int2x(i,k,l,b,spin_int2aa) &
             - g2inp(k,l,c,b,spin_g2aa)*int2x(j,i,l,a,spin_int2aa) &
             - g2inp(k,l,a,c,spin_g2aa)*int2x(j,i,l,b,spin_int2aa)
     end do
     ! diagram 3
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
     ! diagram 8
     do d = norb1+1,act1_ul
        do e = norb1+1,d-1
           i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
                + g3inp(i,j,k,a,d,e,spin_g3aaa)*int2x(d,e,b,c,spin_int2aa) &
                - g3inp(i,j,k,b,d,e,spin_g3aaa)*int2x(d,e,a,c,spin_int2aa) &
                - g3inp(i,j,k,c,d,e,spin_g3aaa)*int2x(d,e,b,a,spin_int2aa)
        end do
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
     ! diagram 1
     i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
          + g2inp(i,j,a,b,1)*fock(k,c,spin_focka) &
          + g2inp(k,j,c,b,4)*fock(i,a,spin_focka) &
          + g2inp(k,j,a,c,6)*fock(i,b,spin_focka) &
          + g2inp(i,k,c,b,5)*fock(j,a,spin_focka) &
          + g2inp(i,k,a,c,3)*fock(j,b,spin_focka)
     ! diagram 2
     do l = 1, norb1
        i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
             - g2inp(i,l,a,b,spin_g2aa)*int2x(j,k,l,c,spin_int2ab) &
             + g2inp(j,l,a,b,spin_g2aa)*int2x(i,k,l,c,spin_int2ab) &
!ZERO        + g2inp(k,l,a,b,spin_g2aa)*int2x(j,i,l,c,spin_int2aa) &
             + g2inp(i,l,b,c,spin_g2ab)*int2x(j,k,a,l,spin_int2ab) &
             - g2inp(i,l,a,c,spin_g2ab)*int2x(j,k,b,l,spin_int2ab) &
             - g2inp(j,l,b,c,spin_g2ab)*int2x(i,k,a,l,spin_int2ab) &
             + g2inp(j,l,a,c,spin_g2ab)*int2x(i,k,b,l,spin_int2ab) &
             - g2inp(k,l,c,b,spin_g2ab)*int2x(j,i,l,a,spin_int2aa) &
             + g2inp(l,k,a,c,spin_g2ab)*int2x(j,i,l,b,spin_int2aa)
     end do
     ! diagram 3
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
     ! diagram 8
     do d = norb1+1,act1_ul
        do e = norb1+1,act1_ul
           i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
                + g3inp(i,j,k,a,d,e,spin_g3aab)*int2x(d,e,b,c,spin_int2ab) &
                - g3inp(i,j,k,b,d,e,spin_g3aab)*int2x(d,e,a,c,spin_int2ab)
        end do
        do e = norb1+1,d-1
           i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
                - g3inp(i,j,k,d,e,c,spin_g3aab)*int2x(d,e,b,a,spin_int2aa)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man08
!##########################################################
