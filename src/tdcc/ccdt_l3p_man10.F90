!##########################################################
subroutine ccdt_l3p_man10(i0,work1,work2,work3)

! i0 ( i j k a b c )_ytv + = -1/2 * P( 9 ) * Sum ( l ) 
!  * i1 ( i j l a )_yt * v ( k l b c )_v 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: &
       work1(1:norb1,1:norb1,1:norb1,(norb1+1):nact),&
       work2(1:norb1,1:norb1,1:norb1,(norb1+1):nact),work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
!  call ccdt_l3p_10_1(1,1,1,1,work1)
!  call ccdt_l3p_10_1(1,2,1,2,work2)
!  work1 = -0.5d0*work1
!  work2 = -0.5d0*work2
  call ccdt_l3p_man10_1(work1,work2)
  call tdcc_fillooovaa(work1)

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
        i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
             + work1(i,j,l,a)*int2x(k,l,b,c,spin_int2aa) &
             - work1(k,j,l,a)*int2x(i,l,b,c,spin_int2aa) &
             - work1(i,k,l,a)*int2x(j,l,b,c,spin_int2aa) &
             - work1(i,j,l,b)*int2x(k,l,a,c,spin_int2aa) &
             - work1(i,j,l,c)*int2x(k,l,b,a,spin_int2aa) &
             + work1(k,j,l,b)*int2x(i,l,a,c,spin_int2aa) &
             + work1(i,k,l,b)*int2x(j,l,a,c,spin_int2aa) &
             + work1(k,j,l,c)*int2x(i,l,b,a,spin_int2aa) &
             + work1(i,k,l,c)*int2x(j,l,b,a,spin_int2aa)
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
        i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
             - work1(i,j,l,a)*int2x(l,k,b,c,spin_int2ab) &
             - work2(k,j,l,a)*int2x(i,l,b,c,spin_int2ab) &
             + work2(k,i,l,a)*int2x(j,l,b,c,spin_int2ab) &
             + work1(i,j,l,b)*int2x(l,k,a,c,spin_int2ab) &
!ZERO        - work1(i,j,l,c)*int2x(k,l,b,a,spin_int2aa) &
             + work2(k,j,l,b)*int2x(i,l,a,c,spin_int2ab) &
             - work2(k,i,l,b)*int2x(j,l,a,c,spin_int2ab) &
             - work2(j,k,l,c)*int2x(i,l,b,a,spin_int2aa) &
             + work2(i,k,l,c)*int2x(j,l,b,a,spin_int2aa)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man10
!##########################################################
subroutine ccdt_l3p_man10_1(i1aa,i1ab)

!     i1 ( i j k a )_yt + = -0.5 * Sum ( l b c ) * t ( b c l k )_t * y ( i j l a b c )_y 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,1:norb1,(norb1+1):nact),&
       i1ab(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooovaa
     i = h1_ooovaa(icc)
     j = h2_ooovaa(icc)
     k = h3_ooovaa(icc)
     a = p4_ooovaa(icc)
     do l = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,b-1
           i1aa(i,j,k,a) = i1aa(i,j,k,a) &
                - t2inp(b,c,l,k,spin_t2aa)*g3inp(i,j,l,a,b,c,spin_g3aaa)
        end do
        do c = norb1+1,act1_ul
           i1aa(i,j,k,a) = i1aa(i,j,k,a) &
                + t2inp(b,c,k,l,spin_t2ab)*g3inp(i,j,l,a,b,c,spin_g3aab)
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
     do l = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,b-1
           i1ab(i,j,k,a) = i1ab(i,j,k,a) &
                - t2inp(b,c,l,k,spin_t2aa)*g3inp(i,l,j,c,b,a,spin_g3aab)
        end do
        do c = norb1+1,act1_ul
           i1ab(i,j,k,a) = i1ab(i,j,k,a) &
                + t2inp(b,c,k,l,spin_t2ab)*g3inp(l,j,i,a,c,b,spin_g3aab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_man10_1
!##########################################################
