!##########################################################
subroutine ccdt_t3p_man01(i0,work1,work2,work3)

! i0 ( a b c i j k )_vt + = -1 * P( 9 ) * Sum ( l ) 
!  * t ( a b i l )_t * i1 ( l c j k )_v 5

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : t2inp,ncc3aaa,ncc3aab,norb1
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa 
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab 
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,(norb1+1):nact,1:norb1,1:norb1), &
       work2(1:norb1,(norb1+1):nact,1:norb1,1:norb1),work3(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_t3p_man01_1(work1,work2)
  call tdcc_fillovooaa(work1)

  !$omp parallel default(shared) private(a,b,c,i,j,k)
  !$omp do
  do icc = 1, ncc3aaa
     !aaaaaa
     a = p1_cc3aaa(icc)
     b = p2_cc3aaa(icc)
     c = p3_cc3aaa(icc)
     i = h1_cc3aaa(icc)
     j = h2_cc3aaa(icc)
     k = h3_cc3aaa(icc)
     do l = 1,norb1
        i0(a,b,c,i,j,k,1) = i0(a,b,c,i,j,k,1) &
             + t2inp(a,b,i,l,spin_t2aa) * work1(l,c,j,k) &
             - t2inp(a,b,j,l,spin_t2aa) * work1(l,c,i,k) &
             - t2inp(a,b,k,l,spin_t2aa) * work1(l,c,j,i) &
             - t2inp(b,c,i,l,spin_t2aa) * work1(l,a,k,j) &
             + t2inp(a,c,i,l,spin_t2aa) * work1(l,b,k,j) &
             + t2inp(b,c,j,l,spin_t2aa) * work1(l,a,k,i) &
             - t2inp(a,c,j,l,spin_t2aa) * work1(l,b,k,i) &
             + t2inp(b,c,l,k,spin_t2aa) * work1(l,a,j,i) &
             - t2inp(a,c,l,k,spin_t2aa) * work1(l,b,j,i)
     end do
  end do
  !$omp end do

  !$omp do
  do icc = 1, ncc3aab
     !aabaab
     a = p1_cc3aab(icc)
     b = p2_cc3aab(icc)
     c = p3_cc3aab(icc)
     i = h1_cc3aab(icc)
     j = h2_cc3aab(icc)
     k = h3_cc3aab(icc)
     do l = 1,norb1
        i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
             + t2inp(a,b,i,l,spin_t2aa) * work2(l,c,j,k) &
             - t2inp(a,b,j,l,spin_t2aa) * work2(l,c,i,k) &
            !- t2inp(a,b,k,l,spin_t2aa) * work1(l,c,j,i) &
             - t2inp(b,c,i,l,spin_t2ab) * work2(l,a,k,j) &
             + t2inp(a,c,i,l,spin_t2ab) * work2(l,b,k,j) &
             + t2inp(b,c,j,l,spin_t2ab) * work2(l,a,k,i) &
             - t2inp(a,c,j,l,spin_t2ab) * work2(l,b,k,i) &
             + t2inp(b,c,l,k,spin_t2ab) * work1(l,a,j,i) &
             - t2inp(a,c,l,k,spin_t2ab) * work1(l,b,j,i)
     end do
  end do
  !$omp end do
  !$omp end parallel

  !debug
  !do icc = 1, ncc3aaa
  !   a = p1_cc3aaa(icc)
  !   b = p2_cc3aaa(icc)
  !   c = p3_cc3aaa(icc)
  !   i = h1_cc3aaa(icc)
  !   j = h2_cc3aaa(icc)
  !   k = h3_cc3aaa(icc)
  !   write(6,"('i0aa: ',6i5,2f20.10)") a,b,c,i,j,k,i0(a,b,c,i,j,k,1)
  !end do
  !debug

end subroutine ccdt_t3p_man01
!##########################################################
subroutine ccdt_t3p_man01_1(i1aa,i1ab)

! i1 ( i a j k )_v + = -1 * v ( i a j k )_v 0
! i1 ( i a j k )_ft + = +1 * Sum ( b ) * t ( a b j k )_t * f ( i b )_f 0
! i1 ( i a j k )_vt + = -1 * P( 2 ) * Sum ( l b ) * t ( a b j l )_t * v ( l i k b )_v 0
! i1 ( i a j k )_vt + = -1/2 * Sum ( b c ) * t ( b c j k )_t * v ( i a b c )_v 0

!!i1 ( k a i j )_vt + = -1/2 * Sum ( l d e ) * t ( a d e i j l )_t * v ( l k d e )_v 0
! i1 ( i a j k )_vt + = -1/2 * Sum ( l b c ) * t ( a b c j k l )_t * v ( l i b c )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,t3inp,fock,int2x
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
     !i1-1
     i1aa(i,a,j,k) = i1aa(i,a,j,k) - int2x(i,a,j,k,spin_int2aa)

     !i1-2
     do b = norb1+1,nact
        i1aa(i,a,j,k) = i1aa(i,a,j,k) + t2inp(a,b,j,k,spin_t2aa)*fock(i,b,spin_focka)
     end do

     !i1-3
     do l = 1,norb1
     do b = norb1+1,nact
        i1aa(i,a,j,k) = i1aa(i,a,j,k) &
             - t2inp(a,b,j,l,spin_t2aa)*int2x(l,i,k,b,spin_int2aa) &
             + t2inp(a,b,j,l,spin_t2ab)*int2x(i,l,k,b,spin_int2ab) &
             + t2inp(a,b,k,l,spin_t2aa)*int2x(l,i,j,b,spin_int2aa) &
             - t2inp(a,b,k,l,spin_t2ab)*int2x(i,l,j,b,spin_int2ab)
     end do
     end do

     !i1-4
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i1aa(i,a,j,k) = i1aa(i,a,j,k) &
                - t2inp(b,c,j,k,spin_t2aa)*int2x(i,a,b,c,spin_int2aa)
        end do
     end do

     !i1-5
     do l = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,b-1
           i1aa(i,a,j,k) = i1aa(i,a,j,k) &
                - t3inp(a,b,c,j,k,l,spin_t3aaa)*int2x(l,i,b,c,spin_int2aa)
        end do
     end do
     end do
     do l = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,act1_ul
           i1aa(i,a,j,k) = i1aa(i,a,j,k) &
                + t3inp(a,b,c,j,k,l,spin_t3aab)*int2x(i,l,b,c,spin_int2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nooovab
     i = h3_ooovab(icc)
     a = p4_ooovab(icc)
     j = h1_ooovab(icc)
     k = h2_ooovab(icc)
     !i1-1
     i1ab(i,a,j,k) = i1ab(i,a,j,k) - int2x(i,a,j,k,spin_int2ab)

     !i1-2
     do b = norb1+1,nact
        i1ab(i,a,j,k) = i1ab(i,a,j,k) - t2inp(b,a,j,k,spin_t2ab)*fock(i,b,spin_focka)
     end do

     !i1-3
     do l = 1,norb1
     do b = norb1+1,nact
        i1ab(i,a,j,k) = i1ab(i,a,j,k) &
             + t2inp(b,a,j,l,spin_t2ab)*int2x(l,i,k,b,spin_int2ab) &
             + t2inp(a,b,k,l,spin_t2ab)*int2x(l,i,j,b,spin_int2aa) &
             - t2inp(a,b,k,l,spin_t2aa)*int2x(i,l,j,b,spin_int2ab)
     end do
     end do

     !i1-4
     do b = norb1+1,nact
        do c = norb1+1,nact
           i1ab(i,a,j,k) = i1ab(i,a,j,k) &
                - t2inp(b,c,j,k,spin_t2ab)*int2x(i,a,b,c,spin_int2ab)
        end do
     end do

     !i1-5
     do l = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,b-1
           i1ab(i,a,j,k) = i1ab(i,a,j,k) &
                - t3inp(c,b,a,j,l,k,spin_t3aab)*int2x(l,i,b,c,spin_int2aa)
        end do
     end do
     end do
     do l = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,act1_ul
           i1ab(i,a,j,k) = i1ab(i,a,j,k) &
                + t3inp(a,c,b,l,k,j,spin_t3aab)*int2x(i,l,b,c,spin_int2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_t3p_man01_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
