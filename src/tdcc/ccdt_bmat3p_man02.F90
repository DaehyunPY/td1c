!##########################################################
subroutine ccdt_bmat3p_man02(i0,work1,work2)

! i0 ( i a )_dytt + = -1/2 * Sum ( j k b ) 
!  * i1 ( j k b i )_dyt * t ( a b j k )_t 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1,(norb1+1):nact,1:norb1),&
       work2(1:norb1,1:norb1,(norb1+1):nact,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_bmat3p_man02_1(work1,work2)
  call tdcc_filloovoaa(work1)

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, nov
     i = h1_ov(icc)
     a = p2_ov(icc)
     do b = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i0(i,a) = i0(i,a) &
                - work1(j,k,b,i)*t2inp(a,b,j,k,spin_t2aa)
        end do
        do k = 1,norb1
           i0(i,a) = i0(i,a) &
                + work2(k,j,b,i)*t2inp(a,b,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_bmat3p_man02
!##########################################################
subroutine ccdt_bmat3p_man02_1(i1aa,i1ab)

! i1 ( i j a k )_dyt + = +1/2 * Sum ( l b c ) 
!  * dy ( l i j b c a )_dy * t ( c b k l )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : norb1,dg3inp,t2inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,(norb1+1):nact,1:norb1),&
       i1ab(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooovaa
     i = h1_ooovaa(icc)
     j = h2_ooovaa(icc)
     a = p4_ooovaa(icc)
     k = h3_ooovaa(icc)
     do l = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,b-1
           i1aa(i,j,a,k) = i1aa(i,j,a,k) &
                + dg3inp(l,i,j,b,c,a,spin_g3aaa)*t2inp(c,b,k,l,spin_t2aa)
        end do
        do c = norb1+1,act1_ul
           i1aa(i,j,a,k) = i1aa(i,j,a,k) &
                + dg3inp(i,j,l,c,a,b,spin_g3aab)*t2inp(c,b,k,l,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nooovab
     i = h1_ooovab(icc)
     j = h2_ooovab(icc)
     a = p4_ooovab(icc)
     k = h3_ooovab(icc)
     do l = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,b-1
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                - dg3inp(l,j,i,b,c,a,spin_g3aab)*t2inp(c,b,k,l,spin_t2aa)
        end do
        do c = norb1+1,act1_ul
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                - dg3inp(l,i,j,b,a,c,spin_g3aab)*t2inp(c,b,k,l,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_bmat3p_man02_1
!##########################################################
