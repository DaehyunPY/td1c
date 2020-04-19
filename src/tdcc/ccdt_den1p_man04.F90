!##########################################################
subroutine ccdt_den1p_man04(i0,work1,work2)

!3: i0 ( a i )_yt + = +1/4 * Sum ( j k b c ) 
!  * y ( j k b c )_y * t ( b c a j k i )_t 0
!4: i0 ( a i )_ytt + = 1/2 * Sum ( j b c ) 
!  * i1 ( j a b c )_yt * t ( b c i j )_t 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       work2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den1p_man04_1(work1,work2)
  call tdcc_fillovvvaa(work1)

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, nov
     a = p2_ov(icc)
     i = h1_ov(icc)
     ! diagram 3
     do j = act1_ll,norb1
     do b = norb1+1,act1_ul
        do k = act1_ll,j-1
        do c = norb1+1,b-1
           i0(a,i) = i0(a,i) + g2inp(j,k,b,c,spin_g2aa)*( &
                t3inp(b,c,a,j,k,i,spin_t3aaa) + &
                t3inp(b,c,a,j,k,i,spin_t3aab))
        end do
        end do
        do k = act1_ll,norb1
        do c = norb1+1,act1_ul
           i0(a,i) = i0(a,i) + g2inp(j,k,b,c,spin_g2ab)* &
                t3inp(a,b,c,i,j,k,spin_t3aab)
        end do
        end do
     end do
     end do

     ! diagram 4
     do j = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,b-1
           i0(a,i) = i0(a,i) + work1(j,a,b,c)*t2inp(b,c,i,j,spin_t2aa)
        end do
        do c = norb1+1,act1_ul
           i0(a,i) = i0(a,i) - work2(j,a,b,c)*t2inp(c,b,i,j,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den1p_man04
!##########################################################
subroutine ccdt_den1p_man04_1(i1aa,i1ab)

! i1 ( i a b c )_yt + = 1/2 * Sum ( j k d ) * y ( i j k b c d )_y * t ( a d j k )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),&
       i1ab(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, novvvaa
     i = h1_ovvvaa(icc)
     a = p2_ovvvaa(icc)
     b = p3_ovvvaa(icc)
     c = p4_ovvvaa(icc)
     do d = norb1+1,act1_ul
     do j = act1_ll,norb1
        do k = act1_ll,j-1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + g3inp(i,j,k,b,c,d,spin_g3aaa)*t2inp(a,d,j,k,spin_t2aa)
        end do
        do k = act1_ll,norb1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + g3inp(i,j,k,b,c,d,spin_g3aab)*t2inp(a,d,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, novvvab
     i = h1_ovvvab(icc)
     a = p2_ovvvab(icc)
     b = p3_ovvvab(icc)
     c = p4_ovvvab(icc)
     do d = norb1+1,act1_ul
     do j = act1_ll,norb1
        do k = act1_ll,j-1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                + g3inp(j,k,i,c,d,b,spin_g3aab)*t2inp(a,d,j,k,spin_t2aa)
        end do
        do k = act1_ll,norb1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                + g3inp(i,k,j,b,d,c,spin_g3aab)*t2inp(a,d,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den1p_man04_1
!##########################################################
