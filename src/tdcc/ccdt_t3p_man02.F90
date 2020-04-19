!##########################################################
subroutine ccdt_t3p_man02(i0,work1,work2,work3)

! i0 ( a b c i j k )_vt + = -1 * P( 9 ) * Sum ( d ) 
!  * t ( a d i j )_t * i1 ( b c k d )_v 4

  use mod_ormas, only : nact
  use mod_cc, only : t2inp,ncc3aaa,ncc3aab,norb1
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa 
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab 
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       work2((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact),work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_t3p_man02_1(work1,work2)
  call tdcc_fillvvovaa(work1)

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
     do d = norb1+1,nact
        i0(a,b,c,i,j,k,1) = i0(a,b,c,i,j,k,1) &
             + t2inp(a,d,i,j,spin_t2aa) * work1(c,b,k,d) &
             - t2inp(a,d,j,k,spin_t2aa) * work1(b,c,i,d) &
             + t2inp(a,d,i,k,spin_t2aa) * work1(b,c,j,d) &
             - t2inp(b,d,i,j,spin_t2aa) * work1(c,a,k,d) &
             + t2inp(c,d,i,j,spin_t2aa) * work1(b,a,k,d) &
             + t2inp(b,d,j,k,spin_t2aa) * work1(a,c,i,d) &
             - t2inp(d,c,j,k,spin_t2aa) * work1(b,a,i,d) &
             - t2inp(b,d,i,k,spin_t2aa) * work1(a,c,j,d) &
             + t2inp(d,c,i,k,spin_t2aa) * work1(b,a,j,d)
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
     do d = norb1+1,nact
        i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
             + t2inp(a,d,i,j,spin_t2aa) * work2(c,b,k,d) &
             - t2inp(a,d,j,k,spin_t2ab) * work2(b,c,i,d) &
             + t2inp(a,d,i,k,spin_t2ab) * work2(b,c,j,d) &
             - t2inp(b,d,i,j,spin_t2aa) * work2(c,a,k,d) &
            !+ t2inp(c,d,i,j,spin_t2aa) * work1(b,a,k,d) &
             + t2inp(b,d,j,k,spin_t2ab) * work2(a,c,i,d) &
             - t2inp(d,c,j,k,spin_t2ab) * work1(b,a,i,d) &
             - t2inp(b,d,i,k,spin_t2ab) * work2(a,c,j,d) &
             + t2inp(d,c,i,k,spin_t2ab) * work1(b,a,j,d)
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man02
!##########################################################
subroutine ccdt_t3p_man02_1(i1aa,i1ab)

! i1 ( a b i c )_v + = 1 * v ( a b i c )_v 0
! i1 ( a b i c )_vt + = 1/2 * Sum ( j k ) * t ( a b j k )_t * v ( j k i c )_v 0
! i1 ( a b i c )_vt + = 1 * P( 2 ) * Sum ( j d ) * t ( a d i j )_t * v ( j b d c )_v 0
! i1 ( a b i c )_vt + = 1/2 * Sum ( l m d ) * t ( a b d i l m )_t * v ( l m d c )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,t3inp,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact), &
       i1ab((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novvvaa
     a = p3_ovvvaa(icc)
     b = p4_ovvvaa(icc)
     i = h1_ovvvaa(icc)
     c = p2_ovvvaa(icc)
     !i1-1
     i1aa(a,b,i,c) = i1aa(a,b,i,c) + int2x(a,b,i,c,spin_int2aa)

     !i1-2
     do j = 1,norb1
        do k = 1,j-1
           i1aa(a,b,i,c) = i1aa(a,b,i,c) &
                + t2inp(a,b,j,k,spin_t2aa)*int2x(j,k,i,c,spin_int2aa)
        end do
     end do

     !i1-3
     do j = 1,norb1
     do d = norb1+1,nact
        i1aa(a,b,i,c) = i1aa(a,b,i,c) &
             + t2inp(a,d,i,j,spin_t2aa)*int2x(j,b,d,c,spin_int2aa) &
             + t2inp(a,d,i,j,spin_t2ab)*int2x(j,b,d,c,spin_int2ab) &
             - t2inp(b,d,i,j,spin_t2aa)*int2x(j,a,d,c,spin_int2aa) &
             - t2inp(b,d,i,j,spin_t2ab)*int2x(j,a,d,c,spin_int2ab)
     end do
     end do

     !i1-4
     do d = norb1+1,nact
     do l = 1,norb1
     do m = 1,l-1
        i1aa(a,b,i,c) = i1aa(a,b,i,c) + t3inp(a,b,d,i,l,m,spin_t3aaa)*int2x(l,m,d,c,spin_int2aa)
     end do
     end do
     end do
     do d = norb1+1,nact
     do l = 1,norb1
     do m = 1,norb1
        i1aa(a,b,i,c) = i1aa(a,b,i,c) - t3inp(a,b,d,i,l,m,spin_t3aab)*int2x(l,m,c,d,spin_int2ab)
     end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, novvvab
     a = p3_ovvvab(icc)
     b = p4_ovvvab(icc)
     i = h1_ovvvab(icc)
     c = p2_ovvvab(icc)
     !i1-1
     i1ab(a,b,i,c) = i1ab(a,b,i,c) + int2x(a,b,i,c,spin_int2ab)

     !i1-2
     do j = 1,norb1
        do k = 1,norb1
           i1ab(a,b,i,c) = i1ab(a,b,i,c) &
                + t2inp(a,b,j,k,spin_t2ab)*int2x(j,k,i,c,spin_int2ab)
        end do
     end do

     !i1-3
     do j = 1,norb1
     do d = norb1+1,nact
        i1ab(a,b,i,c) = i1ab(a,b,i,c) &
             + t2inp(a,d,i,j,spin_t2aa)*int2x(j,b,d,c,spin_int2ab) &
             + t2inp(a,d,i,j,spin_t2ab)*int2x(j,b,d,c,spin_int2aa) &
             - t2inp(d,b,i,j,spin_t2ab)*int2x(a,j,d,c,spin_int2ab)
     end do
     end do

     !i1-4
     do d = norb1+1,nact
     do l = 1,norb1
     do m = 1,l-1
        i1ab(a,b,i,c) = i1ab(a,b,i,c) + t3inp(b,d,a,l,m,i,spin_t3aab)*int2x(l,m,d,c,spin_int2aa)
     end do
     end do
     end do
     do d = norb1+1,nact
     do l = 1,norb1
     do m = 1,norb1
        i1ab(a,b,i,c) = i1ab(a,b,i,c) - t3inp(a,d,b,i,l,m,spin_t3aab)*int2x(l,m,d,c,spin_int2ab)
     end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man02_1
!##########################################################
