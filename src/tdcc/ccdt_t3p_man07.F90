!##########################################################
subroutine ccdt_t3p_man07(i0,work1,work2,work3)

!7:  i0 ( a b c i j k )_vt + = 1/2 * P( 3 ) * Sum ( d e ) 
!     * t ( a d e i j k )_t * v ( b c d e )_v 0
!8:  i0 ( a b c i j k )_vtt + = 1/4 * P( 3 ) * Sum ( d e ) 
!     * t ( c d e i j k )_t * i1 ( a b d e )_vt 1

  use mod_ormas, only : nact
  use mod_cc, only : t3inp,ncc3aaa,ncc3aab,norb1
  use mod_cc, only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa 
  use mod_cc, only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab 
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       work2((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_t3p_man07_1(work1,work2)
  call tdcc_fillvvvvaa(work1)

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
     do e = norb1+1,d-1
        i0(a,b,c,i,j,k,1) = i0(a,b,c,i,j,k,1) &
             + t3inp(a,d,e,i,j,k,spin_t3aaa) * work1(b,c,d,e) &
             + t3inp(d,b,e,i,j,k,spin_t3aaa) * work1(a,c,d,e) &
             + t3inp(d,e,c,i,j,k,spin_t3aaa) * work1(a,b,d,e)
     end do
     end do
  end do
  !$omp end do

  !$omp do
  do icc = 1, ncc3aab
     ! aabaab
     a = p1_cc3aab(icc)
     b = p2_cc3aab(icc)
     c = p3_cc3aab(icc)
     i = h1_cc3aab(icc)
     j = h2_cc3aab(icc)
     k = h3_cc3aab(icc)
     do d = norb1+1,nact
     do e = norb1+1,nact
        i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
             + t3inp(a,d,e,i,j,k,spin_t3aab) * work2(b,c,d,e) &
             + t3inp(d,b,e,i,j,k,spin_t3aab) * work2(a,c,d,e)
     end do
     end do
     do d = norb1+1,nact
     do e = norb1+1,d-1
        i0(a,b,c,i,j,k,2) = i0(a,b,c,i,j,k,2) &
             + t3inp(d,e,c,i,j,k,spin_t3aab) * work1(a,b,d,e)
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man07
!##########################################################
subroutine ccdt_t3p_man07_1(i1aa,i1ab)

! i1 ( a b c d )_vt + = v ( a b c d )_v 0
! i1 ( a b c d )_vt + = 0.5 * Sum ( i j ) 
!  * v ( i j c d )_v * t ( a b i j )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       i1ab((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nvvvvaa
     a = p1_vvvvaa(icc)
     b = p2_vvvvaa(icc)
     c = p3_vvvvaa(icc)
     d = p4_vvvvaa(icc)
     i1aa(a,b,c,d) = i1aa(a,b,c,d) + int2x(a,b,c,d,spin_int2aa)
     do i = 1,norb1
        do j = 1,i-1
           i1aa(a,b,c,d) = i1aa(a,b,c,d) + int2x(i,j,c,d,spin_int2aa)*t2inp(a,b,i,j,spin_t2aa)
        end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nvvvvab
     a = p1_vvvvab(icc)
     b = p2_vvvvab(icc)
     c = p3_vvvvab(icc)
     d = p4_vvvvab(icc)
     i1ab(a,b,c,d) = i1ab(a,b,c,d) + int2x(a,b,c,d,spin_int2ab)
     do i = 1,norb1
        do j = 1,norb1
           i1ab(a,b,c,d) = i1ab(a,b,c,d) + int2x(i,j,c,d,spin_int2ab)*t2inp(a,b,i,j,spin_t2ab)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t3p_man07_1
!##########################################################
