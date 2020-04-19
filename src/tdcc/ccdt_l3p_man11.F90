!##########################################################
subroutine ccdt_l3p_man11(i0,work1,work2,work3)

! i0 ( i j k a b c )_ytv + = 1/4 * P( 3 ) * Sum ( d e ) 
!  * y ( i j k a d e )_y * i1 ( d e b c )_vt 1

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: &
       work1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),&
       work2((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_l3p_man11_1(work1,work2)
  call tdcc_fillvvvvaa(work1)

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
        do e = norb1+1,d-1
           i0(i,j,k,a,b,c,1) = i0(i,j,k,a,b,c,1) &
                + g3inp(i,j,k,a,d,e,spin_g3aaa)*work1(d,e,b,c) &
                - g3inp(i,j,k,b,d,e,spin_g3aaa)*work1(d,e,a,c) &
                - g3inp(i,j,k,c,d,e,spin_g3aaa)*work1(d,e,b,a)
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
     do d = norb1+1,nact
        do e = norb1+1,nact
           i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
                + g3inp(i,j,k,a,d,e,spin_g3aab)*work2(d,e,b,c) &
                - g3inp(i,j,k,b,d,e,spin_g3aab)*work2(d,e,a,c)
        end do
        do e = norb1+1,d-1
           i0(i,j,k,a,b,c,2) = i0(i,j,k,a,b,c,2) &
                - g3inp(i,j,k,d,e,c,spin_g3aab)*work1(d,e,b,a)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man11
!##########################################################
subroutine ccdt_l3p_man11_1(i1aa,i1ab)

! i1 ( a b c d )_vt + = 0.5 * Sum ( i j ) 
!  * t ( a b i j )_t * v ( i j c d )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,int2x,t2inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),&
       i1ab((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, nvvvvaa
     a = p1_vvvvaa(icc)
     b = p2_vvvvaa(icc)
     c = p3_vvvvaa(icc)
     d = p4_vvvvaa(icc)
     do i = 1,norb1
        do j = 1,i-1
           i1aa(a,b,c,d) = i1aa(a,b,c,d) &
                + t2inp(a,b,i,j,spin_t2aa)*int2x(i,j,c,d,spin_int2aa)
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
     do i = 1,norb1
        do j = 1,norb1
           i1ab(a,b,c,d) = i1ab(a,b,c,d) &
                + t2inp(a,b,i,j,spin_t2ab)*int2x(i,j,c,d,spin_int2ab)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_man11_1
!##########################################################
