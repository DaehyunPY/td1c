!##########################################################
subroutine ccdt_l3p_man09(i0,work1,work2,work3)

! i0 ( i j k a b c )_ytv + = -1/2 * P( 9 ) * Sum ( d ) 
!  * i1 ( i d a b )_yt * v ( j k c d )_v 1

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: &
       work1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),&
       work2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
!  call ccdt_l3p_9_1(1,1,1,1,work1)
!  call ccdt_l3p_9_1(1,2,1,2,work2)
!  work1 = -0.5d0*work1
!  work2 = -0.5d0*work2
  call ccdt_l3p_man09_1(work1,work2)
  call tdcc_fillovvvaa(work1)

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
             + work1(i,d,a,b)*int2x(j,k,c,d,spin_int2aa) &
             - work1(j,d,a,b)*int2x(i,k,c,d,spin_int2aa) &
             - work1(k,d,a,b)*int2x(j,i,c,d,spin_int2aa) &
             - work1(i,d,c,b)*int2x(j,k,a,d,spin_int2aa) &
             - work1(i,d,a,c)*int2x(j,k,b,d,spin_int2aa) &
             + work1(j,d,c,b)*int2x(i,k,a,d,spin_int2aa) &
             + work1(k,d,c,b)*int2x(j,i,a,d,spin_int2aa) &
             + work1(j,d,a,c)*int2x(i,k,b,d,spin_int2aa) &
             + work1(k,d,a,c)*int2x(j,i,b,d,spin_int2aa)
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
             - work1(i,d,a,b)*int2x(j,k,d,c,spin_int2ab) &
             + work1(j,d,a,b)*int2x(i,k,d,c,spin_int2ab) &
!ZERO        - work1(k,d,a,b)*int2x(j,i,c,d,spin_int2aa) &
             + work2(i,d,b,c)*int2x(j,k,a,d,spin_int2ab) &
             - work2(i,d,a,c)*int2x(j,k,b,d,spin_int2ab) &
             - work2(j,d,b,c)*int2x(i,k,a,d,spin_int2ab) &
             + work2(k,d,c,b)*int2x(j,i,a,d,spin_int2aa) &
             + work2(j,d,a,c)*int2x(i,k,b,d,spin_int2ab) &
             - work2(k,d,c,a)*int2x(j,i,b,d,spin_int2aa)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man09
!##########################################################
subroutine ccdt_l3p_man09_1(i1aa,i1ab)

! i1 ( i a b c )_yt + = -0.5 * Sum ( j k d ) 
!  * t ( d a j k )_t * y ( i j k b c d )_y 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),&
       i1ab(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novvvaa
     i = h1_ovvvaa(icc)
     a = p2_ovvvaa(icc)
     b = p3_ovvvaa(icc)
     c = p4_ovvvaa(icc)
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                - t2inp(d,a,j,k,spin_t2aa)*g3inp(i,j,k,b,c,d,spin_g3aaa)
        end do
        do k = 1,norb1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + t2inp(a,d,j,k,spin_t2ab)*g3inp(i,j,k,b,c,d,spin_g3aab)
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
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                - t2inp(d,a,j,k,spin_t2aa)*g3inp(j,k,i,c,d,b,spin_g3aab)
        end do
        do k = 1,norb1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                + t2inp(d,a,j,k,spin_t2ab)*g3inp(i,j,k,b,d,c,spin_g3aab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_man09_1
!##########################################################
