!##########################################################
subroutine ccdt_l3p_man05(i0,work1,work2,work3)

! i0 ( i j k a b c )_yf + = 1 * P( 3 ) * Sum ( d ) 
!  * y ( i j k a b d )_y * i1 ( d c )_f 2

  use mod_ormas,only : nact
  use mod_cc,only : int2x,norb1,ncc3aaa,ncc3aab,t2inp,g2inp,t3inp,g3inp
  use mod_cc,only : h1_cc3aaa,h2_cc3aaa,h3_cc3aaa,p1_cc3aaa,p2_cc3aaa,p3_cc3aaa
  use mod_cc,only : h1_cc3aab,h2_cc3aab,h3_cc3aab,p1_cc3aab,p2_cc3aab,p3_cc3aab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:2)
  complex(kind(0d0)) :: work1((norb1+1):nact,(norb1+1):nact),work2(1),work3(1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
!  call ccdt_l3p_5_1(1,1,work1)
!  call ccdt_l3p_5_2(1,1,work1)
  call ccdt_l3p_man05_1(work1)

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
             + g3inp(i,j,k,a,b,d,spin_g3aaa)*work1(d,c) &
             - g3inp(i,j,k,c,b,d,spin_g3aaa)*work1(d,a) &
             - g3inp(i,j,k,a,c,d,spin_g3aaa)*work1(d,b)
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
             + g3inp(i,j,k,a,b,d,spin_g3aab)*work1(d,c) &
             + g3inp(i,j,k,d,b,c,spin_g3aab)*work1(d,a) &
             + g3inp(i,j,k,a,d,c,spin_g3aab)*work1(d,b)
     end do
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine ccdt_l3p_man05
!##########################################################
subroutine ccdt_l3p_man05_1(i1)

! 5-1
! i1 ( a b )_f + = 1 * f ( a b )_f 0
! 5-2
! i1 ( p7 p1 )_vt + = -1/2 * Sum ( h9 h10 p8 ) 
!  * t ( p7 p8 h9 h10 )_t * v ( h9 h10 p1 p8 )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : t2inp,fock,int2x,norb1
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nvv
     a = p1_vv(icc)
     b = p2_vv(icc)
     ! 5-1
     i1(a,b) = i1(a,b) + fock(a,b,spin_focka)

     ! 5-2
     do c = norb1+1,nact
     do i = 1,norb1
        do j = 1,i-1
           i1(a,b) = i1(a,b) &
                - t2inp(a,c,i,j,spin_t2aa)*int2x(i,j,b,c,spin_int2aa)
        end do
        do j = 1,norb1
           i1(a,b) = i1(a,b) &
                - t2inp(a,c,i,j,spin_t2ab)*int2x(i,j,b,c,spin_int2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_man05_1
!##########################################################
