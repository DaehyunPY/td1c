!##########################################################
subroutine ccdt_den2p_man15(i0ab,i0ba,work1,work2,work3)

! i0 ( i a j k )_ytt + = +1/2 * Sum ( b c ) 
!  * i1 ( i a b c )_yt * t ( b c j k )_t 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       work2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den2p_man15_1(work1,work2)
  !call tdcc_filovvvaa(work1) ! not needed

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, nooovab
     i = h3_ooovab(icc)
     a = p4_ooovab(icc)
     j = h1_ooovab(icc)
     k = h2_ooovab(icc)
     tmp = 0d0
     do b = norb1+1,nact
     do c = norb1+1,nact
        tmp = tmp + work2(i,a,b,c)*t2inp(b,c,j,k,spin_t2ab)
     end do
     end do
     i0ab(i,a,j,k) = i0ab(i,a,j,k) + tmp
     i0ba(i,a,k,j) = i0ba(i,a,k,j) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man15
!##########################################################
subroutine ccdt_den2p_man15_1(i1aa,i1ab)

! i1 ( i a b c )_yt + = +1/2 * Sum ( j k d ) 
!  * y ( i j k b c d )_y * t ( a d j k )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       i1ab(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  ! not needed
  i1aa = 0d0

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novvvab
     i = h1_ovvvab(icc)
     a = p2_ovvvab(icc)
     b = p3_ovvvab(icc)
     c = p4_ovvvab(icc)
     do d = norb1+1,act1_ul
     do j = act1_ll,norb1
        do k = act1_ll,j-1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) + &
             g3inp(j,k,i,c,d,b,spin_g3aab)*t2inp(a,d,j,k,spin_t2aa)
        end do
        do k = act1_ll,norb1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) + &
             g3inp(i,k,j,b,d,c,spin_g3aab)*t2inp(a,d,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man15_1
!##########################################################
