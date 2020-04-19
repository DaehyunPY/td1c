!##########################################################
subroutine ccdt_den2p_man21(i0ab,i0ba,work1,work2,work3)

!  i0 ( a b i j )_ytt + = +1/2 * Sum ( c d ) 
!  * i1 ( a b c d )_yt * t ( c d i j )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       work2((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den2p_man21_1(work1,work2)
  ! call tdcc_fillvvvvaa(work1) !i1aa not needed

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, noovvab
     i = h1_oovvab(icc)
     j = h2_oovvab(icc)
     a = p3_oovvab(icc)
     b = p4_oovvab(icc)
     tmp = 0d0
     do c = norb1+1,nact
     do d = norb1+1,nact
        tmp = tmp + work2(a,b,c,d)*t2inp(c,d,i,j,spin_t2ab)
     end do
     end do
     i0ab(a,b,i,j) = i0ab(a,b,i,j) + tmp
     i0ba(a,b,j,i) = i0ba(a,b,j,i) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man21
!##########################################################
subroutine ccdt_den2p_man21_1(i1aa,i1ab)

! i1 ( a b c d )_yt + = +1/6 * Sum ( i j k e ) 
!  * y ( i j k e c d )_y * t ( e a b i j k )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       i1ab((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  ! not needed
  i1aa = 0

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nvvvvab
     a = p1_vvvvab(icc)
     b = p2_vvvvab(icc)
     c = p3_vvvvab(icc)
     d = p4_vvvvab(icc)
     do e = norb1+1,nact
     do i = 1,norb1
     do j = 1,norb1
     do k = 1,j-1
        i1ab(a,b,c,d) = i1ab(a,b,c,d) &
             + g3inp(k,j,i,e,c,d,spin_g3aab)*t3inp(e,a,b,k,j,i,spin_t3aab) &
             + g3inp(k,j,i,e,d,c,spin_g3aab)*t3inp(e,b,a,k,j,i,spin_t3aab)
     end do
     end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man21_1
!##########################################################
