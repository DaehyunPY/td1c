!##########################################################
subroutine ccdt_t2p_man02(i0,work1,work2,work3)

!4: i0 ( a b i j )_vt + = 1/2 * Sum ( k l ) * t ( a b k l )_t * i1 ( k l i j )_v 2

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,t3inp
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,1:norb1,1:norb1,1:norb1), &
       work2(1:norb1,1:norb1,1:norb1,1:norb1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_t2p_man02_1(work1,work2)
  call tdcc_fillooooaa(work1)

  !$omp parallel default(shared) private(a,b,i,j,tmp)
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     tmp = 0d0
     ! diagram 4
     do k = 1,norb1
     do l = 1,k-1
        tmp = tmp + t2inp(a,b,k,l,spin_t2aa)*work1(k,l,i,j)
     end do
     end do
     i0(a,b,i,j,1) = i0(a,b,i,j,1) + tmp
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     tmp = 0d0
     ! diagram 4
     do k = 1,norb1
     do l = 1,norb1
        tmp = tmp + t2inp(a,b,k,l,spin_t2ab)*work2(k,l,i,j)
     end do
     end do
     i0(a,b,i,j,2) = i0(a,b,i,j,2) + tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t2p_man02
!##########################################################
subroutine ccdt_t2p_man02_1(i1aa,i1ab)

! i1 ( i j k l )_v + = 1 * v ( i j k l )_v 0
! i1 ( i j k l )_vt + = 1/2 * Sum ( a b ) 
!  * t ( a b k l )_t * v ( i j a b )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,int2x
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,1:norb1,1:norb1), &
       i1ab(1:norb1,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooooaa
     i = h1_ooooaa(icc)
     j = h2_ooooaa(icc)
     k = h3_ooooaa(icc)
     l = h4_ooooaa(icc)
     i1aa(i,j,k,l) = i1aa(i,j,k,l) + int2x(i,j,k,l,spin_int2aa)
     do a = norb1+1,nact
        do b = norb1+1,a-1
           i1aa(i,j,k,l) = i1aa(i,j,k,l) + t2inp(a,b,k,l,spin_t2aa)*int2x(i,j,a,b,spin_int2aa)
        end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nooooab
     i = h1_ooooab(icc)
     j = h2_ooooab(icc)
     k = h3_ooooab(icc)
     l = h4_ooooab(icc)
     i1ab(i,j,k,l) = i1ab(i,j,k,l) + int2x(i,j,k,l,spin_int2ab)
     do a = norb1+1,nact
        do b = norb1+1,nact
           i1ab(i,j,k,l) = i1ab(i,j,k,l) + t2inp(a,b,k,l,spin_t2ab)*int2x(i,j,a,b,spin_int2ab)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_t2p_man02_1
!##########################################################
