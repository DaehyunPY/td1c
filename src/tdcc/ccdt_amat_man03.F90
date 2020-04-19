!##########################################################
subroutine ccdt_amat_man03(i0,work1,work2)

! i0 ( a b i j )_ytt + = -1/4 * Sum ( k c d ) 
!  * i1 ( k c d a )_yt * t ( b c d j k i )_t 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul,rotaa_mapb
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),&
       work2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: irot,jrot,ai,bj

  work1 = 0d0
  work2 = 0d0
  call ccdt_amat_man03_1(work1,work2)

  !$omp parallel default(shared) private(i,j,k,l,a,b,c,d,ai,bj,irot,jrot)
  !$omp do collapse(2)
  do ai = 1, nXai
  do bj = 1, nXai
     irot = aiX(ai)
     a = rotaa_mapb(irot,1)
     i = rotaa_mapb(irot,2)
     jrot = aiX(bj)
     b = rotaa_mapb(jrot,1)
     j = rotaa_mapb(jrot,2)
     do k = act1_ll,norb1
     do c = norb1+1,act1_ul
        do d = norb1+1,c-1
           i0(a,i,b,j,1) = i0(a,i,b,j,1) &
                - work1(k,a,c,d)*t3inp(b,c,d,j,k,i,spin_t3aaa)
           i0(a,i,b,j,2) = i0(a,i,b,j,2) &
                - work1(k,a,c,d)*t3inp(d,c,b,i,k,j,spin_t3aab)
        end do
        do d = norb1+1,act1_ul
           i0(a,i,b,j,1) = i0(a,i,b,j,1) &
                - work2(k,a,c,d)*t3inp(b,d,c,j,i,k,spin_t3aab)
           i0(a,i,b,j,2) = i0(a,i,b,j,2) &
                - work2(k,a,c,d)*t3inp(b,c,d,j,k,i,spin_t3aab)
        end do
     end do
     end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_amat_man03
!##########################################################
subroutine ccdt_amat_man03_1(i1aa,i1ab)

! i1 ( i a b c )_yt + = +1/2 * Sum ( j k d ) 
!  * y ( i j k b c d )_y * t ( a d j k )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),&
       i1ab(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared)
  !$omp do collapse(4)
  do i = act1_ll,norb1
  do a = norb1+1,nact
  do b = norb1+1,act1_ul
  do c = norb1+1,act1_ul
     do d = norb1+1,act1_ul
     do j = act1_ll,norb1
        do k = act1_ll,j-1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + g3inp(i,j,k,b,c,d,spin_g3aaa)*t2inp(a,d,j,k,spin_t2aa)
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                + g3inp(j,k,i,c,d,b,spin_g3aab)*t2inp(a,d,j,k,spin_t2aa)
        end do
        do k = act1_ll,norb1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + g3inp(i,j,k,b,c,d,spin_g3aab)*t2inp(a,d,j,k,spin_t2ab)
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                + g3inp(i,k,j,b,d,c,spin_g3aab)*t2inp(a,d,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_amat_man03_1
!##########################################################
