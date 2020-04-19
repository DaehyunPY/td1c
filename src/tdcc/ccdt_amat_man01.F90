!##########################################################
subroutine ccdt_amat_man01(i0,work1,work2)

! i0 ( a i b j )_ytt + = -1/2 * Sum ( c ) 
!  * i1 ( c b )_yt * t ( a c i j )_t 1

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1((norb1+1):nact,(norb1+1):nact),work2(1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: irot,jrot,ai,bj

  work1 = 0d0
  call ccdt_amat_man01_1(work1)

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
     do c = norb1+1,nact
        i0(a,i,b,j,1) = i0(a,i,b,j,1) &
             - work1(b,c)*t2inp(a,c,i,j,spin_t2aa)
        i0(a,i,b,j,2) = i0(a,i,b,j,2) &
             - work1(b,c)*t2inp(a,c,i,j,spin_t2ab)
     end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_amat_man01
!##########################################################
subroutine ccdt_amat_man01_1(i1)

! i1 ( a b )_yt + = +1/2 * Sum ( i j c ) 
!  * y ( i j b c )_y * t ( a c i j )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared)
  !$omp do collapse(2)
  do a = norb1+1,nact
  do b = norb1+1,nact
     do c = norb1+1,nact
     do i = 1,norb1
        do j = 1,i-1
           i1(a,b) = i1(a,b) &
                + g2inp(i,j,b,c,spin_g2aa)*t2inp(a,c,i,j,spin_t2aa)
        end do
        do j = 1,norb1
           i1(a,b) = i1(a,b) &
                + g2inp(i,j,b,c,spin_g2ab)*t2inp(a,c,i,j,spin_t2ab)
        end do
     end do
     end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_amat_man01_1
!##########################################################
