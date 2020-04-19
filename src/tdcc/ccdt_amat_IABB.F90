!##########################################################
subroutine ccdt_amat_IABB(i0,work1,work2,work3)

!1 i0 ( a2 i | b1 a2 )_yt <= -1/4 * Sum ( k l m c d e )
!  * y ( k l m c d e )_y * t ( c d k i )_t * t ( b1 e l m )

!1 i0 ( a1 i | b2 a1 )_yt <= -1/4 * Sum ( k l m c d e ) 
!  * y ( k l m c d e )_y * t ( c d k i )_t * t ( b2 e l m )

!2 i0 ( a i | b2 b1 )_yt <= -1/2 * Sum ( k l m c d ) 
!  * y ( k l m b1 c d )_y * t ( c b2 i k )_t * t ( a d l m )

!3 i0 ( a i | b2 b1 )_yt <= -1/4 * Sum ( k l m c d ) 
!  * y ( k l m b1 c d )_y * t ( b2 a k l )_t * t ( c d i m )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1:*),work2(1:*),work3(1:*)

  call ccdt_amat_IABB_1(i0,work1,work2,work3)
  call ccdt_amat_IABB_2(i0,work1,work2,work3)
  call ccdt_amat_IABB_3(i0,work1,work2,work3)

end subroutine ccdt_amat_IABB
!##########################################################
subroutine ccdt_amat_IABB_1(i0,work1,work2,work3)

!1 i0 ( a2 i | b1 a2 )_yt <= -1/4 * Sum ( k l m c d e )
!  * y ( k l m c d e )_y * t ( c d k i )_t * t ( b1 e l m )

!1 i0 ( a1 i | b2 a1 )_yt <= -1/4 * Sum ( k l m c d e ) 
!  * y ( k l m c d e )_y * t ( c d k i )_t * t ( b2 e l m )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1,norb1+1:nact,1:norb1),&
       work2(1:norb1,1:norb1,norb1+1:nact,1:norb1),work3(1)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  work1 = 0d0
  work2 = 0d0
  call ccdt_amat_IABB_1_1(work1,work2)

  !$omp parallel default(shared) private(a2,b1,i,j,k,l,m,a,b,c,d,e,tmp1,tmp2)
  !$omp do collapse(2)
  do b1 = norb1+1,act1_ul
  do i = 1,norb1
     tmp1 = 0d0
     do e = norb1+1,act1_ul
     do l = act1_ll,norb1
        do m = act1_ll,l-1
           tmp1 = tmp1 - work1(l,m,e,i)*t2inp(b1,e,l,m,spin_t2aa)
        end do
        do m = act1_ll,norb1
           tmp1 = tmp1 - work2(l,m,e,i)*t2inp(b1,e,l,m,spin_t2ab)
        end do
     end do
     end do
     do a2 = act1_ul+1,nact
        i0(a2,i,b1,a2,1) = i0(a2,i,b1,a2,1) + tmp1
     end do
  end do
  end do
  !$omp end do

  !$omp do collapse(2)
  do b2 = act1_ul+1,nact
  do i = 1,norb1
     tmp1 = 0d0
     do e = norb1+1,act1_ul
     do l = act1_ll,norb1
        do m = act1_ll,l-1
           tmp1 = tmp1 - work1(l,m,e,i)*t2inp(b2,e,l,m,spin_t2aa)
        end do
        do m = act1_ll,norb1
           tmp1 = tmp1 - work2(l,m,e,i)*t2inp(b2,e,l,m,spin_t2ab)
        end do
     end do
     end do
     do a1 = norb1+1,act1_ul
        i0(a1,i,b2,a1,1) = i0(a1,i,b2,a1,1) + tmp1
     end do
  end do
  end do
  !$omp end do

  !$omp end parallel
end subroutine ccdt_amat_IABB_1
!##########################################################
subroutine ccdt_amat_IABB_1_1(work1,work2)

! i1 ( l m e i )_yt = 1/2 * Sum ( k c d )
!  * y ( k l m c d e )_y * t ( c d k i )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1,norb1+1:nact,1:norb1),&
       work2(1:norb1,1:norb1,norb1+1:nact,1:norb1)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  !$omp parallel default(shared) private(i,j,k,l,m,a,b,c,d,e,tmp1,tmp2)
  !$omp do collapse(4)
  do i = 1,norb1
  do e = norb1+1,act1_ul
  do m = act1_ll,norb1
  do l = act1_ll,norb1
     tmp1 = 0d0
     tmp2 = 0d0
     do k = act1_ll,norb1
     do c = norb1+1,act1_ul
        do d = norb1+1,c-1
           tmp1 = tmp1 + g3inp(k,l,m,c,d,e,spin_g3aaa)*t2inp(c,d,k,i,spin_t2aa)
           tmp2 = tmp2 + g3inp(k,l,m,c,d,e,spin_g3aab)*t2inp(c,d,k,i,spin_t2aa)
        end do
        do d = norb1+1,act1_ul
           tmp1 = tmp1 + g3inp(l,m,k,d,e,c,spin_g3aab)*t2inp(c,d,k,i,spin_t2ab)
           tmp2 = tmp2 + g3inp(k,m,l,c,e,d,spin_g3aab)*t2inp(c,d,k,i,spin_t2ab)
        end do
     end do
     end do
     work1(l,m,e,i) = work1(l,m,e,i) + tmp1
     work2(l,m,e,i) = work2(l,m,e,i) + tmp2
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_amat_IABB_1_1
!##########################################################
subroutine ccdt_amat_IABB_2(i0,work1,work2,work3)

!2 i0 ( a i | b2 b1 )_yt <= -1/2 * Sum ( k l m c d ) 
!  * y ( k l m b1 c d )_y * t ( c b2 i k )_t * t ( a d l m )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,norb1+1:nact,norb1+1:nact,norb1+1:nact),&
       work2(1:norb1,norb1+1:nact,norb1+1:nact,norb1+1:nact),&
       work3(1:norb1,norb1+1:nact,norb1+1:nact,norb1+1:nact)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m,irot,ai
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  work1 = 0d0
  work2 = 0d0
  work3 = 0d0
  call ccdt_amat_IABB_2_1(work1,work2,work3)

  !$omp parallel default(shared) private(ai,irot,j1,j2,i,k,a,b,tmp1,tmp2)  
  !$omp do collapse(3)
  do b1 = norb1+1,act1_ul
  do b2 = act1_ul+1,nact
  do ai = 1,nXai
     irot = aiX(ai)
     a = rotaa_mapb(irot,1)
     i = rotaa_mapb(irot,2)
     tmp1 = 0d0
     tmp2 = 0d0
     do c = norb1+1,act1_ul
     do k = act1_ll,norb1
        tmp1 = tmp1 + t2inp(b2,c,i,k,spin_t2aa)*work1(k,a,b1,c) &
                    + t2inp(b2,c,i,k,spin_t2ab)*work2(k,a,b1,c)
        tmp2 = tmp2 - t2inp(c,b2,i,k,spin_t2ab)*work3(k,a,b1,c)
     end do
     end do
     i0(a,i,b2,b1,1) = i0(a,i,b2,b1,1) + tmp1
     i0(a,i,b2,b1,2) = i0(a,i,b2,b1,2) + tmp2
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_amat_IABB_2
!##########################################################
subroutine ccdt_amat_IABB_2_1(work1,work2,work3)

! i1 ( k a b1 c )_yt = 1/2 * Sum ( d l m )
!  * y ( k l m b1 c d )_y * t ( a d l m )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,norb1+1:nact,norb1+1:nact,norb1+1:nact),&
       work2(1:norb1,norb1+1:nact,norb1+1:nact,norb1+1:nact),&
       work3(1:norb1,norb1+1:nact,norb1+1:nact,norb1+1:nact)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2,tmp3

  !$omp parallel default(shared) private(i,j,k,l,m,a,b,c,d,e,tmp1,tmp2,tmp3)
  !$omp do collapse(4)
  do c = norb1+1,act1_ul
  do b = norb1+1,act1_ul
  do a = norb1+1,nact
  do k = act1_ll,norb1
     tmp1 = 0d0
     tmp2 = 0d0
     tmp3 = 0d0
     do d = norb1+1,act1_ul
     do l = act1_ll,norb1
        do m = act1_ll,l-1
           tmp1 = tmp1 + g3inp(k,l,m,b,c,d,spin_g3aaa)*t2inp(a,d,l,m,spin_t2aa)
           tmp2 = tmp2 + g3inp(m,l,k,b,d,c,spin_g3aab)*t2inp(a,d,l,m,spin_t2aa)
           tmp3 = tmp3 + g3inp(l,m,k,c,d,b,spin_g3aab)*t2inp(a,d,l,m,spin_t2aa)
        end do
        do m = act1_ll,norb1
           tmp1 = tmp1 + g3inp(k,l,m,b,c,d,spin_g3aab)*t2inp(a,d,l,m,spin_t2ab)
           tmp2 = tmp2 + g3inp(k,m,l,d,c,b,spin_g3aab)*t2inp(a,d,l,m,spin_t2ab)
           tmp3 = tmp3 + g3inp(k,m,l,b,d,c,spin_g3aab)*t2inp(a,d,l,m,spin_t2ab)
        end do
     end do
     end do
     work1(k,a,b,c) = work1(k,a,b,c) + tmp1
     work2(k,a,b,c) = work2(k,a,b,c) + tmp2
     work3(k,a,b,c) = work3(k,a,b,c) + tmp3
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_amat_IABB_2_1
!##########################################################
subroutine ccdt_amat_IABB_3(i0,work1,work2,work3)

!2 i0 ( a i | b2 b1 )_yt <= -1/4 * Sum ( k l m c d ) 
!  * y ( k l m b1 c d )_y * t ( b2 a k l )_t * t ( c d i m )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1,1:norb1,norb1+1:nact),&
       work2(1:norb1,1:norb1,1:norb1,norb1+1:nact),work3(1)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m,irot,ai
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  work1 = 0d0
  work2 = 0d0
  call ccdt_amat_IABB_3_1(work1,work2)

  !$omp parallel default(shared) private(ai,irot,j1,j2,i,k,a,b,tmp1,tmp2)  
  !$omp do collapse(3)
  do b1 = norb1+1,act1_ul
  do b2 = act1_ul+1,nact
  do ai = 1,nXai
     irot = aiX(ai)
     a = rotaa_mapb(irot,1)
     i = rotaa_mapb(irot,2)
     tmp1 = 0d0
     tmp2 = 0d0
     do k = act1_ll,norb1
        do l = act1_ll,k-1
           tmp1 = tmp1 + work1(k,l,i,b1)*t2inp(a,b2,k,l,spin_t2aa)
        end do
        do l = act1_ll,norb1
           tmp2 = tmp2 + work2(k,l,i,b1)*t2inp(a,b2,k,l,spin_t2ab)
        end do
     end do
     i0(a,i,b2,b1,1) = i0(a,i,b2,b1,1) + tmp1
     i0(a,i,b2,b1,2) = i0(a,i,b2,b1,2) + tmp2
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_amat_IABB_3
!##########################################################
subroutine ccdt_amat_IABB_3_1(work1,work2)

! i1 ( k l i b1 )_yt = 1/2 * Sum ( m c d )
!  * y ( k l m b1 c d )_y * t ( c d i m )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1,1:norb1,norb1+1:nact),&
       work2(1:norb1,1:norb1,1:norb1,norb1+1:nact)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  !$omp parallel default(shared) private(i,j,k,l,m,a,b,c,d,e,tmp1,tmp2)
  !$omp do collapse(4)
  do b = norb1+1,act1_ul
  do i = 1,norb1
  do l = act1_ll,norb1
  do k = act1_ll,norb1
     tmp1 = 0d0
     tmp2 = 0d0
     do m = act1_ll,norb1
     do c = norb1+1,act1_ul
        do d = norb1+1,c-1
           tmp1 = tmp1 + g3inp(k,l,m,b,c,d,spin_g3aaa)*t2inp(c,d,i,m,spin_t2aa)
           tmp2 = tmp2 + g3inp(k,m,l,d,c,b,spin_g3aab)*t2inp(c,d,i,m,spin_t2aa)
        end do
        do d = norb1+1,act1_ul
           tmp1 = tmp1 + g3inp(k,l,m,b,c,d,spin_g3aab)*t2inp(c,d,i,m,spin_t2ab)
           tmp2 = tmp2 + g3inp(m,l,k,b,d,c,spin_g3aab)*t2inp(c,d,i,m,spin_t2ab)
        end do
     end do
     end do
     work1(k,l,i,b) = work1(k,l,i,b) + tmp1
     work2(k,l,i,b) = work2(k,l,i,b) + tmp2
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_amat_IABB_3_1
!##########################################################
