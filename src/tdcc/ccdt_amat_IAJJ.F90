!##########################################################
subroutine ccdt_amat_IAJJ(i0,work1,work2,work3)

!1 i0 ( a i2 | i2 j1 )_yt + = 1/4 * Sum ( k l m c d e )
!  * y ( k l m c d e )_y * t ( a d k l )_t * t ( c e j1 m )

!1 i0 ( a i1 | i1 j2 )_yt + = 1/4 * Sum ( k l m c d e ) 
!  * y ( k l m c d e )_y * t ( a d k l )_t * t ( c e j2 m )

!2 i0 ( a i | j1 j2 )_yt + = 1/2 * Sum ( k l b c d ) 
!  * y ( j1 k l b c d )_y * t ( a b k j2 )_t * t ( c d i l )

!3 i0 ( a i | j1 j2 )_yt + = 1/4 * Sum ( k l b c d ) 
!  * y ( j1 k l b c d )_y * t ( c b i j2 )_t * t ( a d k l )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1:*),work2(1:*),work3(1:*)

  call ccdt_amat_IAJJ_1(i0,work1,work2,work3)
  call ccdt_amat_IAJJ_2(i0,work1,work2,work3)
  call ccdt_amat_IAJJ_3(i0,work1,work2,work3)

end subroutine ccdt_amat_IAJJ
!##########################################################
subroutine ccdt_amat_IAJJ_1(i0,work1,work2,work3)

! i0 ( a i2 | i2 j1 )_yt + = 1/4 * Sum ( k l m c d e )
!  * y ( k l m c d e )_y * t ( a d k l )_t * t ( c e j1 m )

! i0 ( a i1 | i1 j2 )_yt + = 1/4 * Sum ( k l m c d e ) 
!  * y ( k l m c d e )_y * t ( a d k l )_t * t ( c e j2 m )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &
       work1(1:norb1,1:norb1,1:norb1,norb1+1:nact),&
       work2(1:norb1,1:norb1,1:norb1,norb1+1:nact),work3(1)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  work1 = 0d0
  work2 = 0d0
  call ccdt_amat_IAJJ_1_1(work1,work2)

  !$omp parallel default(shared) private(i1,i2,j1,j2,a1,a2,b1,b2,i,j,k,l,m,a,b,c,d,e,tmp1,tmp2)
  !$omp do collapse(2)
  do a = norb1+1,nact
  do j1 = act1_ll,norb1
     tmp1 = 0d0
     do d = norb1+1,act1_ul
     do k = act1_ll,norb1
        do l = act1_ll,k-1
           tmp1 = tmp1 + work1(k,l,j1,d)*t2inp(a,d,k,l,spin_t2aa)
        end do
        do l = act1_ll,norb1
           tmp1 = tmp1 + work2(k,l,j1,d)*t2inp(a,d,k,l,spin_t2ab)
        end do
     end do
     end do
     do i2 = 1,act1_ll-1
        i0(a,i2,i2,j1,1) = i0(a,i2,i2,j1,1) + tmp1
     end do
  end do
  end do
  !$omp end do

  !$omp do collapse(2)
  do a = norb1+1,nact
  do j2 = 1,act1_ll-1
     tmp1 = 0d0
     do d = norb1+1,act1_ul
     do k = act1_ll,norb1
        do l = act1_ll,k-1
           tmp1 = tmp1 + work1(k,l,j2,d)*t2inp(a,d,k,l,spin_t2aa)
        end do
        do l = act1_ll,norb1
           tmp1 = tmp1 + work2(k,l,j2,d)*t2inp(a,d,k,l,spin_t2ab)
        end do
     end do
     end do
     do i1 = act1_ll,norb1
        i0(a,i1,i1,j2,1) = i0(a,i1,i1,j2,1) + tmp1
     end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_amat_IAJJ_1
!##########################################################
subroutine ccdt_amat_IAJJ_1_1(work1,work2)

! i1 ( k1 l1  j d1 )_yt = 1/2 * Sum ( m1 c1 e1 )
!  * y ( k1 l1 m1 c1 d1 e1 )_y * t ( c1 e1 j m1 )

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

  !$omp parallel default(shared) private(i1,i2,j1,j2,a1,a2,b1,b2,i,j,k,l,m,a,b,c,d,e,tmp1,tmp2)
  !$omp do collapse(4)
  do d = norb1+1,act1_ul
  do j = 1,norb1
  do l = act1_ll,norb1
  do k = act1_ll,norb1
     tmp1 = 0d0
     tmp2 = 0d0
     do m = act1_ll,norb1
     do c = norb1+1,act1_ul
        do e = norb1+1,c-1
           tmp1 = tmp1 - g3inp(k,l,m,c,e,d,spin_g3aaa)*t2inp(c,e,j,m,spin_t2aa)
           tmp2 = tmp2 + g3inp(k,m,l,c,e,d,spin_g3aab)*t2inp(c,e,j,m,spin_t2aa)
        end do
        do e = norb1+1,act1_ul
           tmp1 = tmp1 + g3inp(k,l,m,c,d,e,spin_g3aab)*t2inp(c,e,j,m,spin_t2ab)
           tmp2 = tmp2 - g3inp(l,m,k,e,d,c,spin_g3aab)*t2inp(c,e,j,m,spin_t2ab)
        end do
     end do
     end do
     work1(k,l,j,d) = work1(k,l,j,d) + tmp1
     work2(k,l,j,d) = work2(k,l,j,d) + tmp2
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_amat_IAJJ_1_1
!##########################################################
subroutine ccdt_amat_IAJJ_2(i0,work1,work2,work3)

!2 i0 ( a i | j1 j2 )_yt + = 1/2 * Sum ( k l b c d ) 
!  * y ( j1 k l b c d )_y * t ( a b k j2 )_t * t ( c d i l )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1,norb1+1:nact,1:norb1),&
       work2(1:norb1,1:norb1,norb1+1:nact,1:norb1),&
       work3(1:norb1,1:norb1,norb1+1:nact,1:norb1)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m,irot,ai
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  work1 = 0d0
  work2 = 0d0
  work3 = 0d0
  call ccdt_amat_IAJJ_2_1(work1,work2,work3)

  !$omp parallel default(shared) private(ai,irot,j1,j2,i,k,a,b,tmp1,tmp2)  
  !$omp do collapse(3)
  do j2 = 1,act1_ll-1
  do j1 = act1_ll,norb1
  do ai = 1,nXai
     irot = aiX(ai)
     a = rotaa_mapb(irot,1)
     i = rotaa_mapb(irot,2)
     tmp1 = 0d0
     tmp2 = 0d0
     do b = norb1+1,act1_ul
     do k = act1_ll,norb1
        tmp1 = tmp1 + t2inp(a,b,k,j2,spin_t2aa)*work1(j1,k,b,i) &
                    - t2inp(a,b,j2,k,spin_t2ab)*work2(j1,k,b,i)
        tmp2 = tmp2 + t2inp(a,b,k,j2,spin_t2ab)*work3(j1,k,b,i)
     end do
     end do
     i0(a,i,j1,j2,1) = i0(a,i,j1,j2,1) + tmp1
     i0(a,i,j1,j2,2) = i0(a,i,j1,j2,2) + tmp2

  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_amat_IAJJ_2
!##########################################################
subroutine ccdt_amat_IAJJ_2_1(work1,work2,work3)

! i1 ( j k b i )_yt = 1/2 * Sum ( l c d )
!  * y ( j k l b c d )_y * t ( c d i l )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1,norb1+1:nact,1:norb1),&
       work2(1:norb1,1:norb1,norb1+1:nact,1:norb1),&
       work3(1:norb1,1:norb1,norb1+1:nact,1:norb1)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2,tmp3

  !$omp parallel default(shared) private(i,j,k,l,m,a,b,c,d,e,tmp1,tmp2,tmp3)
  !$omp do collapse(4)
  do i = 1,norb1
  do b = norb1+1,act1_ul
  do k = act1_ll,norb1
  do j = act1_ll,norb1
     tmp1 = 0d0
     tmp2 = 0d0
     tmp3 = 0d0
     do l = act1_ll,norb1
     do c = norb1+1,act1_ul
        do d = norb1+1,c-1
           tmp1 = tmp1 + g3inp(j,k,l,b,c,d,spin_g3aaa)*t2inp(c,d,i,l,spin_t2aa)
           tmp2 = tmp2 + g3inp(j,l,k,d,c,b,spin_g3aab)*t2inp(c,d,i,l,spin_t2aa)
           tmp3 = tmp3 + g3inp(k,l,j,c,d,b,spin_g3aab)*t2inp(c,d,i,l,spin_t2aa)
        end do
        do d = norb1+1,act1_ul
           tmp1 = tmp1 + g3inp(j,k,l,b,c,d,spin_g3aab)*t2inp(c,d,i,l,spin_t2ab)
           tmp2 = tmp2 + g3inp(l,k,j,b,d,c,spin_g3aab)*t2inp(c,d,i,l,spin_t2ab)
           tmp3 = tmp3 + g3inp(j,l,k,b,d,c,spin_g3aab)*t2inp(c,d,i,l,spin_t2ab)
        end do
     end do
     end do
     work1(j,k,b,i) = work1(j,k,b,i) + tmp1
     work2(j,k,b,i) = work2(j,k,b,i) + tmp2
     work3(j,k,b,i) = work3(j,k,b,i) + tmp3
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_amat_IAJJ_2_1
!##########################################################
subroutine ccdt_amat_IAJJ_3(i0,work1,work2,work3)

!2 i0 ( a i | j1 j2 )_yt + = 1/4 * Sum ( k l b c d ) 
!  * y ( j1 k l b c d )_y * t ( d a k l )_t * t ( b c i j2 )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact,1:2)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(norb1+1:nact,1:norb1,norb1+1:nact,norb1+1:nact),&
       work2(norb1+1:nact,1:norb1,norb1+1:nact,norb1+1:nact),work3(1)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m,irot,ai
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  work1 = 0d0
  work2 = 0d0
  call ccdt_amat_IAJJ_3_1(work1,work2)

  !$omp parallel default(shared) private(ai,irot,j1,j2,i,k,a,b,tmp1,tmp2)  
  !$omp do collapse(3)
  do j2 = 1,act1_ll-1
  do j1 = act1_ll,norb1
  do ai = 1,nXai
     irot = aiX(ai)
     a = rotaa_mapb(irot,1)
     i = rotaa_mapb(irot,2)
     tmp1 = 0d0
     tmp2 = 0d0
     do b = norb1+1,act1_ul
        do c = norb1+1,b-1
           tmp1 = tmp1 + work1(a,j1,b,c)*t2inp(b,c,i,j2,spin_t2aa)
        end do
        do c = norb1+1,act1_ul
           tmp2 = tmp2 + work2(a,j1,b,c)*t2inp(b,c,i,j2,spin_t2ab)
        end do
     end do
     i0(a,i,j1,j2,1) = i0(a,i,j1,j2,1) + tmp1
     i0(a,i,j1,j2,2) = i0(a,i,j1,j2,2) + tmp2

  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_amat_IAJJ_3
!##########################################################
subroutine ccdt_amat_IAJJ_3_1(work1,work2)

! i1 ( a j b c )_yt = 1/2 * Sum ( d k l )
!  * y ( j k l b c d )_y * t ( d a k l )

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,rotaa_mapb,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,nXai,aiX
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(norb1+1:nact,1:norb1,norb1+1:nact,norb1+1:nact),&
       work2(norb1+1:nact,1:norb1,norb1+1:nact,norb1+1:nact)

  integer(c_int) :: a,b,c,d,e,i,j,k,l,m
  integer(c_int) :: i1,i2,j1,j2,a1,a2,b1,b2
  complex(kind(0d0)) :: tmp1,tmp2

  !$omp parallel default(shared) private(i,j,k,l,m,a,b,c,d,e,tmp1,tmp2)
  !$omp do collapse(4)
  do c = norb1+1,act1_ul
  do b = norb1+1,act1_ul
  do j = act1_ll,norb1
  do a = norb1+1,nact
     tmp1 = 0d0
     tmp2 = 0d0
     do d = norb1+1,act1_ul
     do k = act1_ll,norb1
        do l = act1_ll,k-1
           tmp1 = tmp1 + g3inp(j,k,l,b,c,d,spin_g3aaa)*t2inp(d,a,k,l,spin_t2aa)
           tmp2 = tmp2 + g3inp(k,l,j,b,d,c,spin_g3aab)*t2inp(a,d,k,l,spin_t2aa)
        end do
        do l = act1_ll,norb1
           tmp1 = tmp1 - g3inp(j,k,l,b,c,d,spin_g3aab)*t2inp(a,d,k,l,spin_t2ab)
           tmp2 = tmp2 + g3inp(j,k,l,c,d,b,spin_g3aab)*t2inp(a,d,l,k,spin_t2ab)
        end do
     end do
     end do
     work1(a,j,b,c) = work1(a,j,b,c) + tmp1
     work2(a,j,b,c) = work2(a,j,b,c) + tmp2
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_amat_IAJJ_3_1
!##########################################################
