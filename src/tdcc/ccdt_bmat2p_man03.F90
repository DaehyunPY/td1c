!##########################################################
subroutine ccdt_bmat2p_man03(i0,work1,work2)

! i0 ( i a )_ytdt + = -1/4 * Sum ( j k b ) 
!  * i1 ( j k b i )_yt * dt ( b a j k )_dt 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,dt2inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1)
       work1(1:norb1,1:norb1,(norb1+1):nact,1:norb1),&
       work2(1:norb1,1:norb1,(norb1+1):nact,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_bmat2p_man03_1(work1,work2)
  call tdcc_filloovoaa(work1)

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d)
  !$omp do
  do icc = 1, nov
     i = h1_ov(icc)
     a = p2_ov(icc)
     do b = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i0(i,a) = i0(i,a) &
                - work1(j,k,b,i)*dt2inp(b,a,j,k,spin_t2aa)
        end do
        do k = 1,norb1
           i0(i,a) = i0(i,a) &
                - work2(j,k,b,i)*dt2inp(b,a,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_bmat2p_man03
!##########################################################
subroutine ccdt_bmat2p_man03_1(i1aa,i1ab)

! i1 ( i j a k )_yt + = +1 * Sum ( l b c ) 
!  * y ( i j l a b c )_y * t ( b c k l )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,g3inp,t2inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,(norb1+1):nact,1:norb1),&
       i1ab(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nooovaa
     i = h1_ooovaa(icc)
     j = h2_ooovaa(icc)
     a = p4_ooovaa(icc)
     k = h3_ooovaa(icc)
     do l = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i1aa(i,j,a,k) = i1aa(i,j,a,k) &
                + g3inp(i,j,l,a,b,c,spin_g3aaa)*t2inp(b,c,k,l,spin_t2aa)
        end do
        do c = norb1+1,nact
           i1aa(i,j,a,k) = i1aa(i,j,a,k) &
                + g3inp(i,j,l,a,b,c,spin_g3aab)*t2inp(b,c,k,l,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, nooovab
     i = h1_ooovab(icc)
     j = h2_ooovab(icc)
     a = p4_ooovab(icc)
     k = h3_ooovab(icc)
     do l = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                + g3inp(j,l,i,b,c,a,spin_g3aab)*t2inp(b,c,k,l,spin_t2aa)
        end do
        do c = norb1+1,nact
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                + g3inp(i,l,j,a,c,b,spin_g3aab)*t2inp(b,c,k,l,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_bmat2p_man03_1
!##########################################################
