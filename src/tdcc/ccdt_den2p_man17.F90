!##########################################################
subroutine ccdt_den2p_man17(i0ab,i0ba,work1,work2,work3)

! i0 ( i a j k )_ytt + = - P( j k ) * Sum ( l b ) 
!  * i1 ( i l b j )_yt * t ( a b k l )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1(1:norb1,1:norb1,(norb1+1):nact,1:norb1), &
       work2(1:norb1,1:norb1,(norb1+1):nact,1:norb1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den2p_man17_1(work1,work2)
  call tdcc_filloovoaa(work1)

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, nooovab
     i = h3_ooovab(icc)
     a = p4_ooovab(icc)
     j = h1_ooovab(icc)
     k = h2_ooovab(icc)
     tmp = 0d0
     do l = 1,norb1
     do b = norb1+1,nact
        tmp = tmp &
             - work1(i,l,b,j)*t2inp(a,b,k,l,spin_t2ab) &
             + work2(l,i,b,j)*t2inp(a,b,k,l,spin_t2aa) &
             - work2(i,l,b,k)*t2inp(a,b,l,j,spin_t2ab)
     end do
     end do
     i0ab(i,a,j,k) = i0ab(i,a,j,k) + tmp
     i0ba(i,a,k,j) = i0ba(i,a,k,j) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man17
!##########################################################
subroutine ccdt_den2p_man17_1(i1aa,i1ab)

! i1 ( i j a k )_yt + = +1/2 * Sum ( l b c ) 
!  * y ( l i j b c a )_y * t ( b c l k )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,1:norb1,(norb1+1):nact,1:norb1), &
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
                + g3inp(l,i,j,b,c,a,spin_g3aaa)*t2inp(b,c,l,k,spin_t2aa)
        end do
        do c = norb1+1,nact
           i1aa(i,j,a,k) = i1aa(i,j,a,k) &
                + g3inp(i,j,l,c,a,b,spin_g3aab)*t2inp(b,c,l,k,spin_t2ab)
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
                - g3inp(l,j,i,b,c,a,spin_g3aab)*t2inp(b,c,l,k,spin_t2aa)
        end do
        do c = norb1+1,nact
           i1ab(i,j,a,k) = i1ab(i,j,a,k) &
                - g3inp(l,i,j,b,a,c,spin_g3aab)*t2inp(b,c,l,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man17_1
!##########################################################
