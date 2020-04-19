!##########################################################
subroutine ccdt_den2p_man25(i0ab,i0ba,work1,work2,work3)

!  i0 ( a b i j )_ytt + = -1/2 * P( i j ) * Sum ( k l c ) 
!  * i1 ( k l c i )_yt * t ( a b c k j l )_t 1

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
  call ccdt_den2p_man25_1(work1,work2)
  call tdcc_filloovoaa(work1)

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, noovvab
     i = h1_oovvab(icc)
     j = h2_oovvab(icc)
     a = p3_oovvab(icc)
     b = p4_oovvab(icc)
     tmp = 0d0
     do c = norb1+1,nact
     do k = 1,norb1
        do l = 1,k-1
           tmp = tmp &
                - work1(k,l,c,i)*t3inp(a,c,b,k,l,j,spin_t3aab) &
                - work1(k,l,c,j)*t3inp(b,c,a,k,l,i,spin_t3aab)
        end do
        do l = 1,norb1
           tmp = tmp &
                + work2(l,k,c,i)*t3inp(b,c,a,j,l,k,spin_t3aab) &
                + work2(k,l,c,j)*t3inp(a,c,b,i,k,l,spin_t3aab)
        end do
     end do
     end do
     i0ab(a,b,i,j) = i0ab(a,b,i,j) + tmp
     i0ba(a,b,j,i) = i0ba(a,b,j,i) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man25
!##########################################################
subroutine ccdt_den2p_man25_1(i1aa,i1ab)

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

end subroutine ccdt_den2p_man25_1
!##########################################################
