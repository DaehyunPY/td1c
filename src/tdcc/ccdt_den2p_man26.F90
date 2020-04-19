!##########################################################
subroutine ccdt_den2p_man26(i0ab,i0ba,work1,work2,work3)

!D##!04: i0 ( a b i j )_ytt + = +1/2 * P( i j ) * P( a b ) * Sum ( k c ) * i1 ( k a c i )_yt * t ( b c j k )_t 1 <== CHECK THIS!
!D##!20: i0 ( a b i j )_ytt + =        P( i j ) * P( a b ) * Sum ( k c ) * i1 ( k a c i )_yt * t ( b c j k )_t 1
!
!B##!05: i0 ( a b i j )_ytt + = -1/2 * P( i j ) * Sum ( k ) * i1 ( k i )_yt * t ( a b k j )_t 1
!B##!23: i0 ( a b i j )_ytt + = - P( i j ) * Sum ( k ) * i1 ( k i )_yt * t ( a b k j )_t 1
!
!C##!06: i0 ( a b i j )_ytt + = P( a b ) * Sum ( c ) * i1 ( a c )_yt * t ( c b i j )_t 1
!C##!24: i0 ( a b i j )_ytt + = - P( a b ) * Sum ( c ) * i1 ( a c )_yt * t ( c b i j )_t 1
!
!A##!07: i0 ( a b i j )_ytt + = +1/2 * Sum ( k l ) * i1 ( k l i j )_yt * t ( a b k l )_t 1
!A##!22: i0 ( a b i j )_ytt + = +1/2 * Sum ( k l ) * i1 ( k l i j )_yt * t ( a b k l )_t 1
!
!###!21: i0 ( a b i j )_ytt + = +1/2 * Sum ( c d ) * i1 ( a b c d )_yt * t ( c d i j )_t 1
!
!###!25: i0 ( a b i j )_ytt + = -1/2 * P( i j ) * Sum ( k l c ) * i1 ( k l c i )_yt * t ( a b c k j l )_t 1

!03: i0 ( a b i j )_t + = +1 * t ( a b i j )_t 0
!26: i0 ( a b i j )_ytt + = -1/4 * P( a b ) * Sum ( k c d ) * i1 ( k c d a )_yt * t ( c b d i j k )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,cc_rank
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       work2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den2p_man26_1(work1,work2)
  call tdcc_fillovvvaa(work1)

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, noovvab
     i = h1_oovvab(icc)
     j = h2_oovvab(icc)
     a = p3_oovvab(icc)
     b = p4_oovvab(icc)
     tmp = 0d0
     ! diagram 3
     tmp = tmp + t2inp(a,b,i,j,spin_t2ab)

     ! diagram 26
     do k = 1,norb1
     do c = norb1+1,nact
        do d = norb1+1,c-1
           tmp = tmp &
                - work1(k,a,c,d)*t3inp(c,d,b,i,k,j,spin_t3aab) &
                + work1(k,b,c,d)*t3inp(c,d,a,k,j,i,spin_t3aab)
        end do
        do d = norb1+1,nact
           tmp = tmp &
                + work2(k,a,d,c)*t3inp(b,d,c,j,k,i,spin_t3aab) &
                - work2(k,b,c,d)*t3inp(c,a,d,i,k,j,spin_t3aab)
        end do
     end do
     end do
     i0ab(a,b,i,j) = i0ab(a,b,i,j) + tmp
     i0ba(a,b,j,i) = i0ba(a,b,j,i) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man26
!##########################################################
subroutine ccdt_den2p_man26_1(i1aa,i1ab)

! i1 ( i a b c )_yt + = +1 * Sum ( j k d ) 
!  * y ( j k i d a b )_y * t ( d c j k )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       i1ab(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novvvaa
     i = h1_ovvvaa(icc)
     a = p2_ovvvaa(icc)
     b = p3_ovvvaa(icc)
     c = p4_ovvvaa(icc)
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + g3inp(j,k,i,d,b,c,spin_g3aaa)*t2inp(d,a,j,k,spin_t2aa)
        end do
        do k = 1,norb1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + g3inp(i,k,j,c,b,d,spin_g3aab)*t2inp(d,a,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, novvvab
     i = h1_ovvvab(icc)
     a = p2_ovvvab(icc)
     b = p3_ovvvab(icc)
     c = p4_ovvvab(icc)
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                - g3inp(j,k,i,d,c,b,spin_g3aab)*t2inp(d,a,j,k,spin_t2aa)
        end do
        do k = 1,norb1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                - g3inp(j,i,k,d,b,c,spin_g3aab)*t2inp(d,a,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man26_1
!##########################################################
