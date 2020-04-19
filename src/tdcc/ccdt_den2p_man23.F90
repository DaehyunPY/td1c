!##########################################################
subroutine ccdt_den2p_man23(i0ab,i0ba,work1,work2,work3)

!5:  i0 ( a b i j )_ytt + = +1/2 * P( i j ) * Sum ( k ) * i1 ( k i )_yt * t ( a b k j )_t 1
!23: i0 ( a b i j )_ytt + = - P( i j ) * Sum ( k ) * i1 ( k i )_yt * t ( a b k j )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1(1:norb1,1:norb1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  call ccdt_den2p_man23_1(work1)

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, noovvab
     i = h1_oovvab(icc)
     j = h2_oovvab(icc)
     a = p3_oovvab(icc)
     b = p4_oovvab(icc)
     tmp = 0d0
     do k = 1,norb1
        tmp = tmp &
             - work1(k,i)*t2inp(a,b,k,j,spin_t2ab) &
             - work1(k,j)*t2inp(a,b,i,k,spin_t2ab)
     end do
     i0ab(a,b,i,j) = i0ab(a,b,i,j) + tmp
     i0ba(a,b,j,i) = i0ba(a,b,j,i) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man23
!##########################################################
subroutine ccdt_den2p_man23_1(i1)

!5-1:  i1 ( i j )_yt + = -1 * Sum ( k a b ) * y ( k i a b )_y * t ( a b k j )_t 0
!23-1: i1 ( i j )_yt + = +1/12 * Sum ( k l a b c ) 
!  * y ( k l i a b c )_y * t ( a b c k l j )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, noo
     i = h1_oo(icc)
     j = h2_oo(icc)
     ! diagram 5-1
     do k = 1,norb1
     do a = norb1+1,nact
        do b = norb1+1,a-1
           i1(i,j) = i1(i,j) + g2inp(k,i,a,b,spin_g2aa)*t2inp(a,b,k,j,spin_t2aa)
        end do
        do b = norb1+1,nact
           i1(i,j) = i1(i,j) + g2inp(k,i,a,b,spin_g2ab)*t2inp(a,b,k,j,spin_t2ab)
        end do
     end do
     end do

     ! diagram 23-1
     do k = 1,norb1
     do a = norb1+1,nact
        do l = 1,k-1
        do b = norb1+1,a-1
        do c = norb1+1,b-1
           i1(i,j) = i1(i,j) + g3inp(k,l,i,a,b,c,spin_g3aaa)*t3inp(a,b,c,k,l,j,spin_t3aaa)
        end do
        end do
        end do

        do l = 1,k-1
        do b = norb1+1,nact
        do c = norb1+1,b-1
           i1(i,j) = i1(i,j) + g3inp(k,l,i,b,c,a,spin_g3aab)*t3inp(b,c,a,k,l,j,spin_t3aab)
        end do
        end do
        end do

        do l = 1,norb1
        do b = norb1+1,nact
        do c = norb1+1,b-1
           i1(i,j) = i1(i,j) + g3inp(k,i,l,c,b,a,spin_g3aab)*t3inp(c,b,a,k,j,l,spin_t3aab)
        end do
        end do
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man23_1
!##########################################################
