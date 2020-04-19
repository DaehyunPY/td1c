!##########################################################
subroutine ccdt_den2p_man24(i0ab,i0ba,work1,work2,work3)

!6:  i0 ( a b i j )_ytt + = - P( a b ) * Sum ( c ) * i1 ( a c )_yt * t ( c b i j )_t 1
!24: i0 ( a b i j )_ytt + = - P( a b ) * Sum ( c ) * i1 ( a c )_yt * t ( c b i j )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1((norb1+1):nact,(norb1+1):nact),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  call ccdt_den2p_man24_1(work1)

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, noovvab
     i = h1_oovvab(icc)
     j = h2_oovvab(icc)
     a = p3_oovvab(icc)
     b = p4_oovvab(icc)
     tmp = 0d0
     do c = norb1+1,nact
        tmp = tmp &
             - work1(a,c)*t2inp(c,b,i,j,spin_t2ab) &
             - work1(b,c)*t2inp(a,c,i,j,spin_t2ab)
     end do
     i0ab(a,b,i,j) = i0ab(a,b,i,j) + tmp
     i0ba(a,b,j,i) = i0ba(a,b,j,i) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man24
!##########################################################
subroutine ccdt_den2p_man24_1(i1)

!5-1:  i1 ( a b )_yt + = +1/2 * Sum ( i j c ) * y ( i j c b )_y * t ( c a i j )_t 0
!24-1: i1 ( a b )_yt + = +1/12 * Sum ( i j k c d ) 
!  * y ( i j k c d b )_y * t ( c d a i j k )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, nvv
     a = p1_vv(icc)
     b = p2_vv(icc)
     ! diagram 6-1
     do c = norb1+1,nact
     do i = 1,norb1
        do j = 1,i-1
           i1(a,b) = i1(a,b) + g2inp(i,j,c,b,spin_g2aa)*t2inp(c,a,i,j,spin_t2aa)
        end do
        do j = 1,norb1
           i1(a,b) = i1(a,b) + g2inp(i,j,c,b,spin_g2ab)*t2inp(c,a,i,j,spin_t2ab)
        end do
     end do
     end do

     ! diagram 24-1
     do c = norb1+1,nact
     do i = 1,norb1
        do d = norb1+1,c-1
        do j = 1,i-1
        do k = 1,j-1
           i1(a,b) = i1(a,b) &
                + g3inp(i,j,k,c,d,b,spin_g3aaa)*t3inp(c,d,a,i,j,k,spin_t3aaa)
        end do
        end do
        end do

        do d = norb1+1,c-1
        do j = 1,norb1
        do k = 1,j-1
           i1(a,b) = i1(a,b) &
                + g3inp(j,k,i,c,d,b,spin_g3aab)*t3inp(c,d,a,j,k,i,spin_t3aab)
        end do
        end do
        end do

        do d = norb1+1,nact
        do j = 1,norb1
        do k = 1,j-1
           i1(a,b) = i1(a,b) &
                + g3inp(k,j,i,c,b,d,spin_g3aab)*t3inp(c,a,d,k,j,i,spin_t3aab)
        end do
        end do
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man24_1
!##########################################################
