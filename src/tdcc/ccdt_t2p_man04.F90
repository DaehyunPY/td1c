!##########################################################
subroutine ccdt_t2p_man04(i0,work1,work2,work3)

!1: i0 ( a b i j )_v + = 1 * v ( a b i j )_v 0
!6: i0 ( a b i j )_vt + = 1/2 * Sum ( c d ) * t ( c d i j )_t * v ( a b c d )_v 0
!7: i0 ( a b i j )_tf + = 1 * Sum ( c k ) * t ( a b c i j k )_t * f ( k c )_f 0
!8: i0 ( a b i j )_vt + = -1/2 * P( 2 ) * Sum ( k l c ) * t ( a b c i k l )_t * v ( k l j c )_v 0
!9: i0 ( a b i j )_vt + = -1/2 * P( 2 ) * Sum ( k c d ) * t ( a c d i j k )_t * v ( k b c d )_v 0

  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc2aa,ncc2ab,t2inp,t3inp,cc_rank
  use mod_cc,only : h1_cc2aa,h2_cc2aa,p1_cc2aa,p2_cc2aa
  use mod_cc,only : h1_cc2ab,h2_cc2ab,p1_cc2ab,p2_cc2ab
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:2)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(a,b,i,j,tmp)
  ! ##### CCD contributions #####
  !$omp do
  do icc = 1,ncc2aa
     a = p1_cc2aa(icc)
     b = p2_cc2aa(icc)
     i = h1_cc2aa(icc)
     j = h2_cc2aa(icc)
     tmp = 0d0
     ! diagram 1
     tmp = tmp + int2x(a,b,i,j,spin_int2aa)
     ! diagram 6
     do c = norb1+1,nact
     do d = norb1+1,c-1
        tmp = tmp + t2inp(c,d,i,j,spin_t2aa)*int2x(a,b,c,d,spin_int2aa)
     end do
     end do
     i0(a,b,i,j,1) = i0(a,b,i,j,1) + tmp
  end do
  !$omp end do
  !$omp do
  do icc = 1,ncc2ab
     a = p1_cc2ab(icc)
     b = p2_cc2ab(icc)
     i = h1_cc2ab(icc)
     j = h2_cc2ab(icc)
     tmp = 0d0
     ! diagram 1
     tmp = tmp + int2x(a,b,i,j,spin_int2ab)
     ! diagram 6
     do c = norb1+1,nact
     do d = norb1+1,nact
        tmp = tmp + t2inp(c,d,i,j,spin_t2ab)*int2x(a,b,c,d,spin_int2ab)
     end do
     end do
     i0(a,b,i,j,2) = i0(a,b,i,j,2) + tmp
  end do
  !$omp end do

  ! ##### CCDT contributions #####
  if (cc_rank >= 3) then
     !$omp do
     do icc = 1,ncc2aa
        a = p1_cc2aa(icc)
        b = p2_cc2aa(icc)
        i = h1_cc2aa(icc)
        j = h2_cc2aa(icc)
        tmp = 0d0
        ! diagram 7
        do c = norb1+1,nact
        do k = 1,norb1
           tmp = tmp + (t3inp(a,b,c,i,j,k,spin_t3aaa) &
                       +t3inp(a,b,c,i,j,k,spin_t3aab))*fock(k,c,spin_focka)
        end do
        end do
        ! diagram 8
        do c = norb1+1,nact
        do k = 1,norb1
           do l = 1,k-1
              tmp = tmp - t3inp(a,b,c,i,k,l,spin_t3aaa)*int2x(k,l,j,c,spin_int2aa) &
                        + t3inp(a,b,c,j,k,l,spin_t3aaa)*int2x(k,l,i,c,spin_int2aa)
           end do
           do l = 1,norb1
              tmp = tmp - t3inp(a,b,c,i,k,l,spin_t3aab)*int2x(k,l,j,c,spin_int2ab) &
                        + t3inp(a,b,c,j,k,l,spin_t3aab)*int2x(k,l,i,c,spin_int2ab)
           end do
        end do
        end do
        ! diagram 9
        do k = 1,norb1
        do c = norb1+1,nact
           do d = norb1+1,c-1
              tmp = tmp - t3inp(a,c,d,i,j,k,spin_t3aaa)*int2x(k,b,c,d,spin_int2aa) &
                        + t3inp(b,c,d,i,j,k,spin_t3aaa)*int2x(k,a,c,d,spin_int2aa)
           end do
           do d = norb1+1,nact
              tmp = tmp + t3inp(a,c,d,i,j,k,spin_t3aab)*int2x(b,k,c,d,spin_int2ab) &
                        - t3inp(b,c,d,i,j,k,spin_t3aab)*int2x(a,k,c,d,spin_int2ab)
           end do
        end do
        end do
        i0(a,b,i,j,1) = i0(a,b,i,j,1) + tmp
     end do
     !$omp end do
     !$omp do
     do icc = 1,ncc2ab
        a = p1_cc2ab(icc)
        b = p2_cc2ab(icc)
        i = h1_cc2ab(icc)
        j = h2_cc2ab(icc)
        tmp = 0d0
        ! diagram 7
        do c = norb1+1,nact
        do k = 1,norb1
           tmp = tmp + (t3inp(a,c,b,i,k,j,spin_t3aab) &
                       +t3inp(b,c,a,j,k,i,spin_t3aab))*fock(k,c,spin_focka)
        end do
        end do
        ! diagram 8
        do c = norb1+1,nact
        do k = 1,norb1
           do l = 1,k-1
              tmp = tmp - t3inp(b,c,a,k,l,i,spin_t3aab)*int2x(k,l,j,c,spin_int2aa) &
                        + t3inp(a,c,b,l,k,j,spin_t3aab)*int2x(k,l,i,c,spin_int2aa)
           end do
           do l = 1,norb1
              tmp = tmp - t3inp(a,c,b,i,l,k,spin_t3aab)*int2x(k,l,j,c,spin_int2ab) &
                        + t3inp(c,b,a,j,l,k,spin_t3aab)*int2x(k,l,i,c,spin_int2ab)
           end do
        end do
        end do
        ! diagram 9
        do k = 1,norb1
        do c = norb1+1,nact
           do d = norb1+1,c-1
              tmp = tmp - t3inp(c,d,a,j,k,i,spin_t3aab)*int2x(k,b,c,d,spin_int2aa) &
                        + t3inp(d,c,b,i,k,j,spin_t3aab)*int2x(k,a,c,d,spin_int2aa)
           end do
           do d = norb1+1,nact
              tmp = tmp + t3inp(a,c,d,i,k,j,spin_t3aab)*int2x(k,b,c,d,spin_int2ab) &
                        - t3inp(b,d,c,k,j,i,spin_t3aab)*int2x(a,k,c,d,spin_int2ab)
           end do
        end do
        end do
        i0(a,b,i,j,2) = i0(a,b,i,j,2) + tmp
     end do
     !$omp end do
  end if
  !$omp end parallel

end subroutine ccdt_t2p_man04
!##########################################################
