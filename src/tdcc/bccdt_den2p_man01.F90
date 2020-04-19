!##########################################################
subroutine bccdt_den2p_man01(i0ab,i0ba,work1,work2,work3)

!3:  i0 ( a b c i )_yt + = +1 * Sum ( j ) * y ( j c )_y * t ( a b j i )_t 0
!4:  i0 ( i a j k )_yt + = -1 * Sum ( b ) * y ( i b )_y * t ( b a j k )_t 0
!22: i0 ( a b i j )_yt + = +1 * Sum ( k c ) * y ( k c )_y * t ( a b c i j k )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,t3inp,g1inp,cc_rank
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,tmp)
  ! diagram 3
  !$omp do
  do icc = 1, novvvab
     a = p3_ovvvab(icc)
     b = p4_ovvvab(icc)
     c = p2_ovvvab(icc)
     i = h1_ovvvab(icc)
     tmp = 0d0
     do j = 1,norb1
        tmp = tmp + g1inp(j,c,spin_g1a)*t2inp(a,b,j,i,spin_t2ab)
     end do
     i0ab(a,b,c,i) = i0ab(a,b,c,i) + tmp
     i0ba(b,a,c,i) = i0ba(b,a,c,i) - tmp
  end do
  !$omp end do

  ! diagram 4
  !$omp do
  do icc = 1, nooovab
     i = h3_ooovab(icc)
     a = p4_ooovab(icc)
     j = h1_ooovab(icc)
     k = h2_ooovab(icc)
     tmp = 0d0
     do b = norb1+1,nact
        tmp = tmp - g1inp(i,b,spin_g1a)*t2inp(b,a,j,k,spin_t2ab)
     end do
     i0ab(i,a,j,k) = i0ab(i,a,j,k) + tmp
     i0ba(i,a,k,j) = i0ba(i,a,k,j) - tmp
  end do
  !$omp end do

  ! diagram 22
  if (cc_rank >= 3) then
     !$omp do
     do icc = 1, noovvab
        a = p3_oovvab(icc)
        b = p4_oovvab(icc)
        i = h1_oovvab(icc)
        j = h2_oovvab(icc)
        tmp = 0d0
        do k = 1,norb1
        do c = norb1+1,nact
           tmp = tmp + g1inp(k,c,spin_g1a) &
                *(t3inp(a,c,b,i,k,j,spin_t3aab) &
                + t3inp(b,c,a,j,k,i,spin_t3aab))
        end do
        end do
        i0ab(a,b,i,j) = i0ab(a,b,i,j) + tmp
        i0ba(a,b,j,i) = i0ba(a,b,j,i) - tmp
     end do
     !$omp end do
  end if
  !$omp end parallel

end subroutine bccdt_den2p_man01
!##########################################################
