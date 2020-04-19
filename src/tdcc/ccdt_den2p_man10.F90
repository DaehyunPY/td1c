!##########################################################
subroutine ccdt_den2p_man10(i0ab,i0ba,work1,work2,work3)

!1:  i0 ( a b c d )_yt + = +1/2 * Sum ( i j ) * y ( i j c d )_y * t ( a b i j )_t 0
!10: i0 ( a b c d )_yt + = +1/6 * Sum ( i j k e ) * y ( i j k c d e )_y * t ( a b e i j k )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp,cc_rank
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, nvvvvab
     a = p1_vvvvab(icc)
     b = p2_vvvvab(icc)
     c = p3_vvvvab(icc)
     d = p4_vvvvab(icc)
     tmp = 0d0
     ! diagram 1
     do i = 1,norb1
     do j = 1,norb1
        tmp = tmp + g2inp(i,j,c,d,spin_g2ab)*t2inp(a,b,i,j,spin_t2ab)
     end do
     end do

     ! diagram 10
     if (cc_rank >= 3) then
        do e = norb1+1,nact
        do i = 1,norb1
        do j = 1,norb1
           do k = 1,j-1
              tmp = tmp &
                   + g3inp(k,j,i,c,e,d,spin_g3aab)*t3inp(a,e,b,k,j,i,spin_t3aab) &
                   + g3inp(j,k,i,d,e,c,spin_g3aab)*t3inp(b,e,a,j,k,i,spin_t3aab)
           end do
        end do
        end do
        end do
     end if
     i0ab(a,b,c,d) = i0ab(a,b,c,d) + tmp
     i0ba(a,b,d,c) = i0ba(a,b,d,c) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man10
!##########################################################
