!##########################################################
subroutine ccdt_den2p_man01(i0ab,i0ba,work1,work2,work3)

!  i0 ( a b c d )_yt + = +1/2 * Sum ( i j ) * y ( i j c d )_y * t ( a b i j )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,tmp)
  !$omp do
  do icc = 1, nvvvvab
     a = p1_vvvvab(icc)
     b = p2_vvvvab(icc)
     c = p3_vvvvab(icc)
     d = p4_vvvvab(icc)
     tmp = 0d0
     do i = 1,norb1
     do j = 1,norb1
        tmp = tmp + g2inp(i,j,c,d,spin_g2ab)*t2inp(a,b,i,j,spin_t2ab)
     end do
     end do
     i0ab(a,b,c,d) = i0ab(a,b,c,d) + tmp
     i0ba(a,b,d,c) = i0ba(a,b,d,c) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man01
!##########################################################
