!##########################################################
subroutine ccdt_den2p_man03(i0ab,i0ba,work1,work2,work3)

!  i0 ( a b i j )_t + = +1 * t ( a b i j )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,a,b,c,d,tmp)
  !$omp do
  do icc = 1, noovvab
     i = h1_oovvab(icc)
     j = h2_oovvab(icc)
     a = p3_oovvab(icc)
     b = p4_oovvab(icc)
     i0ab(a,b,i,j) = i0ab(a,b,i,j) + t2inp(a,b,i,j,spin_t2ab)
     i0ba(a,b,j,i) = i0ba(a,b,j,i) - t2inp(a,b,i,j,spin_t2ab)
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man03
!##########################################################
