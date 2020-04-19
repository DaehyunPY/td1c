!##########################################################
subroutine ccdt_den2p_man02(i0ab,i0ba,work1,work2,work3)

!  i0 ( i j k l )_yt + = +1/2 * Sum ( a b ) * y ( i j a b )_y * t ( a b k l )_t 0

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
  do icc = 1, nooooab
     i = h1_ooooab(icc)
     j = h2_ooooab(icc)
     k = h3_ooooab(icc)
     l = h4_ooooab(icc)
     tmp = 0d0
     do a = norb1+1,nact
     do b = norb1+1,nact
        tmp = tmp + g2inp(i,j,a,b,spin_g2ab)*t2inp(a,b,k,l,spin_t2ab)
     end do
     end do
     i0ab(i,j,k,l) = i0ab(i,j,k,l) + tmp
     i0ba(i,j,l,k) = i0ab(i,j,l,k) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man02
!##########################################################
