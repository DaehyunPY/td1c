!##########################################################
subroutine ccdt_den2p_man13(i0ab,i0ba,work1,work2,work3)

!  i0 ( i a j k )_yt + = -1/2 * Sum ( l b c ) * y ( l i b c )_y * t ( b c a l j k )_t 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,act1_ll,act1_ul
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, nooovab
     i = h3_ooovab(icc)
     a = p4_ooovab(icc)
     j = h1_ooovab(icc)
     k = h2_ooovab(icc)
     tmp = 0d0
     do l = act1_ll,norb1
     do b = norb1+1,act1_ul
        do c = norb1+1,b-1
           tmp = tmp - g2inp(l,i,b,c,spin_g2aa)*t3inp(b,c,a,l,j,k,spin_t3aab)
        end do
        do c = norb1+1,act1_ul
           tmp = tmp - g2inp(l,i,b,c,spin_g2ab)*t3inp(a,b,c,k,l,j,spin_t3aab)
        end do
     end do
     end do
     i0ab(i,a,j,k) = i0ab(i,a,j,k) + tmp
     i0ba(i,a,k,j) = i0ba(i,a,k,j) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man13
!##########################################################
