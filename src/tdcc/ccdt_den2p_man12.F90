!##########################################################
subroutine ccdt_den2p_man12(i0ab,i0ba,work1,work2,work3)

!  i0 ( a b c i )_yt + = +1/2 * Sum ( j k d ) * y ( j k d c )_y * t ( d a b j k i )_t 0

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
  do icc = 1, novvvab
     a = p3_ovvvab(icc)
     b = p4_ovvvab(icc)
     c = p2_ovvvab(icc)
     i = h1_ovvvab(icc)
     tmp = 0d0
     do d = norb1+1,act1_ul
     do j = act1_ll,norb1
        do k = act1_ll,j-1
           tmp = tmp &
                + g2inp(j,k,d,c,spin_g2aa)*t3inp(d,a,b,j,k,i,spin_t3aab)
        end do
        do k = act1_ll,norb1
           tmp = tmp &
                + g2inp(j,k,d,c,spin_g2ab)*t3inp(d,b,a,j,i,k,spin_t3aab)
        end do
     end do
     end do
     i0ab(a,b,c,i) = i0ab(a,b,c,i) + tmp
     i0ba(b,a,c,i) = i0ba(b,a,c,i) - tmp
!ba not needed     do d = norb1+1,nact
!ba not needed     do j = 1,norb1
!ba not needed        do k = 1,j-1
!ba not needed           i0ab(a,b,c,i) = i0ab(a,b,c,i) &
!ba not needed                + g2inp(j,k,d,c,spin_g2aa)*t3inp(d,a,b,j,k,i,spin_t3aab)
!ba not needed           i0ba(a,b,c,i) = i0ba(a,b,c,i) &
!ba not needed                - g2inp(j,k,d,c,spin_g2aa)*t3inp(d,b,a,j,k,i,spin_t3aab)
!ba not needed        end do
!ba not needed        do k = 1,norb1
!ba not needed           i0ab(a,b,c,i) = i0ab(a,b,c,i) &
!ba not needed                + g2inp(j,k,d,c,spin_g2ab)*t3inp(d,b,a,j,i,k,spin_t3aab)
!ba not needed           i0ba(a,b,c,i) = i0ba(a,b,c,i) &
!ba not needed                - g2inp(j,k,d,c,spin_g2ab)*t3inp(d,a,b,j,i,k,spin_t3aab)
!ba not needed        end do
!ba not needed     end do
!ba not needed     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man12
!##########################################################
