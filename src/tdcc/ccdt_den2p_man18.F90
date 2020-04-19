!##########################################################
subroutine ccdt_den2p_man18(i0ab,i0ba,work1,work2,work3)

! i0 ( a i b c )_yt + = +1/2 * Sum ( j k d ) 
!  * y ( j k i d b c )_y * t ( d a j k )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: work1(1),work2(1),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, novvvab
     a = p2_ovvvab(icc)
     i = h1_ovvvab(icc)
     b = p3_ovvvab(icc)
     c = p4_ovvvab(icc)
     tmp = 0d0
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           tmp = tmp + g3inp(j,k,i,d,b,c,spin_g3aab)*t2inp(d,a,j,k,spin_t2aa)
        end do
        do k = 1,norb1
           tmp = tmp + g3inp(j,i,k,d,c,b,spin_g3aab)*t2inp(d,a,j,k,spin_t2ab)
        end do
     end do
     end do
     i0ab(a,i,b,c) = i0ab(a,i,b,c) + tmp
     i0ba(a,i,c,b) = i0ba(a,i,c,b) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man18
!##########################################################
