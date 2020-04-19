!##########################################################
subroutine ccdt_den2p_man19(i0ab,i0ba,work1,work2,work3)

! i0 ( i j k a )_yt + = -1/2 * Sum ( l b c ) 
!  * y ( l i j b c a )_y * t ( b c l k )_t 0

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
  do icc = 1, nooovab
     i = h1_ooovab(icc)
     j = h2_ooovab(icc)
     k = h3_ooovab(icc)
     a = p4_ooovab(icc)
     tmp = 0d0
     do l = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           tmp = tmp - g3inp(l,i,j,b,c,a,spin_g3aab)*t2inp(b,c,l,k,spin_t2aa)
        end do
        do c = norb1+1,nact
           tmp = tmp - g3inp(l,j,i,b,a,c,spin_g3aab)*t2inp(b,c,l,k,spin_t2ab)
        end do
     end do
     end do
     i0ab(i,j,k,a) = i0ab(i,j,k,a) + tmp
     i0ba(j,i,k,a) = i0ba(j,i,k,a) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man19
!##########################################################
