!##########################################################
subroutine ccdt_den2p_man16(i0ab,i0ba,work1,work2,work3)

!  i0 ( a b c i )_ytt + = P( a b ) 
!  * Sum ( j d ) * i1 ( j a c d )_yt * t ( b d i j )_t 1

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: i0ab(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: i0ba(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: &!work1(1),work2(1),work3(1)
       work1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       work2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact),work3(1)
  complex(kind(0d0)) :: tmp

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  work1 = 0d0
  work2 = 0d0
  call ccdt_den2p_man16_1(work1,work2)
  call tdcc_fillovvvaa(work1)

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e,tmp)
  !$omp do
  do icc = 1, novvvab
     a = p3_ovvvab(icc)
     b = p4_ovvvab(icc)
     c = p2_ovvvab(icc)
     i = h1_ovvvab(icc)
     tmp = 0d0
     do j = 1,norb1
        do d = norb1+1,nact
           tmp = tmp &
                + work1(j,a,c,d)*t2inp(b,d,i,j,spin_t2ab) &
                - work2(j,a,d,c)*t2inp(b,d,i,j,spin_t2aa) &
                + work2(j,b,c,d)*t2inp(a,d,j,i,spin_t2ab)
        end do
     end do
     i0ab(a,b,c,i) = i0ab(a,b,c,i) + tmp
     i0ba(b,a,c,i) = i0ba(b,a,c,i) - tmp
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man16
!##########################################################
subroutine ccdt_den2p_man16_1(i1aa,i1ab)

! i1 ( i a b c )_yt + = +1/2 * Sum ( j k d ) 
!  * y ( j k i d b c )_y * t ( d a j k )_t 0

  use mod_ormas,only : nact
  use mod_cc,only : norb1,t2inp,g2inp,t3inp,g3inp
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i1aa(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact), &
       i1ab(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_long) :: icc,a,b,c,d,e,i,j,k,l,m

  !$omp parallel default(shared) private(icc,i,j,k,l,m,a,b,c,d,e)
  !$omp do
  do icc = 1, novvvaa
     i = h1_ovvvaa(icc)
     a = p2_ovvvaa(icc)
     b = p3_ovvvaa(icc)
     c = p4_ovvvaa(icc)
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + g3inp(j,k,i,d,b,c,spin_g3aaa)*t2inp(d,a,j,k,spin_t2aa)
        end do
        do k = 1,norb1
           i1aa(i,a,b,c) = i1aa(i,a,b,c) &
                + g3inp(i,k,j,c,b,d,spin_g3aab)*t2inp(d,a,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp do
  do icc = 1, novvvab
     i = h1_ovvvab(icc)
     a = p2_ovvvab(icc)
     b = p3_ovvvab(icc)
     c = p4_ovvvab(icc)
     do d = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                - g3inp(j,k,i,d,c,b,spin_g3aab)*t2inp(d,a,j,k,spin_t2aa)
        end do
        do k = 1,norb1
           i1ab(i,a,b,c) = i1ab(i,a,b,c) &
                - g3inp(j,i,k,d,b,c,spin_g3aab)*t2inp(d,a,j,k,spin_t2ab)
        end do
     end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine ccdt_den2p_man16_1
!##########################################################
