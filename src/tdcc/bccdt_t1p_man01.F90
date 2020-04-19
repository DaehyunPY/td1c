!##########################################################
subroutine bccdt_t1p_man01(i0)

!1: i0 ( a i )_f + = 1 * f ( a i )_f 0
!2: i0 ( a i )_tf + = 1 * Sum ( b j ) * t ( a b i j )_t * f ( j b )_f 0
!3: i0 ( a i )_vt + = -1/2 * Sum ( j k b ) * t ( a b j k )_t * v ( j k i b )_v 0
!4: i0 ( a i )_vt + = -1/2 * Sum ( j b c ) * t ( b c i j )_t * v ( j a b c )_v 0
!5: i0 ( a i )_vt + = 1/4 * Sum ( j k b c ) * t ( a b c i j k )_t * v ( j k b c )_v 0

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : fock,int2x,norb1,ncc1a,t2inp,t3inp
  use mod_cc,only : h1_cc1a,p1_cc1a,cc_rank
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,1:norb1)

  integer(c_int) :: icc,a,b,c,d,e,i,j,k,l,m

  !##################################################
  !$omp parallel default(shared) private(a,i)
  !$omp do
  do icc = 1,ncc1a
     a = p1_cc1a(icc)
     i = h1_cc1a(icc)
     ! diagram 1
     i0(a,i) = i0(a,i) + fock(a,i,spin_focka)
     ! diagram 2
     do b = norb1+1,nact
     do j = 1,norb1
        i0(a,i) = i0(a,i) + fock(j,b,spin_focka) *( &
             t2inp(a,b,i,j,spin_t2aa) &
           + t2inp(a,b,i,j,spin_t2ab))
     end do
     end do
     ! diagram 3
     do b = norb1+1,nact
     do j = 1,norb1
        do k = 1,j-1
           i0(a,i) = i0(a,i) - t2inp(a,b,j,k,spin_t2aa)*int2x(j,k,i,b,spin_int2aa)
        end do
        do k = 1,norb1
           i0(a,i) = i0(a,i) - t2inp(a,b,j,k,spin_t2ab)*int2x(j,k,i,b,spin_int2ab)
        end do
     end do
     end do
     ! diagram 4
     do j = 1,norb1
     do b = norb1+1,nact
        do c = norb1+1,b-1
           i0(a,i) = i0(a,i) + t2inp(b,c,i,j,spin_t2aa)*int2x(a,j,b,c,spin_int2aa)
        end do
        do c = norb1+1,nact
           i0(a,i) = i0(a,i) + t2inp(b,c,i,j,spin_t2ab)*int2x(a,j,b,c,spin_int2ab)
        end do
     end do
     end do
     ! diagram 5
     if (cc_rank >= 3) then
        do j = 1,norb1
        do b = norb1+1,nact
           do k = 1,j-1
           do c = norb1+1,b-1
              i0(a,i) = i0(a,i) + int2x(j,k,b,c,spin_int2aa) *( &
                   t3inp(a,b,c,i,j,k,spin_t3aaa) &
                 + t3inp(b,c,a,j,k,i,spin_t3aab))
           end do
           end do
        
           do k = 1,norb1
           do c = norb1+1,nact
              i0(a,i) = i0(a,i) + t3inp(a,b,c,i,j,k,spin_t3aab)*int2x(j,k,b,c,spin_int2ab)
           end do
           end do
        end do
        end do
     end if
  end do
  !$omp end do
  !$omp end parallel
  !##################################################

end subroutine bccdt_t1p_man01
!##########################################################
