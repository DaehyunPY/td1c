!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_bmat2p_1(sh1, sp2, i0)

!  i0 ( h1 p2 )_ydt + = +1 * Sum ( h3 p4 ) * y ( h3 p4 )_y * dt ( p2 p4 h1 h3 )_dt 0

  use mod_const, only : czero, runit
  use mod_ormas, only : nact
  use mod_cc, only : norb1,g1inp
  use mod_cc, only : dt2inp

  implicit none
  integer, intent(in) :: sh1, sp2
  complex(kind(0d0)), intent(inout) :: i0(1:nact, 1:nact)
  integer :: h1, p2
  integer :: h3, p4, sh3, sp4
  integer :: spin_g1inp_1
  integer :: spin_dt2inp_2
  integer, external :: tdcc_spin_g1inp
  integer, external :: tdcc_spin_t2inp
  complex(kind(0d0)), parameter :: fact = +1.0d+0 * runit

  do sh3 = 1, 2
  do sp4 = 1, 2
     spin_g1inp_1 = tdcc_spin_g1inp(sh3, sp4)
     spin_dt2inp_2 = tdcc_spin_t2inp(sp2, sp4, sh1, sh3)
     if(spin_g1inp_1 * spin_dt2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1, norb1
     do p2 = norb1+1, nact
     do h3 = 1, norb1
     do p4 = norb1+1, nact
        i0(h1, p2) = i0(h1, p2) + fact * g1inp(h3, p4, spin_g1inp_1) * dt2inp(p2, p4, h1, h3, spin_dt2inp_2)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccd_bmat2p_1
!##########################################################
!##########################################################
!##########################################################
