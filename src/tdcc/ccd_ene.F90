!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_ene_1(i0)

! i0 ( )_vt + = +1/4 * Sum ( h3 h4 p1 p2 ) * t ( p1 p2 h3 h4 )_t * v ( h3 h4 p1 p2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit
  use mod_ormas, only : nact
  use mod_cc, only : norb1,t2inp
  use mod_cc, only : int2x

  implicit none
  complex(kind(0d0)), intent(inout) :: i0
  integer(c_int) :: h3, h4, p1, p2, sh3, sh4, sp1, sp2
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int), external :: tdcc_spin_t2inp
  integer(c_int), external :: tdcc_spin_int2x
  complex(kind(0d0)), parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  do sh3 = 1, 2
  do sh4 = 1, 2
  do sp1 = 1, 2
  do sp2 = 1, 2
     spin_t2inp = tdcc_spin_t2inp(sp1, sp2, sh3, sh4)
     spin_int2x = tdcc_spin_int2x(sh3, sh4, sp1, sp2)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h3 = 1, norb1
     do h4 = 1, norb1
     do p1 = norb1+1, nact
     do p2 = norb1+1, nact
        i0 = i0 + fact * t2inp(p1, p2, h3, h4, spin_t2inp) * int2x(h3, h4, p1, p2, spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  end do
end subroutine ccd_ene_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
complex(kind(0d0)) function ccd_ene_main()

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero

  implicit none
  complex(kind(0d0)) :: cene

  cene = czero
  call ccd_ene_1(cene)
  ccd_ene_main = cene

end function ccd_ene_main
!**********************************************************
