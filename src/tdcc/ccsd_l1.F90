!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_1(sh2,sp1,i0)

! i0 ( h2 p1 )_f + = 1 * f ( h2 p1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh2,sp1)
  do h2 = 1,norb1
  do p1 = norb1+1,nact
     i0(h2,p1) = i0(h2,p1) + fact * fock(h2,p1,spin_fock)
  end do
  end do
end subroutine ccsd_l1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_2(sh2,sp1,i0)

! i0 ( h2 p1 )_yf + = -1 * Sum ( h7 ) * y ( h7 p1 )_y * i1 ( h2 h7 )_f 4

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h7,sh7
  integer(c_int) :: spin_g1inp
  integer(c_int) :: spin_itm_hh
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh7 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh7,sp1)
     spin_itm_hh = tdcc_spin_fock(sh2,sh7)
     if(spin_g1inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_l1_2_1(sh2,sh7,itm_hh)
     call ccsd_l1_2_2(sh2,sh7,itm_hh)
     call ccsd_l1_2_3(sh2,sh7,itm_hh)
     call ccsd_l1_2_4(sh2,sh7,itm_hh)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * g1inp(h7,p1,spin_g1inp) * itm_hh(h2,h7)
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
end subroutine ccsd_l1_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_2_1(sh2,sh7,i1)

!     i1 ( h2 h7 )_f + = 1 * f ( h2 h7 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh2,sh7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h2,h7
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh2,sh7)
  do h2 = 1,norb1
  do h7 = 1,norb1
     i1(h2,h7) = i1(h2,h7) + fact * fock(h2,h7,spin_fock)
  end do
  end do
end subroutine ccsd_l1_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_2_2(sh2,sh7,i1)

!     i1 ( h2 h7 )_ft + = 1 * Sum ( p3 ) * t ( p3 h7 )_t * i2 ( h2 p3 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h2,h7
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh7)
     spin_itm_hp = tdcc_spin_fock(sh2,sp3)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_2_2_1(sh2,sp3,itm_hp)
     call ccsd_l1_2_2_2(sh2,sp3,itm_hp)

     do h2 = 1,norb1
     do h7 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,h7) = i1(h2,h7) + fact * t1inp(p3,h7,spin_t1inp) * itm_hp(h2,p3)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsd_l1_2_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_2_2_1(sh2,sp3,i2)

!         i2 ( h2 p3 )_f + = 1 * f ( h2 p3 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh2,sp3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p3
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh2,sp3)
  do h2 = 1,norb1
  do p3 = norb1+1,nact
     i2(h2,p3) = i2(h2,p3) + fact * fock(h2,p3,spin_fock)
  end do
  end do
end subroutine ccsd_l1_2_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_2_2_2(sh2,sp3,i2)

!         i2 ( h2 p3 )_vt + = 1 * Sum ( h6 p5 ) * t ( p5 h6 )_t * v ( h2 h6 p3 p5 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p3
  integer(c_int) :: h6,p5,sh6,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh6 = 1,2
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_int2x = tdcc_spin_int2x(sh2,sh6,sp3,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do p3 = norb1+1,nact
     do h6 = 1,norb1
     do p5 = norb1+1,nact
        i2(h2,p3) = i2(h2,p3) + fact * t1inp(p5,h6,spin_t1inp) * int2x(h2,h6,p3,p5,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_2_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_2_3(sh2,sh7,i1)

!     i1 ( h2 h7 )_vt + = 1 * Sum ( h4 p3 ) * t ( p3 h4 )_t * v ( h2 h4 h7 p3 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sh7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h2,h7
  integer(c_int) :: h4,p3,sh4,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh4 = 1,2
  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh4)
     spin_int2x = tdcc_spin_int2x(sh2,sh4,sh7,sp3)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do h7 = 1,norb1
     do h4 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,h7) = i1(h2,h7) + fact * t1inp(p3,h4,spin_t1inp) * int2x(h2,h4,h7,p3,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_2_4(sh2,sh7,i1)

!     i1 ( h2 h7 )_vt + = -1/2 * Sum ( h6 p3 p4 ) * t ( p3 p4 h6 h7 )_t * v ( h2 h6 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sh7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h2,h7
  integer(c_int) :: h6,p3,p4,sh6,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh6 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh6,sh7)
     spin_int2x = tdcc_spin_int2x(sh2,sh6,sp3,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do h7 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,h7) = i1(h2,h7) + fact * t2inp(p3,p4,h6,h7,spin_t2inp) * int2x(h2,h6,p3,p4,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_l1_2_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_3(sh2,sp1,i0)

! i0 ( h2 p1 )_yf + = 1 * Sum ( p7 ) * y ( h2 p7 )_y * i1 ( p7 p1 )_f 3

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_g1inp
  integer(c_int) :: spin_itm_pp
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh2,sp7)
     spin_itm_pp = tdcc_spin_fock(sp7,sp1)
     if(spin_g1inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_l1_3_1(sp7,sp1,itm_pp)
     call ccsd_l1_3_2(sp7,sp1,itm_pp)
     call ccsd_l1_3_3(sp7,sp1,itm_pp)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * g1inp(h2,p7,spin_g1inp) * itm_pp(p7,p1)
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
end subroutine ccsd_l1_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_3_1(sp7,sp1,i1)

!     i1 ( p7 p1 )_f + = 1 * f ( p7 p1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sp7,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p7,p1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp7,sp1)
  do p7 = norb1+1,nact
  do p1 = norb1+1,nact
     i1(p7,p1) = i1(p7,p1) + fact * fock(p7,p1,spin_fock)
  end do
  end do
end subroutine ccsd_l1_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_3_2(sp7,sp1,i1)

!     i1 ( p7 p1 )_vt + = -1 * Sum ( h4 p3 ) * t ( p3 h4 )_t * v ( h4 p7 p1 p3 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp7,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p7,p1
  integer(c_int) :: h4,p3,sh4,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh4 = 1,2
  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh4)
     spin_int2x = tdcc_spin_int2x(sh4,sp7,sp1,sp3)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p7 = norb1+1,nact
     do p1 = norb1+1,nact
     do h4 = 1,norb1
     do p3 = norb1+1,nact
        i1(p7,p1) = i1(p7,p1) + fact * t1inp(p3,h4,spin_t1inp) * int2x(h4,p7,p1,p3,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_3_3(sp7,sp1,i1)

!     i1 ( p7 p1 )_vtt + = -1 * Sum ( h4 ) * t ( p7 h4 )_t * i2 ( h4 p1 )_vt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp7,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p7,p1
  integer(c_int) :: h4,sh4
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_dummy1
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sh4 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh4)
     spin_itm_hp = tdcc_spin_dummy1(sh4,sp1)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_3_3_1(sh4,sp1,itm_hp)

     do p7 = norb1+1,nact
     do p1 = norb1+1,nact
     do h4 = 1,norb1
        i1(p7,p1) = i1(p7,p1) + fact * t1inp(p7,h4,spin_t1inp) * itm_hp(h4,p1)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsd_l1_3_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_3_3_1(sh4,sp1,i2)

!         i2 ( h4 p1 )_vt + = 1 * Sum ( h6 p5 ) * t ( p5 h6 )_t * v ( h4 h6 p1 p5 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh4,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h4,p1
  integer(c_int) :: h6,p5,sh6,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh6 = 1,2
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_int2x = tdcc_spin_int2x(sh4,sh6,sp1,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do p5 = norb1+1,nact
        i2(h4,p1) = i2(h4,p1) + fact * t1inp(p5,h6,spin_t1inp) * int2x(h4,h6,p1,p5,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_3_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_4(sh2,sp1,i0)

! i0 ( h2 p1 )_yv + = -1 * Sum ( h4 p3 ) * y ( h4 p3 )_y * v ( h2 p3 h4 p1 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h4,p3,sh4,sp3
  integer(c_int) :: spin_g1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh4 = 1,2
  do sp3 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh4,sp3)
     spin_int2x = tdcc_spin_int2x(sh2,sp3,sh4,sp1)
     if(spin_g1inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h4 = 1,norb1
     do p3 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * g1inp(h4,p3,spin_g1inp) * int2x(h2,p3,h4,p1,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5(sh2,sp1,i0)

! i0 ( h2 p1 )_yv + = -1/2 * Sum ( h11 h12 p9 ) * y ( h11 h12 p1 p9 )_y * i1 ( h2 p9 h11 h12 )_v 6

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h11,h12,p9,sh11,sh12,sp9
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_hphh
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh11 = 1,2
  do sh12 = 1,2
  do sp9 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh11,sh12,sp1,sp9)
     spin_itm_hphh = tdcc_spin_int2x(sh2,sp9,sh11,sh12)
     if(spin_g2inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsd_l1_5_1(sh2,sp9,sh11,sh12,itm_hphh)
     call ccsd_l1_5_2(sh2,sp9,sh11,sh12,itm_hphh)
     call ccsd_l1_5_3(sh2,sp9,sh11,sh12,itm_hphh)
     call ccsd_l1_5_4(sh2,sp9,sh11,sh12,itm_hphh)
     call ccsd_l1_5_5(sh2,sp9,sh11,sh12,itm_hphh)
     call ccsd_l1_5_6(sh2,sp9,sh11,sh12,itm_hphh)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p9 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * g2inp(h11,h12,p1,p9,spin_g2inp) * itm_hphh(h2,p9,h11,h12)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  deallocate(itm_hphh)
end subroutine ccsd_l1_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_1(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_v + = 1 * v ( h2 p9 h11 h12 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh2,sp9,sh11,sh12)
  do h2 = 1,norb1
  do p9 = norb1+1,nact
  do h11 = 1,norb1
  do h12 = 1,norb1
     i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * int2x(h2,p9,h11,h12,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsd_l1_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_2(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_vt + = -1 * Sum ( h7 ) * t ( p9 h7 )_t * i2 ( h2 h7 h11 h12 )_v 3

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: h7,sh7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hhhh
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh7)
     spin_itm_hhhh = tdcc_spin_int2x(sh2,sh7,sh11,sh12)
     if(spin_t1inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsd_l1_5_2_1(sh2,sh7,sh11,sh12,itm_hhhh)
     call ccsd_l1_5_2_2(sh2,sh7,sh11,sh12,itm_hhhh)
     call ccsd_l1_5_2_3(sh2,sh7,sh11,sh12,itm_hhhh)

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do h7 = 1,norb1
        i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * t1inp(p9,h7,spin_t1inp) * itm_hhhh(h2,h7,h11,h12)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsd_l1_5_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_2_1(sh2,sh7,sh11,sh12,i2)

!         i2 ( h2 h7 h11 h12 )_v + = 1 * v ( h2 h7 h11 h12 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sh7,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,h7,h11,h12
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh2,sh7,sh11,sh12)
  do h2 = 1,norb1
  do h7 = 1,norb1
  do h11 = 1,norb1
  do h12 = 1,norb1
     i2(h2,h7,h11,h12) = i2(h2,h7,h11,h12) + fact * int2x(h2,h7,h11,h12,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsd_l1_5_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_2_2(sh2,sh7,sh11,sh12,i2)

!         i2 ( h2 h7 h11 h12 )_vt + = -2 * Sum ( p3 ) * t ( p3 h11 )_t * i3 ( h2 h7 h12 p3 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh7,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,h7,h11,h12
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh11)
     spin_itm_hhhp = tdcc_spin_int2x(sh2,sh7,sh12,sp3)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_5_2_2_1(sh2,sh7,sh12,sp3,itm_hhhp)
     call ccsd_l1_5_2_2_2(sh2,sh7,sh12,sp3,itm_hhhp)

     do h2 = 1,norb1
     do h7 = 1,norb1
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
        i2(h2,h7,h11,h12) = i2(h2,h7,h11,h12) + fact * t1inp(p3,h11,spin_t1inp) * itm_hhhp(h2,h7,h12,p3)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l1_5_2_2
!##########################################################
!##########################################################
subroutine ccsd_l1_5_2_2_1(sh2,sh7,sh12,sp3,i3)

!             i3 ( h2 h7 h12 p3 )_v + = 1 * v ( h2 h7 h12 p3 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sh7,sh12,sp3
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,h7,h12,p3
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh2,sh7,sh12,sp3)
  do h2 = 1,norb1
  do h7 = 1,norb1
  do h12 = 1,norb1
  do p3 = norb1+1,nact
     i3(h2,h7,h12,p3) = i3(h2,h7,h12,p3) + fact * int2x(h2,h7,h12,p3,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsd_l1_5_2_2_1
!##########################################################
!##########################################################
subroutine ccsd_l1_5_2_2_2(sh2,sh7,sh12,sp3,i3)

!             i3 ( h2 h7 h12 p3 )_vt + = -1/2 * Sum ( p5 ) * t ( p5 h12 )_t * v ( h2 h7 p3 p5 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sh7,sh12,sp3
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,h7,h12,p3
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh12)
     spin_int2x = tdcc_spin_int2x(sh2,sh7,sp3,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do h7 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
     do p5 = norb1+1,nact
        i3(h2,h7,h12,p3) = i3(h2,h7,h12,p3) + fact * t1inp(p5,h12,spin_t1inp) * int2x(h2,h7,p3,p5,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_5_2_2_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_2_3(sh2,sh7,sh11,sh12,i2)

!         i2 ( h2 h7 h11 h12 )_vt + = 1/2 * Sum ( p3 p4 ) * t ( p3 p4 h11 h12 )_t * v ( h2 h7 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sh7,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,h7,h11,h12
  integer(c_int) :: p3,p4,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh11,sh12)
     spin_int2x = tdcc_spin_int2x(sh2,sh7,sp3,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do h7 = 1,norb1
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i2(h2,h7,h11,h12) = i2(h2,h7,h11,h12) + fact * t2inp(p3,p4,h11,h12,spin_t2inp) * int2x(h2,h7,p3,p4,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_5_2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_3(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_vt + = -2 * Sum ( p3 ) * t ( p3 h11 )_t * i2 ( h2 p9 h12 p3 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hphp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh11)
     spin_itm_hphp = tdcc_spin_int2x(sh2,sp9,sh12,sp3)
     if(spin_t1inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_5_3_1(sh2,sp9,sh12,sp3,itm_hphp)
     call ccsd_l1_5_3_2(sh2,sp9,sh12,sp3,itm_hphp)

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * t1inp(p3,h11,spin_t1inp) * itm_hphp(h2,p9,h12,p3)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphp)
end subroutine ccsd_l1_5_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_3_1(sh2,sp9,sh12,sp3,i2)

!         i2 ( h2 p9 h12 p3 )_v + = 1 * v ( h2 p9 h12 p3 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh12,sp3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p9,h12,p3
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh2,sp9,sh12,sp3)
  do h2 = 1,norb1
  do p9 = norb1+1,nact
  do h12 = 1,norb1
  do p3 = norb1+1,nact
     i2(h2,p9,h12,p3) = i2(h2,p9,h12,p3) + fact * int2x(h2,p9,h12,p3,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsd_l1_5_3_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_3_2(sh2,sp9,sh12,sp3,i2)

!         i2 ( h2 p9 h12 p3 )_vt + = -1/2 * Sum ( p5 ) * t ( p5 h12 )_t * v ( h2 p9 p3 p5 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh12,sp3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p9,h12,p3
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh12)
     spin_int2x = tdcc_spin_int2x(sh2,sp9,sp3,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h12 = 1,norb1
     do p3 = norb1+1,nact
     do p5 = norb1+1,nact
        i2(h2,p9,h12,p3) = i2(h2,p9,h12,p3) + fact * t1inp(p5,h12,spin_t1inp) * int2x(h2,p9,p3,p5,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_5_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_4(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_ft + = 1 * Sum ( p5 ) * t ( p5 p9 h11 h12 )_t * i2 ( h2 p5 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp9,sh11,sh12)
     spin_itm_hp = tdcc_spin_fock(sh2,sp5)
     if(spin_t2inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_5_4_1(sh2,sp5,itm_hp)
     call ccsd_l1_5_4_2(sh2,sp5,itm_hp)

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p5 = norb1+1,nact
        i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * t2inp(p5,p9,h11,h12,spin_t2inp) * itm_hp(h2,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsd_l1_5_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_4_1(sh2,sp5,i2)

!         i2 ( h2 p5 )_f + = 1 * f ( h2 p5 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh2,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p5
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh2,sp5)
  do h2 = 1,norb1
  do p5 = norb1+1,nact
     i2(h2,p5) = i2(h2,p5) + fact * fock(h2,p5,spin_fock)
  end do
  end do
end subroutine ccsd_l1_5_4_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_4_2(sh2,sp5,i2)

!         i2 ( h2 p5 )_vt + = 1 * Sum ( h8 p7 ) * t ( p7 h8 )_t * v ( h2 h8 p5 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p5
  integer(c_int) :: h8,p7,sh8,sp7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_int2x = tdcc_spin_int2x(sh2,sh8,sp5,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do p5 = norb1+1,nact
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i2(h2,p5) = i2(h2,p5) + fact * t1inp(p7,h8,spin_t1inp) * int2x(h2,h8,p5,p7,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_5_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_5(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_vt + = -2 * Sum ( h6 p4 ) * t ( p4 p9 h6 h11 )_t * i2 ( h2 h6 h12 p4 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: h6,p4,sh6,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh6 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp9,sh6,sh11)
     spin_itm_hhhp = tdcc_spin_int2x(sh2,sh6,sh12,sp4)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_5_5_1(sh2,sh6,sh12,sp4,itm_hhhp)
     call ccsd_l1_5_5_2(sh2,sh6,sh12,sp4,itm_hhhp)

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do h6 = 1,norb1
     do p4 = norb1+1,nact
        i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * t2inp(p4,p9,h6,h11,spin_t2inp) * itm_hhhp(h2,h6,h12,p4)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l1_5_5
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_5_1(sh2,sh6,sh12,sp4,i2)

!         i2 ( h2 h6 h12 p4 )_v + = 1 * v ( h2 h6 h12 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sh6,sh12,sp4
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,h6,h12,p4
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh2,sh6,sh12,sp4)
  do h2 = 1,norb1
  do h6 = 1,norb1
  do h12 = 1,norb1
  do p4 = norb1+1,nact
     i2(h2,h6,h12,p4) = i2(h2,h6,h12,p4) + fact * int2x(h2,h6,h12,p4,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsd_l1_5_5_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_5_2(sh2,sh6,sh12,sp4,i2)

!         i2 ( h2 h6 h12 p4 )_vt + = -1 * Sum ( p7 ) * t ( p7 h12 )_t * v ( h2 h6 p4 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sh6,sh12,sp4
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,h6,h12,p4
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh12)
     spin_int2x = tdcc_spin_int2x(sh2,sh6,sp4,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do h6 = 1,norb1
     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h2,h6,h12,p4) = i2(h2,h6,h12,p4) + fact * t1inp(p7,h12,spin_t1inp) * int2x(h2,h6,p4,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_5_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_5_6(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_vt + = 1/2 * Sum ( p3 p4 ) * t ( p3 p4 h11 h12 )_t * v ( h2 p9 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: p3,p4,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh11,sh12)
     spin_int2x = tdcc_spin_int2x(sh2,sp9,sp3,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * t2inp(p3,p4,h11,h12,spin_t2inp) * int2x(h2,p9,p3,p4,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_5_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_6(sh2,sp1,i0)

! i0 ( h2 p1 )_yv + = 1/2 * Sum ( h7 p8 p5 ) * y ( h2 h7 p5 p8 )_y * i1 ( p5 p8 h7 p1 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h7,p8,p5,sh7,sp8,sp5
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_pphp
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh7 = 1,2
  do sp8 = 1,2
  do sp5 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh2,sh7,sp5,sp8)
     spin_itm_pphp = tdcc_spin_int2x(sp5,sp8,sh7,sp1)
     if(spin_g2inp * spin_itm_pphp == 0) cycle

     itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_6_1(sp5,sp8,sh7,sp1,itm_pphp)
     call ccsd_l1_6_2(sp5,sp8,sh7,sp1,itm_pphp)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do p8 = norb1+1,nact
     do p5 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * g2inp(h2,h7,p5,p8,spin_g2inp) * itm_pphp(p5,p8,h7,p1)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  deallocate(itm_pphp)
end subroutine ccsd_l1_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_6_1(sp5,sp8,sh7,sp1,i1)

!     i1 ( p5 p8 h7 p1 )_v + = -1 * v ( p5 p8 h7 p1 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sp5,sp8,sh7,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: p5,p8,h7,p1
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sp5,sp8,sh7,sp1)
  do p5 = norb1+1,nact
  do p8 = norb1+1,nact
  do h7 = 1,norb1
  do p1 = norb1+1,nact
     i1(p5,p8,h7,p1) = i1(p5,p8,h7,p1) + fact * int2x(p5,p8,h7,p1,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsd_l1_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_6_2(sp5,sp8,sh7,sp1,i1)

!     i1 ( p5 p8 h7 p1 )_vt + = 1 * Sum ( p3 ) * t ( p3 h7 )_t * v ( p5 p8 p1 p3 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp5,sp8,sh7,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: p5,p8,h7,p1
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh7)
     spin_int2x = tdcc_spin_int2x(sp5,sp8,sp1,sp3)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p5 = norb1+1,nact
     do p8 = norb1+1,nact
     do h7 = 1,norb1
     do p1 = norb1+1,nact
     do p3 = norb1+1,nact
        i1(p5,p8,h7,p1) = i1(p5,p8,h7,p1) + fact * t1inp(p3,h7,spin_t1inp) * int2x(p5,p8,p1,p3,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_6_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_7(sh2,sp1,i0)

! i0 ( h2 p1 )_vt + = 1 * Sum ( h10 p9 ) * i1 ( p9 h10 )_t * v ( h2 h10 p1 p9 )_v 4

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h10,p9,sh10,sp9
  integer(c_int) :: spin_itm_ph
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_ph(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_ph((norb1+1):nact,1:norb1))
  do sh10 = 1,2
  do sp9 = 1,2
     spin_itm_ph = tdcc_spin_t1inp(sp9,sh10)
     spin_int2x = tdcc_spin_int2x(sh2,sh10,sp1,sp9)
     if(spin_itm_ph * spin_int2x == 0) cycle

     itm_ph((norb1+1):nact,1:norb1) = czero
     call ccsd_l1_7_1(sp9,sh10,itm_ph)
     call ccsd_l1_7_2(sp9,sh10,itm_ph)
     call ccsd_l1_7_3(sp9,sh10,itm_ph)
     call ccsd_l1_7_4(sp9,sh10,itm_ph)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h10 = 1,norb1
     do p9 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * itm_ph(p9,h10) * int2x(h2,h10,p1,p9,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_ph)
end subroutine ccsd_l1_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_7_1(sp9,sh10,i1)

!     i1 ( p9 h10 )_t + = 1 * t ( p9 h10 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp9,sh10
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer(c_int) :: p9,h10
  integer(c_int) :: sdum
  integer(c_int) :: spin_t1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
  do p9 = norb1+1,nact
  do h10 = 1,norb1
     i1(p9,h10) = i1(p9,h10) + fact * t1inp(p9,h10,spin_t1inp)
  end do
  end do
end subroutine ccsd_l1_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_7_2(sp9,sh10,i1)

!     i1 ( p9 h10 )_yt + = 1 * Sum ( h5 p3 ) * t ( p3 p9 h5 h10 )_t * y ( h5 p3 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g1inp

  implicit none
  integer(c_int),intent(in) :: sp9,sh10
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer(c_int) :: p9,h10
  integer(c_int) :: h5,p3,sh5,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g1inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp9,sh5,sh10)
     spin_g1inp = tdcc_spin_g1inp(sh5,sp3)
     if(spin_t2inp * spin_g1inp == 0) cycle

     do p9 = norb1+1,nact
     do h10 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
        i1(p9,h10) = i1(p9,h10) + fact * t2inp(p3,p9,h5,h10,spin_t2inp) * g1inp(h5,p3,spin_g1inp)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_7_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_7_3(sp9,sh10,i1)

!     i1 ( p9 h10 )_ytt + = -1 * Sum ( h6 ) * t ( p9 h6 )_t * i2 ( h6 h10 )_yt 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp9,sh10
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer(c_int) :: p9,h10
  integer(c_int) :: h6,sh6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hh
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_dummy1
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh6)
     spin_itm_hh = tdcc_spin_dummy1(sh6,sh10)
     if(spin_t1inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_l1_7_3_1(sh6,sh10,itm_hh)
     call ccsd_l1_7_3_2(sh6,sh10,itm_hh)

     do p9 = norb1+1,nact
     do h10 = 1,norb1
     do h6 = 1,norb1
        i1(p9,h10) = i1(p9,h10) + fact * t1inp(p9,h6,spin_t1inp) * itm_hh(h6,h10)
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
end subroutine ccsd_l1_7_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_7_3_1(sh6,sh10,i2)

!         i2 ( h6 h10 )_yt + = 1 * Sum ( p5 ) * t ( p5 h10 )_t * y ( h6 p5 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g1inp

  implicit none
  integer(c_int),intent(in) :: sh6,sh10
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1)
  integer(c_int) :: h6,h10
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh10)
     spin_g1inp = tdcc_spin_g1inp(sh6,sp5)
     if(spin_t1inp * spin_g1inp == 0) cycle

     do h6 = 1,norb1
     do h10 = 1,norb1
     do p5 = norb1+1,nact
        i2(h6,h10) = i2(h6,h10) + fact * t1inp(p5,h10,spin_t1inp) * g1inp(h6,p5,spin_g1inp)
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_7_3_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_7_3_2(sh6,sh10,i2)

!         i2 ( h6 h10 )_yt + = 1/2 * Sum ( h5 p3 p4 ) * t ( p3 p4 h5 h10 )_t * y ( h5 h6 p3 p4 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh6,sh10
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1)
  integer(c_int) :: h6,h10
  integer(c_int) :: h5,p3,p4,sh5,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh10)
     spin_g2inp = tdcc_spin_g2inp(sh5,sh6,sp3,sp4)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h6 = 1,norb1
     do h10 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i2(h6,h10) = i2(h6,h10) + fact * t2inp(p3,p4,h5,h10,spin_t2inp) * g2inp(h5,h6,p3,p4,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_l1_7_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_7_4(sp9,sh10,i1)

!     i1 ( p9 h10 )_ytt + = 1/2 * Sum ( h5 h6 p3 ) * t ( p3 p9 h5 h6 )_t * i2 ( h5 h6 h10 p3 )_yt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp9,sh10
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer(c_int) :: p9,h10
  integer(c_int) :: h5,h6,p3,sh5,sh6,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh5 = 1,2
  do sh6 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp9,sh5,sh6)
     spin_itm_hhhp = tdcc_spin_dummy2(sh5,sh6,sh10,sp3)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_7_4_1(sh5,sh6,sh10,sp3,itm_hhhp)

     do p9 = norb1+1,nact
     do h10 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
        i1(p9,h10) = i1(p9,h10) + fact * t2inp(p3,p9,h5,h6,spin_t2inp) * itm_hhhp(h5,h6,h10,p3)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l1_7_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_7_4_1(sh5,sh6,sh10,sp3,i2)

!         i2 ( h5 h6 h10 p3 )_yt + = -1 * Sum ( p7 ) * t ( p7 h10 )_t * y ( h5 h6 p3 p7 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh5,sh6,sh10,sp3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h5,h6,h10,p3
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh10)
     spin_g2inp = tdcc_spin_g2inp(sh5,sh6,sp3,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h5 = 1,norb1
     do h6 = 1,norb1
     do h10 = 1,norb1
     do p3 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h5,h6,h10,p3) = i2(h5,h6,h10,p3) + fact * t1inp(p7,h10,spin_t1inp) * g2inp(h5,h6,p3,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_7_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_8(sh2,sp1,i0)

! i0 ( h2 p1 )_ytf + = -1 * Sum ( h3 ) * i1 ( h2 h3 )_yt * f ( h3 p1 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h3,sh3
  integer(c_int) :: spin_itm_hh
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh3 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh2,sh3)
     spin_fock = tdcc_spin_fock(sh3,sp1)
     if(spin_itm_hh * spin_fock == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_l1_8_1(sh2,sh3,itm_hh)
     call ccsd_l1_8_2(sh2,sh3,itm_hh)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h3 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * itm_hh(h2,h3) * fock(h3,p1,spin_fock)
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
end subroutine ccsd_l1_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_8_1(sh2,sh3,i1)

!     i1 ( h2 h3 )_yt + = 1 * Sum ( p4 ) * t ( p4 h3 )_t * y ( h2 p4 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h2,h3
  integer(c_int) :: p4,sp4
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp4 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp4,sh3)
     spin_g1inp = tdcc_spin_g1inp(sh2,sp4)
     if(spin_t1inp * spin_g1inp == 0) cycle

     do h2 = 1,norb1
     do h3 = 1,norb1
     do p4 = norb1+1,nact
        i1(h2,h3) = i1(h2,h3) + fact * t1inp(p4,h3,spin_t1inp) * g1inp(h2,p4,spin_g1inp)
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_8_2(sh2,sh3,i1)

!     i1 ( h2 h3 )_yt + = 1/2 * Sum ( h6 p4 p5 ) * t ( p4 p5 h3 h6 )_t * y ( h2 h6 p4 p5 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h2,h3
  integer(c_int) :: h6,p4,p5,sh6,sp4,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh6 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh3,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh6,sp4,sp5)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do h3 = 1,norb1
     do h6 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h2,h3) = i1(h2,h3) + fact * t2inp(p4,p5,h3,h6,spin_t2inp) * g2inp(h2,h6,p4,p5,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_l1_8_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_9(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( h6 h8 ) * i1 ( h6 h8 )_yt * v ( h2 h8 h6 p1 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h6,h8,sh6,sh8
  integer(c_int) :: spin_itm_hh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh6 = 1,2
  do sh8 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh6,sh8)
     spin_int2x = tdcc_spin_int2x(sh2,sh8,sh6,sp1)
     if(spin_itm_hh * spin_int2x == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_l1_9_1(sh6,sh8,itm_hh)
     call ccsd_l1_9_2(sh6,sh8,itm_hh)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do h8 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * itm_hh(h6,h8) * int2x(h2,h8,h6,p1,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hh)
end subroutine ccsd_l1_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_9_1(sh6,sh8,i1)

!     i1 ( h6 h8 )_yt + = 1 * Sum ( p3 ) * t ( p3 h8 )_t * y ( h6 p3 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g1inp

  implicit none
  integer(c_int),intent(in) :: sh6,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h6,h8
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh8)
     spin_g1inp = tdcc_spin_g1inp(sh6,sp3)
     if(spin_t1inp * spin_g1inp == 0) cycle

     do h6 = 1,norb1
     do h8 = 1,norb1
     do p3 = norb1+1,nact
        i1(h6,h8) = i1(h6,h8) + fact * t1inp(p3,h8,spin_t1inp) * g1inp(h6,p3,spin_g1inp)
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_9_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_9_2(sh6,sh8,i1)

!     i1 ( h6 h8 )_yt + = 1/2 * Sum ( h5 p3 p4 ) * t ( p3 p4 h5 h8 )_t * y ( h5 h6 p3 p4 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh6,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h6,h8
  integer(c_int) :: h5,p3,p4,sh5,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh5,sh6,sp3,sp4)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h6 = 1,norb1
     do h8 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h6,h8) = i1(h6,h8) + fact * t2inp(p3,p4,h5,h8,spin_t2inp) * g2inp(h5,h6,p3,p4,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_l1_9_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_10(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( p7 p8 ) * i1 ( p7 p8 )_yt * v ( h2 p8 p1 p7 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_itm_pp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
  do sp8 = 1,2
     spin_itm_pp = tdcc_spin_dummy1(sp7,sp8)
     spin_int2x = tdcc_spin_int2x(sh2,sp8,sp1,sp7)
     if(spin_itm_pp * spin_int2x == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_l1_10_1(sp7,sp8,itm_pp)
     call ccsd_l1_10_2(sp7,sp8,itm_pp)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * itm_pp(p7,p8) * int2x(h2,p8,p1,p7,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_pp)
end subroutine ccsd_l1_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_10_1(sp7,sp8,i1)

!     i1 ( p7 p8 )_yt + = 1 * Sum ( h4 ) * t ( p7 h4 )_t * y ( h4 p8 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g1inp

  implicit none
  integer(c_int),intent(in) :: sp7,sp8
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p7,p8
  integer(c_int) :: h4,sh4
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh4 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh4)
     spin_g1inp = tdcc_spin_g1inp(sh4,sp8)
     if(spin_t1inp * spin_g1inp == 0) cycle

     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
     do h4 = 1,norb1
        i1(p7,p8) = i1(p7,p8) + fact * t1inp(p7,h4,spin_t1inp) * g1inp(h4,p8,spin_g1inp)
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_10_2(sp7,sp8,i1)

!     i1 ( p7 p8 )_yt + = 1/2 * Sum ( h5 h6 p3 ) * t ( p3 p7 h5 h6 )_t * y ( h5 h6 p3 p8 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sp7,sp8
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p7,p8
  integer(c_int) :: h5,h6,p3,sh5,sh6,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp7,sh5,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh5,sh6,sp3,sp8)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
        i1(p7,p8) = i1(p7,p8) + fact * t2inp(p3,p7,h5,h6,spin_t2inp) * g2inp(h5,h6,p3,p8,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_l1_10_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_11(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( h6 h4 p5 ) * i1 ( h2 h6 h4 p5 )_yt * v ( h4 p5 h6 p1 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h6,h4,p5,sh6,sh4,sp5
  integer(c_int) :: spin_itm_hhhp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh6 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh2,sh6,sh4,sp5)
     spin_int2x = tdcc_spin_int2x(sh4,sp5,sh6,sp1)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_11_1(sh2,sh6,sh4,sp5,itm_hhhp)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * itm_hhhp(h2,h6,h4,p5) * int2x(h4,p5,h6,p1,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l1_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_11_1(sh2,sh6,sh4,sp5,i1)

!     i1 ( h2 h6 h4 p5 )_yt + = 1 * Sum ( p3 ) * t ( p3 h4 )_t * y ( h2 h6 p3 p5 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh6,sh4,sp5
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,h6,h4,p5
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh4)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh6,sp3,sp5)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do h6 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p3 = norb1+1,nact
        i1(h2,h6,h4,p5) = i1(h2,h6,h4,p5) + fact * t1inp(p3,h4,spin_t1inp) * g2inp(h2,h6,p3,p5,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_12(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1/2 * Sum ( p9 h12 h6 ) * i1 ( h2 p9 h6 h12 )_yt * v ( h6 h12 p1 p9 )_v 4

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: p9,h12,h6,sp9,sh12,sh6
  integer(c_int) :: spin_itm_hphh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sp9 = 1,2
  do sh12 = 1,2
  do sh6 = 1,2
     spin_itm_hphh = tdcc_spin_dummy2(sh2,sp9,sh6,sh12)
     spin_int2x = tdcc_spin_int2x(sh6,sh12,sp1,sp9)
     if(spin_itm_hphh * spin_int2x == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsd_l1_12_1(sh2,sp9,sh6,sh12,itm_hphh)
     call ccsd_l1_12_2(sh2,sp9,sh6,sh12,itm_hphh)
     call ccsd_l1_12_3(sh2,sp9,sh6,sh12,itm_hphh)
     call ccsd_l1_12_4(sh2,sp9,sh6,sh12,itm_hphh)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p9 = norb1+1,nact
     do h12 = 1,norb1
     do h6 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * itm_hphh(h2,p9,h6,h12) * int2x(h6,h12,p1,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  deallocate(itm_hphh)
end subroutine ccsd_l1_12
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_12_1(sh2,sp9,sh6,sh12,i1)

!     i1 ( h2 p9 h6 h12 )_yt + = -1 * Sum ( p3 ) * t ( p3 p9 h6 h12 )_t * y ( h2 p3 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh6,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h6,h12
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g1inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp9,sh6,sh12)
     spin_g1inp = tdcc_spin_g1inp(sh2,sp3)
     if(spin_t2inp * spin_g1inp == 0) cycle

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h6 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,p9,h6,h12) = i1(h2,p9,h6,h12) + fact * t2inp(p3,p9,h6,h12,spin_t2inp) * g1inp(h2,p3,spin_g1inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_12_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_12_2(sh2,sp9,sh6,sh12,i1)

!     i1 ( h2 p9 h6 h12 )_ytt + = 1/2 * Sum ( h10 ) * t ( p9 h10 )_t * i2 ( h2 h10 h6 h12 )_yt 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh6,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h6,h12
  integer(c_int) :: h10,sh10
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hhhh
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh10 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
     spin_itm_hhhh = tdcc_spin_dummy2(sh2,sh10,sh6,sh12)
     if(spin_t1inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsd_l1_12_2_1(sh2,sh10,sh6,sh12,itm_hhhh)
     call ccsd_l1_12_2_2(sh2,sh10,sh6,sh12,itm_hhhh)

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h6 = 1,norb1
     do h12 = 1,norb1
     do h10 = 1,norb1
        i1(h2,p9,h6,h12) = i1(h2,p9,h6,h12) + fact * t1inp(p9,h10,spin_t1inp) * itm_hhhh(h2,h10,h6,h12)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsd_l1_12_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_12_2_1(sh2,sh10,sh6,sh12,i2)

!         i2 ( h2 h10 h6 h12 )_yt + = 1 * Sum ( p3 p4 ) * t ( p3 p4 h6 h12 )_t * y ( h2 h10 p3 p4 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh10,sh6,sh12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,h10,h6,h12
  integer(c_int) :: p3,p4,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh6,sh12)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh10,sp3,sp4)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do h10 = 1,norb1
     do h6 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i2(h2,h10,h6,h12) = i2(h2,h10,h6,h12) + fact * t2inp(p3,p4,h6,h12,spin_t2inp) * g2inp(h2,h10,p3,p4,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_12_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_12_2_2(sh2,sh10,sh6,sh12,i2)

!         i2 ( h2 h10 h6 h12 )_ytt + = -2 * Sum ( p5 ) * t ( p5 h12 )_t * i3 ( h2 h10 h6 p5 )_yt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh10,sh6,sh12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,h10,h6,h12
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh12)
     spin_itm_hhhp = tdcc_spin_dummy2(sh2,sh10,sh6,sp5)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_12_2_2_1(sh2,sh10,sh6,sp5,itm_hhhp)

     do h2 = 1,norb1
     do h10 = 1,norb1
     do h6 = 1,norb1
     do h12 = 1,norb1
     do p5 = norb1+1,nact
        i2(h2,h10,h6,h12) = i2(h2,h10,h6,h12) + fact * t1inp(p5,h12,spin_t1inp) * itm_hhhp(h2,h10,h6,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l1_12_2_2
!##########################################################
!##########################################################
subroutine ccsd_l1_12_2_2_1(sh2,sh10,sh6,sp5,i3)

!             i3 ( h2 h10 h6 p5 )_yt + = 1 * Sum ( p7 ) * t ( p7 h6 )_t * y ( h2 h10 p5 p7 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh10,sh6,sp5
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,h10,h6,p5
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh10,sp5,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do h10 = 1,norb1
     do h6 = 1,norb1
     do p5 = norb1+1,nact
     do p7 = norb1+1,nact
        i3(h2,h10,h6,p5) = i3(h2,h10,h6,p5) + fact * t1inp(p7,h6,spin_t1inp) * g2inp(h2,h10,p5,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_12_2_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_12_3(sh2,sp9,sh6,sh12,i1)

!     i1 ( h2 p9 h6 h12 )_ytt + = 2 * Sum ( h5 p3 ) * t ( p3 p9 h5 h12 )_t * i2 ( h2 h5 h6 p3 )_yt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh6,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h6,h12
  integer(c_int) :: h5,p3,sh5,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh5 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp9,sh5,sh12)
     spin_itm_hhhp = tdcc_spin_dummy2(sh2,sh5,sh6,sp3)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_12_3_1(sh2,sh5,sh6,sp3,itm_hhhp)

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h6 = 1,norb1
     do h12 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,p9,h6,h12) = i1(h2,p9,h6,h12) + fact * t2inp(p3,p9,h5,h12,spin_t2inp) * itm_hhhp(h2,h5,h6,p3)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l1_12_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_12_3_1(sh2,sh5,sh6,sp3,i2)

!         i2 ( h2 h5 h6 p3 )_yt + = 1 * Sum ( p7 ) * t ( p7 h6 )_t * y ( h2 h5 p3 p7 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh5,sh6,sp3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,h5,h6,p3
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh5,sp3,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h2,h5,h6,p3) = i2(h2,h5,h6,p3) + fact * t1inp(p7,h6,spin_t1inp) * g2inp(h2,h5,p3,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_12_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_12_4(sh2,sp9,sh6,sh12,i1)

!     i1 ( h2 p9 h6 h12 )_ytt + = 1 * t ( p9 h6 )_t * i2 ( h2 h12 )_yt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh6,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h6,h12
  integer(c_int) :: sdum
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hh
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_dummy1
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sdum = 1,1
     spin_t1inp = tdcc_spin_t1inp(sp9,sh6)
     spin_itm_hh = tdcc_spin_dummy1(sh2,sh12)
     if(spin_t1inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_l1_12_4_1(sh2,sh12,itm_hh)

     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h6 = 1,norb1
     do h12 = 1,norb1
        i1(h2,p9,h6,h12) = i1(h2,p9,h6,h12) + fact * t1inp(p9,h6,spin_t1inp) * itm_hh(h2,h12)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
end subroutine ccsd_l1_12_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_12_4_1(sh2,sh12,i2)

!         i2 ( h2 h12 )_yt + = -1 * Sum ( h5 p3 p4 ) * t ( p3 p4 h5 h12 )_t * y ( h2 h5 p3 p4 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1)
  integer(c_int) :: h2,h12
  integer(c_int) :: h5,p3,p4,sh5,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh12)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh5,sp3,sp4)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do h12 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i2(h2,h12) = i2(h2,h12) + fact * t2inp(p3,p4,h5,h12,spin_t2inp) * g2inp(h2,h5,p3,p4,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_l1_12_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_13(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = -1/4 * Sum ( h7 h8 h6 ) * i1 ( h2 h7 h6 h8 )_yt * v ( h6 h8 h7 p1 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h7,h8,h6,sh7,sh8,sh6
  integer(c_int) :: spin_itm_hhhh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh7 = 1,2
  do sh8 = 1,2
  do sh6 = 1,2
     spin_itm_hhhh = tdcc_spin_dummy2(sh2,sh7,sh6,sh8)
     spin_int2x = tdcc_spin_int2x(sh6,sh8,sh7,sp1)
     if(spin_itm_hhhh * spin_int2x == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsd_l1_13_1(sh2,sh7,sh6,sh8,itm_hhhh)
     call ccsd_l1_13_2(sh2,sh7,sh6,sh8,itm_hhhh)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do h6 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * itm_hhhh(h2,h7,h6,h8) * int2x(h6,h8,h7,p1,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsd_l1_13
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_13_1(sh2,sh7,sh6,sh8,i1)

!     i1 ( h2 h7 h6 h8 )_yt + = 1 * Sum ( p3 p4 ) * t ( p3 p4 h6 h8 )_t * y ( h2 h7 p3 p4 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh7,sh6,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,h7,h6,h8
  integer(c_int) :: p3,p4,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh6,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh7,sp3,sp4)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do h7 = 1,norb1
     do h6 = 1,norb1
     do h8 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,h7,h6,h8) = i1(h2,h7,h6,h8) + fact * t2inp(p3,p4,h6,h8,spin_t2inp) * g2inp(h2,h7,p3,p4,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_13_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_13_2(sh2,sh7,sh6,sh8,i1)

!     i1 ( h2 h7 h6 h8 )_ytt + = -2 * Sum ( p3 ) * t ( p3 h8 )_t * i2 ( h2 h7 h6 p3 )_yt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh7,sh6,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,h7,h6,h8
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh8)
     spin_itm_hhhp = tdcc_spin_dummy2(sh2,sh7,sh6,sp3)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_13_2_1(sh2,sh7,sh6,sp3,itm_hhhp)

     do h2 = 1,norb1
     do h7 = 1,norb1
     do h6 = 1,norb1
     do h8 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,h7,h6,h8) = i1(h2,h7,h6,h8) + fact * t1inp(p3,h8,spin_t1inp) * itm_hhhp(h2,h7,h6,p3)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l1_13_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_13_2_1(sh2,sh7,sh6,sp3,i2)

!         i2 ( h2 h7 h6 p3 )_yt + = 1 * Sum ( p5 ) * t ( p5 h6 )_t * y ( h2 h7 p3 p5 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh7,sh6,sp3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,h7,h6,p3
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh7,sp3,sp5)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do h7 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do p5 = norb1+1,nact
        i2(h2,h7,h6,p3) = i2(h2,h7,h6,p3) + fact * t1inp(p5,h6,spin_t1inp) * g2inp(h2,h7,p3,p5,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_13_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_14(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( p8 h6 p7 ) * i1 ( h2 p8 h6 p7 )_yt * v ( h6 p7 p1 p8 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: p8,h6,p7,sp8,sh6,sp7
  integer(c_int) :: spin_itm_hphp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp8 = 1,2
  do sh6 = 1,2
  do sp7 = 1,2
     spin_itm_hphp = tdcc_spin_dummy2(sh2,sp8,sh6,sp7)
     spin_int2x = tdcc_spin_int2x(sh6,sp7,sp1,sp8)
     if(spin_itm_hphp * spin_int2x == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_14_1(sh2,sp8,sh6,sp7,itm_hphp)
     call ccsd_l1_14_2(sh2,sp8,sh6,sp7,itm_hphp)

     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p8 = norb1+1,nact
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * itm_hphp(h2,p8,h6,p7) * int2x(h6,p7,p1,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  deallocate(itm_hphp)
end subroutine ccsd_l1_14
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_14_1(sh2,sp8,sh6,sp7,i1)

!     i1 ( h2 p8 h6 p7 )_yt + = 1 * Sum ( h5 p3 ) * t ( p3 p8 h5 h6 )_t * y ( h2 h5 p3 p7 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp8,sh6,sp7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p8,h6,p7
  integer(c_int) :: h5,p3,sh5,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp8,sh5,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh5,sp3,sp7)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do p8 = norb1+1,nact
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do h5 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,p8,h6,p7) = i1(h2,p8,h6,p7) + fact * t2inp(p3,p8,h5,h6,spin_t2inp) * g2inp(h2,h5,p3,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l1_14_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_14_2(sh2,sp8,sh6,sp7,i1)

!     i1 ( h2 p8 h6 p7 )_ytt + = -1 * Sum ( h4 ) * t ( p8 h4 )_t * i2 ( h2 h4 h6 p7 )_yt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp8,sh6,sp7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p8,h6,p7
  integer(c_int) :: h4,sh4
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh4 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh4)
     spin_itm_hhhp = tdcc_spin_dummy2(sh2,sh4,sh6,sp7)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l1_14_2_1(sh2,sh4,sh6,sp7,itm_hhhp)

     do h2 = 1,norb1
     do p8 = norb1+1,nact
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do h4 = 1,norb1
        i1(h2,p8,h6,p7) = i1(h2,p8,h6,p7) + fact * t1inp(p8,h4,spin_t1inp) * itm_hhhp(h2,h4,h6,p7)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l1_14_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l1_14_2_1(sh2,sh4,sh6,sp7,i2)

!         i2 ( h2 h4 h6 p7 )_yt + = 1 * Sum ( p5 ) * t ( p5 h6 )_t * y ( h2 h4 p5 p7 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh4,sh6,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,h4,h6,p7
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh4,sp5,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h2 = 1,norb1
     do h4 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do p5 = norb1+1,nact
        i2(h2,h4,h6,p7) = i2(h2,h4,h6,p7) + fact * t1inp(p5,h6,spin_t1inp) * g2inp(h2,h4,p5,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l1_14_2_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccsd_l1_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1out

  implicit none

  call ccsd_l1_1(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_2(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_3(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_4(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_5(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_6(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_7(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_8(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_9(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_10(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_11(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_12(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_13(1,1,g1out(1,norb1+1,1))
  call ccsd_l1_14(1,1,g1out(1,norb1+1,1))

end subroutine ccsd_l1_main
!**********************************************************
