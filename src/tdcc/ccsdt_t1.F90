!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_1(sp2,sh1,i0)

! i0 ( p2 h1 )_f + = 1 * f ( p2 h1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp2,sh1)
  do p2 = norb1+1,nact
  do h1 = 1,norb1
     i0(p2,h1) = i0(p2,h1) + fact * fock(p2,h1,spin_fock)
  end do
  end do
end subroutine ccsdt_t1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_2(sp2,sh1,i0)

! i0 ( p2 h1 )_tf + = -1 * Sum ( h7 ) * t ( p2 h7 )_t * i1 ( h7 h1 )_f 4

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: h7,sh7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hh
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp2,sh7)
     spin_itm_hh = tdcc_spin_fock(sh7,sh1)
     if(spin_t1inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsdt_t1_2_1(sh7,sh1,itm_hh)
     call ccsdt_t1_2_2(sh7,sh1,itm_hh)
     call ccsdt_t1_2_3(sh7,sh1,itm_hh)
     call ccsdt_t1_2_4(sh7,sh1,itm_hh)

     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do h7 = 1,norb1
        i0(p2,h1) = i0(p2,h1) + fact * t1inp(p2,h7,spin_t1inp) * itm_hh(h7,h1)
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
end subroutine ccsdt_t1_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_2_1(sh7,sh1,i1)

!     i1 ( h7 h1 )_f + = 1 * f ( h7 h1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh7,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h7,h1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh7,sh1)
  do h7 = 1,norb1
  do h1 = 1,norb1
     i1(h7,h1) = i1(h7,h1) + fact * fock(h7,h1,spin_fock)
  end do
  end do
end subroutine ccsdt_t1_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_2_2(sh7,sh1,i1)

!     i1 ( h7 h1 )_ft + = 1 * Sum ( p8 ) * t ( p8 h1 )_t * i2 ( h7 p8 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh7,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h7,h1
  integer(c_int) :: p8,sp8
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_itm_hp = tdcc_spin_fock(sh7,sp8)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_t1_2_2_1(sh7,sp8,itm_hp)
     call ccsdt_t1_2_2_2(sh7,sp8,itm_hp)

     do h7 = 1,norb1
     do h1 = 1,norb1
     do p8 = norb1+1,nact
        i1(h7,h1) = i1(h7,h1) + fact * t1inp(p8,h1,spin_t1inp) * itm_hp(h7,p8)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_t1_2_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_2_2_1(sh7,sp8,i2)

!         i2 ( h7 p8 )_f + = 1 * f ( h7 p8 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh7,sp8
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h7,p8
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh7,sp8)
  do h7 = 1,norb1
  do p8 = norb1+1,nact
     i2(h7,p8) = i2(h7,p8) + fact * fock(h7,p8,spin_fock)
  end do
  end do
end subroutine ccsdt_t1_2_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_2_2_2(sh7,sp8,i2)

!         i2 ( h7 p8 )_vt + = 1 * Sum ( h6 p5 ) * t ( p5 h6 )_t * v ( h6 h7 p5 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh7,sp8
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h7,p8
  integer(c_int) :: h6,p5,sh6,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh6 = 1,2
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_int2x = tdcc_spin_int2x(sh6,sh7,sp5,sp8)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h7 = 1,norb1
     do p8 = norb1+1,nact
     do h6 = 1,norb1
     do p5 = norb1+1,nact
        i2(h7,p8) = i2(h7,p8) + fact * t1inp(p5,h6,spin_t1inp) * int2x(h6,h7,p5,p8,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t1_2_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_2_3(sh7,sh1,i1)

!     i1 ( h7 h1 )_vt + = -1 * Sum ( h5 p4 ) * t ( p4 h5 )_t * v ( h5 h7 h1 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh7,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h7,h1
  integer(c_int) :: h5,p4,sh5,sp4
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh5 = 1,2
  do sp4 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp4,sh5)
     spin_int2x = tdcc_spin_int2x(sh5,sh7,sh1,sp4)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h7 = 1,norb1
     do h1 = 1,norb1
     do h5 = 1,norb1
     do p4 = norb1+1,nact
        i1(h7,h1) = i1(h7,h1) + fact * t1inp(p4,h5,spin_t1inp) * int2x(h5,h7,h1,p4,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t1_2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_2_4(sh7,sh1,i1)

!     i1 ( h7 h1 )_vt + = -1/2 * Sum ( h5 p3 p4 ) * t ( p3 p4 h1 h5 )_t * v ( h5 h7 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh7,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h7,h1
  integer(c_int) :: h5,p3,p4,sh5,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh1,sh5)
     spin_int2x = tdcc_spin_int2x(sh5,sh7,sp3,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h7 = 1,norb1
     do h1 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h7,h1) = i1(h7,h1) + fact * t2inp(p3,p4,h1,h5,spin_t2inp) * int2x(h5,h7,p3,p4,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_t1_2_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_3(sp2,sh1,i0)

! i0 ( p2 h1 )_tf + = 1 * Sum ( p3 ) * t ( p3 h1 )_t * i1 ( p2 p3 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_pp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh1)
     spin_itm_pp = tdcc_spin_fock(sp2,sp3)
     if(spin_t1inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_t1_3_1(sp2,sp3,itm_pp)
     call ccsdt_t1_3_2(sp2,sp3,itm_pp)

     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do p3 = norb1+1,nact
        i0(p2,h1) = i0(p2,h1) + fact * t1inp(p3,h1,spin_t1inp) * itm_pp(p2,p3)
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
end subroutine ccsdt_t1_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_3_1(sp2,sp3,i1)

!     i1 ( p2 p3 )_f + = 1 * f ( p2 p3 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sp2,sp3
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p2,p3
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp2,sp3)
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i1(p2,p3) = i1(p2,p3) + fact * fock(p2,p3,spin_fock)
  end do
  end do
end subroutine ccsdt_t1_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_3_2(sp2,sp3,i1)

!     i1 ( p2 p3 )_vt + = -1 * Sum ( h5 p4 ) * t ( p4 h5 )_t * v ( h5 p2 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp2,sp3
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p2,p3
  integer(c_int) :: h5,p4,sh5,sp4
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh5 = 1,2
  do sp4 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp4,sh5)
     spin_int2x = tdcc_spin_int2x(sh5,sp2,sp3,sp4)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h5 = 1,norb1
     do p4 = norb1+1,nact
        i1(p2,p3) = i1(p2,p3) + fact * t1inp(p4,h5,spin_t1inp) * int2x(h5,p2,p3,p4,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t1_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_4(sp2,sh1,i0)

! i0 ( p2 h1 )_vt + = -1 * Sum ( h4 p3 ) * t ( p3 h4 )_t * v ( h4 p2 h1 p3 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: h4,p3,sh4,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh4 = 1,2
  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh4)
     spin_int2x = tdcc_spin_int2x(sh4,sp2,sh1,sp3)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do h4 = 1,norb1
     do p3 = norb1+1,nact
        i0(p2,h1) = i0(p2,h1) + fact * t1inp(p3,h4,spin_t1inp) * int2x(h4,p2,h1,p3,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t1_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_5(sp2,sh1,i0)

! i0 ( p2 h1 )_tf + = 1 * Sum ( p7 h8 ) * t ( p2 p7 h1 h8 )_t * i1 ( h8 p7 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: p7,h8,sp7,sh8
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp7 = 1,2
  do sh8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp2,sp7,sh1,sh8)
     spin_itm_hp = tdcc_spin_fock(sh8,sp7)
     if(spin_t2inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_t1_5_1(sh8,sp7,itm_hp)
     call ccsdt_t1_5_2(sh8,sp7,itm_hp)

     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do p7 = norb1+1,nact
     do h8 = 1,norb1
        i0(p2,h1) = i0(p2,h1) + fact * t2inp(p2,p7,h1,h8,spin_t2inp) * itm_hp(h8,p7)
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_t1_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_5_1(sh8,sp7,i1)

!     i1 ( h8 p7 )_f + = 1 * f ( h8 p7 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh8,sp7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer(c_int) :: h8,p7
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh8,sp7)
  do h8 = 1,norb1
  do p7 = norb1+1,nact
     i1(h8,p7) = i1(h8,p7) + fact * fock(h8,p7,spin_fock)
  end do
  end do
end subroutine ccsdt_t1_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_5_2(sh8,sp7,i1)

!     i1 ( h8 p7 )_vt + = 1 * Sum ( h6 p5 ) * t ( p5 h6 )_t * v ( h6 h8 p5 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh8,sp7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer(c_int) :: h8,p7
  integer(c_int) :: h6,p5,sh6,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh6 = 1,2
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_int2x = tdcc_spin_int2x(sh6,sh8,sp5,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h8 = 1,norb1
     do p7 = norb1+1,nact
     do h6 = 1,norb1
     do p5 = norb1+1,nact
        i1(h8,p7) = i1(h8,p7) + fact * t1inp(p5,h6,spin_t1inp) * int2x(h6,h8,p5,p7,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t1_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_6(sp2,sh1,i0)

! i0 ( p2 h1 )_vt + = -1/2 * Sum ( h4 h5 p3 ) * t ( p2 p3 h4 h5 )_t * i1 ( h4 h5 h1 p3 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: h4,h5,p3,sh4,sh5,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh4 = 1,2
  do sh5 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp2,sp3,sh4,sh5)
     spin_itm_hhhp = tdcc_spin_int2x(sh4,sh5,sh1,sp3)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t1_6_1(sh4,sh5,sh1,sp3,itm_hhhp)
     call ccsdt_t1_6_2(sh4,sh5,sh1,sp3,itm_hhhp)

     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
        i0(p2,h1) = i0(p2,h1) + fact * t2inp(p2,p3,h4,h5,spin_t2inp) * itm_hhhp(h4,h5,h1,p3)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_t1_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_6_1(sh4,sh5,sh1,sp3,i1)

!     i1 ( h4 h5 h1 p3 )_v + = 1 * v ( h4 h5 h1 p3 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh4,sh5,sh1,sp3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h4,h5,h1,p3
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh4,sh5,sh1,sp3)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h1 = 1,norb1
  do p3 = norb1+1,nact
     i1(h4,h5,h1,p3) = i1(h4,h5,h1,p3) + fact * int2x(h4,h5,h1,p3,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t1_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_6_2(sh4,sh5,sh1,sp3,i1)

!     i1 ( h4 h5 h1 p3 )_vt + = -1 * Sum ( p6 ) * t ( p6 h1 )_t * v ( h4 h5 p3 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh4,sh5,sh1,sp3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h4,h5,h1,p3
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh1)
     spin_int2x = tdcc_spin_int2x(sh4,sh5,sp3,sp6)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h1 = 1,norb1
     do p3 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h4,h5,h1,p3) = i1(h4,h5,h1,p3) + fact * t1inp(p6,h1,spin_t1inp) * int2x(h4,h5,p3,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t1_6_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_7(sp2,sh1,i0)

! i0 ( p2 h1 )_vt + = -1/2 * Sum ( h5 p3 p4 ) * t ( p3 p4 h1 h5 )_t * v ( h5 p2 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: h5,p3,p4,sh5,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh1,sh5)
     spin_int2x = tdcc_spin_int2x(sh5,sp2,sp3,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i0(p2,h1) = i0(p2,h1) + fact * t2inp(p3,p4,h1,h5,spin_t2inp) * int2x(h5,p2,p3,p4,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_t1_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t1_8(sp2,sh1,i0)

! i0 ( p2 h1 )_vt + = 1/4 * Sum ( h5 h6 p3 p4 ) * t ( p2 p3 p4 h1 h5 h6 )_t * v ( h5 h6 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: h5,h6,p3,p4,sh5,sh6,sp3,sp4
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp2,sp3,sp4,sh1,sh5,sh6)
     spin_int2x = tdcc_spin_int2x(sh5,sh6,sp3,sp4)
     if(spin_t3inp * spin_int2x == 0) cycle

     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i0(p2,h1) = i0(p2,h1) + fact * t3inp(p2,p3,p4,h1,h5,h6,spin_t3inp) * int2x(h5,h6,p3,p4,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  end do
end subroutine ccsdt_t1_8
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccsdt_t1_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1out

  implicit none

  call ccsdt_t1_1(1,1,t1out(norb1+1,1,1))
  call ccsdt_t1_2(1,1,t1out(norb1+1,1,1))
  call ccsdt_t1_3(1,1,t1out(norb1+1,1,1))
  call ccsdt_t1_4(1,1,t1out(norb1+1,1,1))
  call ccsdt_t1_5(1,1,t1out(norb1+1,1,1))
  call ccsdt_t1_6(1,1,t1out(norb1+1,1,1))
  call ccsdt_t1_7(1,1,t1out(norb1+1,1,1))
  call ccsdt_t1_8(1,1,t1out(norb1+1,1,1))

end subroutine ccsdt_t1_main
!**********************************************************
