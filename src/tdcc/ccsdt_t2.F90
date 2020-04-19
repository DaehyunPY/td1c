!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_1(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_v + = 1 * v ( p3 p4 h1 h2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sp3,sp4,sh1,sh2)
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
          int2x(p3,p4,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_vt + = -1 * P( 2 ) * Sum ( h9 ) * t ( p3 h9 )_t * i1 ( h9 p4 h1 h2 )_v 7

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t2_2_perm(sp3,sp4,sh1,sh2,i0_perm)
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_2_perm(sp4,sp3,sh1,sh2,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p4,p3,h1,h2)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_2_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: h9,sh9
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hphh
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh9)
     spin_itm_hphh = tdcc_spin_int2x(sh9,sp4,sh1,sh2)
     if(spin_t1inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_2_1(sh9,sp4,sh1,sh2,itm_hphh)
     call ccsdt_t2_2_2(sh9,sp4,sh1,sh2,itm_hphh)
     call ccsdt_t2_2_3(sh9,sp4,sh1,sh2,itm_hphh)
     call ccsdt_t2_2_4(sh9,sp4,sh1,sh2,itm_hphh)
     call ccsdt_t2_2_5(sh9,sp4,sh1,sh2,itm_hphh)
     call ccsdt_t2_2_6(sh9,sp4,sh1,sh2,itm_hphh)
     call ccsdt_t2_2_7(sh9,sp4,sh1,sh2,itm_hphh)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h9 = 1,norb1
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t1inp(p3,h9,spin_t1inp) * itm_hphh(h9,p4,h1,h2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphh)
  end subroutine ccsdt_t2_2_perm
  !--------------------------------------------
end subroutine ccsdt_t2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_1(sh9,sp3,sh1,sh2,i1)

!     i1 ( h9 p3 h1 h2 )_v + = 1 * v ( h9 p3 h1 h2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h9,p3,h1,h2
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh9,sp3,sh1,sh2)
  do h9 = 1,norb1
  do p3 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact * &
          int2x(h9,p3,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_2(sh9,sp3,sh1,sh2,i1)

!     i1 ( h9 p3 h1 h2 )_vt + = 1/2 * Sum ( h6 ) * t ( p3 h6 )_t * i2 ( h6 h9 h1 h2 )_v 3

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h9,p3,h1,h2
  integer(c_int) :: h6,sh6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hhhh
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh6)
     spin_itm_hhhh = tdcc_spin_int2x(sh6,sh9,sh1,sh2)
     if(spin_t1inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t2_2_2_1(sh6,sh9,sh1,sh2,itm_hhhh)
     call ccsdt_t2_2_2_2(sh6,sh9,sh1,sh2,itm_hhhh)
     call ccsdt_t2_2_2_3(sh6,sh9,sh1,sh2,itm_hhhh)

     do h9 = 1,norb1
     do p3 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h6 = 1,norb1
        i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact * &
             t1inp(p3,h6,spin_t1inp) * itm_hhhh(h6,h9,h1,h2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsdt_t2_2_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_2_1(sh6,sh9,sh1,sh2,i2)

!         i2 ( h6 h9 h1 h2 )_v + = 1 * v ( h6 h9 h1 h2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sh9,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h6,h9,h1,h2
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh6,sh9,sh1,sh2)
  do h6 = 1,norb1
  do h9 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i2(h6,h9,h1,h2) = i2(h6,h9,h1,h2) + fact * &
          int2x(h6,h9,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_2_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_2_2(sh6,sh9,sh1,sh2,i2)

!         i2 ( h6 h9 h1 h2 )_vt + = -1 * P( 2 ) * Sum ( p7 ) * t ( p7 h1 )_t * i3 ( h6 h9 h2 p7 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh6,sh9,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h6,h9,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i2_perm(:,:,:,:)

  allocate(i2_perm(1:norb1,1:norb1,1:norb1,1:norb1))

  i2_perm(1:norb1,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t2_2_2_2_perm(sh6,sh9,sh1,sh2,i2_perm)
  fact_p = +1.0d+0
  do h6 = 1,norb1
  do h9 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i2(h6,h9,h1,h2) = i2(h6,h9,h1,h2) + fact_p * i2_perm(h6,h9,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sh6 * sh9 * sh1 * sh2 == 1)) then
     i2_perm(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t2_2_2_2_perm(sh6,sh9,sh2,sh1,i2_perm)
  end if
  fact_p = -1.0d+0
  do h6 = 1,norb1
  do h9 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i2(h6,h9,h1,h2) = i2(h6,h9,h1,h2) + fact_p * i2_perm(h6,h9,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i2_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_2_2_2_perm(sh6,sh9,sh1,sh2,i2)

  implicit none
  integer(c_int),intent(in) :: sh6,sh9,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: h6,h9,h1,h2
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh1)
     spin_itm_hhhp = tdcc_spin_int2x(sh6,sh9,sh2,sp7)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_2_2_2_1(sh6,sh9,sh2,sp7,itm_hhhp)
     call ccsdt_t2_2_2_2_2(sh6,sh9,sh2,sp7,itm_hhhp)

     do h6 = 1,norb1
     do h9 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p7 = norb1+1,nact
        i2(h6,h9,h1,h2) = i2(h6,h9,h1,h2) + fact * &
             t1inp(p7,h1,spin_t1inp) * itm_hhhp(h6,h9,h2,p7)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_t2_2_2_2_perm
  !--------------------------------------------
end subroutine ccsdt_t2_2_2_2
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_2_2_1(sh6,sh9,sh1,sp7,i3)

!             i3 ( h6 h9 h1 p7 )_v + = 1 * v ( h6 h9 h1 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sh9,sh1,sp7
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,h9,h1,p7
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh6,sh9,sh1,sp7)
  do h6 = 1,norb1
  do h9 = 1,norb1
  do h1 = 1,norb1
  do p7 = norb1+1,nact
     i3(h6,h9,h1,p7) = i3(h6,h9,h1,p7) + fact * &
          int2x(h6,h9,h1,p7,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_2_2_2_1
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_2_2_2(sh6,sh9,sh1,sp7,i3)

!             i3 ( h6 h9 h1 p7 )_vt + = -1/2 * Sum ( p8 ) * t ( p8 h1 )_t * v ( h6 h9 p7 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sh9,sh1,sp7
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,h9,h1,p7
  integer(c_int) :: p8,sp8
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_int2x = tdcc_spin_int2x(sh6,sh9,sp7,sp8)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h6 = 1,norb1
     do h9 = 1,norb1
     do h1 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i3(h6,h9,h1,p7) = i3(h6,h9,h1,p7) + fact * &
             t1inp(p8,h1,spin_t1inp) * int2x(h6,h9,p7,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t2_2_2_2_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_2_3(sh6,sh9,sh1,sh2,i2)

!         i2 ( h6 h9 h1 h2 )_vt + = 1/2 * Sum ( p7 p8 ) * t ( p7 p8 h1 h2 )_t * v ( h6 h9 p7 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sh9,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h6,h9,h1,h2
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sh6,sh9,sp7,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h6 = 1,norb1
     do h9 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i2(h6,h9,h1,h2) = i2(h6,h9,h1,h2) + fact * &
             t2inp(p7,p8,h1,h2,spin_t2inp) * int2x(h6,h9,p7,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_2_2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_3(sh9,sp3,sh1,sh2,i1)

!     i1 ( h9 p3 h1 h2 )_vt + = -1 * P( 2 ) * Sum ( p6 ) * t ( p6 h1 )_t * i2 ( h9 p3 h2 p6 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h9,p3,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:)

  allocate(i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1))

  i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t2_2_3_perm(sh9,sp3,sh1,sh2,i1_perm)
  fact_p = +1.0d+0
  do h9 = 1,norb1
  do p3 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact_p * i1_perm(h9,p3,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sh9 * sp3 * sh1 * sh2 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_2_3_perm(sh9,sp3,sh2,sh1,i1_perm)
  end if
  fact_p = -1.0d+0
  do h9 = 1,norb1
  do p3 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact_p * i1_perm(h9,p3,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_2_3_perm(sh9,sp3,sh1,sh2,i1)

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: h9,p3,h1,h2
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hphp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh1)
     spin_itm_hphp = tdcc_spin_int2x(sh9,sp3,sh2,sp6)
     if(spin_t1inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_2_3_1(sh9,sp3,sh2,sp6,itm_hphp)
     call ccsdt_t2_2_3_2(sh9,sp3,sh2,sp6,itm_hphp)

     do h9 = 1,norb1
     do p3 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p6 = norb1+1,nact
        i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact * &
             t1inp(p6,h1,spin_t1inp) * itm_hphp(h9,p3,h2,p6)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphp)
  end subroutine ccsdt_t2_2_3_perm
  !--------------------------------------------
end subroutine ccsdt_t2_2_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_3_1(sh9,sp3,sh1,sp6,i2)

!         i2 ( h9 p3 h1 p6 )_v + = 1 * v ( h9 p3 h1 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sp6
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h9,p3,h1,p6
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh9,sp3,sh1,sp6)
  do h9 = 1,norb1
  do p3 = norb1+1,nact
  do h1 = 1,norb1
  do p6 = norb1+1,nact
     i2(h9,p3,h1,p6) = i2(h9,p3,h1,p6) + fact * &
          int2x(h9,p3,h1,p6,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_2_3_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_3_2(sh9,sp3,sh1,sp6,i2)

!         i2 ( h9 p3 h1 p6 )_vt + = -1/2 * Sum ( p7 ) * t ( p7 h1 )_t * v ( h9 p3 p6 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sp6
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h9,p3,h1,p6
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh1)
     spin_int2x = tdcc_spin_int2x(sh9,sp3,sp6,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do p3 = norb1+1,nact
     do h1 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h9,p3,h1,p6) = i2(h9,p3,h1,p6) + fact * &
             t1inp(p7,h1,spin_t1inp) * int2x(h9,p3,p6,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t2_2_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_4(sh9,sp3,sh1,sh2,i1)

!     i1 ( h9 p3 h1 h2 )_ft + = -1 * Sum ( p5 ) * t ( p3 p5 h1 h2 )_t * i2 ( h9 p5 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h9,p3,h1,h2
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp5,sh1,sh2)
     spin_itm_hp = tdcc_spin_fock(sh9,sp5)
     if(spin_t2inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_2_4_1(sh9,sp5,itm_hp)
     call ccsdt_t2_2_4_2(sh9,sp5,itm_hp)

     do h9 = 1,norb1
     do p3 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p5 = norb1+1,nact
        i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact * &
             t2inp(p3,p5,h1,h2,spin_t2inp) * itm_hp(h9,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_t2_2_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_4_1(sh9,sp5,i2)

!         i2 ( h9 p5 )_f + = 1 * f ( h9 p5 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh9,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h9,p5
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh9,sp5)
  do h9 = 1,norb1
  do p5 = norb1+1,nact
     i2(h9,p5) = i2(h9,p5) + fact * &
          fock(h9,p5,spin_fock)
  end do
  end do
end subroutine ccsdt_t2_2_4_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_4_2(sh9,sp5,i2)

!         i2 ( h9 p5 )_vt + = -1 * Sum ( h7 p6 ) * t ( p6 h7 )_t * v ( h7 h9 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h9,p5
  integer(c_int) :: h7,p6,sh7,sp6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh7)
     spin_int2x = tdcc_spin_int2x(sh7,sh9,sp5,sp6)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do p5 = norb1+1,nact
     do h7 = 1,norb1
     do p6 = norb1+1,nact
        i2(h9,p5) = i2(h9,p5) + fact * &
             t1inp(p6,h7,spin_t1inp) * int2x(h7,h9,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_2_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_5(sh9,sp3,sh1,sh2,i1)

!     i1 ( h9 p3 h1 h2 )_vt + = 1 * P( 2 ) * Sum ( h5 p10 ) * t ( p3 p10 h1 h5 )_t * i2 ( h5 h9 h2 p10 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h9,p3,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:)

  allocate(i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1))

  i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t2_2_5_perm(sh9,sp3,sh1,sh2,i1_perm)
  fact_p = +1.0d+0
  do h9 = 1,norb1
  do p3 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact_p * i1_perm(h9,p3,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sh9 * sp3 * sh1 * sh2 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_2_5_perm(sh9,sp3,sh2,sh1,i1_perm)
  end if
  fact_p = -1.0d+0
  do h9 = 1,norb1
  do p3 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact_p * i1_perm(h9,p3,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_2_5_perm(sh9,sp3,sh1,sh2,i1)

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: h9,p3,h1,h2
  integer(c_int) :: h5,p10,sh5,sp10
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh5 = 1,2
  do sp10 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp10,sh1,sh5)
     spin_itm_hhhp = tdcc_spin_int2x(sh5,sh9,sh2,sp10)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_2_5_1(sh5,sh9,sh2,sp10,itm_hhhp)
     call ccsdt_t2_2_5_2(sh5,sh9,sh2,sp10,itm_hhhp)

     do h9 = 1,norb1
     do p3 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h5 = 1,norb1
     do p10 = norb1+1,nact
        i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact * &
             t2inp(p3,p10,h1,h5,spin_t2inp) * itm_hhhp(h5,h9,h2,p10)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_t2_2_5_perm
  !--------------------------------------------
end subroutine ccsdt_t2_2_5
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_5_1(sh5,sh9,sh1,sp10,i2)

!         i2 ( h5 h9 h1 p10 )_v + = 1 * v ( h5 h9 h1 p10 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh5,sh9,sh1,sp10
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h5,h9,h1,p10
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh5,sh9,sh1,sp10)
  do h5 = 1,norb1
  do h9 = 1,norb1
  do h1 = 1,norb1
  do p10 = norb1+1,nact
     i2(h5,h9,h1,p10) = i2(h5,h9,h1,p10) + fact * &
          int2x(h5,h9,h1,p10,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_2_5_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_5_2(sh5,sh9,sh1,sp10,i2)

!         i2 ( h5 h9 h1 p10 )_vt + = 1 * Sum ( p6 ) * t ( p6 h1 )_t * v ( h5 h9 p6 p10 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh5,sh9,sh1,sp10
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h5,h9,h1,p10
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh1)
     spin_int2x = tdcc_spin_int2x(sh5,sh9,sp6,sp10)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h5 = 1,norb1
     do h9 = 1,norb1
     do h1 = 1,norb1
     do p10 = norb1+1,nact
     do p6 = norb1+1,nact
        i2(h5,h9,h1,p10) = i2(h5,h9,h1,p10) + fact * &
             t1inp(p6,h1,spin_t1inp) * int2x(h5,h9,p6,p10,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t2_2_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_6(sh9,sp3,sh1,sh2,i1)

!     i1 ( h9 p3 h1 h2 )_vt + = 1/2 * Sum ( p5 p6 ) * t ( p5 p6 h1 h2 )_t * v ( h9 p3 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h9,p3,h1,h2
  integer(c_int) :: p5,p6,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sh9,sp3,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do p3 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact * &
             t2inp(p5,p6,h1,h2,spin_t2inp) * int2x(h9,p3,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_2_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_2_7(sh9,sp3,sh1,sh2,i1)

!     i1 ( h9 p3 h1 h2 )_vt + = 1/2 * Sum ( h7 p5 p6 ) * t ( p3 p5 p6 h1 h2 h7 )_t * v ( h7 h9 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sp3,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h9,p3,h1,h2
  integer(c_int) :: h7,p5,p6,sh7,sp5,sp6
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp5,sp6,sh1,sh2,sh7)
     spin_int2x = tdcc_spin_int2x(sh7,sh9,sp5,sp6)
     if(spin_t3inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do p3 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h9,p3,h1,h2) = i1(h9,p3,h1,h2) + fact * &
             t3inp(p3,p5,p6,h1,h2,h7,spin_t3inp) * int2x(h7,h9,p5,p6,spin_int2x)
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
end subroutine ccsdt_t2_2_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_3(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_vt + = -1 * P( 2 ) * Sum ( p5 ) * t ( p5 h1 )_t * i1 ( p3 p4 h2 p5 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t2_3_perm(sp3,sp4,sh1,sh2,i0_perm)
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_3_perm(sp3,sp4,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_3_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_pphp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh1)
     spin_itm_pphp = tdcc_spin_int2x(sp3,sp4,sh2,sp5)
     if(spin_t1inp * spin_itm_pphp == 0) cycle

     itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_3_1(sp3,sp4,sh2,sp5,itm_pphp)
     call ccsdt_t2_3_2(sp3,sp4,sh2,sp5,itm_pphp)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p5 = norb1+1,nact
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t1inp(p5,h1,spin_t1inp) * itm_pphp(p3,p4,h2,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pphp)
  end subroutine ccsdt_t2_3_perm
  !--------------------------------------------
end subroutine ccsdt_t2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_3_1(sp3,sp4,sh1,sp5,i1)

!     i1 ( p3 p4 h1 p5 )_v + = 1 * v ( p3 p4 h1 p5 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sp5
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: p3,p4,h1,p5
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sp3,sp4,sh1,sp5)
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do p5 = norb1+1,nact
     i1(p3,p4,h1,p5) = i1(p3,p4,h1,p5) + fact * &
          int2x(p3,p4,h1,p5,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_3_2(sp3,sp4,sh1,sp5,i1)

!     i1 ( p3 p4 h1 p5 )_vt + = -1/2 * Sum ( p6 ) * t ( p6 h1 )_t * v ( p3 p4 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sp5
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: p3,p4,h1,p5
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh1)
     spin_int2x = tdcc_spin_int2x(sp3,sp4,sp5,sp6)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(p3,p4,h1,p5) = i1(p3,p4,h1,p5) + fact * &
             t1inp(p6,h1,spin_t1inp) * int2x(p3,p4,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t2_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_4(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_tf + = -1 * P( 2 ) * Sum ( h9 ) * t ( p3 p4 h1 h9 )_t * i1 ( h9 h2 )_f 4

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t2_4_perm(sp3,sp4,sh1,sh2,i0_perm)
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_4_perm(sp3,sp4,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_4_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: h9,sh9
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hh
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh1,sh9)
     spin_itm_hh = tdcc_spin_fock(sh9,sh2)
     if(spin_t2inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsdt_t2_4_1(sh9,sh2,itm_hh)
     call ccsdt_t2_4_2(sh9,sh2,itm_hh)
     call ccsdt_t2_4_3(sh9,sh2,itm_hh)
     call ccsdt_t2_4_4(sh9,sh2,itm_hh)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h9 = 1,norb1
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t2inp(p3,p4,h1,h9,spin_t2inp) * itm_hh(h9,h2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccsdt_t2_4_perm
  !--------------------------------------------
end subroutine ccsdt_t2_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_4_1(sh9,sh1,i1)

!     i1 ( h9 h1 )_f + = 1 * f ( h9 h1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh9,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h9,h1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh9,sh1)
  do h9 = 1,norb1
  do h1 = 1,norb1
     i1(h9,h1) = i1(h9,h1) + fact * &
          fock(h9,h1,spin_fock)
  end do
  end do
end subroutine ccsdt_t2_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_4_2(sh9,sh1,i1)

!     i1 ( h9 h1 )_ft + = 1 * Sum ( p8 ) * t ( p8 h1 )_t * i2 ( h9 p8 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh9,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h9,h1
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
     spin_itm_hp = tdcc_spin_fock(sh9,sp8)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_4_2_1(sh9,sp8,itm_hp)
     call ccsdt_t2_4_2_2(sh9,sp8,itm_hp)

     do h9 = 1,norb1
     do h1 = 1,norb1
     do p8 = norb1+1,nact
        i1(h9,h1) = i1(h9,h1) + fact * &
             t1inp(p8,h1,spin_t1inp) * itm_hp(h9,p8)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_t2_4_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_4_2_1(sh9,sp8,i2)

!         i2 ( h9 p8 )_f + = 1 * f ( h9 p8 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh9,sp8
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h9,p8
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh9,sp8)
  do h9 = 1,norb1
  do p8 = norb1+1,nact
     i2(h9,p8) = i2(h9,p8) + fact * &
          fock(h9,p8,spin_fock)
  end do
  end do
end subroutine ccsdt_t2_4_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_4_2_2(sh9,sp8,i2)

!         i2 ( h9 p8 )_vt + = 1 * Sum ( h7 p6 ) * t ( p6 h7 )_t * v ( h7 h9 p6 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sp8
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer(c_int) :: h9,p8
  integer(c_int) :: h7,p6,sh7,sp6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh7 = 1,2
  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh7)
     spin_int2x = tdcc_spin_int2x(sh7,sh9,sp6,sp8)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do p8 = norb1+1,nact
     do h7 = 1,norb1
     do p6 = norb1+1,nact
        i2(h9,p8) = i2(h9,p8) + fact * &
             t1inp(p6,h7,spin_t1inp) * int2x(h7,h9,p6,p8,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_4_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_4_3(sh9,sh1,i1)

!     i1 ( h9 h1 )_vt + = -1 * Sum ( h7 p6 ) * t ( p6 h7 )_t * v ( h7 h9 h1 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h9,h1
  integer(c_int) :: h7,p6,sh7,sp6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh7)
     spin_int2x = tdcc_spin_int2x(sh7,sh9,sh1,sp6)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do h1 = 1,norb1
     do h7 = 1,norb1
     do p6 = norb1+1,nact
        i1(h9,h1) = i1(h9,h1) + fact * &
             t1inp(p6,h7,spin_t1inp) * int2x(h7,h9,h1,p6,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_4_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_4_4(sh9,sh1,i1)

!     i1 ( h9 h1 )_vt + = -1/2 * Sum ( h8 p6 p7 ) * t ( p6 p7 h1 h8 )_t * v ( h8 h9 p6 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h9,h1
  integer(c_int) :: h8,p6,p7,sh8,sp6,sp7
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp7,sh1,sh8)
     spin_int2x = tdcc_spin_int2x(sh8,sh9,sp6,sp7)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do h1 = 1,norb1
     do h8 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h9,h1) = i1(h9,h1) + fact * &
             t2inp(p6,p7,h1,h8,spin_t2inp) * int2x(h8,h9,p6,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_t2_4_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_5(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_tf + = 1 * P( 2 ) * Sum ( p5 ) * t ( p3 p5 h1 h2 )_t * i1 ( p4 p5 )_f 3

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t2_5_perm(sp3,sp4,sh1,sh2,i0_perm)
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_5_perm(sp4,sp3,sh1,sh2,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p4,p3,h1,h2)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_5_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_pp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp5,sh1,sh2)
     spin_itm_pp = tdcc_spin_fock(sp4,sp5)
     if(spin_t2inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_t2_5_1(sp4,sp5,itm_pp)
     call ccsdt_t2_5_2(sp4,sp5,itm_pp)
     call ccsdt_t2_5_3(sp4,sp5,itm_pp)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p5 = norb1+1,nact
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t2inp(p3,p5,h1,h2,spin_t2inp) * itm_pp(p4,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccsdt_t2_5_perm
  !--------------------------------------------
end subroutine ccsdt_t2_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_5_1(sp3,sp5,i1)

!     i1 ( p3 p5 )_f + = 1 * f ( p3 p5 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sp3,sp5
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p3,p5
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp3,sp5)
  do p3 = norb1+1,nact
  do p5 = norb1+1,nact
     i1(p3,p5) = i1(p3,p5) + fact * &
          fock(p3,p5,spin_fock)
  end do
  end do
end subroutine ccsdt_t2_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_5_2(sp3,sp5,i1)

!     i1 ( p3 p5 )_vt + = -1 * Sum ( h7 p6 ) * t ( p6 h7 )_t * v ( h7 p3 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp3,sp5
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p3,p5
  integer(c_int) :: h7,p6,sh7,sp6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh7)
     spin_int2x = tdcc_spin_int2x(sh7,sp3,sp5,sp6)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p3 = norb1+1,nact
     do p5 = norb1+1,nact
     do h7 = 1,norb1
     do p6 = norb1+1,nact
        i1(p3,p5) = i1(p3,p5) + fact * &
             t1inp(p6,h7,spin_t1inp) * int2x(h7,p3,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_5_3(sp3,sp5,i1)

!     i1 ( p3 p5 )_vt + = -1/2 * Sum ( h7 h8 p6 ) * t ( p3 p6 h7 h8 )_t * v ( h7 h8 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp3,sp5
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p3,p5
  integer(c_int) :: h7,h8,p6,sh7,sh8,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp6,sh7,sh8)
     spin_int2x = tdcc_spin_int2x(sh7,sh8,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do p3 = norb1+1,nact
     do p5 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i1(p3,p5) = i1(p3,p5) + fact * &
             t2inp(p3,p6,h7,h8,spin_t2inp) * int2x(h7,h8,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_t2_5_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_6(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_vt + = 1/2 * Sum ( h9 h10 ) * t ( p3 p4 h9 h10 )_t * i1 ( h9 h10 h1 h2 )_v 3

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: h9,h10,sh9,sh10
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hhhh
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh10 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh9,sh10)
     spin_itm_hhhh = tdcc_spin_int2x(sh9,sh10,sh1,sh2)
     if(spin_t2inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t2_6_1(sh9,sh10,sh1,sh2,itm_hhhh)
     call ccsdt_t2_6_2(sh9,sh10,sh1,sh2,itm_hhhh)
     call ccsdt_t2_6_3(sh9,sh10,sh1,sh2,itm_hhhh)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t2inp(p3,p4,h9,h10,spin_t2inp) * itm_hhhh(h9,h10,h1,h2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsdt_t2_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_6_1(sh9,sh10,sh1,sh2,i1)

!     i1 ( h9 h10 h1 h2 )_v + = 1 * v ( h9 h10 h1 h2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sh10,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h9,h10,h1,h2
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh9,sh10,sh1,sh2)
  do h9 = 1,norb1
  do h10 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h9,h10,h1,h2) = i1(h9,h10,h1,h2) + fact * &
          int2x(h9,h10,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_6_2(sh9,sh10,sh1,sh2,i1)

!     i1 ( h9 h10 h1 h2 )_vt + = -1 * P( 2 ) * Sum ( p8 ) * t ( p8 h1 )_t * i2 ( h9 h10 h2 p8 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh9,sh10,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h9,h10,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:)

  allocate(i1_perm(1:norb1,1:norb1,1:norb1,1:norb1))

  i1_perm(1:norb1,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t2_6_2_perm(sh9,sh10,sh1,sh2,i1_perm)
  fact_p = +1.0d+0
  do h9 = 1,norb1
  do h10 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h9,h10,h1,h2) = i1(h9,h10,h1,h2) + fact_p * i1_perm(h9,h10,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sh9 * sh10 * sh1 * sh2 == 1)) then
     i1_perm(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t2_6_2_perm(sh9,sh10,sh2,sh1,i1_perm)
  end if
  fact_p = -1.0d+0
  do h9 = 1,norb1
  do h10 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h9,h10,h1,h2) = i1(h9,h10,h1,h2) + fact_p * i1_perm(h9,h10,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_6_2_perm(sh9,sh10,sh1,sh2,i1)

  implicit none
  integer(c_int),intent(in) :: sh9,sh10,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: h9,h10,h1,h2
  integer(c_int) :: p8,sp8
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_itm_hhhp = tdcc_spin_int2x(sh9,sh10,sh2,sp8)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_6_2_1(sh9,sh10,sh2,sp8,itm_hhhp)
     call ccsdt_t2_6_2_2(sh9,sh10,sh2,sp8,itm_hhhp)

     do h9 = 1,norb1
     do h10 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p8 = norb1+1,nact
        i1(h9,h10,h1,h2) = i1(h9,h10,h1,h2) + fact * &
             t1inp(p8,h1,spin_t1inp) * itm_hhhp(h9,h10,h2,p8)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_t2_6_2_perm
  !--------------------------------------------
end subroutine ccsdt_t2_6_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_6_2_1(sh9,sh10,sh1,sp8,i2)

!         i2 ( h9 h10 h1 p8 )_v + = 1 * v ( h9 h10 h1 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sh10,sh1,sp8
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h9,h10,h1,p8
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh9,sh10,sh1,sp8)
  do h9 = 1,norb1
  do h10 = 1,norb1
  do h1 = 1,norb1
  do p8 = norb1+1,nact
     i2(h9,h10,h1,p8) = i2(h9,h10,h1,p8) + fact * &
          int2x(h9,h10,h1,p8,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_6_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_6_2_2(sh9,sh10,sh1,sp8,i2)

!         i2 ( h9 h10 h1 p8 )_vt + = 1/2 * Sum ( p6 ) * t ( p6 h1 )_t * v ( h9 h10 p6 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sh10,sh1,sp8
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h9,h10,h1,p8
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh1)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp6,sp8)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do h10 = 1,norb1
     do h1 = 1,norb1
     do p8 = norb1+1,nact
     do p6 = norb1+1,nact
        i2(h9,h10,h1,p8) = i2(h9,h10,h1,p8) + fact * &
             t1inp(p6,h1,spin_t1inp) * int2x(h9,h10,p6,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t2_6_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_6_3(sh9,sh10,sh1,sh2,i1)

!     i1 ( h9 h10 h1 h2 )_vt + = 1/2 * Sum ( p7 p8 ) * t ( p7 p8 h1 h2 )_t * v ( h9 h10 p7 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sh10,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h9,h10,h1,h2
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp7,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do h10 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h9,h10,h1,h2) = i1(h9,h10,h1,h2) + fact * &
             t2inp(p7,p8,h1,h2,spin_t2inp) * int2x(h9,h10,p7,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_6_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_7(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_vt + = -1 * P( 4 ) * Sum ( h6 p5 ) * t ( p3 p5 h1 h6 )_t * i1 ( h6 p4 h2 p5 )_v 3

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t2_7_perm(sp3,sp4,sh1,sh2,i0_perm)
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_7_perm(sp3,sp4,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h2,h1)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_7_perm(sp4,sp3,sh1,sh2,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p4,p3,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_7_perm(sp4,sp3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p4,p3,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_7_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: h6,p5,sh6,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hphp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh6 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp5,sh1,sh6)
     spin_itm_hphp = tdcc_spin_int2x(sh6,sp4,sh2,sp5)
     if(spin_t2inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_7_1(sh6,sp4,sh2,sp5,itm_hphp)
     call ccsdt_t2_7_2(sh6,sp4,sh2,sp5,itm_hphp)
     call ccsdt_t2_7_3(sh6,sp4,sh2,sp5,itm_hphp)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h6 = 1,norb1
     do p5 = norb1+1,nact
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t2inp(p3,p5,h1,h6,spin_t2inp) * itm_hphp(h6,p4,h2,p5)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hphp)
  end subroutine ccsdt_t2_7_perm
  !--------------------------------------------
end subroutine ccsdt_t2_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_7_1(sh6,sp3,sh1,sp5,i1)

!     i1 ( h6 p3 h1 p5 )_v + = 1 * v ( h6 p3 h1 p5 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sp3,sh1,sp5
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,p3,h1,p5
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh6,sp3,sh1,sp5)
  do h6 = 1,norb1
  do p3 = norb1+1,nact
  do h1 = 1,norb1
  do p5 = norb1+1,nact
     i1(h6,p3,h1,p5) = i1(h6,p3,h1,p5) + fact * &
          int2x(h6,p3,h1,p5,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_7_2(sh6,sp3,sh1,sp5,i1)

!     i1 ( h6 p3 h1 p5 )_vt + = -1 * Sum ( p7 ) * t ( p7 h1 )_t * v ( h6 p3 p5 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sp3,sh1,sp5
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,p3,h1,p5
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh1)
     spin_int2x = tdcc_spin_int2x(sh6,sp3,sp5,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do h1 = 1,norb1
     do p5 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h6,p3,h1,p5) = i1(h6,p3,h1,p5) + fact * &
             t1inp(p7,h1,spin_t1inp) * int2x(h6,p3,p5,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t2_7_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_7_3(sh6,sp3,sh1,sp5,i1)

!     i1 ( h6 p3 h1 p5 )_vt + = -1/2 * Sum ( h8 p7 ) * t ( p3 p7 h1 h8 )_t * v ( h6 h8 p5 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sp3,sh1,sp5
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,p3,h1,p5
  integer(c_int) :: h8,p7,sh8,sp7
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp7,sh1,sh8)
     spin_int2x = tdcc_spin_int2x(sh6,sh8,sp5,sp7)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do h1 = 1,norb1
     do p5 = norb1+1,nact
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i1(h6,p3,h1,p5) = i1(h6,p3,h1,p5) + fact * &
             t2inp(p3,p7,h1,h8,spin_t2inp) * int2x(h6,h8,p5,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_7_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_8(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_vt + = 1/2 * Sum ( p5 p6 ) * t ( p5 p6 h1 h2 )_t * v ( p3 p4 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: p5,p6,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sp3,sp4,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t2inp(p5,p6,h1,h2,spin_t2inp) * int2x(p3,p4,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_9(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_tf + = 1 * Sum ( p9 h10 ) * t ( p3 p4 p9 h1 h2 h10 )_t * i1 ( h10 p9 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: p9,h10,sp9,sh10
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_itm_hp
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp9 = 1,2
  do sh10 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp4,sp9,sh1,sh2,sh10)
     spin_itm_hp = tdcc_spin_fock(sh10,sp9)
     if(spin_t3inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_9_1(sh10,sp9,itm_hp)
     call ccsdt_t2_9_2(sh10,sp9,itm_hp)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h10 = 1,norb1
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t3inp(p3,p4,p9,h1,h2,h10,spin_t3inp) * itm_hp(h10,p9)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_t2_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_9_1(sh10,sp9,i1)

!     i1 ( h10 p9 )_f + = 1 * f ( h10 p9 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh10,sp9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer(c_int) :: h10,p9
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh10,sp9)
  do h10 = 1,norb1
  do p9 = norb1+1,nact
     i1(h10,p9) = i1(h10,p9) + fact * &
          fock(h10,p9,spin_fock)
  end do
  end do
end subroutine ccsdt_t2_9_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_9_2(sh10,sp9,i1)

!     i1 ( h10 p9 )_vt + = 1 * Sum ( h8 p7 ) * t ( p7 h8 )_t * v ( h8 h10 p7 p9 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh10,sp9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer(c_int) :: h10,p9
  integer(c_int) :: h8,p7,sh8,sp7
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_int2x = tdcc_spin_int2x(sh8,sh10,sp7,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h10 = 1,norb1
     do p9 = norb1+1,nact
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i1(h10,p9) = i1(h10,p9) + fact * &
             t1inp(p7,h8,spin_t1inp) * int2x(h8,h10,p7,p9,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t2_9_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_10(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_vt + = -1/2 * P( 2 ) * Sum ( h6 h7 p5 ) * t ( p3 p4 p5 h1 h6 h7 )_t * i1 ( h6 h7 h2 p5 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t2_10_perm(sp3,sp4,sh1,sh2,i0_perm)
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_10_perm(sp3,sp4,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_10_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: h6,h7,p5,sh6,sh7,sp5
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_itm_hhhp
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh6 = 1,2
  do sh7 = 1,2
  do sp5 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp4,sp5,sh1,sh6,sh7)
     spin_itm_hhhp = tdcc_spin_int2x(sh6,sh7,sh2,sp5)
     if(spin_t3inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t2_10_1(sh6,sh7,sh2,sp5,itm_hhhp)
     call ccsdt_t2_10_2(sh6,sh7,sh2,sp5,itm_hhhp)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p5 = norb1+1,nact
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t3inp(p3,p4,p5,h1,h6,h7,spin_t3inp) * itm_hhhp(h6,h7,h2,p5)
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
  deallocate(itm_hhhp)
  end subroutine ccsdt_t2_10_perm
  !--------------------------------------------
end subroutine ccsdt_t2_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_10_1(sh6,sh7,sh1,sp5,i1)

!     i1 ( h6 h7 h1 p5 )_v + = 1 * v ( h6 h7 h1 p5 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sh7,sh1,sp5
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,h7,h1,p5
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh6,sh7,sh1,sp5)
  do h6 = 1,norb1
  do h7 = 1,norb1
  do h1 = 1,norb1
  do p5 = norb1+1,nact
     i1(h6,h7,h1,p5) = i1(h6,h7,h1,p5) + fact * &
          int2x(h6,h7,h1,p5,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t2_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_10_2(sh6,sh7,sh1,sp5,i1)

!     i1 ( h6 h7 h1 p5 )_vt + = -1 * Sum ( p8 ) * t ( p8 h1 )_t * v ( h6 h7 p5 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sh7,sh1,sp5
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,h7,h1,p5
  integer(c_int) :: p8,sp8
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_int2x = tdcc_spin_int2x(sh6,sh7,sp5,sp8)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h6 = 1,norb1
     do h7 = 1,norb1
     do h1 = 1,norb1
     do p5 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h6,h7,h1,p5) = i1(h6,h7,h1,p5) + fact * &
             t1inp(p8,h1,spin_t1inp) * int2x(h6,h7,p5,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t2_10_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t2_11(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_vt + = -1/2 * P( 2 ) * Sum ( h7 p5 p6 ) * t ( p3 p5 p6 h1 h2 h7 )_t * v ( h7 p4 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t2_11_perm(sp3,sp4,sh1,sh2,i0_perm)
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t2_11_perm(sp4,sp3,sh1,sh2,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p4,p3,h1,h2)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t2_11_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: h7,p5,p6,sh7,sp5,sp6
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp5,sp6,sh1,sh2,sh7)
     spin_int2x = tdcc_spin_int2x(sh7,sp4,sp5,sp6)
     if(spin_t3inp * spin_int2x == 0) cycle

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * &
             t3inp(p3,p5,p6,h1,h2,h7,spin_t3inp) * int2x(h7,p4,p5,p6,spin_int2x)
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
  end subroutine ccsdt_t2_11_perm
  !--------------------------------------------
end subroutine ccsdt_t2_11
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccsdt_t2_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2out

  implicit none

  call ccsdt_t2_1(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_2(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_3(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_4(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_5(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_6(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_7(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_8(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_9(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_10(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
  call ccsdt_t2_11(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))

  call ccsdt_t2_1(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_2(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_3(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_4(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_5(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_6(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_7(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_8(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_9(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_10(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  call ccsdt_t2_11(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))

end subroutine ccsdt_t2_main
!**********************************************************
