!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_1(sp3,sp4,sh1,sh2,i0)

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
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * int2x(p3,p4,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccd_t2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_2(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_tf + = -1 * P( 2 ) * Sum ( h5 ) * t ( p3 p4 h1 h5 )_t * i1 ( h5 h2 )_f 2

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
  call ccd_t2_2_perm(sp3,sp4,sh1,sh2,i0_perm)
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
     call ccd_t2_2_perm(sp3,sp4,sh2,sh1,i0_perm)
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
  subroutine ccd_t2_2_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hh
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh1,sh5)
     spin_itm_hh = tdcc_spin_fock(sh5,sh2)
     if(spin_t2inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccd_t2_2_1(sh5,sh2,itm_hh)
     call ccd_t2_2_2(sh5,sh2,itm_hh)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h5 = 1,norb1
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * t2inp(p3,p4,h1,h5,spin_t2inp) * itm_hh(h5,h2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccd_t2_2_perm
  !--------------------------------------------
end subroutine ccd_t2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_2_1(sh5,sh1,i1)

!     i1 ( h5 h1 )_f + = 1 * f ( h5 h1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh5,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h5,h1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh5,sh1)
  do h5 = 1,norb1
  do h1 = 1,norb1
     i1(h5,h1) = i1(h5,h1) + fact * fock(h5,h1,spin_fock)
  end do
  end do
end subroutine ccd_t2_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_2_2(sh5,sh1,i1)

!     i1 ( h5 h1 )_vt + = 1/2 * Sum ( h8 p6 p7 ) * t ( p6 p7 h1 h8 )_t * v ( h5 h8 p6 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh5,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h5,h1
  integer(c_int) :: h8,p6,p7,sh8,sp6,sp7
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp7,sh1,sh8)
     spin_int2x = tdcc_spin_int2x(sh5,sh8,sp6,sp7)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h5 = 1,norb1
     do h1 = 1,norb1
     do h8 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h5,h1) = i1(h5,h1) + fact * t2inp(p6,p7,h1,h8,spin_t2inp) * int2x(h5,h8,p6,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccd_t2_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_3(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_tf + = 1 * P( 2 ) * Sum ( p5 ) * t ( p3 p5 h1 h2 )_t * i1 ( p4 p5 )_f 2

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
  call ccd_t2_3_perm(sp3,sp4,sh1,sh2,i0_perm)
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
     call ccd_t2_3_perm(sp4,sp3,sh1,sh2,i0_perm)
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
  subroutine ccd_t2_3_perm(sp3,sp4,sh1,sh2,i0)

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
     call ccd_t2_3_1(sp4,sp5,itm_pp)
     call ccd_t2_3_2(sp4,sp5,itm_pp)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p5 = norb1+1,nact
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * t2inp(p3,p5,h1,h2,spin_t2inp) * itm_pp(p4,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccd_t2_3_perm
  !--------------------------------------------
end subroutine ccd_t2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_3_1(sp3,sp5,i1)

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
     i1(p3,p5) = i1(p3,p5) + fact * fock(p3,p5,spin_fock)
  end do
  end do
end subroutine ccd_t2_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_3_2(sp3,sp5,i1)

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
        i1(p3,p5) = i1(p3,p5) + fact * t2inp(p3,p6,h7,h8,spin_t2inp) * int2x(h7,h8,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccd_t2_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_4(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_vt + = 1/2 * Sum ( h5 h6 ) * t ( p3 p4 h5 h6 )_t * i1 ( h5 h6 h1 h2 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: p3,p4,h1,h2
  integer(c_int) :: h5,h6,sh5,sh6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hhhh
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh5 = 1,2
  do sh6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh6)
     spin_itm_hhhh = tdcc_spin_int2x(sh5,sh6,sh1,sh2)
     if(spin_t2inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccd_t2_4_1(sh5,sh6,sh1,sh2,itm_hhhh)
     call ccd_t2_4_2(sh5,sh6,sh1,sh2,itm_hhhh)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * t2inp(p3,p4,h5,h6,spin_t2inp) * itm_hhhh(h5,h6,h1,h2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccd_t2_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_4_1(sh5,sh6,sh1,sh2,i1)

!     i1 ( h5 h6 h1 h2 )_v + = 1 * v ( h5 h6 h1 h2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh5,sh6,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h5,h6,h1,h2
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh5,sh6,sh1,sh2)
  do h5 = 1,norb1
  do h6 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h5,h6,h1,h2) = i1(h5,h6,h1,h2) + fact * int2x(h5,h6,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccd_t2_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_4_2(sh5,sh6,sh1,sh2,i1)

!     i1 ( h5 h6 h1 h2 )_vt + = 1/2 * Sum ( p7 p8 ) * t ( p7 p8 h1 h2 )_t * v ( h5 h6 p7 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh5,sh6,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h5,h6,h1,h2
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sh5,sh6,sp7,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h5 = 1,norb1
     do h6 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h5,h6,h1,h2) = i1(h5,h6,h1,h2) + fact * t2inp(p7,p8,h1,h2,spin_t2inp) * int2x(h5,h6,p7,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_t2_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_5(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_vt + = -1 * P( 4 ) * Sum ( h6 p5 ) * t ( p3 p5 h1 h6 )_t * i1 ( h6 p4 h2 p5 )_v 2

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
  call ccd_t2_5_perm(sp3,sp4,sh1,sh2,i0_perm)
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
     call ccd_t2_5_perm(sp3,sp4,sh2,sh1,i0_perm)
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
     call ccd_t2_5_perm(sp4,sp3,sh1,sh2,i0_perm)
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
     call ccd_t2_5_perm(sp4,sp3,sh2,sh1,i0_perm)
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
  subroutine ccd_t2_5_perm(sp3,sp4,sh1,sh2,i0)

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
     call ccd_t2_5_1(sh6,sp4,sh2,sp5,itm_hphp)
     call ccd_t2_5_2(sh6,sp4,sh2,sp5,itm_hphp)

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h6 = 1,norb1
     do p5 = norb1+1,nact
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * t2inp(p3,p5,h1,h6,spin_t2inp) * itm_hphp(h6,p4,h2,p5)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hphp)
  end subroutine ccd_t2_5_perm
  !--------------------------------------------
end subroutine ccd_t2_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_5_1(sh6,sp3,sh1,sp5,i1)

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
     i1(h6,p3,h1,p5) = i1(h6,p3,h1,p5) + fact * int2x(h6,p3,h1,p5,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccd_t2_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_5_2(sh6,sp3,sh1,sp5,i1)

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
        i1(h6,p3,h1,p5) = i1(h6,p3,h1,p5) + fact * t2inp(p3,p7,h1,h8,spin_t2inp) * int2x(h6,h8,p5,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_t2_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_t2_6(sp3,sp4,sh1,sh2,i0)

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
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * t2inp(p5,p6,h1,h2,spin_t2inp) * int2x(p3,p4,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_t2_6
!##########################################################
!##########################################################
!##########################################################
