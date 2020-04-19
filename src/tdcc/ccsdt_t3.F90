!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( h11 ) * t ( p4 p5 h1 h11 )_t * i1 ( h11 p6 h2 h3 )_v 7

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_1_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_perm(sp4,sp5,sp6,sh2,sh1,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h2,h1,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_perm(sp6,sp5,sp4,sh2,sh1,sh3,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h2,h1,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_perm(sp4,sp6,sp5,sh2,sh1,sh3,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h2,h1,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_perm(sp6,sp5,sp4,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_perm(sp4,sp6,sp5,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_1_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h11,sh11
  integer :: spin_t2inp
  integer :: spin_itm_hphh
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh11 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh1,sh11)
     spin_itm_hphh = tdcc_spin_int2x(sh11,sp6,sh2,sh3)
     if(spin_t2inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_1(sh11,sp6,sh2,sh3,itm_hphh)
     call ccsdt_t3_1_2(sh11,sp6,sh2,sh3,itm_hphh)
     call ccsdt_t3_1_3(sh11,sp6,sh2,sh3,itm_hphh)
     call ccsdt_t3_1_4(sh11,sp6,sh2,sh3,itm_hphh)
     call ccsdt_t3_1_5(sh11,sp6,sh2,sh3,itm_hphh)
     call ccsdt_t3_1_6(sh11,sp6,sh2,sh3,itm_hphh)
     call ccsdt_t3_1_7(sh11,sp6,sh2,sh3,itm_hphh)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h11 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t2inp(p4,p5,h1,h11,spin_t2inp) * itm_hphh(h11,p6,h2,h3)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphh)
  end subroutine ccsdt_t3_1_perm
  !--------------------------------------------
end subroutine ccsdt_t3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_1(sh11,sp4,sh1,sh2,i1)

!     i1 ( h11 p4 h1 h2 )_v + = 1 * v ( h11 p4 h1 h2 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h11,p4,h1,h2
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh11,sp4,sh1,sh2)
  do h11 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact * &
          int2x(h11,p4,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_2(sh11,sp4,sh1,sh2,i1)

!     i1 ( h11 p4 h1 h2 )_vt + = 1 * Sum ( h7 ) * t ( p4 h7 )_t * i2 ( h7 h11 h1 h2 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h11,p4,h1,h2
  integer :: h7,sh7
  integer :: spin_t1inp
  integer :: spin_itm_hhhh
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp4,sh7)
     spin_itm_hhhh = tdcc_spin_int2x(sh7,sh11,sh1,sh2)
     if(spin_t1inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_2_1(sh7,sh11,sh1,sh2,itm_hhhh)
     call ccsdt_t3_1_2_2(sh7,sh11,sh1,sh2,itm_hhhh)
     call ccsdt_t3_1_2_3(sh7,sh11,sh1,sh2,itm_hhhh)

     do h11 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h7 = 1,norb1
        i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact * &
             t1inp(p4,h7,spin_t1inp) * itm_hhhh(h7,h11,h1,h2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsdt_t3_1_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_2_1(sh7,sh11,sh1,sh2,i2)

!         i2 ( h7 h11 h1 h2 )_v + = 1 * v ( h7 h11 h1 h2 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh7,sh11,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h7,h11,h1,h2
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh7,sh11,sh1,sh2)
  do h7 = 1,norb1
  do h11 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i2(h7,h11,h1,h2) = i2(h7,h11,h1,h2) + fact * &
          int2x(h7,h11,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_1_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_2_2(sh7,sh11,sh1,sh2,i2)

!         i2 ( h7 h11 h1 h2 )_vt + = -1 * P( 2 ) * Sum ( p8 ) * t ( p8 h1 )_t * i3 ( h7 h11 h2 p8 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh7,sh11,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h7,h11,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i2_perm(:,:,:,:)

  allocate(i2_perm(1:norb1,1:norb1,1:norb1,1:norb1))

  i2_perm(1:norb1,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_1_2_2_perm(sh7,sh11,sh1,sh2,i2_perm)
  fact_p = +1.0d+0
  do h7 = 1,norb1
  do h11 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i2(h7,h11,h1,h2) = i2(h7,h11,h1,h2) + fact_p * i2_perm(h7,h11,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sh7 * sh11 * sh1 * sh2 == 1)) then
     i2_perm(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_2_2_perm(sh7,sh11,sh2,sh1,i2_perm)
  end if
  fact_p = -1.0d+0
  do h7 = 1,norb1
  do h11 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i2(h7,h11,h1,h2) = i2(h7,h11,h1,h2) + fact_p * i2_perm(h7,h11,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i2_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_1_2_2_perm(sh7,sh11,sh1,sh2,i2)

  implicit none
  integer,intent(in) :: sh7,sh11,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)

  integer :: h7,h11,h1,h2
  integer :: p8,sp8
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_itm_hhhp = tdcc_spin_int2x(sh7,sh11,sh2,sp8)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_1_2_2_1(sh7,sh11,sh2,sp8,itm_hhhp)
     call ccsdt_t3_1_2_2_2(sh7,sh11,sh2,sp8,itm_hhhp)

     do h7 = 1,norb1
     do h11 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p8 = norb1+1,nact
        i2(h7,h11,h1,h2) = i2(h7,h11,h1,h2) + fact * &
             t1inp(p8,h1,spin_t1inp) * itm_hhhp(h7,h11,h2,p8)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_t3_1_2_2_perm
  !--------------------------------------------
end subroutine ccsdt_t3_1_2_2
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_2_2_1(sh7,sh11,sh1,sp8,i3)

!             i3 ( h7 h11 h1 p8 )_v + = 1 * v ( h7 h11 h1 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh7,sh11,sh1,sp8
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h7,h11,h1,p8
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh7,sh11,sh1,sp8)
  do h7 = 1,norb1
  do h11 = 1,norb1
  do h1 = 1,norb1
  do p8 = norb1+1,nact
     i3(h7,h11,h1,p8) = i3(h7,h11,h1,p8) + fact * &
          int2x(h7,h11,h1,p8,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_1_2_2_1
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_2_2_2(sh7,sh11,sh1,sp8,i3)

!             i3 ( h7 h11 h1 p8 )_vt + = -1/2 * Sum ( p9 ) * t ( p9 h1 )_t * v ( h7 h11 p8 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh7,sh11,sh1,sp8
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h7,h11,h1,p8
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh1)
     spin_int2x = tdcc_spin_int2x(sh7,sh11,sp8,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h7 = 1,norb1
     do h11 = 1,norb1
     do h1 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
        i3(h7,h11,h1,p8) = i3(h7,h11,h1,p8) + fact * &
             t1inp(p9,h1,spin_t1inp) * int2x(h7,h11,p8,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_1_2_2_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_2_3(sh7,sh11,sh1,sh2,i2)

!         i2 ( h7 h11 h1 h2 )_vt + = 1/2 * Sum ( p8 p9 ) * t ( p8 p9 h1 h2 )_t * v ( h7 h11 p8 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh7,sh11,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h7,h11,h1,h2
  integer :: p8,p9,sp8,sp9
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp8 = 1,2
  do sp9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp8,sp9,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sh7,sh11,sp8,sp9)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h7 = 1,norb1
     do h11 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
        i2(h7,h11,h1,h2) = i2(h7,h11,h1,h2) + fact * &
             t2inp(p8,p9,h1,h2,spin_t2inp) * int2x(h7,h11,p8,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_1_2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_3(sh11,sp4,sh1,sh2,i1)

!     i1 ( h11 p4 h1 h2 )_vt + = -1 * P( 2 ) * Sum ( p9 ) * t ( p9 h1 )_t * i2 ( h11 p4 h2 p9 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h11,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:)

  allocate(i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1))

  i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t3_1_3_perm(sh11,sp4,sh1,sh2,i1_perm)
  fact_p = +1.0d+0
  do h11 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact_p * i1_perm(h11,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sh11 * sp4 * sh1 * sh2 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_3_perm(sh11,sp4,sh2,sh1,i1_perm)
  end if
  fact_p = -1.0d+0
  do h11 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact_p * i1_perm(h11,p4,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_1_3_perm(sh11,sp4,sh1,sh2,i1)

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)

  integer :: h11,p4,h1,h2
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh1)
     spin_itm_hphp = tdcc_spin_int2x(sh11,sp4,sh2,sp9)
     if(spin_t1inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_1_3_1(sh11,sp4,sh2,sp9,itm_hphp)
     call ccsdt_t3_1_3_2(sh11,sp4,sh2,sp9,itm_hphp)

     do h11 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p9 = norb1+1,nact
        i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact * &
             t1inp(p9,h1,spin_t1inp) * itm_hphp(h11,p4,h2,p9)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphp)
  end subroutine ccsdt_t3_1_3_perm
  !--------------------------------------------
end subroutine ccsdt_t3_1_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_3_1(sh11,sp4,sh1,sp9,i2)

!         i2 ( h11 p4 h1 p9 )_v + = 1 * v ( h11 p4 h1 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sp9
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h11,p4,h1,p9
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh11,sp4,sh1,sp9)
  do h11 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do p9 = norb1+1,nact
     i2(h11,p4,h1,p9) = i2(h11,p4,h1,p9) + fact * &
          int2x(h11,p4,h1,p9,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_1_3_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_3_2(sh11,sp4,sh1,sp9,i2)

!         i2 ( h11 p4 h1 p9 )_vt + = 1/2 * Sum ( p8 ) * t ( p8 h1 )_t * v ( h11 p4 p8 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sp9
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h11,p4,h1,p9
  integer :: p8,sp8
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_int2x = tdcc_spin_int2x(sh11,sp4,sp8,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h11 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do p9 = norb1+1,nact
     do p8 = norb1+1,nact
        i2(h11,p4,h1,p9) = i2(h11,p4,h1,p9) + fact * &
             t1inp(p8,h1,spin_t1inp) * int2x(h11,p4,p8,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_1_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_4(sh11,sp4,sh1,sh2,i1)

!     i1 ( h11 p4 h1 h2 )_ft + = -1 * Sum ( p12 ) * t ( p4 p12 h1 h2 )_t * i2 ( h11 p12 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h11,p4,h1,h2
  integer :: p12,sp12
  integer :: spin_t2inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp12 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp12,sh1,sh2)
     spin_itm_hp = tdcc_spin_fock(sh11,sp12)
     if(spin_t2inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_1_4_1(sh11,sp12,itm_hp)
     call ccsdt_t3_1_4_2(sh11,sp12,itm_hp)

     do h11 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p12 = norb1+1,nact
        i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact * &
             t2inp(p4,p12,h1,h2,spin_t2inp) * itm_hp(h11,p12)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_t3_1_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_4_1(sh11,sp12,i2)

!         i2 ( h11 p12 )_f + = 1 * f ( h11 p12 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh11,sp12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h11,p12
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh11,sp12)
  do h11 = 1,norb1
  do p12 = norb1+1,nact
     i2(h11,p12) = i2(h11,p12) + fact * &
          fock(h11,p12,spin_fock)
  end do
  end do
end subroutine ccsdt_t3_1_4_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_4_2(sh11,sp12,i2)

!         i2 ( h11 p12 )_vt + = 1 * Sum ( h10 p9 ) * t ( p9 h10 )_t * v ( h10 h11 p9 p12 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh11,sp12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h11,p12
  integer :: h10,p9,sh10,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh10 = 1,2
  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
     spin_int2x = tdcc_spin_int2x(sh10,sh11,sp9,sp12)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h11 = 1,norb1
     do p12 = norb1+1,nact
     do h10 = 1,norb1
     do p9 = norb1+1,nact
        i2(h11,p12) = i2(h11,p12) + fact * &
             t1inp(p9,h10,spin_t1inp) * int2x(h10,h11,p9,p12,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_1_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_5(sh11,sp4,sh1,sh2,i1)

!     i1 ( h11 p4 h1 h2 )_vt + = 1 * P( 2 ) * Sum ( h9 p8 ) * t ( p4 p8 h1 h9 )_t * i2 ( h9 h11 h2 p8 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h11,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:)

  allocate(i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1))

  i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccsdt_t3_1_5_perm(sh11,sp4,sh1,sh2,i1_perm)
  fact_p = +1.0d+0
  do h11 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact_p * i1_perm(h11,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sh11 * sp4 * sh1 * sh2 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_t3_1_5_perm(sh11,sp4,sh2,sh1,i1_perm)
  end if
  fact_p = -1.0d+0
  do h11 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact_p * i1_perm(h11,p4,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_1_5_perm(sh11,sp4,sh1,sh2,i1)

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)

  integer :: h11,p4,h1,h2
  integer :: h9,p8,sh9,sp8
  integer :: spin_t2inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp8,sh1,sh9)
     spin_itm_hhhp = tdcc_spin_int2x(sh9,sh11,sh2,sp8)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_1_5_1(sh9,sh11,sh2,sp8,itm_hhhp)
     call ccsdt_t3_1_5_2(sh9,sh11,sh2,sp8,itm_hhhp)

     do h11 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h9 = 1,norb1
     do p8 = norb1+1,nact
        i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact * &
             t2inp(p4,p8,h1,h9,spin_t2inp) * itm_hhhp(h9,h11,h2,p8)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_t3_1_5_perm
  !--------------------------------------------
end subroutine ccsdt_t3_1_5
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_5_1(sh9,sh11,sh1,sp8,i2)

!         i2 ( h9 h11 h1 p8 )_v + = 1 * v ( h9 h11 h1 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh9,sh11,sh1,sp8
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h9,h11,h1,p8
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh9,sh11,sh1,sp8)
  do h9 = 1,norb1
  do h11 = 1,norb1
  do h1 = 1,norb1
  do p8 = norb1+1,nact
     i2(h9,h11,h1,p8) = i2(h9,h11,h1,p8) + fact * &
          int2x(h9,h11,h1,p8,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_1_5_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_5_2(sh9,sh11,sh1,sp8,i2)

!         i2 ( h9 h11 h1 p8 )_vt + = -1 * Sum ( p10 ) * t ( p10 h1 )_t * v ( h9 h11 p8 p10 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh9,sh11,sh1,sp8
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h9,h11,h1,p8
  integer :: p10,sp10
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp10 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp10,sh1)
     spin_int2x = tdcc_spin_int2x(sh9,sh11,sp8,sp10)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do h11 = 1,norb1
     do h1 = 1,norb1
     do p8 = norb1+1,nact
     do p10 = norb1+1,nact
        i2(h9,h11,h1,p8) = i2(h9,h11,h1,p8) + fact * &
             t1inp(p10,h1,spin_t1inp) * int2x(h9,h11,p8,p10,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_1_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_6(sh11,sp4,sh1,sh2,i1)

!     i1 ( h11 p4 h1 h2 )_vt + = 1/2 * Sum ( p8 p9 ) * t ( p8 p9 h1 h2 )_t * v ( h11 p4 p8 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h11,p4,h1,h2
  integer :: p8,p9,sp8,sp9
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp8 = 1,2
  do sp9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp8,sp9,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sh11,sp4,sp8,sp9)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h11 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
        i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact * &
             t2inp(p8,p9,h1,h2,spin_t2inp) * int2x(h11,p4,p8,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_1_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_1_7(sh11,sp4,sh1,sh2,i1)

!     i1 ( h11 p4 h1 h2 )_vt + = 1/2 * Sum ( h9 p7 p8 ) * t ( p4 p7 p8 h1 h2 h9 )_t * v ( h9 h11 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh11,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h11,p4,h1,h2
  integer :: h9,p7,p8,sh9,sp7,sp8
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh9 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp7,sp8,sh1,sh2,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sh11,sp7,sp8)
     if(spin_t3inp * spin_int2x == 0) cycle

     do h11 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h9 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h11,p4,h1,h2) = i1(h11,p4,h1,h2) + fact * &
             t3inp(p4,p7,p8,h1,h2,h9,spin_t3inp) * int2x(h9,h11,p7,p8,spin_int2x)
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
end subroutine ccsdt_t3_1_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_2(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( p12 ) * t ( p4 p12 h1 h2 )_t * i1 ( p5 p6 h3 p12 )_v 5

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_2_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_2_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_2_perm(sp4,sp5,sp6,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_2_perm(sp5,sp4,sp6,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p5,p4,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_2_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_2_perm(sp5,sp4,sp6,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p5,p4,p6,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_2_perm(sp6,sp5,sp4,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_2_perm(sp5,sp4,sp6,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p5,p4,p6,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_2_perm(sp6,sp5,sp4,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_2_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: p12,sp12
  integer :: spin_t2inp
  integer :: spin_itm_pphp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp12 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp12,sh1,sh2)
     spin_itm_pphp = tdcc_spin_int2x(sp5,sp6,sh3,sp12)
     if(spin_t2inp * spin_itm_pphp == 0) cycle

     itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_2_1(sp5,sp6,sh3,sp12,itm_pphp)
     call ccsdt_t3_2_2(sp5,sp6,sh3,sp12,itm_pphp)
     call ccsdt_t3_2_3(sp5,sp6,sh3,sp12,itm_pphp)
     call ccsdt_t3_2_4(sp5,sp6,sh3,sp12,itm_pphp)
     call ccsdt_t3_2_5(sp5,sp6,sh3,sp12,itm_pphp)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p12 = norb1+1,nact
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t2inp(p4,p12,h1,h2,spin_t2inp) * itm_pphp(p5,p6,h3,p12)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pphp)
  end subroutine ccsdt_t3_2_perm
  !--------------------------------------------
end subroutine ccsdt_t3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_2_1(sp4,sp5,sh1,sp12,i1)

!     i1 ( p4 p5 h1 p12 )_v + = 1 * v ( p4 p5 h1 p12 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp12
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p4,p5,h1,p12
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sp4,sp5,sh1,sp12)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do p12 = norb1+1,nact
     i1(p4,p5,h1,p12) = i1(p4,p5,h1,p12) + fact * &
          int2x(p4,p5,h1,p12,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_2_2(sp4,sp5,sh1,sp12,i1)

!     i1 ( p4 p5 h1 p12 )_vt + = 1 * Sum ( p8 ) * t ( p8 h1 )_t * v ( p4 p5 p8 p12 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp12
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p4,p5,h1,p12
  integer :: p8,sp8
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_int2x = tdcc_spin_int2x(sp4,sp5,sp8,sp12)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do p12 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(p4,p5,h1,p12) = i1(p4,p5,h1,p12) + fact * &
             t1inp(p8,h1,spin_t1inp) * int2x(p4,p5,p8,p12,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_2_3(sp4,sp5,sh1,sp12,i1)

!     i1 ( p4 p5 h1 p12 )_vt + = 1/2 * Sum ( h8 h9 ) * t ( p4 p5 h8 h9 )_t * i2 ( h8 h9 h1 p12 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp12
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p4,p5,h1,p12
  integer :: h8,h9,sh8,sh9
  integer :: spin_t2inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh8 = 1,2
  do sh9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh8,sh9)
     spin_itm_hhhp = tdcc_spin_int2x(sh8,sh9,sh1,sp12)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_2_3_1(sh8,sh9,sh1,sp12,itm_hhhp)
     call ccsdt_t3_2_3_2(sh8,sh9,sh1,sp12,itm_hhhp)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do p12 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
        i1(p4,p5,h1,p12) = i1(p4,p5,h1,p12) + fact * &
             t2inp(p4,p5,h8,h9,spin_t2inp) * itm_hhhp(h8,h9,h1,p12)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_t3_2_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_2_3_1(sh8,sh9,sh1,sp12,i2)

!         i2 ( h8 h9 h1 p12 )_v + = 1 * v ( h8 h9 h1 p12 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh8,sh9,sh1,sp12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h8,h9,h1,p12
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh8,sh9,sh1,sp12)
  do h8 = 1,norb1
  do h9 = 1,norb1
  do h1 = 1,norb1
  do p12 = norb1+1,nact
     i2(h8,h9,h1,p12) = i2(h8,h9,h1,p12) + fact * &
          int2x(h8,h9,h1,p12,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_2_3_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_2_3_2(sh8,sh9,sh1,sp12,i2)

!         i2 ( h8 h9 h1 p12 )_vt + = 1 * Sum ( p10 ) * t ( p10 h1 )_t * v ( h8 h9 p10 p12 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh8,sh9,sh1,sp12
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h8,h9,h1,p12
  integer :: p10,sp10
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp10 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp10,sh1)
     spin_int2x = tdcc_spin_int2x(sh8,sh9,sp10,sp12)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h8 = 1,norb1
     do h9 = 1,norb1
     do h1 = 1,norb1
     do p12 = norb1+1,nact
     do p10 = norb1+1,nact
        i2(h8,h9,h1,p12) = i2(h8,h9,h1,p12) + fact * &
             t1inp(p10,h1,spin_t1inp) * int2x(h8,h9,p10,p12,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_2_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_2_4(sp4,sp5,sh1,sp12,i1)

!     i1 ( p4 p5 h1 p12 )_vt + = 1 * P( 2 ) * Sum ( h9 p8 ) * t ( p4 p8 h1 h9 )_t * v ( h9 p5 p8 p12 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp12
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p4,p5,h1,p12
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:)

  allocate(i1_perm((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact))

  i1_perm((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
  call ccsdt_t3_2_4_perm(sp4,sp5,sh1,sp12,i1_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do p12 = norb1+1,nact
     i1(p4,p5,h1,p12) = i1(p4,p5,h1,p12) + fact_p * i1_perm(p4,p5,h1,p12)
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sh1 * sp12 == 1)) then
     i1_perm((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_2_4_perm(sp5,sp4,sh1,sp12,i1_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do p12 = norb1+1,nact
     i1(p4,p5,h1,p12) = i1(p4,p5,h1,p12) + fact_p * i1_perm(p5,p4,h1,p12)
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_2_4_perm(sp4,sp5,sh1,sp12,i1)

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp12
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer :: p4,p5,h1,p12
  integer :: h9,p8,sh9,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh9 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp8,sh1,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sp5,sp8,sp12)
     if(spin_t2inp * spin_int2x == 0) cycle

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do p12 = norb1+1,nact
     do h9 = 1,norb1
     do p8 = norb1+1,nact
        i1(p4,p5,h1,p12) = i1(p4,p5,h1,p12) + fact * &
             t2inp(p4,p8,h1,h9,spin_t2inp) * int2x(h9,p5,p8,p12,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end subroutine ccsdt_t3_2_4_perm
  !--------------------------------------------
end subroutine ccsdt_t3_2_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_2_5(sp4,sp5,sh1,sp12,i1)

!     i1 ( p4 p5 h1 p12 )_vt + = 1/2 * Sum ( h8 h9 p7 ) * t ( p4 p5 p7 h1 h8 h9 )_t * v ( h8 h9 p7 p12 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp12
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p4,p5,h1,p12
  integer :: h8,h9,p7,sh8,sh9,sp7
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sh9 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp7,sh1,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh8,sh9,sp7,sp12)
     if(spin_t3inp * spin_int2x == 0) cycle

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do p12 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p7 = norb1+1,nact
        i1(p4,p5,h1,p12) = i1(p4,p5,h1,p12) + fact * &
             t3inp(p4,p5,p7,h1,h8,h9,spin_t3inp) * int2x(h8,h9,p7,p12,spin_int2x)
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
end subroutine ccsdt_t3_2_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_3(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_tf + = -1 * P( 3 ) * Sum ( h11 ) * t ( p4 p5 p6 h1 h2 h11 )_t * i1 ( h11 h3 )_f 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_3_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_3_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_3_perm(sp4,sp5,sp6,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_3_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h11,sh11
  integer :: spin_t3inp
  integer :: spin_itm_hh
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh11 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp6,sh1,sh2,sh11)
     spin_itm_hh = tdcc_spin_fock(sh11,sh3)
     if(spin_t3inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsdt_t3_3_1(sh11,sh3,itm_hh)
     call ccsdt_t3_3_2(sh11,sh3,itm_hh)
     call ccsdt_t3_3_3(sh11,sh3,itm_hh)
     call ccsdt_t3_3_4(sh11,sh3,itm_hh)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h11 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p6,h1,h2,h11,spin_t3inp) * itm_hh(h11,h3)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccsdt_t3_3_perm
  !--------------------------------------------
end subroutine ccsdt_t3_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_3_1(sh11,sh1,i1)

!     i1 ( h11 h1 )_f + = 1 * f ( h11 h1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh11,sh1
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1)
  integer :: h11,h1
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh11,sh1)
  do h11 = 1,norb1
  do h1 = 1,norb1
     i1(h11,h1) = i1(h11,h1) + fact * &
          fock(h11,h1,spin_fock)
  end do
  end do
end subroutine ccsdt_t3_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_3_2(sh11,sh1,i1)

!     i1 ( h11 h1 )_ft + = 1 * Sum ( p10 ) * t ( p10 h1 )_t * i2 ( h11 p10 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh11,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h11,h1
  integer :: p10,sp10
  integer :: spin_t1inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp10 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp10,sh1)
     spin_itm_hp = tdcc_spin_fock(sh11,sp10)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_3_2_1(sh11,sp10,itm_hp)
     call ccsdt_t3_3_2_2(sh11,sp10,itm_hp)

     do h11 = 1,norb1
     do h1 = 1,norb1
     do p10 = norb1+1,nact
        i1(h11,h1) = i1(h11,h1) + fact * &
             t1inp(p10,h1,spin_t1inp) * itm_hp(h11,p10)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_t3_3_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_3_2_1(sh11,sp10,i2)

!         i2 ( h11 p10 )_f + = 1 * f ( h11 p10 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh11,sp10
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h11,p10
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh11,sp10)
  do h11 = 1,norb1
  do p10 = norb1+1,nact
     i2(h11,p10) = i2(h11,p10) + fact * &
          fock(h11,p10,spin_fock)
  end do
  end do
end subroutine ccsdt_t3_3_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_3_2_2(sh11,sp10,i2)

!         i2 ( h11 p10 )_vt + = 1 * Sum ( h9 p8 ) * t ( p8 h9 )_t * v ( h9 h11 p8 p10 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh11,sp10
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h11,p10
  integer :: h9,p8,sh9,sp8
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh9 = 1,2
  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sh11,sp8,sp10)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h11 = 1,norb1
     do p10 = norb1+1,nact
     do h9 = 1,norb1
     do p8 = norb1+1,nact
        i2(h11,p10) = i2(h11,p10) + fact * &
             t1inp(p8,h9,spin_t1inp) * int2x(h9,h11,p8,p10,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_3_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_3_3(sh11,sh1,i1)

!     i1 ( h11 h1 )_vt + = -1 * Sum ( h9 p8 ) * t ( p8 h9 )_t * v ( h9 h11 h1 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh11,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h11,h1
  integer :: h9,p8,sh9,sp8
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh9 = 1,2
  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sh11,sh1,sp8)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h11 = 1,norb1
     do h1 = 1,norb1
     do h9 = 1,norb1
     do p8 = norb1+1,nact
        i1(h11,h1) = i1(h11,h1) + fact * &
             t1inp(p8,h9,spin_t1inp) * int2x(h9,h11,h1,p8,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_3_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_3_4(sh11,sh1,i1)

!     i1 ( h11 h1 )_vt + = -1/2 * Sum ( h10 p8 p9 ) * t ( p8 p9 h1 h10 )_t * v ( h10 h11 p8 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh11,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h11,h1
  integer :: h10,p8,p9,sh10,sp8,sp9
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh10 = 1,2
  do sp8 = 1,2
  do sp9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp8,sp9,sh1,sh10)
     spin_int2x = tdcc_spin_int2x(sh10,sh11,sp8,sp9)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h11 = 1,norb1
     do h1 = 1,norb1
     do h10 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
        i1(h11,h1) = i1(h11,h1) + fact * &
             t2inp(p8,p9,h1,h10,spin_t2inp) * int2x(h10,h11,p8,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_t3_3_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_4(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_tf + = 1 * P( 3 ) * Sum ( p7 ) * t ( p4 p5 p7 h1 h2 h3 )_t * i1 ( p6 p7 )_f 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_4_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_4_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_4_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_4_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: p7,sp7
  integer :: spin_t3inp
  integer :: spin_itm_pp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp7,sh1,sh2,sh3)
     spin_itm_pp = tdcc_spin_fock(sp6,sp7)
     if(spin_t3inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_t3_4_1(sp6,sp7,itm_pp)
     call ccsdt_t3_4_2(sp6,sp7,itm_pp)
     call ccsdt_t3_4_3(sp6,sp7,itm_pp)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p7 = norb1+1,nact
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p7,h1,h2,h3,spin_t3inp) * itm_pp(p6,p7)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccsdt_t3_4_perm
  !--------------------------------------------
end subroutine ccsdt_t3_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_4_1(sp4,sp7,i1)

!     i1 ( p4 p7 )_f + = 1 * f ( p4 p7 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sp4,sp7
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p4,p7
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp4,sp7)
  do p4 = norb1+1,nact
  do p7 = norb1+1,nact
     i1(p4,p7) = i1(p4,p7) + fact * &
          fock(p4,p7,spin_fock)
  end do
  end do
end subroutine ccsdt_t3_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_4_2(sp4,sp7,i1)

!     i1 ( p4 p7 )_vt + = -1 * Sum ( h9 p8 ) * t ( p8 h9 )_t * v ( h9 p4 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp7
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact)
  integer :: p4,p7
  integer :: h9,p8,sh9,sp8
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh9 = 1,2
  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sp4,sp7,sp8)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p4 = norb1+1,nact
     do p7 = norb1+1,nact
     do h9 = 1,norb1
     do p8 = norb1+1,nact
        i1(p4,p7) = i1(p4,p7) + fact * &
             t1inp(p8,h9,spin_t1inp) * int2x(h9,p4,p7,p8,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_4_3(sp4,sp7,i1)

!     i1 ( p4 p7 )_vt + = -1/2 * Sum ( h9 h10 p8 ) * t ( p4 p8 h9 h10 )_t * v ( h9 h10 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp7
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact)
  integer :: p4,p7
  integer :: h9,h10,p8,sh9,sh10,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh9 = 1,2
  do sh10 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp8,sh9,sh10)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp7,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     do p4 = norb1+1,nact
     do p7 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p8 = norb1+1,nact
        i1(p4,p7) = i1(p4,p7) + fact * &
             t2inp(p4,p8,h9,h10,spin_t2inp) * int2x(h9,h10,p7,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_t3_4_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_5(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = 1/2 * P( 3 ) * Sum ( h11 h12 ) * t ( p4 p5 p6 h1 h11 h12 )_t * i1 ( h11 h12 h2 h3 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_5_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_5_perm(sp4,sp5,sp6,sh2,sh1,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h2,h1,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_5_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_5_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h11,h12,sh11,sh12
  integer :: spin_t3inp
  integer :: spin_itm_hhhh
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh11 = 1,2
  do sh12 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp6,sh1,sh11,sh12)
     spin_itm_hhhh = tdcc_spin_int2x(sh11,sh12,sh2,sh3)
     if(spin_t3inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_5_1(sh11,sh12,sh2,sh3,itm_hhhh)
     call ccsdt_t3_5_2(sh11,sh12,sh2,sh3,itm_hhhh)
     call ccsdt_t3_5_3(sh11,sh12,sh2,sh3,itm_hhhh)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h11 = 1,norb1
     do h12 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p6,h1,h11,h12,spin_t3inp) * itm_hhhh(h11,h12,h2,h3)
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
  deallocate(itm_hhhh)
  end subroutine ccsdt_t3_5_perm
  !--------------------------------------------
end subroutine ccsdt_t3_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_5_1(sh11,sh12,sh1,sh2,i1)

!     i1 ( h11 h12 h1 h2 )_v + = 1 * v ( h11 h12 h1 h2 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh11,sh12,sh1,sh2
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h11,h12,h1,h2
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh11,sh12,sh1,sh2)
  do h11 = 1,norb1
  do h12 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h11,h12,h1,h2) = i1(h11,h12,h1,h2) + fact * &
          int2x(h11,h12,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_5_2(sh11,sh12,sh1,sh2,i1)

!     i1 ( h11 h12 h1 h2 )_vt + = -1 * P( 2 ) * Sum ( p10 ) * t ( p10 h1 )_t * i2 ( h11 h12 h2 p10 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh11,sh12,sh1,sh2
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h11,h12,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:)

  allocate(i1_perm(1:norb1,1:norb1,1:norb1,1:norb1))

  i1_perm(1:norb1,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_5_2_perm(sh11,sh12,sh1,sh2,i1_perm)
  fact_p = +1.0d+0
  do h11 = 1,norb1
  do h12 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h11,h12,h1,h2) = i1(h11,h12,h1,h2) + fact_p * i1_perm(h11,h12,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sh11 * sh12 * sh1 * sh2 == 1)) then
     i1_perm(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_5_2_perm(sh11,sh12,sh2,sh1,i1_perm)
  end if
  fact_p = -1.0d+0
  do h11 = 1,norb1
  do h12 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h11,h12,h1,h2) = i1(h11,h12,h1,h2) + fact_p * i1_perm(h11,h12,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_5_2_perm(sh11,sh12,sh1,sh2,i1)

  implicit none
  integer,intent(in) :: sh11,sh12,sh1,sh2
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1)

  integer :: h11,h12,h1,h2
  integer :: p10,sp10
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp10 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp10,sh1)
     spin_itm_hhhp = tdcc_spin_int2x(sh11,sh12,sh2,sp10)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_5_2_1(sh11,sh12,sh2,sp10,itm_hhhp)
     call ccsdt_t3_5_2_2(sh11,sh12,sh2,sp10,itm_hhhp)

     do h11 = 1,norb1
     do h12 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p10 = norb1+1,nact
        i1(h11,h12,h1,h2) = i1(h11,h12,h1,h2) + fact * &
             t1inp(p10,h1,spin_t1inp) * itm_hhhp(h11,h12,h2,p10)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_t3_5_2_perm
  !--------------------------------------------
end subroutine ccsdt_t3_5_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_5_2_1(sh11,sh12,sh1,sp10,i2)

!         i2 ( h11 h12 h1 p10 )_v + = 1 * v ( h11 h12 h1 p10 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh11,sh12,sh1,sp10
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h11,h12,h1,p10
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh11,sh12,sh1,sp10)
  do h11 = 1,norb1
  do h12 = 1,norb1
  do h1 = 1,norb1
  do p10 = norb1+1,nact
     i2(h11,h12,h1,p10) = i2(h11,h12,h1,p10) + fact * &
          int2x(h11,h12,h1,p10,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_5_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_5_2_2(sh11,sh12,sh1,sp10,i2)

!         i2 ( h11 h12 h1 p10 )_vt + = 1/2 * Sum ( p8 ) * t ( p8 h1 )_t * v ( h11 h12 p8 p10 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh11,sh12,sh1,sp10
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h11,h12,h1,p10
  integer :: p8,sp8
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_int2x = tdcc_spin_int2x(sh11,sh12,sp8,sp10)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h11 = 1,norb1
     do h12 = 1,norb1
     do h1 = 1,norb1
     do p10 = norb1+1,nact
     do p8 = norb1+1,nact
        i2(h11,h12,h1,p10) = i2(h11,h12,h1,p10) + fact * &
             t1inp(p8,h1,spin_t1inp) * int2x(h11,h12,p8,p10,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_5_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_5_3(sh11,sh12,sh1,sh2,i1)

!     i1 ( h11 h12 h1 h2 )_vt + = 1/2 * Sum ( p9 p10 ) * t ( p9 p10 h1 h2 )_t * v ( h11 h12 p9 p10 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh11,sh12,sh1,sh2
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h11,h12,h1,h2
  integer :: p9,p10,sp9,sp10
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp9 = 1,2
  do sp10 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp9,sp10,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sh11,sh12,sp9,sp10)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h11 = 1,norb1
     do h12 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
        i1(h11,h12,h1,h2) = i1(h11,h12,h1,h2) + fact * &
             t2inp(p9,p10,h1,h2,spin_t2inp) * int2x(h11,h12,p9,p10,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_5_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_6(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( h8 p7 ) * t ( p4 p5 p7 h1 h2 h8 )_t * i1 ( h8 p6 h3 p7 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_6_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_6_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_6_perm(sp4,sp5,sp6,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_6_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_6_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_6_perm(sp6,sp5,sp4,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_6_perm(sp4,sp6,sp5,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_6_perm(sp6,sp5,sp4,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_6_perm(sp4,sp6,sp5,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_6_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h8,p7,sh8,sp7
  integer :: spin_t3inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh8 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp7,sh1,sh2,sh8)
     spin_itm_hphp = tdcc_spin_int2x(sh8,sp6,sh3,sp7)
     if(spin_t3inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_6_1(sh8,sp6,sh3,sp7,itm_hphp)
     call ccsdt_t3_6_2(sh8,sp6,sh3,sp7,itm_hphp)
     call ccsdt_t3_6_3(sh8,sp6,sh3,sp7,itm_hphp)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p7,h1,h2,h8,spin_t3inp) * itm_hphp(h8,p6,h3,p7)
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
  deallocate(itm_hphp)
  end subroutine ccsdt_t3_6_perm
  !--------------------------------------------
end subroutine ccsdt_t3_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_6_1(sh8,sp4,sh1,sp7,i1)

!     i1 ( h8 p4 h1 p7 )_v + = 1 * v ( h8 p4 h1 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh8,sp4,sh1,sp7
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h8,p4,h1,p7
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh8,sp4,sh1,sp7)
  do h8 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do p7 = norb1+1,nact
     i1(h8,p4,h1,p7) = i1(h8,p4,h1,p7) + fact * &
          int2x(h8,p4,h1,p7,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_6_2(sh8,sp4,sh1,sp7,i1)

!     i1 ( h8 p4 h1 p7 )_vt + = -1 * Sum ( p9 ) * t ( p9 h1 )_t * v ( h8 p4 p7 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh8,sp4,sh1,sp7
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h8,p4,h1,p7
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh1)
     spin_int2x = tdcc_spin_int2x(sh8,sp4,sp7,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h8 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do p7 = norb1+1,nact
     do p9 = norb1+1,nact
        i1(h8,p4,h1,p7) = i1(h8,p4,h1,p7) + fact * &
             t1inp(p9,h1,spin_t1inp) * int2x(h8,p4,p7,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_6_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_6_3(sh8,sp4,sh1,sp7,i1)

!     i1 ( h8 p4 h1 p7 )_vt + = -1 * Sum ( h10 p9 ) * t ( p4 p9 h1 h10 )_t * v ( h8 h10 p7 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh8,sp4,sh1,sp7
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h8,p4,h1,p7
  integer :: h10,p9,sh10,sp9
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh10 = 1,2
  do sp9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp9,sh1,sh10)
     spin_int2x = tdcc_spin_int2x(sh8,sh10,sp7,sp9)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h8 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do p7 = norb1+1,nact
     do h10 = 1,norb1
     do p9 = norb1+1,nact
        i1(h8,p4,h1,p7) = i1(h8,p4,h1,p7) + fact * &
             t2inp(p4,p9,h1,h10,spin_t2inp) * int2x(h8,h10,p7,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_6_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_7(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = 1/2 * P( 3 ) * Sum ( p7 p8 ) * t ( p4 p7 p8 h1 h2 h3 )_t * v ( p5 p6 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_7_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_7_perm(sp5,sp4,sp6,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p5,p4,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_7_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_7_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: p7,p8,sp7,sp8
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp7,sp8,sh1,sh2,sh3)
     spin_int2x = tdcc_spin_int2x(sp5,sp6,sp7,sp8)
     if(spin_t3inp * spin_int2x == 0) cycle

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p7,p8,h1,h2,h3,spin_t3inp) * int2x(p5,p6,p7,p8,spin_int2x)
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
  end subroutine ccsdt_t3_7_perm
  !--------------------------------------------
end subroutine ccsdt_t3_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vtt + = 6 * P( 3 ) * Sum ( h12 ) * t ( p4 h12 )_t * i1 ( h12 p5 p6 h1 h2 h3 )_vt 5

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_8_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_perm(sp5,sp4,sp6,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p5,p4,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_8_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h12,sh12
  integer :: spin_t1inp
  integer :: spin_itm_hpphhh
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hpphhh(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 6.0d+0 * runit

  allocate(itm_hpphhh(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))
  do sh12 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp4,sh12)
     spin_itm_hpphhh = tdcc_spin_dummy3(sh12,sp5,sp6,sh1,sh2,sh3)
     if(spin_t1inp * spin_itm_hpphhh == 0) cycle

     itm_hpphhh(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_1(sh12,sp5,sp6,sh1,sh2,sh3,itm_hpphhh)
     call ccsdt_t3_8_2(sh12,sp5,sp6,sh1,sh2,sh3,itm_hpphhh)
     call ccsdt_t3_8_3(sh12,sp5,sp6,sh1,sh2,sh3,itm_hpphhh)
     call ccsdt_t3_8_4(sh12,sp5,sp6,sh1,sh2,sh3,itm_hpphhh)
     call ccsdt_t3_8_5(sh12,sp5,sp6,sh1,sh2,sh3,itm_hpphhh)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h12 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t1inp(p4,h12,spin_t1inp) * itm_hpphhh(h12,p5,p6,h1,h2,h3)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hpphhh)
  end subroutine ccsdt_t3_8_perm
  !--------------------------------------------
end subroutine ccsdt_t3_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_1(sh12,sp4,sp5,sh1,sh2,sh3,i1)

!     i1 ( h12 p4 p5 h1 h2 h3 )_vt + = -1/6 * P( 6 ) * Sum ( p11 ) * t ( p4 p11 h1 h2 )_t * i2 ( h12 p5 h3 p11 )_v 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh12,sp4,sp5,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h12,p4,p5,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:,:,:)

  allocate(i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_8_1_perm(sh12,sp4,sp5,sh1,sh2,sh3,i1_perm)
  fact_p = +1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p4,p5,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh12 * sp4 * sp5 * sh1 * sh2 * sh3 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_1_perm(sh12,sp4,sp5,sh3,sh2,sh1,i1_perm)
  end if
  fact_p = -1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p4,p5,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh12 * sp4 * sp5 * sh1 * sh2 * sh3 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_1_perm(sh12,sp4,sp5,sh1,sh3,sh2,i1_perm)
  end if
  fact_p = -1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p4,p5,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh12 * sp4 * sp5 * sh1 * sh2 * sh3 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_1_perm(sh12,sp5,sp4,sh1,sh2,sh3,i1_perm)
  end if
  fact_p = -1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh12 * sp4 * sp5 * sh1 * sh2 * sh3 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_1_perm(sh12,sp5,sp4,sh3,sh2,sh1,i1_perm)
  end if
  fact_p = +1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p5,p4,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh12 * sp4 * sp5 * sh1 * sh2 * sh3 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_1_perm(sh12,sp5,sp4,sh1,sh3,sh2,i1_perm)
  end if
  fact_p = +1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p5,p4,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_8_1_perm(sh12,sp4,sp5,sh1,sh2,sh3,i1)

  implicit none
  integer,intent(in) :: sh12,sp4,sp5,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: h12,p4,p5,h1,h2,h3
  integer :: p11,sp11
  integer :: spin_t2inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 6.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp11 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp11,sh1,sh2)
     spin_itm_hphp = tdcc_spin_int2x(sh12,sp5,sh3,sp11)
     if(spin_t2inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_8_1_1(sh12,sp5,sh3,sp11,itm_hphp)
     call ccsdt_t3_8_1_2(sh12,sp5,sh3,sp11,itm_hphp)
     call ccsdt_t3_8_1_3(sh12,sp5,sh3,sp11,itm_hphp)
     call ccsdt_t3_8_1_4(sh12,sp5,sh3,sp11,itm_hphp)

     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p11 = norb1+1,nact
        i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact * &
             t2inp(p4,p11,h1,h2,spin_t2inp) * itm_hphp(h12,p5,h3,p11)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphp)
  end subroutine ccsdt_t3_8_1_perm
  !--------------------------------------------
end subroutine ccsdt_t3_8_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_1_1(sh12,sp4,sh1,sp11,i2)

!         i2 ( h12 p4 h1 p11 )_v + = 1 * v ( h12 p4 h1 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h12,p4,h1,p11
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh12,sp4,sh1,sp11)
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do p11 = norb1+1,nact
     i2(h12,p4,h1,p11) = i2(h12,p4,h1,p11) + fact * &
          int2x(h12,p4,h1,p11,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_8_1_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_1_2(sh12,sp4,sh1,sp11,i2)

!         i2 ( h12 p4 h1 p11 )_vt + = 1/2 * Sum ( h8 ) * t ( p4 h8 )_t * i3 ( h8 h12 h1 p11 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h12,p4,h1,p11
  integer :: h8,sh8
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp4,sh8)
     spin_itm_hhhp = tdcc_spin_int2x(sh8,sh12,sh1,sp11)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_8_1_2_1(sh8,sh12,sh1,sp11,itm_hhhp)
     call ccsdt_t3_8_1_2_2(sh8,sh12,sh1,sp11,itm_hhhp)

     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do p11 = norb1+1,nact
     do h8 = 1,norb1
        i2(h12,p4,h1,p11) = i2(h12,p4,h1,p11) + fact * &
             t1inp(p4,h8,spin_t1inp) * itm_hhhp(h8,h12,h1,p11)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_t3_8_1_2
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_1_2_1(sh8,sh12,sh1,sp11,i3)

!             i3 ( h8 h12 h1 p11 )_v + = 1 * v ( h8 h12 h1 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh8,sh12,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h8,h12,h1,p11
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh8,sh12,sh1,sp11)
  do h8 = 1,norb1
  do h12 = 1,norb1
  do h1 = 1,norb1
  do p11 = norb1+1,nact
     i3(h8,h12,h1,p11) = i3(h8,h12,h1,p11) + fact * &
          int2x(h8,h12,h1,p11,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_8_1_2_1
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_1_2_2(sh8,sh12,sh1,sp11,i3)

!             i3 ( h8 h12 h1 p11 )_vt + = 1 * Sum ( p9 ) * t ( p9 h1 )_t * v ( h8 h12 p9 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh8,sh12,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h8,h12,h1,p11
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh1)
     spin_int2x = tdcc_spin_int2x(sh8,sh12,sp9,sp11)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h8 = 1,norb1
     do h12 = 1,norb1
     do h1 = 1,norb1
     do p11 = norb1+1,nact
     do p9 = norb1+1,nact
        i3(h8,h12,h1,p11) = i3(h8,h12,h1,p11) + fact * &
             t1inp(p9,h1,spin_t1inp) * int2x(h8,h12,p9,p11,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_8_1_2_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_1_3(sh12,sp4,sh1,sp11,i2)

!         i2 ( h12 p4 h1 p11 )_vt + = 1 * Sum ( p8 ) * t ( p8 h1 )_t * v ( h12 p4 p8 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h12,p4,h1,p11
  integer :: p8,sp8
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_int2x = tdcc_spin_int2x(sh12,sp4,sp8,sp11)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do p11 = norb1+1,nact
     do p8 = norb1+1,nact
        i2(h12,p4,h1,p11) = i2(h12,p4,h1,p11) + fact * &
             t1inp(p8,h1,spin_t1inp) * int2x(h12,p4,p8,p11,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_8_1_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_1_4(sh12,sp4,sh1,sp11,i2)

!         i2 ( h12 p4 h1 p11 )_vt + = -1 * Sum ( h9 p8 ) * t ( p4 p8 h1 h9 )_t * v ( h9 h12 p8 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h12,p4,h1,p11
  integer :: h9,p8,sh9,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh9 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp8,sh1,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sh12,sp8,sp11)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do p11 = norb1+1,nact
     do h9 = 1,norb1
     do p8 = norb1+1,nact
        i2(h12,p4,h1,p11) = i2(h12,p4,h1,p11) + fact * &
             t2inp(p4,p8,h1,h9,spin_t2inp) * int2x(h9,h12,p8,p11,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_8_1_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_2(sh12,sp4,sp5,sh1,sh2,sh3,i1)

!     i1 ( h12 p4 p5 h1 h2 h3 )_ft + = -1/6 * Sum ( p7 ) * t ( p4 p5 p7 h1 h2 h3 )_t * i2 ( h12 p7 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sh12,sp4,sp5,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h12,p4,p5,h1,h2,h3
  integer :: p7,sp7
  integer :: spin_t3inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 6.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp7,sh1,sh2,sh3)
     spin_itm_hp = tdcc_spin_fock(sh12,sp7)
     if(spin_t3inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_8_2_1(sh12,sp7,itm_hp)
     call ccsdt_t3_8_2_2(sh12,sp7,itm_hp)

     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p7 = norb1+1,nact
        i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact * &
             t3inp(p4,p5,p7,h1,h2,h3,spin_t3inp) * itm_hp(h12,p7)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_t3_8_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_2_1(sh12,sp7,i2)

!         i2 ( h12 p7 )_f + = 1 * f ( h12 p7 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh12,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h12,p7
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh12,sp7)
  do h12 = 1,norb1
  do p7 = norb1+1,nact
     i2(h12,p7) = i2(h12,p7) + fact * &
          fock(h12,p7,spin_fock)
  end do
  end do
end subroutine ccsdt_t3_8_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_2_2(sh12,sp7,i2)

!         i2 ( h12 p7 )_vt + = -1 * Sum ( h9 p8 ) * t ( p8 h9 )_t * v ( h9 h12 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh12,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h12,p7
  integer :: h9,p8,sh9,sp8
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh9 = 1,2
  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sh12,sp7,sp8)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h12 = 1,norb1
     do p7 = norb1+1,nact
     do h9 = 1,norb1
     do p8 = norb1+1,nact
        i2(h12,p7) = i2(h12,p7) + fact * &
             t1inp(p8,h9,spin_t1inp) * int2x(h9,h12,p7,p8,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_t3_8_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_3(sh12,sp4,sp5,sh1,sh2,sh3,i1)

!     i1 ( h12 p4 p5 h1 h2 h3 )_vt + = 1/6 * P( 3 ) * Sum ( h7 p11 ) * t ( p4 p5 p11 h1 h2 h7 )_t * i2 ( h7 h12 h3 p11 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sh12,sp4,sp5,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h12,p4,p5,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:,:,:)

  allocate(i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_8_3_perm(sh12,sp4,sp5,sh1,sh2,sh3,i1_perm)
  fact_p = +1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p4,p5,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh12 * sp4 * sp5 * sh1 * sh2 * sh3 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_3_perm(sh12,sp4,sp5,sh3,sh2,sh1,i1_perm)
  end if
  fact_p = -1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p4,p5,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh12 * sp4 * sp5 * sh1 * sh2 * sh3 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_3_perm(sh12,sp4,sp5,sh1,sh3,sh2,i1_perm)
  end if
  fact_p = -1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p4,p5,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_8_3_perm(sh12,sp4,sp5,sh1,sh2,sh3,i1)

  implicit none
  integer,intent(in) :: sh12,sp4,sp5,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: h12,p4,p5,h1,h2,h3
  integer :: h7,p11,sh7,sp11
  integer :: spin_t3inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 6.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh7 = 1,2
  do sp11 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp11,sh1,sh2,sh7)
     spin_itm_hhhp = tdcc_spin_int2x(sh7,sh12,sh3,sp11)
     if(spin_t3inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_t3_8_3_1(sh7,sh12,sh3,sp11,itm_hhhp)
     call ccsdt_t3_8_3_2(sh7,sh12,sh3,sp11,itm_hhhp)

     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h7 = 1,norb1
     do p11 = norb1+1,nact
        i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact * &
             t3inp(p4,p5,p11,h1,h2,h7,spin_t3inp) * itm_hhhp(h7,h12,h3,p11)
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
  end subroutine ccsdt_t3_8_3_perm
  !--------------------------------------------
end subroutine ccsdt_t3_8_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_3_1(sh7,sh12,sh1,sp11,i2)

!         i2 ( h7 h12 h1 p11 )_v + = 1 * v ( h7 h12 h1 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh7,sh12,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h7,h12,h1,p11
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh7,sh12,sh1,sp11)
  do h7 = 1,norb1
  do h12 = 1,norb1
  do h1 = 1,norb1
  do p11 = norb1+1,nact
     i2(h7,h12,h1,p11) = i2(h7,h12,h1,p11) + fact * &
          int2x(h7,h12,h1,p11,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_t3_8_3_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_3_2(sh7,sh12,sh1,sp11,i2)

!         i2 ( h7 h12 h1 p11 )_vt + = 1 * Sum ( p8 ) * t ( p8 h1 )_t * v ( h7 h12 p8 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh7,sh12,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h7,h12,h1,p11
  integer :: p8,sp8
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh1)
     spin_int2x = tdcc_spin_int2x(sh7,sh12,sp8,sp11)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h7 = 1,norb1
     do h12 = 1,norb1
     do h1 = 1,norb1
     do p11 = norb1+1,nact
     do p8 = norb1+1,nact
        i2(h7,h12,h1,p11) = i2(h7,h12,h1,p11) + fact * &
             t1inp(p8,h1,spin_t1inp) * int2x(h7,h12,p8,p11,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_t3_8_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_4(sh12,sp4,sp5,sh1,sh2,sh3,i1)

!     i1 ( h12 p4 p5 h1 h2 h3 )_vt + = 1/12 * P( 2 ) * Sum ( p7 p8 ) * t ( p4 p7 p8 h1 h2 h3 )_t * v ( h12 p5 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh12,sp4,sp5,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h12,p4,p5,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:,:,:)

  allocate(i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_8_4_perm(sh12,sp4,sp5,sh1,sh2,sh3,i1_perm)
  fact_p = +1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p4,p5,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh12 * sp4 * sp5 * sh1 * sh2 * sh3 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_4_perm(sh12,sp5,sp4,sh1,sh2,sh3,i1_perm)
  end if
  fact_p = -1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_8_4_perm(sh12,sp4,sp5,sh1,sh2,sh3,i1)

  implicit none
  integer,intent(in) :: sh12,sp4,sp5,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: h12,p4,p5,h1,h2,h3
  integer :: p7,p8,sp7,sp8
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 12.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp7,sp8,sh1,sh2,sh3)
     spin_int2x = tdcc_spin_int2x(sh12,sp5,sp7,sp8)
     if(spin_t3inp * spin_int2x == 0) cycle

     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact * &
             t3inp(p4,p7,p8,h1,h2,h3,spin_t3inp) * int2x(h12,p5,p7,p8,spin_int2x)
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
  end subroutine ccsdt_t3_8_4_perm
  !--------------------------------------------
end subroutine ccsdt_t3_8_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_5(sh12,sp4,sp5,sh1,sh2,sh3,i1)

!     i1 ( h12 p4 p5 h1 h2 h3 )_vtt + = -1/24 * P( 2 ) * Sum ( h8 ) * t ( p4 h8 )_t * i2 ( h8 h12 p5 h1 h2 h3 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh12,sp4,sp5,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h12,p4,p5,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:,:,:)

  allocate(i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_8_5_perm(sh12,sp4,sp5,sh1,sh2,sh3,i1_perm)
  fact_p = +1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p4,p5,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh12 * sp4 * sp5 * sh1 * sh2 * sh3 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_5_perm(sh12,sp5,sp4,sh1,sh2,sh3,i1_perm)
  end if
  fact_p = -1.0d+0
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact_p * i1_perm(h12,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_8_5_perm(sh12,sp4,sp5,sh1,sh2,sh3,i1)

  implicit none
  integer,intent(in) :: sh12,sp4,sp5,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: h12,p4,p5,h1,h2,h3
  integer :: h8,sh8
  integer :: spin_t1inp
  integer :: spin_itm_hhphhh
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhphhh(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 24.0d+0 * runit

  allocate(itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1))
  do sh8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp4,sh8)
     spin_itm_hhphhh = tdcc_spin_dummy3(sh8,sh12,sp5,sh1,sh2,sh3)
     if(spin_t1inp * spin_itm_hhphhh == 0) cycle

     itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_8_5_1(sh8,sh12,sp5,sh1,sh2,sh3,itm_hhphhh)

     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h8 = 1,norb1
        i1(h12,p4,p5,h1,h2,h3) = i1(h12,p4,p5,h1,h2,h3) + fact * &
             t1inp(p4,h8,spin_t1inp) * itm_hhphhh(h8,h12,p5,h1,h2,h3)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhphhh)
  end subroutine ccsdt_t3_8_5_perm
  !--------------------------------------------
end subroutine ccsdt_t3_8_5
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_8_5_1(sh8,sh12,sp4,sh1,sh2,sh3,i2)

!         i2 ( h8 h12 p4 h1 h2 h3 )_vt + = 1 * Sum ( p9 p10 ) * t ( p4 p9 p10 h1 h2 h3 )_t * v ( h8 h12 p9 p10 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh8,sh12,sp4,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i2(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h8,h12,p4,h1,h2,h3
  integer :: p9,p10,sp9,sp10
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp9 = 1,2
  do sp10 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp9,sp10,sh1,sh2,sh3)
     spin_int2x = tdcc_spin_int2x(sh8,sh12,sp9,sp10)
     if(spin_t3inp * spin_int2x == 0) cycle

     do h8 = 1,norb1
     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
        i2(h8,h12,p4,h1,h2,h3) = i2(h8,h12,p4,h1,h2,h3) + fact * &
             t3inp(p4,p9,p10,h1,h2,h3,spin_t3inp) * int2x(h8,h12,p9,p10,spin_int2x)
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
end subroutine ccsdt_t3_8_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_9(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vtt + = 1/4 * P( 3 ) * Sum ( h9 h10 ) * t ( p4 p5 h9 h10 )_t * i1 ( h9 h10 p6 h1 h2 h3 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccsdt_t3_9_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_9_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_9_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_t3_9_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h9,h10,sh9,sh10
  integer :: spin_t2inp
  integer :: spin_itm_hhphhh
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhphhh(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh10 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh9,sh10)
     spin_itm_hhphhh = tdcc_spin_dummy3(sh9,sh10,sp6,sh1,sh2,sh3)
     if(spin_t2inp * spin_itm_hhphhh == 0) cycle

     itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_t3_9_1(sh9,sh10,sp6,sh1,sh2,sh3,itm_hhphhh)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t2inp(p4,p5,h9,h10,spin_t2inp) * itm_hhphhh(h9,h10,p6,h1,h2,h3)
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
  deallocate(itm_hhphhh)
  end subroutine ccsdt_t3_9_perm
  !--------------------------------------------
end subroutine ccsdt_t3_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_t3_9_1(sh9,sh10,sp4,sh1,sh2,sh3,i1)

!     i1 ( h9 h10 p4 h1 h2 h3 )_vt + = 1 * Sum ( p7 p8 ) * t ( p4 p7 p8 h1 h2 h3 )_t * v ( h9 h10 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh9,sh10,sp4,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h9,h10,p4,h1,h2,h3
  integer :: p7,p8,sp7,sp8
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp7,sp8,sh1,sh2,sh3)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp7,sp8)
     if(spin_t3inp * spin_int2x == 0) cycle

     do h9 = 1,norb1
     do h10 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h9,h10,p4,h1,h2,h3) = i1(h9,h10,p4,h1,h2,h3) + fact * &
             t3inp(p4,p7,p8,h1,h2,h3,spin_t3inp) * int2x(h9,h10,p7,p8,spin_int2x)
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
end subroutine ccsdt_t3_9_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccsdt_t3_main()

  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3out

  implicit none

  call ccsdt_t3_1(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
  call ccsdt_t3_2(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
  call ccsdt_t3_3(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
  call ccsdt_t3_4(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
  call ccsdt_t3_5(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
  call ccsdt_t3_6(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
  call ccsdt_t3_7(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
  call ccsdt_t3_8(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
  call ccsdt_t3_9(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))

  call ccsdt_t3_1(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
  call ccsdt_t3_2(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
  call ccsdt_t3_3(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
  call ccsdt_t3_4(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
  call ccsdt_t3_5(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
  call ccsdt_t3_6(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
  call ccsdt_t3_7(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
  call ccsdt_t3_8(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
  call ccsdt_t3_9(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))

end subroutine ccsdt_t3_main
!**********************************************************
