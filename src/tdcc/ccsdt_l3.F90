!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_1(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = 1 * P( 9 ) * y ( h4 p1 )_y * v ( h5 h6 p2 p3 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_1_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_1_perm(sh5,sh4,sh6,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_1_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_1_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_1_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_1_perm(sh5,sh4,sh6,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_1_perm(sh5,sh4,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_1_perm(sh6,sh5,sh4,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_1_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_1_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: sdum
  integer :: spin_g1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sdum = 1,1
     spin_g1inp = tdcc_spin_g1inp(sh4,sp1)
     spin_int2x = tdcc_spin_int2x(sh5,sh6,sp2,sp3)
     if(spin_g1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g1inp(h4,p1,spin_g1inp) * int2x(h5,h6,p2,p3,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccsdt_l3_1_perm
  !--------------------------------------------
end subroutine ccsdt_l3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_2(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yf + = 1 * P( 9 ) * y ( h4 h5 p1 p2 )_y * i1 ( h6 p3 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_2_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_2_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_2_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_2_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_2_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_2_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_2_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_2_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_2_perm(sh4,sh6,sh5,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_2_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: sdum
  integer :: spin_g2inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sdum = 1,1
     spin_g2inp = tdcc_spin_g2inp(sh4,sh5,sp1,sp2)
     spin_itm_hp = tdcc_spin_fock(sh6,sp3)
     if(spin_g2inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_l3_2_1(sh6,sp3,itm_hp)
     call ccsdt_l3_2_2(sh6,sp3,itm_hp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g2inp(h4,h5,p1,p2,spin_g2inp) * itm_hp(h6,p3)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
  end subroutine ccsdt_l3_2_perm
  !--------------------------------------------
end subroutine ccsdt_l3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_2_1(sh4,sp1,i1)

!     i1 ( h4 p1 )_f + = 1 * f ( h4 p1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh4,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer :: h4,p1
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh4,sp1)
  do h4 = 1,norb1
  do p1 = norb1+1,nact
     i1(h4,p1) = i1(h4,p1) + fact * fock(h4,p1,spin_fock)
  end do
  end do
end subroutine ccsdt_l3_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_2_2(sh4,sp1,i1)

!     i1 ( h4 p1 )_vt + = 1 * Sum ( h8 p7 ) * t ( p7 h8 )_t * v ( h4 h8 p1 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer :: h4,p1
  integer :: h8,p7,sh8,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_int2x = tdcc_spin_int2x(sh4,sh8,sp1,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i1(h4,p1) = i1(h4,p1) + fact * t1inp(p7,h8,spin_t1inp) * int2x(h4,h8,p1,p7,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l3_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_3(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = -1 * P( 9 ) * Sum ( h9 ) * y ( h4 h9 p1 p2 )_y * i1 ( h5 h6 h9 p3 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_3_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_3_perm(sh5,sh4,sh6,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_3_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_3_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_3_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_3_perm(sh5,sh4,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_3_perm(sh5,sh4,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_3_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_3_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_3_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h9,sh9
  integer :: spin_g2inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh4,sh9,sp1,sp2)
     spin_itm_hhhp = tdcc_spin_int2x(sh5,sh6,sh9,sp3)
     if(spin_g2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l3_3_1(sh5,sh6,sh9,sp3,itm_hhhp)
     call ccsdt_l3_3_2(sh5,sh6,sh9,sp3,itm_hhhp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h9 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g2inp(h4,h9,p1,p2,spin_g2inp) * itm_hhhp(h5,h6,h9,p3)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_l3_3_perm
  !--------------------------------------------
end subroutine ccsdt_l3_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_3_1(sh4,sh5,sh9,sp1,i1)

!     i1 ( h4 h5 h9 p1 )_v + = 1 * v ( h4 h5 h9 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h9,p1
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh4,sh5,sh9,sp1)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h9 = 1,norb1
  do p1 = norb1+1,nact
     i1(h4,h5,h9,p1) = i1(h4,h5,h9,p1) + fact * int2x(h4,h5,h9,p1,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l3_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_3_2(sh4,sh5,sh9,sp1,i1)

!     i1 ( h4 h5 h9 p1 )_vt + = -1 * Sum ( p7 ) * t ( p7 h9 )_t * v ( h4 h5 p1 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h9,p1
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh9)
     spin_int2x = tdcc_spin_int2x(sh4,sh5,sp1,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h9 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h4,h5,h9,p1) = i1(h4,h5,h9,p1) + fact * &
             t1inp(p7,h9,spin_t1inp) * int2x(h4,h5,p1,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l3_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_4(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = -1 * P( 9 ) * Sum ( p7 ) * y ( h4 h5 p1 p7 )_y * v ( h6 p7 p2 p3 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_4_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_4_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_4_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_4_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_4_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_4_perm(sh6,sh5,sh4,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_4_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_4_perm(sh4,sh6,sh5,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_4_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_4_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: p7,sp7
  integer :: spin_g2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh4,sh5,sp1,sp7)
     spin_int2x = tdcc_spin_int2x(sh6,sp7,sp2,sp3)
     if(spin_g2inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g2inp(h4,h5,p1,p7,spin_g2inp) * int2x(h6,p7,p2,p3,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccsdt_l3_4_perm
  !--------------------------------------------
end subroutine ccsdt_l3_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_5(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yf + = -1 * P( 3 ) * Sum ( h11 ) * y ( h4 h5 h11 p1 p2 p3 )_y * i1 ( h6 h11 )_f 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_5_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_5_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_5_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_5_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h11,sh11
  integer :: spin_g3inp
  integer :: spin_itm_hh
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh11 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh11,sp1,sp2,sp3)
     spin_itm_hh = tdcc_spin_fock(sh6,sh11)
     if(spin_g3inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsdt_l3_5_1(sh6,sh11,itm_hh)
     call ccsdt_l3_5_2(sh6,sh11,itm_hh)
     call ccsdt_l3_5_3(sh6,sh11,itm_hh)
     call ccsdt_l3_5_4(sh6,sh11,itm_hh)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h11 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h11,p1,p2,p3,spin_g3inp) * itm_hh(h6,h11)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccsdt_l3_5_perm
  !--------------------------------------------
end subroutine ccsdt_l3_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_5_1(sh4,sh11,i1)

!     i1 ( h4 h11 )_f + = 1 * f ( h4 h11 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh4,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h4,h11
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh4,sh11)
  do h4 = 1,norb1
  do h11 = 1,norb1
     i1(h4,h11) = i1(h4,h11) + fact * fock(h4,h11,spin_fock)
  end do
  end do
end subroutine ccsdt_l3_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_5_2(sh4,sh11,i1)

!     i1 ( h4 h11 )_ft + = 1 * Sum ( p7 ) * t ( p7 h11 )_t * i2 ( h4 p7 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh4,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h4,h11
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh11)
     spin_itm_hp = tdcc_spin_fock(sh4,sp7)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_l3_5_2_1(sh4,sp7,itm_hp)
     call ccsdt_l3_5_2_2(sh4,sp7,itm_hp)

     do h4 = 1,norb1
     do h11 = 1,norb1
     do p7 = norb1+1,nact
        i1(h4,h11) = i1(h4,h11) + fact * t1inp(p7,h11,spin_t1inp) * itm_hp(h4,p7)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_l3_5_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_5_2_1(sh4,sp7,i2)

!         i2 ( h4 p7 )_f + = 1 * f ( h4 p7 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh4,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h4,p7
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh4,sp7)
  do h4 = 1,norb1
  do p7 = norb1+1,nact
     i2(h4,p7) = i2(h4,p7) + fact * fock(h4,p7,spin_fock)
  end do
  end do
end subroutine ccsdt_l3_5_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_5_2_2(sh4,sp7,i2)

!         i2 ( h4 p7 )_vt + = 1 * Sum ( h10 p9 ) * t ( p9 h10 )_t * v ( h4 h10 p7 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h4,p7
  integer :: h10,p9,sh10,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh10 = 1,2
  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
     spin_int2x = tdcc_spin_int2x(sh4,sh10,sp7,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do p7 = norb1+1,nact
     do h10 = 1,norb1
     do p9 = norb1+1,nact
        i2(h4,p7) = i2(h4,p7) + fact * &
             t1inp(p9,h10,spin_t1inp) * int2x(h4,h10,p7,p9,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l3_5_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_5_3(sh4,sh11,i1)

!     i1 ( h4 h11 )_vt + = 1 * Sum ( h8 p7 ) * t ( p7 h8 )_t * v ( h4 h8 h11 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h4,h11
  integer :: h8,p7,sh8,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_int2x = tdcc_spin_int2x(sh4,sh8,sh11,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h11 = 1,norb1
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i1(h4,h11) = i1(h4,h11) + fact * &
             t1inp(p7,h8,spin_t1inp) * int2x(h4,h8,h11,p7,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l3_5_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_5_4(sh4,sh11,i1)

!     i1 ( h4 h11 )_vt + = -1/2 * Sum ( h10 p7 p8 ) * t ( p7 p8 h10 h11 )_t * v ( h4 h10 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h4,h11
  integer :: h10,p7,p8,sh10,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh10 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh10,sh11)
     spin_int2x = tdcc_spin_int2x(sh4,sh10,sp7,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h11 = 1,norb1
     do h10 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h4,h11) = i1(h4,h11) + fact * t2inp(p7,p8,h10,h11,spin_t2inp) * int2x(h4,h10,p7,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_l3_5_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_6(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yf + = 1 * P( 3 ) * Sum ( p12 ) * y ( h4 h5 h6 p1 p2 p12 )_y * i1 ( p12 p3 )_f 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_6_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_6_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_6_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_6_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: p12,sp12
  integer :: spin_g3inp
  integer :: spin_itm_pp
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp12 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp2,sp12)
     spin_itm_pp = tdcc_spin_fock(sp12,sp3)
     if(spin_g3inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_6_1(sp12,sp3,itm_pp)
     call ccsdt_l3_6_2(sp12,sp3,itm_pp)
     call ccsdt_l3_6_3(sp12,sp3,itm_pp)
     call ccsdt_l3_6_4(sp12,sp3,itm_pp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p12 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h6,p1,p2,p12,spin_g3inp) * itm_pp(p12,p3)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccsdt_l3_6_perm
  !--------------------------------------------
end subroutine ccsdt_l3_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_6_1(sp12,sp1,i1)

!     i1 ( p12 p1 )_f + = 1 * f ( p12 p1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sp12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p12,p1
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp12,sp1)
  do p12 = norb1+1,nact
  do p1 = norb1+1,nact
     i1(p12,p1) = i1(p12,p1) + fact * fock(p12,p1,spin_fock)
  end do
  end do
end subroutine ccsdt_l3_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_6_2(sp12,sp1,i1)

!     i1 ( p12 p1 )_vt + = -1 * Sum ( h8 p7 ) * t ( p7 h8 )_t * v ( h8 p12 p1 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p12,p1
  integer :: h8,p7,sh8,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh8 = 1,2
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_int2x = tdcc_spin_int2x(sh8,sp12,sp1,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p12 = norb1+1,nact
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i1(p12,p1) = i1(p12,p1) + fact * t1inp(p7,h8,spin_t1inp) * int2x(h8,p12,p1,p7,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l3_6_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_6_3(sp12,sp1,i1)

!     i1 ( p12 p1 )_vt + = 1/2 * Sum ( h9 h10 p8 ) * t ( p8 p12 h9 h10 )_t * v ( h9 h10 p1 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p12,p1
  integer :: h9,h10,p8,sh9,sh10,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh9 = 1,2
  do sh10 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp8,sp12,sh9,sh10)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp1,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     do p12 = norb1+1,nact
     do p1 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p8 = norb1+1,nact
        i1(p12,p1) = i1(p12,p1) + fact * t2inp(p8,p12,h9,h10,spin_t2inp) * int2x(h9,h10,p1,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_l3_6_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_6_4(sp12,sp1,i1)

!     i1 ( p12 p1 )_vtt + = -1 * Sum ( h8 ) * t ( p12 h8 )_t * i2 ( h8 p1 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sp12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p12,p1
  integer :: h8,sh8
  integer :: spin_t1inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy1
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sh8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp12,sh8)
     spin_itm_hp = tdcc_spin_dummy1(sh8,sp1)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_l3_6_4_1(sh8,sp1,itm_hp)

     do p12 = norb1+1,nact
     do p1 = norb1+1,nact
     do h8 = 1,norb1
        i1(p12,p1) = i1(p12,p1) + fact * t1inp(p12,h8,spin_t1inp) * itm_hp(h8,p1)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_l3_6_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_6_4_1(sh8,sp1,i2)

!         i2 ( h8 p1 )_vt + = 1 * Sum ( h10 p9 ) * t ( p9 h10 )_t * v ( h8 h10 p1 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh8,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h8,p1
  integer :: h10,p9,sh10,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh10 = 1,2
  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
     spin_int2x = tdcc_spin_int2x(sh8,sh10,sp1,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h8 = 1,norb1
     do p1 = norb1+1,nact
     do h10 = 1,norb1
     do p9 = norb1+1,nact
        i2(h8,p1) = i2(h8,p1) + fact * t1inp(p9,h10,spin_t1inp) * int2x(h8,h10,p1,p9,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l3_6_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_7(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = 1/2 * P( 3 ) * Sum ( h11 h12 ) * y ( h4 h11 h12 p1 p2 p3 )_y * i1 ( h5 h6 h11 h12 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_7_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_7_perm(sh5,sh4,sh6,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_7_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_7_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h11,h12,sh11,sh12
  integer :: spin_g3inp
  integer :: spin_itm_hhhh
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh11 = 1,2
  do sh12 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh11,sh12,sp1,sp2,sp3)
     spin_itm_hhhh = tdcc_spin_int2x(sh5,sh6,sh11,sh12)
     if(spin_g3inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_l3_7_1(sh5,sh6,sh11,sh12,itm_hhhh)
     call ccsdt_l3_7_2(sh5,sh6,sh11,sh12,itm_hhhh)
     call ccsdt_l3_7_3(sh5,sh6,sh11,sh12,itm_hhhh)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h11,h12,p1,p2,p3,spin_g3inp) * itm_hhhh(h5,h6,h11,h12)
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
  end subroutine ccsdt_l3_7_perm
  !--------------------------------------------
end subroutine ccsdt_l3_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_7_1(sh4,sh5,sh11,sh12,i1)

!     i1 ( h4 h5 h11 h12 )_v + = 1 * v ( h4 h5 h11 h12 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h4,h5,h11,h12
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh4,sh5,sh11,sh12)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h11 = 1,norb1
  do h12 = 1,norb1
     i1(h4,h5,h11,h12) = i1(h4,h5,h11,h12) + fact * int2x(h4,h5,h11,h12,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l3_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_7_2(sh4,sh5,sh11,sh12,i1)

!     i1 ( h4 h5 h11 h12 )_vt + = -2 * Sum ( p7 ) * t ( p7 h11 )_t * i2 ( h4 h5 h12 p7 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h4,h5,h11,h12
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh11)
     spin_itm_hhhp = tdcc_spin_int2x(sh4,sh5,sh12,sp7)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l3_7_2_1(sh4,sh5,sh12,sp7,itm_hhhp)
     call ccsdt_l3_7_2_2(sh4,sh5,sh12,sp7,itm_hhhp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p7 = norb1+1,nact
        i1(h4,h5,h11,h12) = i1(h4,h5,h11,h12) + fact * t1inp(p7,h11,spin_t1inp) * itm_hhhp(h4,h5,h12,p7)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l3_7_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_7_2_1(sh4,sh5,sh12,sp7,i2)

!         i2 ( h4 h5 h12 p7 )_v + = 1 * v ( h4 h5 h12 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh12,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h12,p7
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh4,sh5,sh12,sp7)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h12 = 1,norb1
  do p7 = norb1+1,nact
     i2(h4,h5,h12,p7) = i2(h4,h5,h12,p7) + fact * int2x(h4,h5,h12,p7,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l3_7_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_7_2_2(sh4,sh5,sh12,sp7,i2)

!         i2 ( h4 h5 h12 p7 )_vt + = -1/2 * Sum ( p9 ) * t ( p9 h12 )_t * v ( h4 h5 p7 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh12,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h12,p7
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh12)
     spin_int2x = tdcc_spin_int2x(sh4,sh5,sp7,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h12 = 1,norb1
     do p7 = norb1+1,nact
     do p9 = norb1+1,nact
        i2(h4,h5,h12,p7) = i2(h4,h5,h12,p7) + fact * &
             t1inp(p9,h12,spin_t1inp) * int2x(h4,h5,p7,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l3_7_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_7_3(sh4,sh5,sh11,sh12,i1)

!     i1 ( h4 h5 h11 h12 )_vt + = 1/2 * Sum ( p7 p8 ) * t ( p7 p8 h11 h12 )_t * v ( h4 h5 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h4,h5,h11,h12
  integer :: p7,p8,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh11,sh12)
     spin_int2x = tdcc_spin_int2x(sh4,sh5,sp7,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h4,h5,h11,h12) = i1(h4,h5,h11,h12) + fact * &
             t2inp(p7,p8,h11,h12,spin_t2inp) * int2x(h4,h5,p7,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l3_7_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_8(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = -1 * P( 9 ) * Sum ( h12 p9 ) * y ( h4 h5 h12 p1 p2 p9 )_y * i1 ( h6 p9 h12 p3 )_v 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_8_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_8_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_8_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_8_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_8_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_8_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_8_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_8_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_8_perm(sh4,sh6,sh5,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_8_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h12,p9,sh12,sp9
  integer :: spin_g3inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh12 = 1,2
  do sp9 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh12,sp1,sp2,sp9)
     spin_itm_hphp = tdcc_spin_int2x(sh6,sp9,sh12,sp3)
     if(spin_g3inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l3_8_1(sh6,sp9,sh12,sp3,itm_hphp)
     call ccsdt_l3_8_2(sh6,sp9,sh12,sp3,itm_hphp)
     call ccsdt_l3_8_3(sh6,sp9,sh12,sp3,itm_hphp)
     call ccsdt_l3_8_4(sh6,sp9,sh12,sp3,itm_hphp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h12 = 1,norb1
     do p9 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h12,p1,p2,p9,spin_g3inp) * itm_hphp(h6,p9,h12,p3)
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
  end subroutine ccsdt_l3_8_perm
  !--------------------------------------------
end subroutine ccsdt_l3_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_8_1(sh4,sp9,sh12,sp1,i1)

!     i1 ( h4 p9 h12 p1 )_v + = 1 * v ( h4 p9 h12 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sp9,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h4,p9,h12,p1
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh4,sp9,sh12,sp1)
  do h4 = 1,norb1
  do p9 = norb1+1,nact
  do h12 = 1,norb1
  do p1 = norb1+1,nact
     i1(h4,p9,h12,p1) = i1(h4,p9,h12,p1) + fact * int2x(h4,p9,h12,p1,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l3_8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_8_2(sh4,sp9,sh12,sp1,i1)

!     i1 ( h4 p9 h12 p1 )_vt + = -1 * Sum ( p7 ) * t ( p7 h12 )_t * v ( h4 p9 p1 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sp9,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h4,p9,h12,p1
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh12)
     spin_int2x = tdcc_spin_int2x(sh4,sp9,sp1,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do p9 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h4,p9,h12,p1) = i1(h4,p9,h12,p1) + fact * t1inp(p7,h12,spin_t1inp) * int2x(h4,p9,p1,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l3_8_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_8_3(sh4,sp9,sh12,sp1,i1)

!     i1 ( h4 p9 h12 p1 )_vt + = -1 * Sum ( h10 p8 ) * t ( p8 p9 h10 h12 )_t * v ( h4 h10 p1 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sp9,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h4,p9,h12,p1
  integer :: h10,p8,sh10,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh10 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp8,sp9,sh10,sh12)
     spin_int2x = tdcc_spin_int2x(sh4,sh10,sp1,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do p9 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h10 = 1,norb1
     do p8 = norb1+1,nact
        i1(h4,p9,h12,p1) = i1(h4,p9,h12,p1) + fact * &
             t2inp(p8,p9,h10,h12,spin_t2inp) * int2x(h4,h10,p1,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l3_8_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_8_4(sh4,sp9,sh12,sp1,i1)

!     i1 ( h4 p9 h12 p1 )_vtt + = 1 * Sum ( h10 ) * t ( p9 h10 )_t * i2 ( h4 h10 h12 p1 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh4,sp9,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h4,p9,h12,p1
  integer :: h10,sh10
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh10 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
     spin_itm_hhhp = tdcc_spin_dummy2(sh4,sh10,sh12,sp1)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l3_8_4_1(sh4,sh10,sh12,sp1,itm_hhhp)

     do h4 = 1,norb1
     do p9 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h10 = 1,norb1
        i1(h4,p9,h12,p1) = i1(h4,p9,h12,p1) + fact * &
             t1inp(p9,h10,spin_t1inp) * itm_hhhp(h4,h10,h12,p1)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l3_8_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_8_4_1(sh4,sh10,sh12,sp1,i2)

!         i2 ( h4 h10 h12 p1 )_vt + = 1 * Sum ( p7 ) * t ( p7 h12 )_t * v ( h4 h10 p1 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh10,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h10,h12,p1
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh12)
     spin_int2x = tdcc_spin_int2x(sh4,sh10,sp1,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h10 = 1,norb1
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h4,h10,h12,p1) = i2(h4,h10,h12,p1) + fact * &
             t1inp(p7,h12,spin_t1inp) * int2x(h4,h10,p1,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l3_8_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_9(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = 1/2 * P( 3 ) * Sum ( p7 p8 ) * y ( h4 h5 h6 p1 p7 p8 )_y * v ( p7 p8 p2 p3 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_9_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_9_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_9_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_9_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: p7,p8,sp7,sp8
  integer :: spin_g3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp7,sp8)
     spin_int2x = tdcc_spin_int2x(sp7,sp8,sp2,sp3)
     if(spin_g3inp * spin_int2x == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h6,p1,p7,p8,spin_g3inp) * int2x(p7,p8,p2,p3,spin_int2x)
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
  end subroutine ccsdt_l3_9_perm
  !--------------------------------------------
end subroutine ccsdt_l3_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_10(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = -1 * P( 9 ) * Sum ( h11 ) * i1 ( h4 h5 h11 p1 )_yt * v ( h6 h11 p2 p3 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_10_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_10_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_10_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_10_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_10_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_10_perm(sh6,sh5,sh4,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_10_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_10_perm(sh4,sh6,sh5,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_10_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_10_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h11,sh11
  integer :: spin_itm_hhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh11 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh4,sh5,sh11,sp1)
     spin_int2x = tdcc_spin_int2x(sh6,sh11,sp2,sp3)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l3_10_1(sh4,sh5,sh11,sp1,itm_hhhp)
     call ccsdt_l3_10_2(sh4,sh5,sh11,sp1,itm_hhhp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h11 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             itm_hhhp(h4,h5,h11,p1) * int2x(h6,h11,p2,p3,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_l3_10_perm
  !--------------------------------------------
end subroutine ccsdt_l3_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_10_1(sh4,sh5,sh11,sp1,i1)

!     i1 ( h4 h5 h11 p1 )_yt + = -1 * Sum ( p7 ) * t ( p7 h11 )_t * y ( h4 h5 p1 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh11,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h11,p1
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh11)
     spin_g2inp = tdcc_spin_g2inp(sh4,sh5,sp1,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h11 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h4,h5,h11,p1) = i1(h4,h5,h11,p1) + fact * &
             t1inp(p7,h11,spin_t1inp) * g2inp(h4,h5,p1,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l3_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_10_2(sh4,sh5,sh11,sp1,i1)

!     i1 ( h4 h5 h11 p1 )_yt + = 1/2 * Sum ( h9 p7 p8 ) * t ( p7 p8 h9 h11 )_t * y ( h4 h5 h9 p1 p7 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh11,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h11,p1
  integer :: h9,p7,p8,sh9,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh9 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh11)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh9,sp1,sp7,sp8)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h11 = 1,norb1
     do p1 = norb1+1,nact
     do h9 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h4,h5,h11,p1) = i1(h4,h5,h11,p1) + fact * &
             t2inp(p7,p8,h9,h11,spin_t2inp) * g3inp(h4,h5,h9,p1,p7,p8,spin_g3inp)
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
end subroutine ccsdt_l3_10_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_11(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytf + = -1 * P( 3 ) * Sum ( h7 ) * i1 ( h4 h5 h6 h7 p1 p2 )_yt * f ( h7 p3 )_f 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_11_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_11_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_11_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_11_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h7,sh7
  integer :: spin_itm_hhhhpp
  integer :: spin_fock
  integer,external :: tdcc_spin_dummy3
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hhhhpp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sh7 = 1,2
     spin_itm_hhhhpp = tdcc_spin_dummy3(sh4,sh5,sh6,sh7,sp1,sp2)
     spin_fock = tdcc_spin_fock(sh7,sp3)
     if(spin_itm_hhhhpp * spin_fock == 0) cycle

     itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_11_1(sh4,sh5,sh6,sh7,sp1,sp2,itm_hhhhpp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h7 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             itm_hhhhpp(h4,h5,h6,h7,p1,p2) * fock(h7,p3,spin_fock)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhhpp)
  end subroutine ccsdt_l3_11_perm
  !--------------------------------------------
end subroutine ccsdt_l3_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_11_1(sh4,sh5,sh6,sh7,sp1,sp2,i1)

!     i1 ( h4 h5 h6 h7 p1 p2 )_yt + = 1 * Sum ( p8 ) * t ( p8 h7 )_t * y ( h4 h5 h6 p1 p2 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sh7,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,h7,p1,p2
  integer :: p8,sp8
  integer :: spin_t1inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp8,sh7)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp2,sp8)
     if(spin_t1inp * spin_g3inp == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h4,h5,h6,h7,p1,p2) = i1(h4,h5,h6,h7,p1,p2) + fact * &
             t1inp(p8,h7,spin_t1inp) * g3inp(h4,h5,h6,p1,p2,p8,spin_g3inp)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l3_11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_12(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = 1 * P( 9 ) * Sum ( h9 h8 ) * i1 ( h4 h5 h9 h8 p1 p2 )_yt * v ( h6 h8 h9 p3 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_12_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_12_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_12_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_12_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_12_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_12_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_12_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_12_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_12_perm(sh4,sh6,sh5,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_12_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h9,h8,sh9,sh8
  integer :: spin_itm_hhhhpp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy3
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhhpp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sh9 = 1,2
  do sh8 = 1,2
     spin_itm_hhhhpp = tdcc_spin_dummy3(sh4,sh5,sh9,sh8,sp1,sp2)
     spin_int2x = tdcc_spin_int2x(sh6,sh8,sh9,sp3)
     if(spin_itm_hhhhpp * spin_int2x == 0) cycle

     itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_12_1(sh4,sh5,sh9,sh8,sp1,sp2,itm_hhhhpp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h9 = 1,norb1
     do h8 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             itm_hhhhpp(h4,h5,h9,h8,p1,p2) * int2x(h6,h8,h9,p3,spin_int2x)
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
  deallocate(itm_hhhhpp)
  end subroutine ccsdt_l3_12_perm
  !--------------------------------------------
end subroutine ccsdt_l3_12
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_12_1(sh4,sh5,sh9,sh8,sp1,sp2,i1)

!     i1 ( h4 h5 h9 h8 p1 p2 )_yt + = 1 * Sum ( p7 ) * t ( p7 h8 )_t * y ( h4 h5 h9 p1 p2 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh9,sh8,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h9,h8,p1,p2
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh9,sp1,sp2,sp7)
     if(spin_t1inp * spin_g3inp == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h9 = 1,norb1
     do h8 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h4,h5,h9,h8,p1,p2) = i1(h4,h5,h9,h8,p1,p2) + fact * &
             t1inp(p7,h8,spin_t1inp) * g3inp(h4,h5,h9,p1,p2,p7,spin_g3inp)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l3_12_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_13(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = 1 * P( 3 ) * Sum ( h8 p9 ) * i1 ( h4 h5 h6 h8 p1 p9 )_yt * v ( h8 p9 p2 p3 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_13_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_13_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_13_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_13_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h8,p9,sh8,sp9
  integer :: spin_itm_hhhhpp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy3
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhhpp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sh8 = 1,2
  do sp9 = 1,2
     spin_itm_hhhhpp = tdcc_spin_dummy3(sh4,sh5,sh6,sh8,sp1,sp9)
     spin_int2x = tdcc_spin_int2x(sh8,sp9,sp2,sp3)
     if(spin_itm_hhhhpp * spin_int2x == 0) cycle

     itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_13_1(sh4,sh5,sh6,sh8,sp1,sp9,itm_hhhhpp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h8 = 1,norb1
     do p9 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             itm_hhhhpp(h4,h5,h6,h8,p1,p9) * int2x(h8,p9,p2,p3,spin_int2x)
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
  deallocate(itm_hhhhpp)
  end subroutine ccsdt_l3_13_perm
  !--------------------------------------------
end subroutine ccsdt_l3_13
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_13_1(sh4,sh5,sh6,sh8,sp1,sp9,i1)

!     i1 ( h4 h5 h6 h8 p1 p9 )_yt + = -1 * Sum ( p7 ) * t ( p7 h8 )_t * y ( h4 h5 h6 p1 p7 p9 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sh8,sp1,sp9
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,h8,p1,p9
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp7,sp9)
     if(spin_t1inp * spin_g3inp == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h8 = 1,norb1
     do p1 = norb1+1,nact
     do p9 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h4,h5,h6,h8,p1,p9) = i1(h4,h5,h6,h8,p1,p9) + fact * &
             t1inp(p7,h8,spin_t1inp) * g3inp(h4,h5,h6,p1,p7,p9,spin_g3inp)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l3_13_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_14(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = -1/2 * P( 9 ) * Sum ( p8 ) * i1 ( h4 p8 p1 p2 )_yt * v ( h5 h6 p3 p8 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_14_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_14_perm(sh5,sh4,sh6,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_14_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_14_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_14_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_14_perm(sh5,sh4,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_14_perm(sh5,sh4,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_14_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_14_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_14_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: p8,sp8
  integer :: spin_itm_hppp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp8 = 1,2
     spin_itm_hppp = tdcc_spin_dummy2(sh4,sp8,sp1,sp2)
     spin_int2x = tdcc_spin_int2x(sh5,sh6,sp3,sp8)
     if(spin_itm_hppp * spin_int2x == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_14_1(sh4,sp8,sp1,sp2,itm_hppp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             itm_hppp(h4,p8,p1,p2) * int2x(h5,h6,p3,p8,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hppp)
  end subroutine ccsdt_l3_14_perm
  !--------------------------------------------
end subroutine ccsdt_l3_14
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_14_1(sh4,sp8,sp1,sp2,i1)

!     i1 ( h4 p8 p1 p2 )_yt + = 1 * Sum ( h9 h10 p7 ) * t ( p7 p8 h9 h10 )_t * y ( h4 h9 h10 p1 p2 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sp8,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,p8,p1,p2
  integer :: h9,h10,p7,sh9,sh10,sp7
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh9 = 1,2
  do sh10 = 1,2
  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh9,sh10,sp1,sp2,sp7)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h4 = 1,norb1
     do p8 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p7 = norb1+1,nact
        i1(h4,p8,p1,p2) = i1(h4,p8,p1,p2) + fact * &
             t2inp(p7,p8,h9,h10,spin_t2inp) * g3inp(h4,h9,h10,p1,p2,p7,spin_g3inp)
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
end subroutine ccsdt_l3_14_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_15(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = -1/4 * P( 3 ) * Sum ( h11 h10 ) * i1 ( h4 h5 h6 h10 h11 p1 )_yt * v ( h10 h11 p2 p3 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l3_15_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_15_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_15_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l3_15_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h11,h10,sh11,sh10
  integer :: spin_itm_hhhhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy3
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhhhp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh11 = 1,2
  do sh10 = 1,2
     spin_itm_hhhhhp = tdcc_spin_dummy3(sh4,sh5,sh6,sh10,sh11,sp1)
     spin_int2x = tdcc_spin_int2x(sh10,sh11,sp2,sp3)
     if(spin_itm_hhhhhp * spin_int2x == 0) cycle

     itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l3_15_1(sh4,sh5,sh6,sh10,sh11,sp1,itm_hhhhhp)
     call ccsdt_l3_15_2(sh4,sh5,sh6,sh10,sh11,sp1,itm_hhhhhp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h11 = 1,norb1
     do h10 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             itm_hhhhhp(h4,h5,h6,h10,h11,p1) * int2x(h10,h11,p2,p3,spin_int2x)
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
  deallocate(itm_hhhhhp)
  end subroutine ccsdt_l3_15_perm
  !--------------------------------------------
end subroutine ccsdt_l3_15
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_15_1(sh4,sh5,sh6,sh10,sh11,sp1,i1)

!     i1 ( h4 h5 h6 h10 h11 p1 )_yt + = -1 * Sum ( p7 p8 ) * t ( p7 p8 h10 h11 )_t * y ( h4 h5 h6 p1 p7 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sh10,sh11,sp1
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h6,h10,h11,p1
  integer :: p7,p8,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh10,sh11)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp7,sp8)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h10 = 1,norb1
     do h11 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h4,h5,h6,h10,h11,p1) = i1(h4,h5,h6,h10,h11,p1) + fact * &
             t2inp(p7,p8,h10,h11,spin_t2inp) * g3inp(h4,h5,h6,p1,p7,p8,spin_g3inp)
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
end subroutine ccsdt_l3_15_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_15_2(sh4,sh5,sh6,sh10,sh11,sp1,i1)

!     i1 ( h4 h5 h6 h10 h11 p1 )_ytt + = -2 * Sum ( p7 ) * t ( p7 h11 )_t * i2 ( h4 h5 h6 h10 p1 p7 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sh10,sh11,sp1
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h6,h10,h11,p1
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_itm_hhhhpp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhhhpp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh11)
     spin_itm_hhhhpp = tdcc_spin_dummy3(sh4,sh5,sh6,sh10,sp1,sp7)
     if(spin_t1inp * spin_itm_hhhhpp == 0) cycle

     itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l3_15_2_1(sh4,sh5,sh6,sh10,sp1,sp7,itm_hhhhpp)

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h10 = 1,norb1
     do h11 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h4,h5,h6,h10,h11,p1) = i1(h4,h5,h6,h10,h11,p1) + fact * &
             t1inp(p7,h11,spin_t1inp) * itm_hhhhpp(h4,h5,h6,h10,p1,p7)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhhpp)
end subroutine ccsdt_l3_15_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l3_15_2_1(sh4,sh5,sh6,sh10,sp1,sp7,i2)

!         i2 ( h4 h5 h6 h10 p1 p7 )_yt + = -1 * Sum ( p9 ) * t ( p9 h10 )_t * y ( h4 h5 h6 p1 p7 p9 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sh10,sp1,sp7
  complex(kind(0d0)),intent(inout) :: &
       i2(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,h10,p1,p7
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp7,sp9)
     if(spin_t1inp * spin_g3inp == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h10 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
     do p9 = norb1+1,nact
        i2(h4,h5,h6,h10,p1,p7) = i2(h4,h5,h6,h10,p1,p7) + fact * &
             t1inp(p9,h10,spin_t1inp) * g3inp(h4,h5,h6,p1,p7,p9,spin_g3inp)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l3_15_2_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccsdt_l3_main()

  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3out

  implicit none

  call ccsdt_l3_1(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_2(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_3(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_4(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_5(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_6(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_7(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_8(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_9(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_10(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_11(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_12(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_13(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_14(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccsdt_l3_15(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))

  call ccsdt_l3_1(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_2(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_3(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_4(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_5(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_6(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_7(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_8(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_9(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_10(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_11(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_12(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_13(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_14(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccsdt_l3_15(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))

end subroutine ccsdt_l3_main
!**********************************************************
