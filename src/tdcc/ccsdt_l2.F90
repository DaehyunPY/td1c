!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_1(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_v + = 1 * v ( h3 h4 p1 p2 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh4,sp1,sp2)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
          int2x(h3,h4,p1,p2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_2(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = 1 * P( 4 ) * y ( h3 p1 )_y * i1 ( h4 p2 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_2_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_2_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_2_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_2_perm(sh4,sh3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_2_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: sdum
  integer :: spin_g1inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_g1inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sdum = 1,1
     spin_g1inp = tdcc_spin_g1inp(sh3,sp1)
     spin_itm_hp = tdcc_spin_fock(sh4,sp2)
     if(spin_g1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_2_1(sh4,sp2,itm_hp)
     call ccsdt_l2_2_2(sh4,sp2,itm_hp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g1inp(h3,p1,spin_g1inp) * itm_hp(h4,p2)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
  end subroutine ccsdt_l2_2_perm
  !--------------------------------------------
end subroutine ccsdt_l2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_2_1(sh3,sp1,i1)

!     i1 ( h3 p1 )_f + = 1 * f ( h3 p1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh3,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer :: h3,p1
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh3,sp1)
  do h3 = 1,norb1
  do p1 = norb1+1,nact
     i1(h3,p1) = i1(h3,p1) + fact * fock(h3,p1,spin_fock)
  end do
  end do
end subroutine ccsdt_l2_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_2_2(sh3,sp1,i1)

!     i1 ( h3 p1 )_vt + = 1 * Sum ( h6 p5 ) * t ( p5 h6 )_t * v ( h3 h6 p1 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer :: h3,p1
  integer :: h6,p5,sh6,sp5
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh6 = 1,2
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_int2x = tdcc_spin_int2x(sh3,sh6,sp1,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,p1) = i1(h3,p1) + fact * t1inp(p5,h6,spin_t1inp) * int2x(h3,h6,p1,p5,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_3(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1 * P( 2 ) * Sum ( h7 ) * y ( h7 p1 )_y * i1 ( h3 h4 h7 p2 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_3_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_3_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_3_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h7,sh7
  integer :: spin_g1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_g1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh7 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh7,sp1)
     spin_itm_hhhp = tdcc_spin_int2x(sh3,sh4,sh7,sp2)
     if(spin_g1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_3_1(sh3,sh4,sh7,sp2,itm_hhhp)
     call ccsdt_l2_3_2(sh3,sh4,sh7,sp2,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h7 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g1inp(h7,p1,spin_g1inp) * itm_hhhp(h3,h4,h7,p2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_l2_3_perm
  !--------------------------------------------
end subroutine ccsdt_l2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_3_1(sh3,sh4,sh7,sp1,i1)

!     i1 ( h3 h4 h7 p1 )_v + = 1 * v ( h3 h4 h7 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sh7,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h7,p1
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh4,sh7,sp1)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do h7 = 1,norb1
  do p1 = norb1+1,nact
     i1(h3,h4,h7,p1) = i1(h3,h4,h7,p1) + fact * int2x(h3,h4,h7,p1,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_3_2(sh3,sh4,sh7,sp1,i1)

!     i1 ( h3 h4 h7 p1 )_vt + = -1 * Sum ( p5 ) * t ( p5 h7 )_t * v ( h3 h4 p1 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sh7,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h7,p1
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh7)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp1,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h7 = 1,norb1
     do p1 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h3,h4,h7,p1) = i1(h3,h4,h7,p1) + fact * t1inp(p5,h7,spin_t1inp) * int2x(h3,h4,p1,p5,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_4(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1 * P( 2 ) * Sum ( p5 ) * y ( h3 p5 )_y * v ( h4 p5 p1 p2 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_4_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_4_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_4_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p5,sp5
  integer :: spin_g1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh3,sp5)
     spin_int2x = tdcc_spin_int2x(sh4,sp5,sp1,sp2)
     if(spin_g1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p5 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g1inp(h3,p5,spin_g1inp) * int2x(h4,p5,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccsdt_l2_4_perm
  !--------------------------------------------
end subroutine ccsdt_l2_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_5(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = -1 * P( 2 ) * Sum ( h9 ) * y ( h3 h9 p1 p2 )_y * i1 ( h4 h9 )_f 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_5_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_5_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_5_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h9,sh9
  integer :: spin_g2inp
  integer :: spin_itm_hh
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh9 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh9,sp1,sp2)
     spin_itm_hh = tdcc_spin_fock(sh4,sh9)
     if(spin_g2inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsdt_l2_5_1(sh4,sh9,itm_hh)
     call ccsdt_l2_5_2(sh4,sh9,itm_hh)
     call ccsdt_l2_5_3(sh4,sh9,itm_hh)
     call ccsdt_l2_5_4(sh4,sh9,itm_hh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g2inp(h3,h9,p1,p2,spin_g2inp) * itm_hh(h4,h9)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccsdt_l2_5_perm
  !--------------------------------------------
end subroutine ccsdt_l2_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_5_1(sh3,sh9,i1)

!     i1 ( h3 h9 )_f + = 1 * f ( h3 h9 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh3,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h9
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh3,sh9)
  do h3 = 1,norb1
  do h9 = 1,norb1
     i1(h3,h9) = i1(h3,h9) + fact * fock(h3,h9,spin_fock)
  end do
  end do
end subroutine ccsdt_l2_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_5_2(sh3,sh9,i1)

!     i1 ( h3 h9 )_ft + = 1 * Sum ( p5 ) * t ( p5 h9 )_t * i2 ( h3 p5 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h9
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh9)
     spin_itm_hp = tdcc_spin_fock(sh3,sp5)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_5_2_1(sh3,sp5,itm_hp)
     call ccsdt_l2_5_2_2(sh3,sp5,itm_hp)

     do h3 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,h9) = i1(h3,h9) + fact * t1inp(p5,h9,spin_t1inp) * itm_hp(h3,p5)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_l2_5_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_5_2_1(sh3,sp5,i2)

!         i2 ( h3 p5 )_f + = 1 * f ( h3 p5 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh3,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h3,p5
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh3,sp5)
  do h3 = 1,norb1
  do p5 = norb1+1,nact
     i2(h3,p5) = i2(h3,p5) + fact * fock(h3,p5,spin_fock)
  end do
  end do
end subroutine ccsdt_l2_5_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_5_2_2(sh3,sp5,i2)

!         i2 ( h3 p5 )_vt + = 1 * Sum ( h8 p7 ) * t ( p7 h8 )_t * v ( h3 h8 p5 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h3,p5
  integer :: h8,p7,sh8,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_int2x = tdcc_spin_int2x(sh3,sh8,sp5,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do p5 = norb1+1,nact
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i2(h3,p5) = i2(h3,p5) + fact * t1inp(p7,h8,spin_t1inp) * int2x(h3,h8,p5,p7,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_5_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_5_3(sh3,sh9,i1)

!     i1 ( h3 h9 )_vt + = 1 * Sum ( h6 p5 ) * t ( p5 h6 )_t * v ( h3 h6 h9 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h9
  integer :: h6,p5,sh6,sp5
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh6 = 1,2
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_int2x = tdcc_spin_int2x(sh3,sh6,sh9,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h9 = 1,norb1
     do h6 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,h9) = i1(h3,h9) + fact * t1inp(p5,h6,spin_t1inp) * int2x(h3,h6,h9,p5,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_5_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_5_4(sh3,sh9,i1)

!     i1 ( h3 h9 )_vt + = -1/2 * Sum ( h8 p5 p6 ) * t ( p5 p6 h8 h9 )_t * v ( h3 h8 p5 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h9
  integer :: h8,p5,p6,sh8,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh3,sh8,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h9 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h9) = i1(h3,h9) + fact * t2inp(p5,p6,h8,h9,spin_t2inp) * int2x(h3,h8,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_l2_5_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_6(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = 1 * P( 2 ) * Sum ( p10 ) * y ( h3 h4 p1 p10 )_y * i1 ( p10 p2 )_f 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_6_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_6_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_6_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p10,sp10
  integer :: spin_g2inp
  integer :: spin_itm_pp
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp10 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp1,sp10)
     spin_itm_pp = tdcc_spin_fock(sp10,sp2)
     if(spin_g2inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_6_1(sp10,sp2,itm_pp)
     call ccsdt_l2_6_2(sp10,sp2,itm_pp)
     call ccsdt_l2_6_3(sp10,sp2,itm_pp)
     call ccsdt_l2_6_4(sp10,sp2,itm_pp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p10 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g2inp(h3,h4,p1,p10,spin_g2inp) * itm_pp(p10,p2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccsdt_l2_6_perm
  !--------------------------------------------
end subroutine ccsdt_l2_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_6_1(sp10,sp1,i1)

!     i1 ( p10 p1 )_f + = 1 * f ( p10 p1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sp10,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p10,p1
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp10,sp1)
  do p10 = norb1+1,nact
  do p1 = norb1+1,nact
     i1(p10,p1) = i1(p10,p1) + fact * fock(p10,p1,spin_fock)
  end do
  end do
end subroutine ccsdt_l2_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_6_2(sp10,sp1,i1)

!     i1 ( p10 p1 )_vt + = -1 * Sum ( h6 p5 ) * t ( p5 h6 )_t * v ( h6 p10 p1 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp10,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p10,p1
  integer :: h6,p5,sh6,sp5
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh6 = 1,2
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_int2x = tdcc_spin_int2x(sh6,sp10,sp1,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p10 = norb1+1,nact
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do p5 = norb1+1,nact
        i1(p10,p1) = i1(p10,p1) + fact * t1inp(p5,h6,spin_t1inp) * int2x(h6,p10,p1,p5,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_6_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_6_3(sp10,sp1,i1)

!     i1 ( p10 p1 )_vt + = 1/2 * Sum ( h7 h8 p6 ) * t ( p6 p10 h7 h8 )_t * v ( h7 h8 p1 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp10,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p10,p1
  integer :: h7,h8,p6,sh7,sh8,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp10,sh7,sh8)
     spin_int2x = tdcc_spin_int2x(sh7,sh8,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do p10 = norb1+1,nact
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i1(p10,p1) = i1(p10,p1) + fact * t2inp(p6,p10,h7,h8,spin_t2inp) * int2x(h7,h8,p1,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_l2_6_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_6_4(sp10,sp1,i1)

!     i1 ( p10 p1 )_vtt + = -1 * Sum ( h6 ) * t ( p10 h6 )_t * i2 ( h6 p1 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sp10,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p10,p1
  integer :: h6,sh6
  integer :: spin_t1inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy1
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sh6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp10,sh6)
     spin_itm_hp = tdcc_spin_dummy1(sh6,sp1)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_6_4_1(sh6,sp1,itm_hp)

     do p10 = norb1+1,nact
     do p1 = norb1+1,nact
     do h6 = 1,norb1
        i1(p10,p1) = i1(p10,p1) + fact * t1inp(p10,h6,spin_t1inp) * itm_hp(h6,p1)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_l2_6_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_6_4_1(sh6,sp1,i2)

!         i2 ( h6 p1 )_vt + = 1 * Sum ( h8 p7 ) * t ( p7 h8 )_t * v ( h6 h8 p1 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh6,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h6,p1
  integer :: h8,p7,sh8,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_int2x = tdcc_spin_int2x(sh6,sh8,sp1,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i2(h6,p1) = i2(h6,p1) + fact * t1inp(p7,h8,spin_t1inp) * int2x(h6,h8,p1,p7,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_6_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_7(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = 1/2 * Sum ( h9 h10 ) * y ( h9 h10 p1 p2 )_y * i1 ( h3 h4 h9 h10 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: h9,h10,sh9,sh10
  integer :: spin_g2inp
  integer :: spin_itm_hhhh
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh10 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh9,sh10,sp1,sp2)
     spin_itm_hhhh = tdcc_spin_int2x(sh3,sh4,sh9,sh10)
     if(spin_g2inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_l2_7_1(sh3,sh4,sh9,sh10,itm_hhhh)
     call ccsdt_l2_7_2(sh3,sh4,sh9,sh10,itm_hhhh)
     call ccsdt_l2_7_3(sh3,sh4,sh9,sh10,itm_hhhh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g2inp(h9,h10,p1,p2,spin_g2inp) * itm_hhhh(h3,h4,h9,h10)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsdt_l2_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_7_1(sh3,sh4,sh9,sh10,i1)

!     i1 ( h3 h4 h9 h10 )_v + = 1 * v ( h3 h4 h9 h10 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sh9,sh10
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h9,h10
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh4,sh9,sh10)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do h9 = 1,norb1
  do h10 = 1,norb1
     i1(h3,h4,h9,h10) = i1(h3,h4,h9,h10) + fact * int2x(h3,h4,h9,h10,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_7_2(sh3,sh4,sh9,sh10,i1)

!     i1 ( h3 h4 h9 h10 )_vt + = -2 * Sum ( p5 ) * t ( p5 h9 )_t * i2 ( h3 h4 h10 p5 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh9,sh10
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h9,h10
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh9)
     spin_itm_hhhp = tdcc_spin_int2x(sh3,sh4,sh10,sp5)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_7_2_1(sh3,sh4,sh10,sp5,itm_hhhp)
     call ccsdt_l2_7_2_2(sh3,sh4,sh10,sp5,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,h4,h9,h10) = i1(h3,h4,h9,h10) + fact * t1inp(p5,h9,spin_t1inp) * itm_hhhp(h3,h4,h10,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l2_7_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_7_2_1(sh3,sh4,sh10,sp5,i2)

!         i2 ( h3 h4 h10 p5 )_v + = 1 * v ( h3 h4 h10 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sh10,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h10,p5
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh4,sh10,sp5)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do h10 = 1,norb1
  do p5 = norb1+1,nact
     i2(h3,h4,h10,p5) = i2(h3,h4,h10,p5) + fact * int2x(h3,h4,h10,p5,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_7_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_7_2_2(sh3,sh4,sh10,sp5,i2)

!         i2 ( h3 h4 h10 p5 )_vt + = -1/2 * Sum ( p7 ) * t ( p7 h10 )_t * v ( h3 h4 p5 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sh10,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h10,p5
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh10)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp5,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h10 = 1,norb1
     do p5 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h3,h4,h10,p5) = i2(h3,h4,h10,p5) + fact * t1inp(p7,h10,spin_t1inp) * int2x(h3,h4,p5,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_7_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_7_3(sh3,sh4,sh9,sh10,i1)

!     i1 ( h3 h4 h9 h10 )_vt + = 1/2 * Sum ( p5 p6 ) * t ( p5 p6 h9 h10 )_t * v ( h3 h4 p5 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sh9,sh10
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h9,h10
  integer :: p5,p6,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh9,sh10)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h4,h9,h10) = i1(h3,h4,h9,h10) + fact * t2inp(p5,p6,h9,h10,spin_t2inp) * int2x(h3,h4,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_7_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_8(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1 * P( 4 ) * Sum ( h9 p7 ) * y ( h3 h9 p1 p7 )_y * i1 ( h4 p7 h9 p2 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_8_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_8_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_8_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_8_perm(sh4,sh3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_8_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h9,p7,sh9,sp7
  integer :: spin_g2inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sp7 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh9,sp1,sp7)
     spin_itm_hphp = tdcc_spin_int2x(sh4,sp7,sh9,sp2)
     if(spin_g2inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_8_1(sh4,sp7,sh9,sp2,itm_hphp)
     call ccsdt_l2_8_2(sh4,sp7,sh9,sp2,itm_hphp)
     call ccsdt_l2_8_3(sh4,sp7,sh9,sp2,itm_hphp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do p7 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g2inp(h3,h9,p1,p7,spin_g2inp) * itm_hphp(h4,p7,h9,p2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hphp)
  end subroutine ccsdt_l2_8_perm
  !--------------------------------------------
end subroutine ccsdt_l2_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_8_1(sh3,sp7,sh9,sp1,i1)

!     i1 ( h3 p7 h9 p1 )_v + = 1 * v ( h3 p7 h9 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sp7,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p7,h9,p1
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sp7,sh9,sp1)
  do h3 = 1,norb1
  do p7 = norb1+1,nact
  do h9 = 1,norb1
  do p1 = norb1+1,nact
     i1(h3,p7,h9,p1) = i1(h3,p7,h9,p1) + fact * int2x(h3,p7,h9,p1,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_8_2(sh3,sp7,sh9,sp1,i1)

!     i1 ( h3 p7 h9 p1 )_vt + = -1 * Sum ( p5 ) * t ( p5 h9 )_t * v ( h3 p7 p1 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sp7,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p7,h9,p1
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh9)
     spin_int2x = tdcc_spin_int2x(sh3,sp7,sp1,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do p7 = norb1+1,nact
     do h9 = 1,norb1
     do p1 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h3,p7,h9,p1) = i1(h3,p7,h9,p1) + fact * t1inp(p5,h9,spin_t1inp) * int2x(h3,p7,p1,p5,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_8_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_8_3(sh3,sp7,sh9,sp1,i1)

!     i1 ( h3 p7 h9 p1 )_vt + = -1 * Sum ( h8 p6 ) * t ( p6 p7 h8 h9 )_t * v ( h3 h8 p1 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sp7,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p7,h9,p1
  integer :: h8,p6,sh8,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp7,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh3,sh8,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do p7 = norb1+1,nact
     do h9 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i1(h3,p7,h9,p1) = i1(h3,p7,h9,p1) + fact * t2inp(p6,p7,h8,h9,spin_t2inp) * int2x(h3,h8,p1,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_8_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_9(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = 1/2 * Sum ( p5 p6 ) * y ( h3 h4 p5 p6 )_y * v ( p5 p6 p1 p2 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: p5,p6,sp5,sp6
  integer :: spin_g2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     spin_int2x = tdcc_spin_int2x(sp5,sp6,sp1,sp2)
     if(spin_g2inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g2inp(h3,h4,p5,p6,spin_g2inp) * int2x(p5,p6,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = 1/2 * P( 2 ) * Sum ( h14 h11 p13 ) * y ( h3 h11 h14 p1 p2 p13 )_y * i1 ( h4 p13 h11 h14 )_v 7

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_10_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_10_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_10_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h14,h11,p13,sh14,sh11,sp13
  integer :: spin_g3inp
  integer :: spin_itm_hphh
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh14 = 1,2
  do sh11 = 1,2
  do sp13 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh3,sh11,sh14,sp1,sp2,sp13)
     spin_itm_hphh = tdcc_spin_int2x(sh4,sp13,sh11,sh14)
     if(spin_g3inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsdt_l2_10_1(sh4,sp13,sh11,sh14,itm_hphh)
     call ccsdt_l2_10_2(sh4,sp13,sh11,sh14,itm_hphh)
     call ccsdt_l2_10_3(sh4,sp13,sh11,sh14,itm_hphh)
     call ccsdt_l2_10_4(sh4,sp13,sh11,sh14,itm_hphh)
     call ccsdt_l2_10_5(sh4,sp13,sh11,sh14,itm_hphh)
     call ccsdt_l2_10_6(sh4,sp13,sh11,sh14,itm_hphh)
     call ccsdt_l2_10_7(sh4,sp13,sh11,sh14,itm_hphh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h14 = 1,norb1
     do h11 = 1,norb1
     do p13 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g3inp(h3,h11,h14,p1,p2,p13,spin_g3inp) * itm_hphh(h4,p13,h11,h14)
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
  deallocate(itm_hphh)
  end subroutine ccsdt_l2_10_perm
  !--------------------------------------------
end subroutine ccsdt_l2_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_1(sh3,sp13,sh11,sh14,i1)

!     i1 ( h3 p13 h11 h14 )_v + = -1 * v ( h3 p13 h11 h14 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sp13,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h3,p13,h11,h14
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sp13,sh11,sh14)
  do h3 = 1,norb1
  do p13 = norb1+1,nact
  do h11 = 1,norb1
  do h14 = 1,norb1
     i1(h3,p13,h11,h14) = i1(h3,p13,h11,h14) + fact * int2x(h3,p13,h11,h14,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_2(sh3,sp13,sh11,sh14,i1)

!     i1 ( h3 p13 h11 h14 )_vt + = 1 * Sum ( h9 ) * t ( p13 h9 )_t * i2 ( h3 h9 h11 h14 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sp13,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h3,p13,h11,h14
  integer :: h9,sh9
  integer :: spin_t1inp
  integer :: spin_itm_hhhh
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp13,sh9)
     spin_itm_hhhh = tdcc_spin_int2x(sh3,sh9,sh11,sh14)
     if(spin_t1inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_l2_10_2_1(sh3,sh9,sh11,sh14,itm_hhhh)
     call ccsdt_l2_10_2_2(sh3,sh9,sh11,sh14,itm_hhhh)
     call ccsdt_l2_10_2_3(sh3,sh9,sh11,sh14,itm_hhhh)

     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h14 = 1,norb1
     do h9 = 1,norb1
        i1(h3,p13,h11,h14) = i1(h3,p13,h11,h14) + fact * t1inp(p13,h9,spin_t1inp) * itm_hhhh(h3,h9,h11,h14)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsdt_l2_10_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_2_1(sh3,sh9,sh11,sh14,i2)

!         i2 ( h3 h9 h11 h14 )_v + = 1 * v ( h3 h9 h11 h14 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh9,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h9,h11,h14
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh9,sh11,sh14)
  do h3 = 1,norb1
  do h9 = 1,norb1
  do h11 = 1,norb1
  do h14 = 1,norb1
     i2(h3,h9,h11,h14) = i2(h3,h9,h11,h14) + fact * int2x(h3,h9,h11,h14,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_10_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_2_2(sh3,sh9,sh11,sh14,i2)

!         i2 ( h3 h9 h11 h14 )_vt + = 2 * Sum ( p5 ) * t ( p5 h14 )_t * i3 ( h3 h9 h11 p5 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh9,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h9,h11,h14
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh14)
     spin_itm_hhhp = tdcc_spin_int2x(sh3,sh9,sh11,sp5)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_10_2_2_1(sh3,sh9,sh11,sp5,itm_hhhp)
     call ccsdt_l2_10_2_2_2(sh3,sh9,sh11,sp5,itm_hhhp)

     do h3 = 1,norb1
     do h9 = 1,norb1
     do h11 = 1,norb1
     do h14 = 1,norb1
     do p5 = norb1+1,nact
        i2(h3,h9,h11,h14) = i2(h3,h9,h11,h14) + fact * t1inp(p5,h14,spin_t1inp) * itm_hhhp(h3,h9,h11,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l2_10_2_2
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_2_2_1(sh3,sh9,sh11,sp5,i3)

!             i3 ( h3 h9 h11 p5 )_v + = 1 * v ( h3 h9 h11 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh9,sh11,sp5
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h9,h11,p5
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh9,sh11,sp5)
  do h3 = 1,norb1
  do h9 = 1,norb1
  do h11 = 1,norb1
  do p5 = norb1+1,nact
     i3(h3,h9,h11,p5) = i3(h3,h9,h11,p5) + fact * int2x(h3,h9,h11,p5,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_10_2_2_1
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_2_2_2(sh3,sh9,sh11,sp5,i3)

!             i3 ( h3 h9 h11 p5 )_vt + = -1/2 * Sum ( p7 ) * t ( p7 h11 )_t * v ( h3 h9 p5 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh9,sh11,sp5
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h9,h11,p5
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh11)
     spin_int2x = tdcc_spin_int2x(sh3,sh9,sp5,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h9 = 1,norb1
     do h11 = 1,norb1
     do p5 = norb1+1,nact
     do p7 = norb1+1,nact
        i3(h3,h9,h11,p5) = i3(h3,h9,h11,p5) + fact * t1inp(p7,h11,spin_t1inp) * int2x(h3,h9,p5,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_10_2_2_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_2_3(sh3,sh9,sh11,sh14,i2)

!         i2 ( h3 h9 h11 h14 )_vt + = 1/2 * Sum ( p5 p6 ) * t ( p5 p6 h11 h14 )_t * v ( h3 h9 p5 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh9,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h9,h11,h14
  integer :: p5,p6,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh11,sh14)
     spin_int2x = tdcc_spin_int2x(sh3,sh9,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h9 = 1,norb1
     do h11 = 1,norb1
     do h14 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i2(h3,h9,h11,h14) = i2(h3,h9,h11,h14) + fact * t2inp(p5,p6,h11,h14,spin_t2inp) * int2x(h3,h9,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_10_2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_3(sh3,sp13,sh11,sh14,i1)

!     i1 ( h3 p13 h11 h14 )_vt + = -2 * Sum ( p5 ) * t ( p5 h14 )_t * i2 ( h3 p13 h11 p5 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sp13,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h3,p13,h11,h14
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh14)
     spin_itm_hphp = tdcc_spin_int2x(sh3,sp13,sh11,sp5)
     if(spin_t1inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_10_3_1(sh3,sp13,sh11,sp5,itm_hphp)
     call ccsdt_l2_10_3_2(sh3,sp13,sh11,sp5,itm_hphp)

     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h14 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,p13,h11,h14) = i1(h3,p13,h11,h14) + fact * t1inp(p5,h14,spin_t1inp) * itm_hphp(h3,p13,h11,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphp)
end subroutine ccsdt_l2_10_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_3_1(sh3,sp13,sh11,sp5,i2)

!         i2 ( h3 p13 h11 p5 )_v + = 1 * v ( h3 p13 h11 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sp13,sh11,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p13,h11,p5
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sp13,sh11,sp5)
  do h3 = 1,norb1
  do p13 = norb1+1,nact
  do h11 = 1,norb1
  do p5 = norb1+1,nact
     i2(h3,p13,h11,p5) = i2(h3,p13,h11,p5) + fact * int2x(h3,p13,h11,p5,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_10_3_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_3_2(sh3,sp13,sh11,sp5,i2)

!         i2 ( h3 p13 h11 p5 )_vt + = -1/2 * Sum ( p7 ) * t ( p7 h11 )_t * v ( h3 p13 p5 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sp13,sh11,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p13,h11,p5
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh11)
     spin_int2x = tdcc_spin_int2x(sh3,sp13,sp5,sp7)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do p5 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h3,p13,h11,p5) = i2(h3,p13,h11,p5) + fact * t1inp(p7,h11,spin_t1inp) * int2x(h3,p13,p5,p7,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_10_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_4(sh3,sp13,sh11,sh14,i1)

!     i1 ( h3 p13 h11 h14 )_ft + = -1 * Sum ( p7 ) * t ( p7 p13 h11 h14 )_t * i2 ( h3 p7 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh3,sp13,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h3,p13,h11,h14
  integer :: p7,sp7
  integer :: spin_t2inp
  integer :: spin_itm_hp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp13,sh11,sh14)
     spin_itm_hp = tdcc_spin_fock(sh3,sp7)
     if(spin_t2inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_10_4_1(sh3,sp7,itm_hp)
     call ccsdt_l2_10_4_2(sh3,sp7,itm_hp)

     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h14 = 1,norb1
     do p7 = norb1+1,nact
        i1(h3,p13,h11,h14) = i1(h3,p13,h11,h14) + fact * t2inp(p7,p13,h11,h14,spin_t2inp) * itm_hp(h3,p7)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsdt_l2_10_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_4_1(sh3,sp7,i2)

!         i2 ( h3 p7 )_f + = 1 * f ( h3 p7 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh3,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h3,p7
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh3,sp7)
  do h3 = 1,norb1
  do p7 = norb1+1,nact
     i2(h3,p7) = i2(h3,p7) + fact * fock(h3,p7,spin_fock)
  end do
  end do
end subroutine ccsdt_l2_10_4_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_4_2(sh3,sp7,i2)

!         i2 ( h3 p7 )_vt + = 1 * Sum ( h10 p9 ) * t ( p9 h10 )_t * v ( h3 h10 p7 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact)
  integer :: h3,p7
  integer :: h10,p9,sh10,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh10 = 1,2
  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
     spin_int2x = tdcc_spin_int2x(sh3,sh10,sp7,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do p7 = norb1+1,nact
     do h10 = 1,norb1
     do p9 = norb1+1,nact
        i2(h3,p7) = i2(h3,p7) + fact * t1inp(p9,h10,spin_t1inp) * int2x(h3,h10,p7,p9,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_10_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_5(sh3,sp13,sh11,sh14,i1)

!     i1 ( h3 p13 h11 h14 )_vt + = -2 * Sum ( h8 p6 ) * t ( p6 p13 h8 h14 )_t * i2 ( h3 h8 h11 p6 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh3,sp13,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h3,p13,h11,h14
  integer :: h8,p6,sh8,sp6
  integer :: spin_t2inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp13,sh8,sh14)
     spin_itm_hhhp = tdcc_spin_int2x(sh3,sh8,sh11,sp6)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_10_5_1(sh3,sh8,sh11,sp6,itm_hhhp)
     call ccsdt_l2_10_5_2(sh3,sh8,sh11,sp6,itm_hhhp)

     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h14 = 1,norb1
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i1(h3,p13,h11,h14) = i1(h3,p13,h11,h14) + fact * t2inp(p6,p13,h8,h14,spin_t2inp) * itm_hhhp(h3,h8,h11,p6)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l2_10_5
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_5_1(sh3,sh8,sh11,sp6,i2)

!         i2 ( h3 h8 h11 p6 )_v + = 1 * v ( h3 h8 h11 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh8,sh11,sp6
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h8,h11,p6
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh8,sh11,sp6)
  do h3 = 1,norb1
  do h8 = 1,norb1
  do h11 = 1,norb1
  do p6 = norb1+1,nact
     i2(h3,h8,h11,p6) = i2(h3,h8,h11,p6) + fact * int2x(h3,h8,h11,p6,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_10_5_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_5_2(sh3,sh8,sh11,sp6,i2)

!         i2 ( h3 h8 h11 p6 )_vt + = -1 * Sum ( p9 ) * t ( p9 h11 )_t * v ( h3 h8 p6 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh8,sh11,sp6
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h8,h11,p6
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh11)
     spin_int2x = tdcc_spin_int2x(sh3,sh8,sp6,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h8 = 1,norb1
     do h11 = 1,norb1
     do p6 = norb1+1,nact
     do p9 = norb1+1,nact
        i2(h3,h8,h11,p6) = i2(h3,h8,h11,p6) + fact * t1inp(p9,h11,spin_t1inp) * int2x(h3,h8,p6,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_10_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_6(sh3,sp13,sh11,sh14,i1)

!     i1 ( h3 p13 h11 h14 )_vt + = -1/2 * Sum ( p5 p6 ) * t ( p5 p6 h11 h14 )_t * v ( h3 p13 p5 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sp13,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h3,p13,h11,h14
  integer :: p5,p6,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh11,sh14)
     spin_int2x = tdcc_spin_int2x(sh3,sp13,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h14 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,p13,h11,h14) = i1(h3,p13,h11,h14) + fact * &
             t2inp(p5,p6,h11,h14,spin_t2inp) * int2x(h3,p13,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_10_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_10_7(sh3,sp13,sh11,sh14,i1)

!     i1 ( h3 p13 h11 h14 )_vt + = 1/2 * Sum ( h10 p6 p7 ) * t ( p6 p7 p13 h10 h11 h14 )_t * v ( h3 h10 p6 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sp13,sh11,sh14
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h3,p13,h11,h14
  integer :: h10,p6,p7,sh10,sp6,sp7
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh10 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp6,sp7,sp13,sh10,sh11,sh14)
     spin_int2x = tdcc_spin_int2x(sh3,sh10,sp6,sp7)
     if(spin_t3inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h14 = 1,norb1
     do h10 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h3,p13,h11,h14) = i1(h3,p13,h11,h14) + fact * &
             t3inp(p6,p7,p13,h10,h11,h14,spin_t3inp) * int2x(h3,h10,p6,p7,spin_int2x)
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
end subroutine ccsdt_l2_10_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = 1/2 * P( 2 ) * Sum ( h12 p13 p11 ) * y ( h3 h4 h12 p1 p11 p13 )_y * i1 ( p11 p13 h12 p2 )_v 7

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_11_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_11_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_11_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h12,p13,p11,sh12,sp13,sp11
  integer :: spin_g3inp
  integer :: spin_itm_pphp
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh12 = 1,2
  do sp13 = 1,2
  do sp11 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh12,sp1,sp11,sp13)
     spin_itm_pphp = tdcc_spin_int2x(sp11,sp13,sh12,sp2)
     if(spin_g3inp * spin_itm_pphp == 0) cycle

     itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_11_1(sp11,sp13,sh12,sp2,itm_pphp)
     call ccsdt_l2_11_2(sp11,sp13,sh12,sp2,itm_pphp)
     call ccsdt_l2_11_3(sp11,sp13,sh12,sp2,itm_pphp)
     call ccsdt_l2_11_4(sp11,sp13,sh12,sp2,itm_pphp)
     call ccsdt_l2_11_5(sp11,sp13,sh12,sp2,itm_pphp)
     call ccsdt_l2_11_6(sp11,sp13,sh12,sp2,itm_pphp)
     call ccsdt_l2_11_7(sp11,sp13,sh12,sp2,itm_pphp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h12 = 1,norb1
     do p13 = norb1+1,nact
     do p11 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             g3inp(h3,h4,h12,p1,p11,p13,spin_g3inp) * itm_pphp(p11,p13,h12,p2)
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
  deallocate(itm_pphp)
  end subroutine ccsdt_l2_11_perm
  !--------------------------------------------
end subroutine ccsdt_l2_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_1(sp11,sp13,sh12,sp1,i1)

!     i1 ( p11 p13 h12 p1 )_v + = -1 * v ( p11 p13 h12 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sp11,sp13,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p11,p13,h12,p1
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sp11,sp13,sh12,sp1)
  do p11 = norb1+1,nact
  do p13 = norb1+1,nact
  do h12 = 1,norb1
  do p1 = norb1+1,nact
     i1(p11,p13,h12,p1) = i1(p11,p13,h12,p1) + fact * int2x(p11,p13,h12,p1,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_2(sp11,sp13,sh12,sp1,i1)

!     i1 ( p11 p13 h12 p1 )_vt + = 1 * Sum ( p5 ) * t ( p5 h12 )_t * v ( p11 p13 p1 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp11,sp13,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p11,p13,h12,p1
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh12)
     spin_int2x = tdcc_spin_int2x(sp11,sp13,sp1,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do p11 = norb1+1,nact
     do p13 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(p11,p13,h12,p1) = i1(p11,p13,h12,p1) + fact * t1inp(p5,h12,spin_t1inp) * int2x(p11,p13,p1,p5,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_11_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_3(sp11,sp13,sh12,sp1,i1)

!     i1 ( p11 p13 h12 p1 )_vt + = -1/2 * Sum ( h7 h8 ) * t ( p11 p13 h7 h8 )_t * i2 ( h7 h8 h12 p1 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp11,sp13,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p11,p13,h12,p1
  integer :: h7,h8,sh7,sh8
  integer :: spin_t2inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh7 = 1,2
  do sh8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp11,sp13,sh7,sh8)
     spin_itm_hhhp = tdcc_spin_int2x(sh7,sh8,sh12,sp1)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_11_3_1(sh7,sh8,sh12,sp1,itm_hhhp)
     call ccsdt_l2_11_3_2(sh7,sh8,sh12,sp1,itm_hhhp)

     do p11 = norb1+1,nact
     do p13 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
        i1(p11,p13,h12,p1) = i1(p11,p13,h12,p1) + fact * t2inp(p11,p13,h7,h8,spin_t2inp) * itm_hhhp(h7,h8,h12,p1)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l2_11_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_3_1(sh7,sh8,sh12,sp1,i2)

!         i2 ( h7 h8 h12 p1 )_v + = 1 * v ( h7 h8 h12 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh7,sh8,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h7,h8,h12,p1
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh7,sh8,sh12,sp1)
  do h7 = 1,norb1
  do h8 = 1,norb1
  do h12 = 1,norb1
  do p1 = norb1+1,nact
     i2(h7,h8,h12,p1) = i2(h7,h8,h12,p1) + fact * int2x(h7,h8,h12,p1,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsdt_l2_11_3_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_3_2(sh7,sh8,sh12,sp1,i2)

!         i2 ( h7 h8 h12 p1 )_vt + = -1 * Sum ( p9 ) * t ( p9 h12 )_t * v ( h7 h8 p1 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh7,sh8,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h7,h8,h12,p1
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh12)
     spin_int2x = tdcc_spin_int2x(sh7,sh8,sp1,sp9)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h7 = 1,norb1
     do h8 = 1,norb1
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do p9 = norb1+1,nact
        i2(h7,h8,h12,p1) = i2(h7,h8,h12,p1) + fact * t1inp(p9,h12,spin_t1inp) * int2x(h7,h8,p1,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_11_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_4(sp11,sp13,sh12,sp1,i1)

!     i1 ( p11 p13 h12 p1 )_vt + = -2 * Sum ( h8 p6 ) * t ( p6 p13 h8 h12 )_t * v ( h8 p11 p1 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp11,sp13,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p11,p13,h12,p1
  integer :: h8,p6,sh8,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp13,sh8,sh12)
     spin_int2x = tdcc_spin_int2x(sh8,sp11,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do p11 = norb1+1,nact
     do p13 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i1(p11,p13,h12,p1) = i1(p11,p13,h12,p1) + fact * &
             t2inp(p6,p13,h8,h12,spin_t2inp) * int2x(h8,p11,p1,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_11_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_5(sp11,sp13,sh12,sp1,i1)

!     i1 ( p11 p13 h12 p1 )_vt + = 1/2 * Sum ( h9 h10 p7 ) * t ( p7 p11 p13 h9 h10 h12 )_t * v ( h9 h10 p1 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp11,sp13,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p11,p13,h12,p1
  integer :: h9,h10,p7,sh9,sh10,sp7
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh9 = 1,2
  do sh10 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp7,sp11,sp13,sh9,sh10,sh12)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp1,sp7)
     if(spin_t3inp * spin_int2x == 0) cycle

     do p11 = norb1+1,nact
     do p13 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p7 = norb1+1,nact
        i1(p11,p13,h12,p1) = i1(p11,p13,h12,p1) + fact * &
             t3inp(p7,p11,p13,h9,h10,h12,spin_t3inp) * int2x(h9,h10,p1,p7,spin_int2x)
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
end subroutine ccsdt_l2_11_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_6(sp11,sp13,sh12,sp1,i1)

!     i1 ( p11 p13 h12 p1 )_vtt + = -2 * Sum ( h8 ) * t ( p13 h8 )_t * i2 ( h8 p11 h12 p1 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sp11,sp13,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p11,p13,h12,p1
  integer :: h8,sh8
  integer :: spin_t1inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp13,sh8)
     spin_itm_hphp = tdcc_spin_dummy2(sh8,sp11,sh12,sp1)
     if(spin_t1inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_11_6_1(sh8,sp11,sh12,sp1,itm_hphp)

     do p11 = norb1+1,nact
     do p13 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
        i1(p11,p13,h12,p1) = i1(p11,p13,h12,p1) + fact * &
             t1inp(p13,h8,spin_t1inp) * itm_hphp(h8,p11,h12,p1)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphp)
end subroutine ccsdt_l2_11_6
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_6_1(sh8,sp11,sh12,sp1,i2)

!         i2 ( h8 p11 h12 p1 )_vt + = -1 * Sum ( p5 ) * t ( p5 h12 )_t * v ( h8 p11 p1 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh8,sp11,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h8,p11,h12,p1
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh12)
     spin_int2x = tdcc_spin_int2x(sh8,sp11,sp1,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h8 = 1,norb1
     do p11 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do p5 = norb1+1,nact
        i2(h8,p11,h12,p1) = i2(h8,p11,h12,p1) + fact * &
             t1inp(p5,h12,spin_t1inp) * int2x(h8,p11,p1,p5,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_11_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_7(sp11,sp13,sh12,sp1,i1)

!     i1 ( p11 p13 h12 p1 )_vtt + = 2 * Sum ( h10 ) * t ( p11 h10 )_t * i2 ( h10 p13 h12 p1 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sp11,sp13,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p11,p13,h12,p1
  integer :: h10,sh10
  integer :: spin_t1inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 2.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh10 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp11,sh10)
     spin_itm_hphp = tdcc_spin_dummy2(sh10,sp13,sh12,sp1)
     if(spin_t1inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_11_7_1(sh10,sp13,sh12,sp1,itm_hphp)

     do p11 = norb1+1,nact
     do p13 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h10 = 1,norb1
        i1(p11,p13,h12,p1) = i1(p11,p13,h12,p1) + fact * &
             t1inp(p11,h10,spin_t1inp) * itm_hphp(h10,p13,h12,p1)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphp)
end subroutine ccsdt_l2_11_7
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_11_7_1(sh10,sp13,sh12,sp1,i2)

!         i2 ( h10 p13 h12 p1 )_vt + = 1 * Sum ( h8 p6 ) * t ( p6 p13 h8 h12 )_t * v ( h8 h10 p1 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh10,sp13,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h10,p13,h12,p1
  integer :: h8,p6,sh8,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp13,sh8,sh12)
     spin_int2x = tdcc_spin_int2x(sh8,sh10,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h10 = 1,norb1
     do p13 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i2(h10,p13,h12,p1) = i2(h10,p13,h12,p1) + fact * &
             t2inp(p6,p13,h8,h12,spin_t2inp) * int2x(h8,h10,p1,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_11_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_12(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1 * P( 2 ) * Sum ( h11 ) * i1 ( h3 h11 )_yt * v ( h4 h11 p1 p2 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_12_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_12_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_12_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h11,sh11
  integer :: spin_itm_hh
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh11 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh3,sh11)
     spin_int2x = tdcc_spin_int2x(sh4,sh11,sp1,sp2)
     if(spin_itm_hh * spin_int2x == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsdt_l2_12_1(sh3,sh11,itm_hh)
     call ccsdt_l2_12_2(sh3,sh11,itm_hh)
     call ccsdt_l2_12_3(sh3,sh11,itm_hh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h11 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hh(h3,h11) * int2x(h4,h11,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccsdt_l2_12_perm
  !--------------------------------------------
end subroutine ccsdt_l2_12
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_12_1(sh3,sh11,i1)

!     i1 ( h3 h11 )_yt + = 1 * Sum ( p5 ) * t ( p5 h11 )_t * y ( h3 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g1inp

  implicit none
  integer,intent(in) :: sh3,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h11
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_g1inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh11)
     spin_g1inp = tdcc_spin_g1inp(sh3,sp5)
     if(spin_t1inp * spin_g1inp == 0) cycle

     do h3 = 1,norb1
     do h11 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,h11) = i1(h3,h11) + fact * t1inp(p5,h11,spin_t1inp) * g1inp(h3,p5,spin_g1inp)
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_12_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_12_2(sh3,sh11,i1)

!     i1 ( h3 h11 )_yt + = -1/2 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h11 )_t * y ( h3 h7 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h11
  integer :: h7,p5,p6,sh7,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh11)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh7,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h11 = 1,norb1
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h11) = i1(h3,h11) + fact * &
             t2inp(p5,p6,h7,h11,spin_t2inp) * g2inp(h3,h7,p5,p6,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_l2_12_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_12_3(sh3,sh11,i1)

!     i1 ( h3 h11 )_yt + = 1/12 * Sum ( h8 h9 p5 p6 p7 ) * t ( p5 p6 p7 h8 h9 h11 )_t * y ( h3 h8 h9 p5 p6 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h11
  integer :: h8,h9,p5,p6,p7,sh8,sh9,sp5,sp6,sp7
  integer :: spin_t3inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 12.0d+0 * runit

  do sh8 = 1,2
  do sh9 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp7,sh8,sh9,sh11)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh8,sh9,sp5,sp6,sp7)
     if(spin_t3inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h11 = 1,norb1
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h3,h11) = i1(h3,h11) + fact * &
             t3inp(p5,p6,p7,h8,h9,h11,spin_t3inp) * g3inp(h3,h8,h9,p5,p6,p7,spin_g3inp)
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
  end do
  end do
end subroutine ccsdt_l2_12_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_13(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytf + = 1 * P( 2 ) * Sum ( h5 ) * i1 ( h3 h4 h5 p1 )_yt * f ( h5 p2 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_13_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_13_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_13_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h5,sh5
  integer :: spin_itm_hhhp
  integer :: spin_fock
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh5 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh5,sp1)
     spin_fock = tdcc_spin_fock(sh5,sp2)
     if(spin_itm_hhhp * spin_fock == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_13_1(sh3,sh4,sh5,sp1,itm_hhhp)
     call ccsdt_l2_13_2(sh3,sh4,sh5,sp1,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h5 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hhhp(h3,h4,h5,p1) * fock(h5,p2,spin_fock)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_l2_13_perm
  !--------------------------------------------
end subroutine ccsdt_l2_13
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_13_1(sh3,sh4,sh5,sp1,i1)

!     i1 ( h3 h4 h5 p1 )_yt + = -1 * Sum ( p6 ) * t ( p6 h5 )_t * y ( h3 h4 p1 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh5,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h5,p1
  integer :: p6,sp6
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp6,sh5)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp1,sp6)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p1 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h4,h5,p1) = i1(h3,h4,h5,p1) + fact * &
             t1inp(p6,h5,spin_t1inp) * g2inp(h3,h4,p1,p6,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_13_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_13_2(sh3,sh4,sh5,sp1,i1)

!     i1 ( h3 h4 h5 p1 )_yt + = -1/2 * Sum ( h8 p6 p7 ) * t ( p6 p7 h5 h8 )_t * y ( h3 h4 h8 p1 p6 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh5,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h5,p1
  integer :: h8,p6,p7,sh8,sp6,sp7
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp7,sh5,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh8,sp1,sp6,sp7)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h3,h4,h5,p1) = i1(h3,h4,h5,p1) + fact * &
             t2inp(p6,p7,h5,h8,spin_t2inp) * g3inp(h3,h4,h8,p1,p6,p7,spin_g3inp)
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
end subroutine ccsdt_l2_13_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_14(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1 * P( 4 ) * Sum ( h8 h10 ) * i1 ( h3 h8 h10 p1 )_yt * v ( h4 h10 h8 p2 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_14_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_14_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_14_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_14_perm(sh4,sh3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_14_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h8,h10,sh8,sh10
  integer :: spin_itm_hhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh8 = 1,2
  do sh10 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh8,sh10,sp1)
     spin_int2x = tdcc_spin_int2x(sh4,sh10,sh8,sp2)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_14_1(sh3,sh8,sh10,sp1,itm_hhhp)
     call ccsdt_l2_14_2(sh3,sh8,sh10,sp1,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h8 = 1,norb1
     do h10 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hhhp(h3,h8,h10,p1) * int2x(h4,h10,h8,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsdt_l2_14_perm
  !--------------------------------------------
end subroutine ccsdt_l2_14
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_14_1(sh3,sh8,sh10,sp1,i1)

!     i1 ( h3 h8 h10 p1 )_yt + = 1 * Sum ( p5 ) * t ( p5 h10 )_t * y ( h3 h8 p1 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh8,sh10,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h8,h10,p1
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh10)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh8,sp1,sp5)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h8 = 1,norb1
     do h10 = 1,norb1
     do p1 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h3,h8,h10,p1) = i1(h3,h8,h10,p1) + fact * &
             t1inp(p5,h10,spin_t1inp) * g2inp(h3,h8,p1,p5,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_14_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_14_2(sh3,sh8,sh10,sp1,i1)

!     i1 ( h3 h8 h10 p1 )_yt + = 1/2 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h10 )_t * y ( h3 h7 h8 p1 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh8,sh10,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h8,h10,p1
  integer :: h7,p5,p6,sh7,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh7,sh8,sp1,sp5,sp6)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h8 = 1,norb1
     do h10 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h8,h10,p1) = i1(h3,h8,h10,p1) + fact * &
             t2inp(p5,p6,h7,h10,spin_t2inp) * g3inp(h3,h7,h8,p1,p5,p6,spin_g3inp)
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
end subroutine ccsdt_l2_14_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_15(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1 * Sum ( h10 p8 ) * i1 ( h3 h4 h10 p8 )_yt * v ( h10 p8 p1 p2 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: h10,p8,sh10,sp8
  integer :: spin_itm_hhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh10 = 1,2
  do sp8 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh10,sp8)
     spin_int2x = tdcc_spin_int2x(sh10,sp8,sp1,sp2)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_15_1(sh3,sh4,sh10,sp8,itm_hhhp)
     call ccsdt_l2_15_2(sh3,sh4,sh10,sp8,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h10 = 1,norb1
     do p8 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hhhp(h3,h4,h10,p8) * int2x(h10,p8,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l2_15
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_15_1(sh3,sh4,sh10,sp8,i1)

!     i1 ( h3 h4 h10 p8 )_yt + = -1 * Sum ( p5 ) * t ( p5 h10 )_t * y ( h3 h4 p5 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh10,sp8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h10,p8
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh10)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp8)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h10 = 1,norb1
     do p8 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h3,h4,h10,p8) = i1(h3,h4,h10,p8) + fact * t1inp(p5,h10,spin_t1inp) * g2inp(h3,h4,p5,p8,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_15_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_15_2(sh3,sh4,sh10,sp8,i1)

!     i1 ( h3 h4 h10 p8 )_yt + = -1/2 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h10 )_t * y ( h3 h4 h7 p5 p6 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh10,sp8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h10,p8
  integer :: h7,p5,p6,sh7,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh7,sp5,sp6,sp8)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h10 = 1,norb1
     do p8 = norb1+1,nact
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h4,h10,p8) = i1(h3,h4,h10,p8) + fact * &
             t2inp(p5,p6,h7,h10,spin_t2inp) * g3inp(h3,h4,h7,p5,p6,p8,spin_g3inp)
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
end subroutine ccsdt_l2_15_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_16(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1 * P( 2 ) * Sum ( h8 h6 p7 ) * i1 ( h3 h4 h8 h6 p1 p7 )_yt * v ( h6 p7 h8 p2 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_16_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_16_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_16_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h8,h6,p7,sh8,sh6,sp7
  integer :: spin_itm_hhhhpp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy3
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhhpp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sh8 = 1,2
  do sh6 = 1,2
  do sp7 = 1,2
     spin_itm_hhhhpp = tdcc_spin_dummy3(sh3,sh4,sh8,sh6,sp1,sp7)
     spin_int2x = tdcc_spin_int2x(sh6,sp7,sh8,sp2)
     if(spin_itm_hhhhpp * spin_int2x == 0) cycle

     itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_16_1(sh3,sh4,sh8,sh6,sp1,sp7,itm_hhhhpp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h8 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hhhhpp(h3,h4,h8,h6,p1,p7) * int2x(h6,p7,h8,p2,spin_int2x)
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
  end subroutine ccsdt_l2_16_perm
  !--------------------------------------------
end subroutine ccsdt_l2_16
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_16_1(sh3,sh4,sh8,sh6,sp1,sp7,i1)

!     i1 ( h3 h4 h8 h6 p1 p7 )_yt + = -1 * Sum ( p5 ) * t ( p5 h6 )_t * y ( h3 h4 h8 p1 p5 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh8,sh6,sp1,sp7
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,h8,h6,p1,p7
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh8,sp1,sp5,sp7)
     if(spin_t1inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h8 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h3,h4,h8,h6,p1,p7) = i1(h3,h4,h8,h6,p1,p7) + fact * &
             t1inp(p5,h6,spin_t1inp) * g3inp(h3,h4,h8,p1,p5,p7,spin_g3inp)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_16_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_17(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/2 * P( 2 ) * Sum ( p11 ) * i1 ( p11 p1 )_yt * v ( h3 h4 p2 p11 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_17_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_17_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_17_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p11,sp11
  integer :: spin_itm_pp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp11 = 1,2
     spin_itm_pp = tdcc_spin_dummy1(sp11,sp1)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp2,sp11)
     if(spin_itm_pp * spin_int2x == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_17_1(sp11,sp1,itm_pp)
     call ccsdt_l2_17_2(sp11,sp1,itm_pp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p11 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_pp(p11,p1) * int2x(h3,h4,p2,p11,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccsdt_l2_17_perm
  !--------------------------------------------
end subroutine ccsdt_l2_17
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_17_1(sp11,sp1,i1)

!     i1 ( p11 p1 )_yt + = -1 * Sum ( h7 h8 p5 ) * t ( p5 p11 h7 h8 )_t * y ( h7 h8 p1 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sp11,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p11,p1
  integer :: h7,h8,p5,sh7,sh8,sp5
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp11,sh7,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh7,sh8,sp1,sp5)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do p11 = norb1+1,nact
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(p11,p1) = i1(p11,p1) + fact * t2inp(p5,p11,h7,h8,spin_t2inp) * g2inp(h7,h8,p1,p5,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsdt_l2_17_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_17_2(sp11,sp1,i1)

!     i1 ( p11 p1 )_yt + = 1/6 * Sum ( h8 h9 h10 p5 p6 ) * t ( p5 p6 p11 h8 h9 h10 )_t * y ( h8 h9 h10 p1 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sp11,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p11,p1
  integer :: h8,h9,h10,p5,p6,sh8,sh9,sh10,sp5,sp6
  integer :: spin_t3inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 6.0d+0 * runit

  do sh8 = 1,2
  do sh9 = 1,2
  do sh10 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp11,sh8,sh9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh8,sh9,sh10,sp1,sp5,sp6)
     if(spin_t3inp * spin_g3inp == 0) cycle

     do p11 = norb1+1,nact
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(p11,p1) = i1(p11,p1) + fact * t3inp(p5,p6,p11,h8,h9,h10,spin_t3inp) * g3inp(h8,h9,h10,p1,p5,p6,spin_g3inp)
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
  end do
  end do
end subroutine ccsdt_l2_17_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_18(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/4 * Sum ( h12 h13 ) * i1 ( h3 h4 h12 h13 )_yt * v ( h12 h13 p1 p2 )_v 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: h12,h13,sh12,sh13
  integer :: spin_itm_hhhh
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh12 = 1,2
  do sh13 = 1,2
     spin_itm_hhhh = tdcc_spin_dummy2(sh3,sh4,sh12,sh13)
     spin_int2x = tdcc_spin_int2x(sh12,sh13,sp1,sp2)
     if(spin_itm_hhhh * spin_int2x == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_l2_18_1(sh3,sh4,sh12,sh13,itm_hhhh)
     call ccsdt_l2_18_2(sh3,sh4,sh12,sh13,itm_hhhh)
     call ccsdt_l2_18_3(sh3,sh4,sh12,sh13,itm_hhhh)
     call ccsdt_l2_18_4(sh3,sh4,sh12,sh13,itm_hhhh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h12 = 1,norb1
     do h13 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hhhh(h3,h4,h12,h13) * int2x(h12,h13,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsdt_l2_18
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_18_1(sh3,sh4,sh12,sh13,i1)

!     i1 ( h3 h4 h12 h13 )_yt + = 1 * Sum ( p5 p6 ) * t ( p5 p6 h12 h13 )_t * y ( h3 h4 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh12,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h12,h13
  integer :: p5,p6,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh12,sh13)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h12 = 1,norb1
     do h13 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h4,h12,h13) = i1(h3,h4,h12,h13) + fact * t2inp(p5,p6,h12,h13,spin_t2inp) * g2inp(h3,h4,p5,p6,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsdt_l2_18_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_18_2(sh3,sh4,sh12,sh13,i1)

!     i1 ( h3 h4 h12 h13 )_yt + = 1/3 * Sum ( h8 p5 p6 p7 ) * t ( p5 p6 p7 h8 h12 h13 )_t * y ( h3 h4 h8 p5 p6 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh12,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h12,h13
  integer :: h8,p5,p6,p7,sh8,sp5,sp6,sp7
  integer :: spin_t3inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 3.0d+0 * runit

  do sh8 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp7,sh8,sh12,sh13)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh8,sp5,sp6,sp7)
     if(spin_t3inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h12 = 1,norb1
     do h13 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h3,h4,h12,h13) = i1(h3,h4,h12,h13) + fact * &
             t3inp(p5,p6,p7,h8,h12,h13,spin_t3inp) * g3inp(h3,h4,h8,p5,p6,p7,spin_g3inp)
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
  end do
  end do
end subroutine ccsdt_l2_18_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_18_3(sh3,sh4,sh12,sh13,i1)

!     i1 ( h3 h4 h12 h13 )_ytt + = -2 * Sum ( p5 ) * t ( p5 h12 )_t * i2 ( h3 h4 h13 p5 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh12,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h12,h13
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh12)
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh13,sp5)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_18_3_1(sh3,sh4,sh13,sp5,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h12 = 1,norb1
     do h13 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,h4,h12,h13) = i1(h3,h4,h12,h13) + fact * t1inp(p5,h12,spin_t1inp) * itm_hhhp(h3,h4,h13,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l2_18_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_18_3_1(sh3,sh4,sh13,sp5,i2)

!         i2 ( h3 h4 h13 p5 )_yt + = -1 * Sum ( p7 ) * t ( p7 h13 )_t * y ( h3 h4 p5 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh13,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h13,p5
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh13)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h13 = 1,norb1
     do p5 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h3,h4,h13,p5) = i2(h3,h4,h13,p5) + fact * t1inp(p7,h13,spin_t1inp) * g2inp(h3,h4,p5,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_18_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_18_4(sh3,sh4,sh12,sh13,i1)

!     i1 ( h3 h4 h12 h13 )_ytt + = 2 * Sum ( p9 ) * t ( p9 h13 )_t * i2 ( h3 h4 h12 p9 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh12,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h12,h13
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh13)
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh12,sp9)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_18_4_1(sh3,sh4,sh12,sp9,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h12 = 1,norb1
     do h13 = 1,norb1
     do p9 = norb1+1,nact
        i1(h3,h4,h12,h13) = i1(h3,h4,h12,h13) + fact * t1inp(p9,h13,spin_t1inp) * itm_hhhp(h3,h4,h12,p9)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l2_18_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_18_4_1(sh3,sh4,sh12,sp9,i2)

!         i2 ( h3 h4 h12 p9 )_yt + = 1 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h12 )_t * y ( h3 h4 h7 p5 p6 p9 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh12,sp9
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h12,p9
  integer :: h7,p5,p6,sh7,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh12)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh7,sp5,sp6,sp9)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h12 = 1,norb1
     do p9 = norb1+1,nact
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i2(h3,h4,h12,p9) = i2(h3,h4,h12,p9) + fact * &
             t2inp(p5,p6,h7,h12,spin_t2inp) * g3inp(h3,h4,h7,p5,p6,p9,spin_g3inp)
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
end subroutine ccsdt_l2_18_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_19(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1/2 * Sum ( h9 p6 ) * i1 ( h9 p6 p1 p2 )_yt * v ( h3 h4 h9 p6 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: h9,p6,sh9,sp6
  integer :: spin_itm_hppp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sh9 = 1,2
  do sp6 = 1,2
     spin_itm_hppp = tdcc_spin_dummy2(sh9,sp6,sp1,sp2)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sh9,sp6)
     if(spin_itm_hppp * spin_int2x == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_19_1(sh9,sp6,sp1,sp2,itm_hppp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do p6 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hppp(h9,p6,p1,p2) * int2x(h3,h4,h9,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hppp)
end subroutine ccsdt_l2_19
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_19_1(sh9,sp6,sp1,sp2,i1)

!     i1 ( h9 p6 p1 p2 )_yt + = 1 * Sum ( h7 h8 p5 ) * t ( p5 p6 h7 h8 )_t * y ( h7 h8 h9 p1 p2 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh9,sp6,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h9,p6,p1,p2
  integer :: h7,h8,p5,sh7,sh8,sp5
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh7,sh8,sh9,sp1,sp2,sp5)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h9 = 1,norb1
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(h9,p6,p1,p2) = i1(h9,p6,p1,p2) + fact * &
             t2inp(p5,p6,h7,h8,spin_t2inp) * g3inp(h7,h8,h9,p1,p2,p5,spin_g3inp)
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
end subroutine ccsdt_l2_19_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_20(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1/2 * P( 4 ) * Sum ( p6 p9 ) * i1 ( h3 p6 p1 p9 )_yt * v ( h4 p9 p2 p6 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_20_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_20_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_20_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_20_perm(sh4,sh3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_20_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p6,p9,sp6,sp9
  integer :: spin_itm_hppp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp6 = 1,2
  do sp9 = 1,2
     spin_itm_hppp = tdcc_spin_dummy2(sh3,sp6,sp1,sp9)
     spin_int2x = tdcc_spin_int2x(sh4,sp9,sp2,sp6)
     if(spin_itm_hppp * spin_int2x == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_20_1(sh3,sp6,sp1,sp9,itm_hppp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p6 = norb1+1,nact
     do p9 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hppp(h3,p6,p1,p9) * int2x(h4,p9,p2,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hppp)
  end subroutine ccsdt_l2_20_perm
  !--------------------------------------------
end subroutine ccsdt_l2_20
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_20_1(sh3,sp6,sp1,sp9,i1)

!     i1 ( h3 p6 p1 p9 )_yt + = -1 * Sum ( h7 h8 p5 ) * t ( p5 p6 h7 h8 )_t * y ( h3 h7 h8 p1 p5 p9 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sp6,sp1,sp9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,p6,p1,p9
  integer :: h7,h8,p5,sh7,sh8,sp5
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh7,sh8,sp1,sp5,sp9)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do p9 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,p6,p1,p9) = i1(h3,p6,p1,p9) + fact * &
             t2inp(p5,p6,h7,h8,spin_t2inp) * g3inp(h3,h7,h8,p1,p5,p9,spin_g3inp)
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
end subroutine ccsdt_l2_20_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_21(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_vty + = 1/12 * Sum ( p5 h8 h9 h10 ) * y ( h8 h9 h10 p1 p2 p5 )_y * i1 ( h3 h4 p5 h8 h9 h10 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: p5,h8,h9,h10,sp5,sh8,sh9,sh10
  integer :: spin_g3inp
  integer :: spin_itm_hhphhh
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhphhh(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 12.0d+0 * runit

  allocate(itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1))
  do sp5 = 1,2
  do sh8 = 1,2
  do sh9 = 1,2
  do sh10 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh8,sh9,sh10,sp1,sp2,sp5)
     spin_itm_hhphhh = tdcc_spin_dummy3(sh3,sh4,sp5,sh8,sh9,sh10)
     if(spin_g3inp * spin_itm_hhphhh == 0) cycle

     itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccsdt_l2_21_1(sh3,sh4,sp5,sh8,sh9,sh10,itm_hhphhh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p5 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             g3inp(h8,h9,h10,p1,p2,p5,spin_g3inp) * itm_hhphhh(h3,h4,p5,h8,h9,h10)
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
  end do
  end do
  deallocate(itm_hhphhh)
end subroutine ccsdt_l2_21
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_21_1(sh3,sh4,sp5,sh8,sh9,sh10,i1)

!     i1 ( h3 h4 p5 h8 h9 h10 )_vt + = 1 * Sum ( p6 p7 ) * t ( p5 p6 p7 h8 h9 h10 )_t * v ( h3 h4 p6 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp5,sh8,sh9,sh10
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,p5,h8,h9,h10
  integer :: p6,p7,sp6,sp7
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp6 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp7,sh8,sh9,sh10)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp6,sp7)
     if(spin_t3inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h3,h4,p5,h8,h9,h10) = i1(h3,h4,p5,h8,h9,h10) + fact * &
             t3inp(p5,p6,p7,h8,h9,h10,spin_t3inp) * int2x(h3,h4,p6,p7,spin_int2x)
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
end subroutine ccsdt_l2_21_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_22(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/4 * P( 4 ) * Sum ( p12 h13 ) * i1 ( h3 p12 h13 p1 )_yt * v ( h4 h13 p2 p12 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_22_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_22_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_22_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_22_perm(sh4,sh3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_22_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p12,h13,sp12,sh13
  integer :: spin_itm_hphp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp12 = 1,2
  do sh13 = 1,2
     spin_itm_hphp = tdcc_spin_dummy2(sh3,sp12,sh13,sp1)
     spin_int2x = tdcc_spin_int2x(sh4,sh13,sp2,sp12)
     if(spin_itm_hphp * spin_int2x == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_22_1(sh3,sp12,sh13,sp1,itm_hphp)
     call ccsdt_l2_22_2(sh3,sp12,sh13,sp1,itm_hphp)
     call ccsdt_l2_22_3(sh3,sp12,sh13,sp1,itm_hphp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p12 = norb1+1,nact
     do h13 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hphp(h3,p12,h13,p1) * int2x(h4,h13,p2,p12,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hphp)
  end subroutine ccsdt_l2_22_perm
  !--------------------------------------------
end subroutine ccsdt_l2_22
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_22_1(sh3,sp12,sh13,sp1,i1)

!     i1 ( h3 p12 h13 p1 )_yt + = 1 * Sum ( h8 h9 p5 p6 ) * t ( p5 p6 p12 h8 h9 h13 )_t * y ( h3 h8 h9 p1 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sp12,sh13,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p12,h13,p1
  integer :: h8,h9,p5,p6,sh8,sh9,sp5,sp6
  integer :: spin_t3inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sh9 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp12,sh8,sh9,sh13)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh8,sh9,sp1,sp5,sp6)
     if(spin_t3inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do p12 = norb1+1,nact
     do h13 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,p12,h13,p1) = i1(h3,p12,h13,p1) + fact * &
             t3inp(p5,p6,p12,h8,h9,h13,spin_t3inp) * g3inp(h3,h8,h9,p1,p5,p6,spin_g3inp)
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
  end do
  end do
end subroutine ccsdt_l2_22_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_22_2(sh3,sp12,sh13,sp1,i1)

!     i1 ( h3 p12 h13 p1 )_ytt + = -4 * Sum ( h8 ) * t ( p12 h8 )_t * i2 ( h3 h8 h13 p1 )_yt 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sp12,sh13,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p12,h13,p1
  integer :: h8,sh8
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -4.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh8 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp12,sh8)
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh8,sh13,sp1)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_22_2_1(sh3,sh8,sh13,sp1,itm_hhhp)
     call ccsdt_l2_22_2_2(sh3,sh8,sh13,sp1,itm_hhhp)

     do h3 = 1,norb1
     do p12 = norb1+1,nact
     do h13 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
        i1(h3,p12,h13,p1) = i1(h3,p12,h13,p1) + fact * &
             t1inp(p12,h8,spin_t1inp) * itm_hhhp(h3,h8,h13,p1)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l2_22_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_22_2_1(sh3,sh8,sh13,sp1,i2)

!         i2 ( h3 h8 h13 p1 )_yt + = 1 * Sum ( p7 ) * t ( p7 h13 )_t * y ( h3 h8 p1 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh8,sh13,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h8,h13,p1
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh13)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh8,sp1,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h8 = 1,norb1
     do h13 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h3,h8,h13,p1) = i2(h3,h8,h13,p1) + fact * &
             t1inp(p7,h13,spin_t1inp) * g2inp(h3,h8,p1,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_22_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_22_2_2(sh3,sh8,sh13,sp1,i2)

!         i2 ( h3 h8 h13 p1 )_yt + = 1/2 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h13 )_t * y ( h3 h7 h8 p1 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh8,sh13,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h8,h13,p1
  integer :: h7,p5,p6,sh7,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh13)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh7,sh8,sp1,sp5,sp6)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h8 = 1,norb1
     do h13 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i2(h3,h8,h13,p1) = i2(h3,h8,h13,p1) + fact * &
             t2inp(p5,p6,h7,h13,spin_t2inp) * g3inp(h3,h7,h8,p1,p5,p6,spin_g3inp)
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
end subroutine ccsdt_l2_22_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_22_3(sh3,sp12,sh13,sp1,i1)

!     i1 ( h3 p12 h13 p1 )_ytt + = 2 * Sum ( h7 h8 p5 ) * t ( p5 p12 h7 h8 )_t * i2 ( h3 h7 h8 h13 p1 p5 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh3,sp12,sh13,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p12,h13,p1
  integer :: h7,h8,p5,sh7,sh8,sp5
  integer :: spin_t2inp
  integer :: spin_itm_hhhhpp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhhhpp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 2.0d+0 * runit

  allocate(itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp12,sh7,sh8)
     spin_itm_hhhhpp = tdcc_spin_dummy3(sh3,sh7,sh8,sh13,sp1,sp5)
     if(spin_t2inp * spin_itm_hhhhpp == 0) cycle

     itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_22_3_1(sh3,sh7,sh8,sh13,sp1,sp5,itm_hhhhpp)

     do h3 = 1,norb1
     do p12 = norb1+1,nact
     do h13 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,p12,h13,p1) = i1(h3,p12,h13,p1) + fact * &
             t2inp(p5,p12,h7,h8,spin_t2inp) * itm_hhhhpp(h3,h7,h8,h13,p1,p5)
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
end subroutine ccsdt_l2_22_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_22_3_1(sh3,sh7,sh8,sh13,sp1,sp5,i2)

!         i2 ( h3 h7 h8 h13 p1 p5 )_yt + = -1 * Sum ( p9 ) * t ( p9 h13 )_t * y ( h3 h7 h8 p1 p5 p9 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh7,sh8,sh13,sp1,sp5
  complex(kind(0d0)),intent(inout) :: &
       i2(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h7,h8,h13,p1,p5
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh13)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh7,sh8,sp1,sp5,sp9)
     if(spin_t1inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h7 = 1,norb1
     do h8 = 1,norb1
     do h13 = 1,norb1
     do p1 = norb1+1,nact
     do p5 = norb1+1,nact
     do p9 = norb1+1,nact
        i2(h3,h7,h8,h13,p1,p5) = i2(h3,h7,h8,h13,p1,p5) + fact * &
             t1inp(p9,h13,spin_t1inp) * g3inp(h3,h7,h8,p1,p5,p9,spin_g3inp)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_22_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_23(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yttv + = 1/2 * P( 2 ) * Sum ( h9 h6 h8 ) * i1 ( h3 h4 h9 h6 h8 p1 )_ytt * v ( h6 h8 h9 p2 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_23_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_23_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_23_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h9,h6,h8,sh9,sh6,sh8
  integer :: spin_itm_hhhhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy3
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhhhp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sh6 = 1,2
  do sh8 = 1,2
     spin_itm_hhhhhp = tdcc_spin_dummy3(sh3,sh4,sh9,sh6,sh8,sp1)
     spin_int2x = tdcc_spin_int2x(sh6,sh8,sh9,sp2)
     if(spin_itm_hhhhhp * spin_int2x == 0) cycle

     itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_23_1(sh3,sh4,sh9,sh6,sh8,sp1,itm_hhhhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do h6 = 1,norb1
     do h8 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hhhhhp(h3,h4,h9,h6,h8,p1) * int2x(h6,h8,h9,p2,spin_int2x)
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
  end subroutine ccsdt_l2_23_perm
  !--------------------------------------------
end subroutine ccsdt_l2_23
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_23_1(sh3,sh4,sh9,sh6,sh8,sp1,i1)

!     i1 ( h3 h4 h9 h6 h8 p1 )_ytt + = 1 * Sum ( p5 ) * t ( p5 h6 )_t * i2 ( h3 h4 h9 h8 p1 p5 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh9,sh6,sh8,sp1
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h9,h6,h8,p1
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_itm_hhhhpp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhhhpp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_itm_hhhhpp = tdcc_spin_dummy3(sh3,sh4,sh9,sh8,sp1,sp5)
     if(spin_t1inp * spin_itm_hhhhpp == 0) cycle

     itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_23_1_1(sh3,sh4,sh9,sh8,sp1,sp5,itm_hhhhpp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h9 = 1,norb1
     do h6 = 1,norb1
     do h8 = 1,norb1
     do p1 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h3,h4,h9,h6,h8,p1) = i1(h3,h4,h9,h6,h8,p1) + fact * &
             t1inp(p5,h6,spin_t1inp) * itm_hhhhpp(h3,h4,h9,h8,p1,p5)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhhpp)
end subroutine ccsdt_l2_23_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_23_1_1(sh3,sh4,sh9,sh8,sp1,sp5,i2)

!         i2 ( h3 h4 h9 h8 p1 p5 )_yt + = -1 * Sum ( p7 ) * t ( p7 h8 )_t * y ( h3 h4 h9 p1 p5 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh9,sh8,sp1,sp5
  complex(kind(0d0)),intent(inout) :: &
       i2(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,h9,h8,p1,p5
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh9,sp1,sp5,sp7)
     if(spin_t1inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h9 = 1,norb1
     do h8 = 1,norb1
     do p1 = norb1+1,nact
     do p5 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h3,h4,h9,h8,p1,p5) = i2(h3,h4,h9,h8,p1,p5) + fact * &
             t1inp(p7,h8,spin_t1inp) * g3inp(h3,h4,h9,p1,p5,p7,spin_g3inp)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_23_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_24(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yttv + = 1/2 * Sum ( p6 p9 ) * i1 ( p6 p9 p1 p2 )_ytt * v ( h3 h4 p6 p9 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: p6,p9,sp6,sp9
  integer :: spin_itm_pppp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp6 = 1,2
  do sp9 = 1,2
     spin_itm_pppp = tdcc_spin_dummy2(sp6,sp9,sp1,sp2)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp6,sp9)
     if(spin_itm_pppp * spin_int2x == 0) cycle

     itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_24_1(sp6,sp9,sp1,sp2,itm_pppp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p6 = norb1+1,nact
     do p9 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_pppp(p6,p9,p1,p2) * int2x(h3,h4,p6,p9,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_pppp)
end subroutine ccsdt_l2_24
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_24_1(sp6,sp9,sp1,sp2,i1)

!     i1 ( p6 p9 p1 p2 )_ytt + = 1 * Sum ( h10 ) * t ( p9 h10 )_t * i2 ( h10 p6 p1 p2 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sp6,sp9,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: p6,p9,p1,p2
  integer :: h10,sh10
  integer :: spin_t1inp
  integer :: spin_itm_hppp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sh10 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
     spin_itm_hppp = tdcc_spin_dummy2(sh10,sp6,sp1,sp2)
     if(spin_t1inp * spin_itm_hppp == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_24_1_1(sh10,sp6,sp1,sp2,itm_hppp)

     do p6 = norb1+1,nact
     do p9 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h10 = 1,norb1
        i1(p6,p9,p1,p2) = i1(p6,p9,p1,p2) + fact * &
             t1inp(p9,h10,spin_t1inp) * itm_hppp(h10,p6,p1,p2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hppp)
end subroutine ccsdt_l2_24_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_24_1_1(sh10,sp6,sp1,sp2,i2)

!         i2 ( h10 p6 p1 p2 )_yt + = 1 * Sum ( h7 h8 p5 ) * t ( p5 p6 h7 h8 )_t * y ( h7 h8 h10 p1 p2 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh10,sp6,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h10,p6,p1,p2
  integer :: h7,h8,p5,sh7,sh8,sp5
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh7,sh8,sh10,sp1,sp2,sp5)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h10 = 1,norb1
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i2(h10,p6,p1,p2) = i2(h10,p6,p1,p2) + fact * &
             t2inp(p5,p6,h7,h8,spin_t2inp) * g3inp(h7,h8,h10,p1,p2,p5,spin_g3inp)
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
end subroutine ccsdt_l2_24_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_25(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yttv + = 1/2 * P( 2 ) * Sum ( p11 h8 h10 ) * i1 ( h3 h4 p11 h8 h10 p1 )_ytt * v ( h8 h10 p2 p11 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccsdt_l2_25_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_25_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsdt_l2_25_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p11,h8,h10,sp11,sh8,sh10
  integer :: spin_itm_hhphhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy3
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhphhp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhphhp(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,(norb1+1):nact))
  do sp11 = 1,2
  do sh8 = 1,2
  do sh10 = 1,2
     spin_itm_hhphhp = tdcc_spin_dummy3(sh3,sh4,sp11,sh8,sh10,sp1)
     spin_int2x = tdcc_spin_int2x(sh8,sh10,sp2,sp11)
     if(spin_itm_hhphhp * spin_int2x == 0) cycle

     itm_hhphhp(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_25_1(sh3,sh4,sp11,sh8,sh10,sp1,itm_hhphhp)
     call ccsdt_l2_25_2(sh3,sh4,sp11,sh8,sh10,sp1,itm_hhphhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p11 = norb1+1,nact
     do h8 = 1,norb1
     do h10 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hhphhp(h3,h4,p11,h8,h10,p1) * int2x(h8,h10,p2,p11,spin_int2x)
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
  deallocate(itm_hhphhp)
  end subroutine ccsdt_l2_25_perm
  !--------------------------------------------
end subroutine ccsdt_l2_25
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_25_1(sh3,sh4,sp11,sh8,sh10,sp1,i1)

!     i1 ( h3 h4 p11 h8 h10 p1 )_ytt + = 1 * t ( p11 h10 )_t * i2 ( h3 h4 h8 p1 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp11,sh8,sh10,sp1
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,p11,h8,h10,p1
  integer :: sdum
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sdum = 1,1
     spin_t1inp = tdcc_spin_t1inp(sp11,sh10)
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh8,sp1)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_25_1_1(sh3,sh4,sh8,sp1,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p11 = norb1+1,nact
     do h8 = 1,norb1
     do h10 = 1,norb1
     do p1 = norb1+1,nact
        i1(h3,h4,p11,h8,h10,p1) = i1(h3,h4,p11,h8,h10,p1) + fact * &
             t1inp(p11,h10,spin_t1inp) * itm_hhhp(h3,h4,h8,p1)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsdt_l2_25_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_25_1_1(sh3,sh4,sh8,sp1,i2)

!         i2 ( h3 h4 h8 p1 )_yt + = 1 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h8 )_t * y ( h3 h4 h7 p1 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh8,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h8,p1
  integer :: h7,p5,p6,sh7,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh7,sp1,sp5,sp6)
     if(spin_t2inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h8 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i2(h3,h4,h8,p1) = i2(h3,h4,h8,p1) + fact * &
             t2inp(p5,p6,h7,h8,spin_t2inp) * g3inp(h3,h4,h7,p1,p5,p6,spin_g3inp)
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
end subroutine ccsdt_l2_25_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_25_2(sh3,sh4,sp11,sh8,sh10,sp1,i1)

!     i1 ( h3 h4 p11 h8 h10 p1 )_yttt + = -1 * Sum ( h6 ) * t ( p11 h6 )_t * i2 ( h3 h4 h6 h8 h10 p1 )_ytt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp11,sh8,sh10,sp1
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,p11,h8,h10,p1
  integer :: h6,sh6
  integer :: spin_t1inp
  integer :: spin_itm_hhhhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhhhhp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp11,sh6)
     spin_itm_hhhhhp = tdcc_spin_dummy3(sh3,sh4,sh6,sh8,sh10,sp1)
     if(spin_t1inp * spin_itm_hhhhhp == 0) cycle

     itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsdt_l2_25_2_1(sh3,sh4,sh6,sh8,sh10,sp1,itm_hhhhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p11 = norb1+1,nact
     do h8 = 1,norb1
     do h10 = 1,norb1
     do p1 = norb1+1,nact
     do h6 = 1,norb1
        i1(h3,h4,p11,h8,h10,p1) = i1(h3,h4,p11,h8,h10,p1) + fact * &
             t1inp(p11,h6,spin_t1inp) * itm_hhhhhp(h3,h4,h6,h8,h10,p1)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhhhp)
end subroutine ccsdt_l2_25_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsdt_l2_25_2_1(sh3,sh4,sh6,sh8,sh10,sp1,i2)

!         i2 ( h3 h4 h6 h8 h10 p1 )_ytt + = 1 * Sum ( p7 ) * t ( p7 h8 )_t * i3 ( h3 h4 h6 h10 p1 p7 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh6,sh8,sh10,sp1
  complex(kind(0d0)),intent(inout) :: &
       i2(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h6,h8,h10,p1
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_itm_hhhhpp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhhhpp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_itm_hhhhpp = tdcc_spin_dummy3(sh3,sh4,sh6,sh10,sp1,sp7)
     if(spin_t1inp * spin_itm_hhhhpp == 0) cycle

     itm_hhhhpp(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsdt_l2_25_2_1_1(sh3,sh4,sh6,sh10,sp1,sp7,itm_hhhhpp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h6 = 1,norb1
     do h8 = 1,norb1
     do h10 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h3,h4,h6,h8,h10,p1) = i2(h3,h4,h6,h8,h10,p1) + fact * &
             t1inp(p7,h8,spin_t1inp) * itm_hhhhpp(h3,h4,h6,h10,p1,p7)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhhpp)
end subroutine ccsdt_l2_25_2_1
!##########################################################
!##########################################################
subroutine ccsdt_l2_25_2_1_1(sh3,sh4,sh6,sh10,sp1,sp7,i3)

!             i3 ( h3 h4 h6 h10 p1 p7 )_yt + = -1 * Sum ( p9 ) * t ( p9 h10 )_t * y ( h3 h4 h6 p1 p7 p9 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh6,sh10,sp1,sp7
  complex(kind(0d0)),intent(inout) :: &
       i3(1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,h6,h10,p1,p7
  integer :: p9,sp9
  integer :: spin_t1inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp9 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh6,sp1,sp7,sp9)
     if(spin_t1inp * spin_g3inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h6 = 1,norb1
     do h10 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
     do p9 = norb1+1,nact
        i3(h3,h4,h6,h10,p1,p7) = i3(h3,h4,h6,h10,p1,p7) + fact * &
             t1inp(p9,h10,spin_t1inp) * g3inp(h3,h4,h6,p1,p7,p9,spin_g3inp)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsdt_l2_25_2_1_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccsdt_l2_main()

  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2out

  implicit none

  call ccsdt_l2_1(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_2(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_3(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_4(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_5(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_6(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_7(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_8(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_9(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_10(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_11(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_12(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_13(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_14(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_15(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_16(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_17(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_18(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_19(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_20(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_21(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_22(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_23(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_24(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsdt_l2_25(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))

  call ccsdt_l2_1(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_2(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_3(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_4(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_5(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_6(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_7(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_8(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_9(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_10(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_11(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_12(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_13(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_14(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_15(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_16(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_17(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_18(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_19(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_20(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_21(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_22(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_23(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_24(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsdt_l2_25(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))

end subroutine ccsdt_l2_main
!**********************************************************
