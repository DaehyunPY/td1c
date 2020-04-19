!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_1(sh3,sh4,sp1,sp2,i0)

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
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * int2x(h3,h4,p1,p2,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccsd_l2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_2(sh3,sh4,sp1,sp2,i0)

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
  call ccsd_l2_2_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_2_perm(sh4,sh3,sp1,sp2,i0_perm)
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
     call ccsd_l2_2_perm(sh3,sh4,sp2,sp1,i0_perm)
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
     call ccsd_l2_2_perm(sh4,sh3,sp2,sp1,i0_perm)
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
  subroutine ccsd_l2_2_perm(sh3,sh4,sp1,sp2,i0)

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
     call ccsd_l2_2_1(sh4,sp2,itm_hp)
     call ccsd_l2_2_2(sh4,sp2,itm_hp)

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
  end subroutine ccsd_l2_2_perm
  !--------------------------------------------
end subroutine ccsd_l2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_2_1(sh3,sp1,i1)

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
end subroutine ccsd_l2_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_2_2(sh3,sp1,i1)

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
end subroutine ccsd_l2_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_3(sh3,sh4,sp1,sp2,i0)

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
  call ccsd_l2_3_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_3_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccsd_l2_3_perm(sh3,sh4,sp1,sp2,i0)

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
     call ccsd_l2_3_1(sh3,sh4,sh7,sp2,itm_hhhp)
     call ccsd_l2_3_2(sh3,sh4,sh7,sp2,itm_hhhp)

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
  end subroutine ccsd_l2_3_perm
  !--------------------------------------------
end subroutine ccsd_l2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_3_1(sh3,sh4,sh7,sp1,i1)

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
end subroutine ccsd_l2_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_3_2(sh3,sh4,sh7,sp1,i1)

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
end subroutine ccsd_l2_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_4(sh3,sh4,sp1,sp2,i0)

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
  call ccsd_l2_4_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_4_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccsd_l2_4_perm(sh3,sh4,sp1,sp2,i0)

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
  end subroutine ccsd_l2_4_perm
  !--------------------------------------------
end subroutine ccsd_l2_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_5(sh3,sh4,sp1,sp2,i0)

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
  call ccsd_l2_5_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_5_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccsd_l2_5_perm(sh3,sh4,sp1,sp2,i0)

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
     call ccsd_l2_5_1(sh4,sh9,itm_hh)
     call ccsd_l2_5_2(sh4,sh9,itm_hh)
     call ccsd_l2_5_3(sh4,sh9,itm_hh)
     call ccsd_l2_5_4(sh4,sh9,itm_hh)

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
  end subroutine ccsd_l2_5_perm
  !--------------------------------------------
end subroutine ccsd_l2_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_5_1(sh3,sh9,i1)

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
end subroutine ccsd_l2_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_5_2(sh3,sh9,i1)

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
     call ccsd_l2_5_2_1(sh3,sp5,itm_hp)
     call ccsd_l2_5_2_2(sh3,sp5,itm_hp)

     do h3 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,h9) = i1(h3,h9) + fact * t1inp(p5,h9,spin_t1inp) * itm_hp(h3,p5)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsd_l2_5_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_5_2_1(sh3,sp5,i2)

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
end subroutine ccsd_l2_5_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_5_2_2(sh3,sp5,i2)

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
end subroutine ccsd_l2_5_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_5_3(sh3,sh9,i1)

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
end subroutine ccsd_l2_5_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_5_4(sh3,sh9,i1)

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
end subroutine ccsd_l2_5_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_6(sh3,sh4,sp1,sp2,i0)

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
  call ccsd_l2_6_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_6_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccsd_l2_6_perm(sh3,sh4,sp1,sp2,i0)

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
     call ccsd_l2_6_1(sp10,sp2,itm_pp)
     call ccsd_l2_6_2(sp10,sp2,itm_pp)
     call ccsd_l2_6_3(sp10,sp2,itm_pp)
     call ccsd_l2_6_4(sp10,sp2,itm_pp)

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
  end subroutine ccsd_l2_6_perm
  !--------------------------------------------
end subroutine ccsd_l2_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_6_1(sp10,sp1,i1)

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
end subroutine ccsd_l2_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_6_2(sp10,sp1,i1)

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
end subroutine ccsd_l2_6_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_6_3(sp10,sp1,i1)

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
end subroutine ccsd_l2_6_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_6_4(sp10,sp1,i1)

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
     call ccsd_l2_6_4_1(sh6,sp1,itm_hp)

     do p10 = norb1+1,nact
     do p1 = norb1+1,nact
     do h6 = 1,norb1
        i1(p10,p1) = i1(p10,p1) + fact * t1inp(p10,h6,spin_t1inp) * itm_hp(h6,p1)
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
end subroutine ccsd_l2_6_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_6_4_1(sh6,sp1,i2)

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
end subroutine ccsd_l2_6_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_7(sh3,sh4,sp1,sp2,i0)

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
     call ccsd_l2_7_1(sh3,sh4,sh9,sh10,itm_hhhh)
     call ccsd_l2_7_2(sh3,sh4,sh9,sh10,itm_hhhh)
     call ccsd_l2_7_3(sh3,sh4,sh9,sh10,itm_hhhh)

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
end subroutine ccsd_l2_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_7_1(sh3,sh4,sh9,sh10,i1)

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
end subroutine ccsd_l2_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_7_2(sh3,sh4,sh9,sh10,i1)

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
     call ccsd_l2_7_2_1(sh3,sh4,sh10,sp5,itm_hhhp)
     call ccsd_l2_7_2_2(sh3,sh4,sh10,sp5,itm_hhhp)

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
end subroutine ccsd_l2_7_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_7_2_1(sh3,sh4,sh10,sp5,i2)

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
end subroutine ccsd_l2_7_2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_7_2_2(sh3,sh4,sh10,sp5,i2)

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
end subroutine ccsd_l2_7_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_7_3(sh3,sh4,sh9,sh10,i1)

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
end subroutine ccsd_l2_7_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_8(sh3,sh4,sp1,sp2,i0)

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
  call ccsd_l2_8_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_8_perm(sh4,sh3,sp1,sp2,i0_perm)
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
     call ccsd_l2_8_perm(sh3,sh4,sp2,sp1,i0_perm)
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
     call ccsd_l2_8_perm(sh4,sh3,sp2,sp1,i0_perm)
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
  subroutine ccsd_l2_8_perm(sh3,sh4,sp1,sp2,i0)

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
     call ccsd_l2_8_1(sh4,sp7,sh9,sp2,itm_hphp)
     call ccsd_l2_8_2(sh4,sp7,sh9,sp2,itm_hphp)
     call ccsd_l2_8_3(sh4,sp7,sh9,sp2,itm_hphp)

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
  end subroutine ccsd_l2_8_perm
  !--------------------------------------------
end subroutine ccsd_l2_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_8_1(sh3,sp7,sh9,sp1,i1)

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
end subroutine ccsd_l2_8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_8_2(sh3,sp7,sh9,sp1,i1)

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
end subroutine ccsd_l2_8_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_8_3(sh3,sp7,sh9,sp1,i1)

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
end subroutine ccsd_l2_8_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_9(sh3,sh4,sp1,sp2,i0)

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
end subroutine ccsd_l2_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_10(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1 * P( 2 ) * Sum ( h9 ) * i1 ( h3 h9 )_yt * v ( h4 h9 p1 p2 )_v 2

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
  call ccsd_l2_10_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_10_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccsd_l2_10_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h9,sh9
  integer :: spin_itm_hh
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh9 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh3,sh9)
     spin_int2x = tdcc_spin_int2x(sh4,sh9,sp1,sp2)
     if(spin_itm_hh * spin_int2x == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_l2_10_1(sh3,sh9,itm_hh)
     call ccsd_l2_10_2(sh3,sh9,itm_hh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hh(h3,h9) * int2x(h4,h9,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccsd_l2_10_perm
  !--------------------------------------------
end subroutine ccsd_l2_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_10_1(sh3,sh9,i1)

!     i1 ( h3 h9 )_yt + = 1 * Sum ( p5 ) * t ( p5 h9 )_t * y ( h3 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g1inp

  implicit none
  integer,intent(in) :: sh3,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h9
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_g1inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh9)
     spin_g1inp = tdcc_spin_g1inp(sh3,sp5)
     if(spin_t1inp * spin_g1inp == 0) cycle

     do h3 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,h9) = i1(h3,h9) + fact * t1inp(p5,h9,spin_t1inp) * g1inp(h3,p5,spin_g1inp)
     end do
     end do
     end do
  end do
end subroutine ccsd_l2_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_10_2(sh3,sh9,i1)

!     i1 ( h3 h9 )_yt + = -1/2 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h9 )_t * y ( h3 h7 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h9
  integer :: h7,p5,p6,sh7,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh9)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh7,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h9 = 1,norb1
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h9) = i1(h3,h9) + fact * t2inp(p5,p6,h7,h9,spin_t2inp) * g2inp(h3,h7,p5,p6,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_l2_10_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_11(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytf + = 1 * P( 2 ) * Sum ( h5 ) * i1 ( h3 h4 h5 p1 )_yt * f ( h5 p2 )_f 1

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
  call ccsd_l2_11_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_11_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccsd_l2_11_perm(sh3,sh4,sp1,sp2,i0)

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
     call ccsd_l2_11_1(sh3,sh4,sh5,sp1,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h5 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hhhp(h3,h4,h5,p1) * fock(h5,p2,spin_fock)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsd_l2_11_perm
  !--------------------------------------------
end subroutine ccsd_l2_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_11_1(sh3,sh4,sh5,sp1,i1)

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
        i1(h3,h4,h5,p1) = i1(h3,h4,h5,p1) + fact * t1inp(p6,h5,spin_t1inp) * g2inp(h3,h4,p1,p6,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l2_11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_12(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1 * P( 4 ) * Sum ( h7 h6 ) * i1 ( h3 h7 h6 p1 )_yt * v ( h4 h6 h7 p2 )_v 1

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
  call ccsd_l2_12_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_12_perm(sh4,sh3,sp1,sp2,i0_perm)
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
     call ccsd_l2_12_perm(sh3,sh4,sp2,sp1,i0_perm)
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
     call ccsd_l2_12_perm(sh4,sh3,sp2,sp1,i0_perm)
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
  subroutine ccsd_l2_12_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h7,h6,sh7,sh6
  integer :: spin_itm_hhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh7 = 1,2
  do sh6 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh7,sh6,sp1)
     spin_int2x = tdcc_spin_int2x(sh4,sh6,sh7,sp2)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l2_12_1(sh3,sh7,sh6,sp1,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h7 = 1,norb1
     do h6 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hhhp(h3,h7,h6,p1) * int2x(h4,h6,h7,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
  end subroutine ccsd_l2_12_perm
  !--------------------------------------------
end subroutine ccsd_l2_12
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_12_1(sh3,sh7,sh6,sp1,i1)

!     i1 ( h3 h7 h6 p1 )_yt + = 1 * Sum ( p5 ) * t ( p5 h6 )_t * y ( h3 h7 p1 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh7,sh6,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h7,h6,p1
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh7,sp1,sp5)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h7 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h3,h7,h6,p1) = i1(h3,h7,h6,p1) + fact * t1inp(p5,h6,spin_t1inp) * g2inp(h3,h7,p1,p5,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l2_12_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_13(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1 * Sum ( h6 p7 ) * i1 ( h3 h4 h6 p7 )_yt * v ( h6 p7 p1 p2 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: h6,p7,sh6,sp7
  integer :: spin_itm_hhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh6 = 1,2
  do sp7 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh6,sp7)
     spin_int2x = tdcc_spin_int2x(sh6,sp7,sp1,sp2)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l2_13_1(sh3,sh4,sh6,sp7,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hhhp(h3,h4,h6,p7) * int2x(h6,p7,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l2_13
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_13_1(sh3,sh4,sh6,sp7,i1)

!     i1 ( h3 h4 h6 p7 )_yt + = -1 * Sum ( p5 ) * t ( p5 h6 )_t * y ( h3 h4 p5 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh6,sp7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h6,p7
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h3,h4,h6,p7) = i1(h3,h4,h6,p7) + fact * t1inp(p5,h6,spin_t1inp) * g2inp(h3,h4,p5,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l2_13_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_14(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/2 * P( 2 ) * Sum ( p6 ) * i1 ( p6 p1 )_yt * v ( h3 h4 p2 p6 )_v 1

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
  call ccsd_l2_14_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_14_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccsd_l2_14_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p6,sp6
  integer :: spin_itm_pp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp6 = 1,2
     spin_itm_pp = tdcc_spin_dummy1(sp6,sp1)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp2,sp6)
     if(spin_itm_pp * spin_int2x == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_l2_14_1(sp6,sp1,itm_pp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_pp(p6,p1) * int2x(h3,h4,p2,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccsd_l2_14_perm
  !--------------------------------------------
end subroutine ccsd_l2_14
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_14_1(sp6,sp1,i1)

!     i1 ( p6 p1 )_yt + = -1 * Sum ( h7 h8 p5 ) * t ( p5 p6 h7 h8 )_t * y ( h7 h8 p1 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sp6,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p6,p1
  integer :: h7,h8,p5,sh7,sh8,sp5
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh7,sh8,sp1,sp5)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(p6,p1) = i1(p6,p1) + fact * t2inp(p5,p6,h7,h8,spin_t2inp) * g2inp(h7,h8,p1,p5,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_l2_14_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_15(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1/4 * Sum ( h9 h8 ) * i1 ( h3 h4 h8 h9 )_yt * v ( h8 h9 p1 p2 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: h9,h8,sh9,sh8
  integer :: spin_itm_hhhh
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh8 = 1,2
     spin_itm_hhhh = tdcc_spin_dummy2(sh3,sh4,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh8,sh9,sp1,sp2)
     if(spin_itm_hhhh * spin_int2x == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsd_l2_15_1(sh3,sh4,sh8,sh9,itm_hhhh)
     call ccsd_l2_15_2(sh3,sh4,sh8,sh9,itm_hhhh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do h8 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hhhh(h3,h4,h8,h9) * int2x(h8,h9,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsd_l2_15
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_15_1(sh3,sh4,sh8,sh9,i1)

!     i1 ( h3 h4 h8 h9 )_yt + = -1 * Sum ( p5 p6 ) * t ( p5 p6 h8 h9 )_t * y ( h3 h4 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh8,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h8,h9
  integer :: p5,p6,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh8,sh9)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h4,h8,h9) = i1(h3,h4,h8,h9) + fact * t2inp(p5,p6,h8,h9,spin_t2inp) * g2inp(h3,h4,p5,p6,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_l2_15_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_15_2(sh3,sh4,sh8,sh9,i1)

!     i1 ( h3 h4 h8 h9 )_ytt + = -2 * Sum ( p5 ) * t ( p5 h9 )_t * i2 ( h3 h4 h8 p5 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh8,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h8,h9
  integer :: p5,sp5
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sp5 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh9)
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh8,sp5)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l2_15_2_1(sh3,sh4,sh8,sp5,itm_hhhp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,h4,h8,h9) = i1(h3,h4,h8,h9) + fact * t1inp(p5,h9,spin_t1inp) * itm_hhhp(h3,h4,h8,p5)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l2_15_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_15_2_1(sh3,sh4,sh8,sp5,i2)

!         i2 ( h3 h4 h8 p5 )_yt + = -1 * Sum ( p7 ) * t ( p7 h8 )_t * y ( h3 h4 p5 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh8,sp5
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h4,h8,p5
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h3,h4,h8,p5) = i2(h3,h4,h8,p5) + fact * t1inp(p7,h8,spin_t1inp) * g2inp(h3,h4,p5,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l2_15_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_16(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yttv + = -1 * P( 4 ) * Sum ( p5 h8 ) * i1 ( h3 p5 h8 p1 )_ytt * v ( h4 h8 p2 p5 )_v 1

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
  call ccsd_l2_16_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccsd_l2_16_perm(sh4,sh3,sp1,sp2,i0_perm)
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
     call ccsd_l2_16_perm(sh3,sh4,sp2,sp1,i0_perm)
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
     call ccsd_l2_16_perm(sh4,sh3,sp2,sp1,i0_perm)
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
  subroutine ccsd_l2_16_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p5,h8,sp5,sh8
  integer :: spin_itm_hphp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp5 = 1,2
  do sh8 = 1,2
     spin_itm_hphp = tdcc_spin_dummy2(sh3,sp5,sh8,sp1)
     spin_int2x = tdcc_spin_int2x(sh4,sh8,sp2,sp5)
     if(spin_itm_hphp * spin_int2x == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsd_l2_16_1(sh3,sp5,sh8,sp1,itm_hphp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p5 = norb1+1,nact
     do h8 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hphp(h3,p5,h8,p1) * int2x(h4,h8,p2,p5,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hphp)
  end subroutine ccsd_l2_16_perm
  !--------------------------------------------
end subroutine ccsd_l2_16
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_16_1(sh3,sp5,sh8,sp1,i1)

!     i1 ( h3 p5 h8 p1 )_ytt + = 1 * Sum ( h6 ) * t ( p5 h6 )_t * i2 ( h3 h6 h8 p1 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer,intent(in) :: sh3,sp5,sh8,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p5,h8,p1
  integer :: h6,sh6
  integer :: spin_t1inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh6,sh8,sp1)
     if(spin_t1inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccsd_l2_16_1_1(sh3,sh6,sh8,sp1,itm_hhhp)

     do h3 = 1,norb1
     do p5 = norb1+1,nact
     do h8 = 1,norb1
     do p1 = norb1+1,nact
     do h6 = 1,norb1
        i1(h3,p5,h8,p1) = i1(h3,p5,h8,p1) + fact * t1inp(p5,h6,spin_t1inp) * itm_hhhp(h3,h6,h8,p1)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhp)
end subroutine ccsd_l2_16_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_l2_16_1_1(sh3,sh6,sh8,sp1,i2)

!         i2 ( h3 h6 h8 p1 )_yt + = 1 * Sum ( p7 ) * t ( p7 h8 )_t * y ( h3 h6 p1 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh6,sh8,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h3,h6,h8,p1
  integer :: p7,sp7
  integer :: spin_t1inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp7,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh6,sp1,sp7)
     if(spin_t1inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h6 = 1,norb1
     do h8 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i2(h3,h6,h8,p1) = i2(h3,h6,h8,p1) + fact * t1inp(p7,h8,spin_t1inp) * g2inp(h3,h6,p1,p7,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_l2_16_1_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccsd_l2_main()

  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2out

  implicit none

  call ccsd_l2_1(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_2(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_3(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_4(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_5(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_6(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_7(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_8(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_9(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_10(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_11(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_12(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_13(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_14(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_15(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccsd_l2_16(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))

  call ccsd_l2_1(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_2(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_3(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_4(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_5(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_6(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_7(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_8(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_9(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_10(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_11(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_12(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_13(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_14(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_15(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccsd_l2_16(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))

end subroutine ccsd_l2_main
!**********************************************************
