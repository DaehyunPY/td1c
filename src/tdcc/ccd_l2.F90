!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_1(sh3,sh4,sp1,sp2,i0)

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
end subroutine ccd_l2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_2(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = -1 * P( 2 ) * Sum ( h9 ) * y ( h3 h9 p1 p2 )_y * i1 ( h4 h9 )_f 2

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
  call ccd_l2_2_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccd_l2_2_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccd_l2_2_perm(sh3,sh4,sp1,sp2,i0)

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
     call ccd_l2_2_1(sh4,sh9,itm_hh)
     call ccd_l2_2_2(sh4,sh9,itm_hh)

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
  end subroutine ccd_l2_2_perm
  !--------------------------------------------
end subroutine ccd_l2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_2_1(sh3,sh9,i1)

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
end subroutine ccd_l2_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_2_2(sh3,sh9,i1)

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
end subroutine ccd_l2_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_3(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = 1 * P( 2 ) * Sum ( p5 ) * y ( h3 h4 p1 p5 )_y * i1 ( p5 p2 )_f 2

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
  call ccd_l2_3_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccd_l2_3_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccd_l2_3_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p5,sp5
  integer :: spin_g2inp
  integer :: spin_itm_pp
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp1,sp5)
     spin_itm_pp = tdcc_spin_fock(sp5,sp2)
     if(spin_g2inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccd_l2_3_1(sp5,sp2,itm_pp)
     call ccd_l2_3_2(sp5,sp2,itm_pp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p5 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g2inp(h3,h4,p1,p5,spin_g2inp) * itm_pp(p5,p2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccd_l2_3_perm
  !--------------------------------------------
end subroutine ccd_l2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_3_1(sp5,sp1,i1)

!     i1 ( p5 p1 )_f + = 1 * f ( p5 p1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sp5,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p5,p1
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp5,sp1)
  do p5 = norb1+1,nact
  do p1 = norb1+1,nact
     i1(p5,p1) = i1(p5,p1) + fact * fock(p5,p1,spin_fock)
  end do
  end do
end subroutine ccd_l2_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_3_2(sp5,sp1,i1)

!     i1 ( p5 p1 )_vt + = -1/2 * Sum ( h7 h8 p6 ) * t ( p5 p6 h7 h8 )_t * v ( h7 h8 p1 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp5,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p5,p1
  integer :: h7,h8,p6,sh7,sh8,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_int2x = tdcc_spin_int2x(sh7,sh8,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do p5 = norb1+1,nact
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i1(p5,p1) = i1(p5,p1) + fact * t2inp(p5,p6,h7,h8,spin_t2inp) * int2x(h7,h8,p1,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccd_l2_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_4(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1/2 * Sum ( h9 h7 ) * y ( h7 h9 p1 p2 )_y * i1 ( h3 h4 h7 h9 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: h9,h7,sh9,sh7
  integer :: spin_g2inp
  integer :: spin_itm_hhhh
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh7 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh7,sh9,sp1,sp2)
     spin_itm_hhhh = tdcc_spin_int2x(sh3,sh4,sh7,sh9)
     if(spin_g2inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccd_l2_4_1(sh3,sh4,sh7,sh9,itm_hhhh)
     call ccd_l2_4_2(sh3,sh4,sh7,sh9,itm_hhhh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do h7 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g2inp(h7,h9,p1,p2,spin_g2inp) * itm_hhhh(h3,h4,h7,h9)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccd_l2_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_4_1(sh3,sh4,sh7,sh9,i1)

!     i1 ( h3 h4 h7 h9 )_v + = -1 * v ( h3 h4 h7 h9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sh7,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h7,h9
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh4,sh7,sh9)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do h7 = 1,norb1
  do h9 = 1,norb1
     i1(h3,h4,h7,h9) = i1(h3,h4,h7,h9) + fact * int2x(h3,h4,h7,h9,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccd_l2_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_4_2(sh3,sh4,sh7,sh9,i1)

!     i1 ( h3 h4 h7 h9 )_vt + = -1/2 * Sum ( p5 p6 ) * t ( p5 p6 h7 h9 )_t * v ( h3 h4 p5 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sh7,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h7,h9
  integer :: p5,p6,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh9)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h7 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h4,h7,h9) = i1(h3,h4,h7,h9) + fact * t2inp(p5,p6,h7,h9,spin_t2inp) * int2x(h3,h4,p5,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_l2_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_5(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1 * P( 4 ) * Sum ( h9 p5 ) * y ( h3 h9 p1 p5 )_y * i1 ( h4 p5 h9 p2 )_v 2

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
  call ccd_l2_5_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccd_l2_5_perm(sh4,sh3,sp1,sp2,i0_perm)
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
     call ccd_l2_5_perm(sh3,sh4,sp2,sp1,i0_perm)
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
     call ccd_l2_5_perm(sh4,sh3,sp2,sp1,i0_perm)
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
  subroutine ccd_l2_5_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h9,p5,sh9,sp5
  integer :: spin_g2inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sp5 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh9,sp1,sp5)
     spin_itm_hphp = tdcc_spin_int2x(sh4,sp5,sh9,sp2)
     if(spin_g2inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccd_l2_5_1(sh4,sp5,sh9,sp2,itm_hphp)
     call ccd_l2_5_2(sh4,sp5,sh9,sp2,itm_hphp)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do p5 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g2inp(h3,h9,p1,p5,spin_g2inp) * itm_hphp(h4,p5,h9,p2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hphp)
  end subroutine ccd_l2_5_perm
  !--------------------------------------------
end subroutine ccd_l2_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_5_1(sh3,sp5,sh9,sp1,i1)

!     i1 ( h3 p5 h9 p1 )_v + = 1 * v ( h3 p5 h9 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sp5,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p5,h9,p1
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sp5,sh9,sp1)
  do h3 = 1,norb1
  do p5 = norb1+1,nact
  do h9 = 1,norb1
  do p1 = norb1+1,nact
     i1(h3,p5,h9,p1) = i1(h3,p5,h9,p1) + fact * int2x(h3,p5,h9,p1,spin_int2x)
  end do
  end do
  end do
  end do
end subroutine ccd_l2_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_5_2(sh3,sp5,sh9,sp1,i1)

!     i1 ( h3 p5 h9 p1 )_vt + = 1 * Sum ( h8 p6 ) * t ( p5 p6 h8 h9 )_t * v ( h3 h8 p1 p6 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh3,sp5,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h3,p5,h9,p1
  integer :: h8,p6,sh8,sp6
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh3,sh8,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do p5 = norb1+1,nact
     do h9 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i1(h3,p5,h9,p1) = i1(h3,p5,h9,p1) + fact * t2inp(p5,p6,h8,h9,spin_t2inp) * int2x(h3,h8,p1,p6,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_l2_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_6(sh3,sh4,sp1,sp2,i0)

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
end subroutine ccd_l2_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_7(sh3,sh4,sp1,sp2,i0)

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
  call ccd_l2_7_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccd_l2_7_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccd_l2_7_perm(sh3,sh4,sp1,sp2,i0)

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
     call ccd_l2_7_1(sp6,sp1,itm_pp)

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
  end subroutine ccd_l2_7_perm
  !--------------------------------------------
end subroutine ccd_l2_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_7_1(sp6,sp1,i1)

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
end subroutine ccd_l2_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_8(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/2 * P( 2 ) * Sum ( h8 ) * i1 ( h3 h8 )_yt * v ( h4 h8 p1 p2 )_v 1

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
  call ccd_l2_8_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccd_l2_8_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccd_l2_8_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h8,sh8
  integer :: spin_itm_hh
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh8 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh3,sh8)
     spin_int2x = tdcc_spin_int2x(sh4,sh8,sp1,sp2)
     if(spin_itm_hh * spin_int2x == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccd_l2_8_1(sh3,sh8,itm_hh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h8 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hh(h3,h8) * int2x(h4,h8,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccd_l2_8_perm
  !--------------------------------------------
end subroutine ccd_l2_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_8_1(sh3,sh8,i1)

!     i1 ( h3 h8 )_yt + = -1 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h8 )_t * y ( h3 h7 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h3,h8
  integer :: h7,p5,p6,sh7,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh7,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h8 = 1,norb1
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h8) = i1(h3,h8) + fact * t2inp(p5,p6,h7,h8,spin_t2inp) * g2inp(h3,h7,p5,p6,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccd_l2_8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_9(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/4 * Sum ( h7 h8 ) * i1 ( h3 h4 h7 h8 )_yt * v ( h7 h8 p1 p2 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: h7,h8,sh7,sh8
  integer :: spin_itm_hhhh
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh7 = 1,2
  do sh8 = 1,2
     spin_itm_hhhh = tdcc_spin_dummy2(sh3,sh4,sh7,sh8)
     spin_int2x = tdcc_spin_int2x(sh7,sh8,sp1,sp2)
     if(spin_itm_hhhh * spin_int2x == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccd_l2_9_1(sh3,sh4,sh7,sh8,itm_hhhh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * itm_hhhh(h3,h4,h7,h8) * int2x(h7,h8,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccd_l2_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_l2_9_1(sh3,sh4,sh7,sh8,i1)

!     i1 ( h3 h4 h7 h8 )_yt + = 1 * Sum ( p5 p6 ) * t ( p5 p6 h7 h8 )_t * y ( h3 h4 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sh7,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h3,h4,h7,h8
  integer :: p5,p6,sp5,sp6
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h4,h7,h8) = i1(h3,h4,h7,h8) + fact * t2inp(p5,p6,h7,h8,spin_t2inp) * g2inp(h3,h4,p5,p6,spin_g2inp)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_l2_9_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccd_l2_main()

  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2out

  implicit none

  call ccd_l2_1(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccd_l2_2(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccd_l2_3(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccd_l2_4(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccd_l2_5(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccd_l2_6(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccd_l2_7(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccd_l2_8(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccd_l2_9(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))

  call ccd_l2_1(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccd_l2_2(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccd_l2_3(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccd_l2_4(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccd_l2_5(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccd_l2_6(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccd_l2_7(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccd_l2_8(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccd_l2_9(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))

end subroutine ccd_l2_main
!**********************************************************
