!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2_1e_1(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = -1 * P( 2 ) * Sum ( h5 ) * y ( h3 h5 p1 p2 )_y * f ( h4 h5 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l2_1e_1_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2_1e_1_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccdt_l2_1e_1_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h5,sh5
  integer :: spin_g2inp
  integer :: spin_fock
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh5 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh5,sp1,sp2)
     spin_fock = tdcc_spin_fock(sh4,sh5)
     if(spin_g2inp * spin_fock == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h5 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             g2inp(h3,h5,p1,p2,spin_g2inp) * fock(h4,h5,spin_fock)
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccdt_l2_1e_1_perm
  !--------------------------------------------
end subroutine ccdt_l2_1e_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2_1e_2(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = 1 * P( 2 ) * Sum ( p5 ) * y ( h3 h4 p1 p5 )_y * f ( p5 p2 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l2_1e_2_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2_1e_2_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2_1e_2_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p5,sp5
  integer :: spin_g2inp
  integer :: spin_fock
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp1,sp5)
     spin_fock = tdcc_spin_fock(sp5,sp2)
     if(spin_g2inp * spin_fock == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p5 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             g2inp(h3,h4,p1,p5,spin_g2inp) * fock(p5,p2,spin_fock)
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccdt_l2_1e_2_perm
  !--------------------------------------------
end subroutine ccdt_l2_1e_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2_1e_3(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_fty + = 1/2 * P( 2 ) * Sum ( p6 h7 h8 ) * y ( h3 h7 h8 p1 p2 p6 )_y * i1 ( h4 p6 h7 h8 )_ft 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l2_1e_3_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2_1e_3_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccdt_l2_1e_3_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: p6,h7,h8,sp6,sh7,sh8
  integer :: spin_g3inp
  integer :: spin_itm_hphh
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sp6 = 1,2
  do sh7 = 1,2
  do sh8 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh3,sh7,sh8,sp1,sp2,sp6)
     spin_itm_hphh = tdcc_spin_dummy2(sh4,sp6,sh7,sh8)
     if(spin_g3inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccdt_l2_1e_3_1(sh4,sp6,sh7,sh8,itm_hphh)

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p6 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             g3inp(h3,h7,h8,p1,p2,p6,spin_g3inp) * itm_hphh(h4,p6,h7,h8)
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
  end subroutine ccdt_l2_1e_3_perm
  !--------------------------------------------
end subroutine ccdt_l2_1e_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2_1e_3_1(sh3,sp6,sh7,sh8,i1)

!     i1 ( h3 p6 h7 h8 )_ft + = -1 * Sum ( p5 ) * t ( p5 p6 h7 h8 )_t * f ( h3 p5 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh3,sp6,sh7,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h3,p6,h7,h8
  integer :: p5,sp5
  integer :: spin_t2inp
  integer :: spin_fock
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_fock = tdcc_spin_fock(sh3,sp5)
     if(spin_t2inp * spin_fock == 0) cycle

     do h3 = 1,norb1
     do p6 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,p6,h7,h8) = i1(h3,p6,h7,h8) + fact * &
             t2inp(p5,p6,h7,h8,spin_t2inp) * fock(h3,p5,spin_fock)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccdt_l2_1e_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2_1e_4(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytf + = 1/2 * P( 2 ) * Sum ( h5 ) * i1 ( h3 h4 h5 p1 )_yt * f ( h5 p2 )_f 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l2_1e_4_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2_1e_4_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2_1e_4_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer :: h3,h4,p1,p2
  integer :: h5,sh5
  integer :: spin_itm_hhhp
  integer :: spin_fock
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh5 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh5,sp1)
     spin_fock = tdcc_spin_fock(sh5,sp2)
     if(spin_itm_hhhp * spin_fock == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccdt_l2_1e_4_1(sh3,sh4,sh5,sp1,itm_hhhp)

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
  end subroutine ccdt_l2_1e_4_perm
  !--------------------------------------------
end subroutine ccdt_l2_1e_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2_1e_4_1(sh3,sh4,sh5,sp1,i1)

!     i1 ( h3 h4 h5 p1 )_yt + = -1 * Sum ( h8 p6 p7 ) * t ( p6 p7 h5 h8 )_t * y ( h3 h4 h8 p1 p6 p7 )_y 0

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
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

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
end subroutine ccdt_l2_1e_4_1
!##########################################################
!##########################################################
!##########################################################
