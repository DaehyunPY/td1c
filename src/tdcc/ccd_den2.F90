!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_1(sp1,sp2,sp3,sp4,i0)

!  i0 ( p1 p2 p3 p4 )_yt + = +1/2 * Sum ( h5 h6 ) * y ( h5 h6 p3 p4 )_y * t ( p1 p2 h5 h6 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sp1,sp2,sp3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,p3,p4
  integer :: h5,h6,sh5,sh6
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh6,sp3,sp4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh6)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
        i0(p1,p2,p3,p4) = i0(p1,p2,p3,p4) + fact * g2inp(h5,h6,p3,p4,spin_g2inp_1) * t2inp(p1,p2,h5,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_den2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_2(sh1,sh2,sh3,sh4,i0)

!  i0 ( h1 h2 h3 h4 )_yt + = +1/2 * Sum ( p5 p6 ) * y ( h1 h2 p5 p6 )_y * t ( p5 p6 h3 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh1,sh2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: h1,h2,h3,h4
  integer :: p5,p6,sp5,sp6
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh1,sh2,sp5,sp6)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp6,sh3,sh4)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(h1,h2,h3,h4) = i0(h1,h2,h3,h4) + fact * g2inp(h1,h2,p5,p6,spin_g2inp_1) * t2inp(p5,p6,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_den2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_3(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_t + = +1 * t ( p1 p2 h3 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,h3,h4
  integer :: sdum
  integer :: spin_t2inp_1
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  spin_t2inp_1 = tdcc_spin_t2inp(sp1,sp2,sh3,sh4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * t2inp(p1,p2,h3,h4,spin_t2inp_1)
  end do
  end do
  end do
  end do
end subroutine ccd_den2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_4(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = +1/2 * P( h3 h4 ) * P( p1 p2 ) * Sum ( h5 p6 ) * i1 ( h5 p6 p1 h3 )_yt * t ( p2 p6 h4 h5 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccd_den2_4_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccd_den2_4_perm(sp1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h4,h3)
  end do
  end do
  end do
  end do

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccd_den2_4_perm(sp2,sp1,sh3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccd_den2_4_perm(sp2,sp1,sh4,sh3,i0_perm)
  end if
  fact_p = +1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h4,h3)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccd_den2_4_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer :: p1,p2,h3,h4
  integer :: h5,p6,sh5,sp6
  integer :: spin_itm_hpph_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hpph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  allocate(itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1))
  do sh5 = 1,2
  do sp6 = 1,2
     spin_itm_hpph_1 = tdcc_spin_dummy2(sh5,sp6,sp1,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp6,sh4,sh5)
     if(spin_itm_hpph_1 * spin_t2inp_2 == 0) cycle

     itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1) = czero
     call ccd_den2_4_1(sh5,sp6,sp1,sh3,itm_hpph)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_hpph(h5,p6,p1,h3) * t2inp(p2,p6,h4,h5,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hpph)
  end subroutine ccd_den2_4_perm
  !--------------------------------------------
end subroutine ccd_den2_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_4_1(sh5,sp6,sp1,sh3,i1)

!      i1 ( h5 p6 p1 h3 )_yt + = +1 * Sum ( h7 p8 ) * y ( h7 h5 p8 p6 )_y * t ( p8 p1 h7 h3 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh5,sp6,sp1,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)
  integer :: h5,p6,p1,h3
  integer :: h7,p8,sh7,sp8
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh7 = 1,2
  do sp8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh7,sh5,sp8,sp6)
     spin_t2inp_2 = tdcc_spin_t2inp(sp8,sp1,sh7,sh3)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do h3 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
        i1(h5,p6,p1,h3) = i1(h5,p6,p1,h3) + fact * g2inp(h7,h5,p8,p6,spin_g2inp_1) * t2inp(p8,p1,h7,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_den2_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_5(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/2 * P( h3 h4 ) * Sum ( h5 ) * i1 ( h5 h3 )_yt * t ( p1 p2 h5 h4 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccd_den2_5_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccd_den2_5_perm(sp1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h4,h3)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccd_den2_5_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer :: p1,p2,h3,h4
  integer :: h5,sh5
  integer :: spin_itm_hh_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh5 = 1,2
     spin_itm_hh_1 = tdcc_spin_dummy1(sh5,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh4)
     if(spin_itm_hh_1 * spin_t2inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccd_den2_5_1(sh5,sh3,itm_hh)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_hh(h5,h3) * t2inp(p1,p2,h5,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccd_den2_5_perm
  !--------------------------------------------
end subroutine ccd_den2_5
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_5_1(sh5,sh3,i1)

!      i1 ( h5 h3 )_yt + = +1 * Sum ( h6 p7 p8 ) * y ( h6 h5 p7 p8 )_y * t ( p7 p8 h6 h3 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh5,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h5,h3
  integer :: h6,p7,p8,sh6,sp7,sp8
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh6,sh5,sp7,sp8)
     spin_t2inp_2 = tdcc_spin_t2inp(sp7,sp8,sh6,sh3)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h5 = 1,norb1
     do h3 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h5,h3) = i1(h5,h3) + fact * g2inp(h6,h5,p7,p8,spin_g2inp_1) * t2inp(p7,p8,h6,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccd_den2_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_6(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/2 * P( p1 p2 ) * Sum ( p5 ) * i1 ( p5 p1 )_yt * t ( p5 p2 h3 h4 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccd_den2_6_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccd_den2_6_perm(sp2,sp1,sh3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h3,h4)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccd_den2_6_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer :: p1,p2,h3,h4
  integer :: p5,sp5
  integer :: spin_itm_pp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_pp_1 = tdcc_spin_dummy1(sp5,sp1)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp2,sh3,sh4)
     if(spin_itm_pp_1 * spin_t2inp_2 == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccd_den2_6_1(sp5,sp1,itm_pp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_pp(p5,p1) * t2inp(p5,p2,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccd_den2_6_perm
  !--------------------------------------------
end subroutine ccd_den2_6
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_6_1(sp5,sp1,i1)

!      i1 ( p5 p1 )_yt + = +1 * Sum ( h6 h7 p8 ) * y ( h6 h7 p8 p5 )_y * t ( p8 p1 h6 h7 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sp5,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p5,p1
  integer :: h6,h7,p8,sh6,sh7,sp8
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
  do sp8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh6,sh7,sp8,sp5)
     spin_t2inp_2 = tdcc_spin_t2inp(sp8,sp1,sh6,sh7)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do p5 = norb1+1,nact
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
        i1(p5,p1) = i1(p5,p1) + fact * g2inp(h6,h7,p8,p5,spin_g2inp_1) * t2inp(p8,p1,h6,h7,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccd_den2_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_7(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = +1/4 * Sum ( h5 h6 ) * i1 ( h5 h6 h3 h4 )_yt * t ( p1 p2 h5 h6 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,h3,h4
  integer :: h5,h6,sh5,sh6
  integer :: spin_itm_hhhh_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh5 = 1,2
  do sh6 = 1,2
     spin_itm_hhhh_1 = tdcc_spin_dummy2(sh5,sh6,sh3,sh4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh6)
     if(spin_itm_hhhh_1 * spin_t2inp_2 == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccd_den2_7_1(sh5,sh6,sh3,sh4,itm_hhhh)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_hhhh(h5,h6,h3,h4) * t2inp(p1,p2,h5,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccd_den2_7
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_7_1(sh5,sh6,sh3,sh4,i1)

!      i1 ( h5 h6 h3 h4 )_yt + = +1 * Sum ( p7 p8 ) * y ( h5 h6 p7 p8 )_y * t ( p7 p8 h3 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh5,sh6,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h5,h6,h3,h4
  integer :: p7,p8,sp7,sp8
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh6,sp7,sp8)
     spin_t2inp_2 = tdcc_spin_t2inp(sp7,sp8,sh3,sh4)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h5 = 1,norb1
     do h6 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h5,h6,h3,h4) = i1(h5,h6,h3,h4) + fact * g2inp(h5,h6,p7,p8,spin_g2inp_1) * t2inp(p7,p8,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_den2_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_8(sh1,sp2,sp3,sh4,i0)

!  i0 ( h1 p2 p3 h4 )_yt + = +1 * Sum ( h5 p6 ) * y ( h5 h1 p6 p3 )_y * t ( p6 p2 h5 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: h1,p2,p3,h4
  integer :: h5,p6,sh5,sp6
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh5 = 1,2
  do sp6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh1,sp6,sp3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp6,sp2,sh5,sh4)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
        i0(h1,p2,p3,h4) = i0(h1,p2,p3,h4) + fact * g2inp(h5,h1,p6,p3,spin_g2inp_1) * t2inp(p6,p2,h5,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccd_den2_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den2_9(sh1,sh2,sp3,sp4,i0)

!  i0 ( h1 h2 p3 p4 )_y + = +1 * y ( h1 h2 p3 p4 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer,intent(in) :: sh1,sh2,sp3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: h1,h2,p3,p4
  integer :: sdum
  integer :: spin_g2inp_1
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  spin_g2inp_1 = tdcc_spin_g2inp(sh1,sh2,sp3,sp4)
  do h1 = 1,norb1
  do h2 = 1,norb1
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
     i0(h1,h2,p3,p4) = i0(h1,h2,p3,p4) + fact * g2inp(h1,h2,p3,p4,spin_g2inp_1)
  end do
  end do
  end do
  end do
end subroutine ccd_den2_9
!##########################################################
!##########################################################
!##########################################################
