!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_1(sp1,sp2,sp3,sh4,i0)

!  i0 ( p1 p2 p3 h4 )_ytt + = +1 * P( p1 p2 ) * i1 ( p3 p1 )_yt * t ( p2 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,p3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_1_perm(sp1,sp2,sp3,sh4,i0_perm)
  fact_p = +1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
  do h4 = 1,norb1
     i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact_p * i0_perm(p1,p2,p3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sp1 * sp2 * sp3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccsd_den2_1_perm(sp2,sp1,sp3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
  do h4 = 1,norb1
     i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact_p * i0_perm(p2,p1,p3,h4)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsd_den2_1_perm(sp1,sp2,sp3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,p3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_itm_pp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sdum = 1,1
     spin_itm_pp_1 = tdcc_spin_dummy1(sp3,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_itm_pp_1 * spin_t1inp_2 == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_1_1(sp3,sp1,itm_pp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
        i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact * itm_pp(p3,p1) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccsd_den2_1_perm
  !--------------------------------------------
end subroutine ccsd_den2_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_1_1(sp3,sp1,i1)

!      i1 ( p3 p1 )_yt + = +1 * Sum ( h5 ) * y ( h5 p3 )_y * t ( p1 h5 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p3,p1
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh5 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh5,sp3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh5)
     if(spin_g1inp_1 * spin_t1inp_2 == 0) cycle

     do p3 = norb1+1,nact
     do p1 = norb1+1,nact
     do h5 = 1,norb1
        i1(p3,p1) = i1(p3,p1) + fact * g1inp(h5,p3,spin_g1inp_1) * t1inp(p1,h5,spin_t1inp_2)
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_2(sh1,sp2,sh3,sh4,i0)

!  i0 ( h1 p2 h3 h4 )_ytt + = -1 * P( h3 h4 ) * i1 ( h1 h3 )_yt * t ( p2 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_2_perm(sh1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  do h1 = 1,norb1
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact_p * i0_perm(h1,p2,h3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sh1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccsd_den2_2_perm(sh1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h1 = 1,norb1
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact_p * i0_perm(h1,p2,h4,h3)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsd_den2_2_perm(sh1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: h1,p2,h3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_itm_hh_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sdum = 1,1
     spin_itm_hh_1 = tdcc_spin_dummy1(sh1,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_itm_hh_1 * spin_t1inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_den2_2_1(sh1,sh3,itm_hh)

     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
        i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact * itm_hh(h1,h3) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccsd_den2_2_perm
  !--------------------------------------------
end subroutine ccsd_den2_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_2_1(sh1,sh3,i1)

!      i1 ( h1 h3 )_yt + = +1 * Sum ( p5 ) * y ( h1 p5 )_y * t ( p5 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h1,h3
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp5 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh1,sp5)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh3)
     if(spin_g1inp_1 * spin_t1inp_2 == 0) cycle

     do h1 = 1,norb1
     do h3 = 1,norb1
     do p5 = norb1+1,nact
        i1(h1,h3) = i1(h1,h3) + fact * g1inp(h1,p5,spin_g1inp_1) * t1inp(p5,h3,spin_t1inp_2)
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_3(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_tt + = +1/2 * P( h3 h4 ) * P( p1 p2 ) * t ( p1 h3 )_t * t ( p2 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_3_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_3_perm(sp1,sp2,sh4,sh3,i0_perm)
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
     call ccsd_den2_3_perm(sp2,sp1,sh3,sh4,i0_perm)
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
     call ccsd_den2_3_perm(sp2,sp1,sh4,sh3,i0_perm)
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
  subroutine ccsd_den2_3_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_t1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  do sdum = 1,1
     spin_t1inp_1 = tdcc_spin_t1inp(sp1,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_t1inp_1 * spin_t1inp_2 == 0) cycle

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * t1inp(p1,h3,spin_t1inp_1) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
  end subroutine ccsd_den2_3_perm
  !--------------------------------------------
end subroutine ccsd_den2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_4(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_yttt + = -1 * P( h3 h4 ) * P( p1 p2 ) * i1 ( h3 p1 )_ytt * t ( p2 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_4_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_4_perm(sp1,sp2,sh4,sh3,i0_perm)
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
     call ccsd_den2_4_perm(sp2,sp1,sh3,sh4,i0_perm)
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
     call ccsd_den2_4_perm(sp2,sp1,sh4,sh3,i0_perm)
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
  subroutine ccsd_den2_4_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_itm_hp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sdum = 1,1
     spin_itm_hp_1 = tdcc_spin_dummy1(sh3,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_itm_hp_1 * spin_t1inp_2 == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsd_den2_4_1(sh3,sp1,itm_hp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_hp(h3,p1) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
  end subroutine ccsd_den2_4_perm
  !--------------------------------------------
end subroutine ccsd_den2_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_4_1(sh3,sp1,i1)

!      i1 ( h3 p1 )_ytt + = +1 * Sum ( h5 ) * i2 ( h5 h3 )_yt * t ( p1 h5 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh3,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,p1
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hh_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh5 = 1,2
     spin_itm_hh_1 = tdcc_spin_dummy1(sh5,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh5)
     if(spin_itm_hh_1 * spin_t1inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_den2_4_1_1(sh5,sh3,itm_hh)

     do h3 = 1,norb1
     do p1 = norb1+1,nact
     do h5 = 1,norb1
        i1(h3,p1) = i1(h3,p1) + fact * itm_hh(h5,h3) * t1inp(p1,h5,spin_t1inp_2)
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
end subroutine ccsd_den2_4_1
!##########################################################
!##########################################################
subroutine ccsd_den2_4_1_1(sh5,sh3,i2)

!          i2 ( h5 h3 )_yt + = +1 * Sum ( p6 ) * y ( h5 p6 )_y * t ( p6 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh5,sh3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1)
  integer(c_int) :: h5,h3
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp6 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh5,sp6)
     spin_t1inp_2 = tdcc_spin_t1inp(sp6,sh3)
     if(spin_g1inp_1 * spin_t1inp_2 == 0) cycle

     do h5 = 1,norb1
     do h3 = 1,norb1
     do p6 = norb1+1,nact
        i2(h5,h3) = i2(h5,h3) + fact * g1inp(h5,p6,spin_g1inp_1) * t1inp(p6,h3,spin_t1inp_2)
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_4_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_5(sh1,sp2,sp3,sh4,i0)

!  i0 ( h1 p2 p3 h4 )_yt + = +1 * y ( h1 p3 )_y * t ( p2 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,p2,p3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sdum = 1,1
     spin_g1inp_1 = tdcc_spin_g1inp(sh1,sp3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_g1inp_1 * spin_t1inp_2 == 0) cycle

     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
        i0(h1,p2,p3,h4) = i0(h1,p2,p3,h4) + fact * g1inp(h1,p3,spin_g1inp_1) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_6(sp1,sp2,sp3,sp4,i0)

!  i0 ( p1 p2 p3 p4 )_yt + = +1/2 * Sum ( h5 h6 ) * y ( h5 h6 p3 p4 )_y * t ( p1 p2 h5 h6 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,p3,p4
  integer(c_int) :: h5,h6,sh5,sh6
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
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
end subroutine ccsd_den2_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_7(sh1,sh2,sh3,sh4,i0)

!  i0 ( h1 h2 h3 h4 )_yt + = +1/2 * Sum ( p5 p6 ) * y ( h1 h2 p5 p6 )_y * t ( p5 p6 h3 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,h2,h3,h4
  integer(c_int) :: p5,p6,sp5,sp6
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
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
end subroutine ccsd_den2_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_8(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_t + = +1 * t ( p1 p2 h3 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_t2inp_1
  integer(c_int),external :: tdcc_spin_t2inp
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
end subroutine ccsd_den2_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_9(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = +1/2 * P( h3 h4 ) * P( p1 p2 ) * Sum ( h5 p6 ) * i1 ( h5 p6 p1 h3 )_yt * t ( p2 p6 h4 h5 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_9_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_9_perm(sp1,sp2,sh4,sh3,i0_perm)
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
     call ccsd_den2_9_perm(sp2,sp1,sh3,sh4,i0_perm)
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
     call ccsd_den2_9_perm(sp2,sp1,sh4,sh3,i0_perm)
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
  subroutine ccsd_den2_9_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: h5,p6,sh5,sp6
  integer(c_int) :: spin_itm_hpph_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hpph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  allocate(itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1))
  do sh5 = 1,2
  do sp6 = 1,2
     spin_itm_hpph_1 = tdcc_spin_dummy2(sh5,sp6,sp1,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp6,sh4,sh5)
     if(spin_itm_hpph_1 * spin_t2inp_2 == 0) cycle

     itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1) = czero
     call ccsd_den2_9_1(sh5,sp6,sp1,sh3,itm_hpph)

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
  end subroutine ccsd_den2_9_perm
  !--------------------------------------------
end subroutine ccsd_den2_9
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_9_1(sh5,sp6,sp1,sh3,i1)

!      i1 ( h5 p6 p1 h3 )_yt + = +1 * Sum ( h7 p8 ) * y ( h7 h5 p8 p6 )_y * t ( p8 p1 h7 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh5,sp6,sp1,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)
  integer(c_int) :: h5,p6,p1,h3
  integer(c_int) :: h7,p8,sh7,sp8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
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
end subroutine ccsd_den2_9_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_10(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/2 * P( h3 h4 ) * Sum ( h5 ) * i1 ( h5 h3 )_yt * t ( p1 p2 h5 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_10_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_10_perm(sp1,sp2,sh4,sh3,i0_perm)
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
  subroutine ccsd_den2_10_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hh_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh5 = 1,2
     spin_itm_hh_1 = tdcc_spin_dummy1(sh5,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh4)
     if(spin_itm_hh_1 * spin_t2inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_den2_10_1(sh5,sh3,itm_hh)

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
  end subroutine ccsd_den2_10_perm
  !--------------------------------------------
end subroutine ccsd_den2_10
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_10_1(sh5,sh3,i1)

!      i1 ( h5 h3 )_yt + = +1 * Sum ( h6 p7 p8 ) * y ( h6 h5 p7 p8 )_y * t ( p7 p8 h6 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh5,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h5,h3
  integer(c_int) :: h6,p7,p8,sh6,sp7,sp8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
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
end subroutine ccsd_den2_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_11(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/2 * P( p1 p2 ) * Sum ( p5 ) * i1 ( p5 p1 )_yt * t ( p5 p2 h3 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_11_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_11_perm(sp2,sp1,sh3,sh4,i0_perm)
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
  subroutine ccsd_den2_11_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_pp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_pp_1 = tdcc_spin_dummy1(sp5,sp1)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp2,sh3,sh4)
     if(spin_itm_pp_1 * spin_t2inp_2 == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_11_1(sp5,sp1,itm_pp)

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
  end subroutine ccsd_den2_11_perm
  !--------------------------------------------
end subroutine ccsd_den2_11
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_11_1(sp5,sp1,i1)

!      i1 ( p5 p1 )_yt + = +1 * Sum ( h6 h7 p8 ) * y ( h6 h7 p8 p5 )_y * t ( p8 p1 h6 h7 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp5,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p5,p1
  integer(c_int) :: h6,h7,p8,sh6,sh7,sp8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
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
end subroutine ccsd_den2_11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_12(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = +1/4 * Sum ( h5 h6 ) * i1 ( h5 h6 h3 h4 )_yt * t ( p1 p2 h5 h6 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: h5,h6,sh5,sh6
  integer(c_int) :: spin_itm_hhhh_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh5 = 1,2
  do sh6 = 1,2
     spin_itm_hhhh_1 = tdcc_spin_dummy2(sh5,sh6,sh3,sh4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh6)
     if(spin_itm_hhhh_1 * spin_t2inp_2 == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsd_den2_12_1(sh5,sh6,sh3,sh4,itm_hhhh)

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
end subroutine ccsd_den2_12
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_12_1(sh5,sh6,sh3,sh4,i1)

!      i1 ( h5 h6 h3 h4 )_yt + = +1 * Sum ( p7 p8 ) * y ( h5 h6 p7 p8 )_y * t ( p7 p8 h3 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh5,sh6,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h5,h6,h3,h4
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
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
end subroutine ccsd_den2_12_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_13(sh1,sp2,sp3,sh4,i0)

!  i0 ( h1 p2 p3 h4 )_yt + = +1 * Sum ( h5 p6 ) * y ( h5 h1 p6 p3 )_y * t ( p6 p2 h5 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,p2,p3,h4
  integer(c_int) :: h5,p6,sh5,sp6
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
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
end subroutine ccsd_den2_13
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_14(sh1,sh2,sp3,sp4,i0)

!  i0 ( h1 h2 p3 p4 )_y + = +1 * y ( h1 h2 p3 p4 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh2,sp3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,h2,p3,p4
  integer(c_int) :: sdum
  integer(c_int) :: spin_g2inp_1
  integer(c_int),external :: tdcc_spin_g2inp
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
end subroutine ccsd_den2_14
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_15(sp1,sp2,sp3,sp4,i0)

!  i0 ( p1 p2 p3 p4 )_ytt + = +1 * Sum ( h5 ) * i1 ( h5 p3 p4 p1 )_yt * t ( p2 h5 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,p3,p4
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hppp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sh5 = 1,2
     spin_itm_hppp_1 = tdcc_spin_dummy2(sh5,sp3,sp4,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh5)
     if(spin_itm_hppp_1 * spin_t1inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_15_1(sh5,sp3,sp4,sp1,itm_hppp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h5 = 1,norb1
        i0(p1,p2,p3,p4) = i0(p1,p2,p3,p4) + fact * itm_hppp(h5,p3,p4,p1) * t1inp(p2,h5,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hppp)
end subroutine ccsd_den2_15
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_15_1(sh5,sp3,sp4,sp1,i1)

!      i1 ( h5 p3 p4 p1 )_yt + = +1 * Sum ( h6 ) * y ( h6 h5 p3 p4 )_y * t ( p1 h6 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh5,sp3,sp4,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h5,p3,p4,p1
  integer(c_int) :: h6,sh6
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh6,sh5,sp3,sp4)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh6)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do p1 = norb1+1,nact
     do h6 = 1,norb1
        i1(h5,p3,p4,p1) = i1(h5,p3,p4,p1) + fact * g2inp(h6,h5,p3,p4,spin_g2inp_1) * t1inp(p1,h6,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_15_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_16(sh1,sh2,sh3,sh4,i0)

!  i0 ( h1 h2 h3 h4 )_ytt + = +1 * Sum ( p5 ) * i1 ( h1 h2 p5 h3 )_yt * t ( p5 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,h2,h3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_hhph_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sp5 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh1,sh2,sp5,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh4)
     if(spin_itm_hhph_1 * spin_t1inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call ccsd_den2_16_1(sh1,sh2,sp5,sh3,itm_hhph)

     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(h1,h2,h3,h4) = i0(h1,h2,h3,h4) + fact * itm_hhph(h1,h2,p5,h3) * t1inp(p5,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhph)
end subroutine ccsd_den2_16
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_16_1(sh1,sh2,sp5,sh3,i1)

!      i1 ( h1 h2 p5 h3 )_yt + = +1 * Sum ( p6 ) * y ( h1 h2 p6 p5 )_y * t ( p6 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh2,sp5,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer(c_int) :: h1,h2,p5,h3
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh1,sh2,sp6,sp5)
     spin_t1inp_2 = tdcc_spin_t1inp(sp6,sh3)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h1 = 1,norb1
     do h2 = 1,norb1
     do p5 = norb1+1,nact
     do h3 = 1,norb1
     do p6 = norb1+1,nact
        i1(h1,h2,p5,h3) = i1(h1,h2,p5,h3) + fact * g2inp(h1,h2,p6,p5,spin_g2inp_1) * t1inp(p6,h3,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_16_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_17(sp1,sp2,sp3,sh4,i0)

!  i0 ( p1 p2 p3 h4 )_yt + = +1 * Sum ( h5 ) * y ( h5 p3 )_y * t ( p1 p2 h5 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,p3,h4
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh5 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh5,sp3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh4)
     if(spin_g1inp_1 * spin_t2inp_2 == 0) cycle

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do h5 = 1,norb1
        i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact * g1inp(h5,p3,spin_g1inp_1) * t2inp(p1,p2,h5,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_17
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_18(sh1,sp2,sh3,sh4,i0)

!  i0 ( h1 p2 h3 h4 )_yt + = -1 * Sum ( p5 ) * y ( h1 p5 )_y * t ( p5 p2 h3 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,p2,h3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh1,sp5)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp2,sh3,sh4)
     if(spin_g1inp_1 * spin_t2inp_2 == 0) cycle

     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact * g1inp(h1,p5,spin_g1inp_1) * t2inp(p5,p2,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_18
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_19(sp1,sp2,sp3,sh4,i0)

!  i0 ( p1 p2 p3 h4 )_ytt + = +1/2 * P( p1 p2 ) * i1 ( p3 p1 )_yt * t ( p2 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,p3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_19_perm(sp1,sp2,sp3,sh4,i0_perm)
  fact_p = +1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
  do h4 = 1,norb1
     i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact_p * i0_perm(p1,p2,p3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sp1 * sp2 * sp3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccsd_den2_19_perm(sp2,sp1,sp3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
  do h4 = 1,norb1
     i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact_p * i0_perm(p2,p1,p3,h4)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsd_den2_19_perm(sp1,sp2,sp3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,p3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_itm_pp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sdum = 1,1
     spin_itm_pp_1 = tdcc_spin_dummy1(sp3,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_itm_pp_1 * spin_t1inp_2 == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_19_1(sp3,sp1,itm_pp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
        i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact * itm_pp(p3,p1) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
  end subroutine ccsd_den2_19_perm
  !--------------------------------------------
end subroutine ccsd_den2_19
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_19_1(sp3,sp1,i1)

!      i1 ( p3 p1 )_yt + = +1 * Sum ( h5 h6 p7 ) * y ( h5 h6 p7 p3 )_y * t ( p7 p1 h5 h6 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p3,p1
  integer(c_int) :: h5,h6,p7,sh5,sh6,sp7
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh6,sp7,sp3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp7,sp1,sh5,sh6)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do p3 = norb1+1,nact
     do p1 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i1(p3,p1) = i1(p3,p1) + fact * g2inp(h5,h6,p7,p3,spin_g2inp_1) * t2inp(p7,p1,h5,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_den2_19_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_20(sh1,sp2,sh3,sh4,i0)

!  i0 ( h1 p2 h3 h4 )_ytt + = -1/2 * P( h3 h4 ) * i1 ( h1 h3 )_yt * t ( p2 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_20_perm(sh1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  do h1 = 1,norb1
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact_p * i0_perm(h1,p2,h3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sh1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccsd_den2_20_perm(sh1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h1 = 1,norb1
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact_p * i0_perm(h1,p2,h4,h3)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsd_den2_20_perm(sh1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: h1,p2,h3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_itm_hh_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sdum = 1,1
     spin_itm_hh_1 = tdcc_spin_dummy1(sh1,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_itm_hh_1 * spin_t1inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_den2_20_1(sh1,sh3,itm_hh)

     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
        i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact * itm_hh(h1,h3) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
  end subroutine ccsd_den2_20_perm
  !--------------------------------------------
end subroutine ccsd_den2_20
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_20_1(sh1,sh3,i1)

!      i1 ( h1 h3 )_yt + = +1 * Sum ( h5 p6 p7 ) * y ( h5 h1 p6 p7 )_y * t ( p6 p7 h5 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h1,h3
  integer(c_int) :: h5,p6,p7,sh5,sp6,sp7
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh1,sp6,sp7)
     spin_t2inp_2 = tdcc_spin_t2inp(sp6,sp7,sh5,sh3)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h1 = 1,norb1
     do h3 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h1,h3) = i1(h1,h3) + fact * g2inp(h5,h1,p6,p7,spin_g2inp_1) * t2inp(p6,p7,h5,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_den2_20_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_21(sp1,sp2,sp3,sh4,i0)

!  i0 ( p1 p2 p3 h4 )_ytt + = +1 * P( p1 p2 ) * Sum ( h5 ) * i1 ( h5 p3 p2 h4 )_yt * t ( p1 h5 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,p3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_21_perm(sp1,sp2,sp3,sh4,i0_perm)
  fact_p = +1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
  do h4 = 1,norb1
     i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact_p * i0_perm(p1,p2,p3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sp1 * sp2 * sp3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccsd_den2_21_perm(sp2,sp1,sp3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
  do h4 = 1,norb1
     i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact_p * i0_perm(p2,p1,p3,h4)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsd_den2_21_perm(sp1,sp2,sp3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,p3,h4
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hpph_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hpph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1))
  do sh5 = 1,2
     spin_itm_hpph_1 = tdcc_spin_dummy2(sh5,sp3,sp2,sh4)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh5)
     if(spin_itm_hpph_1 * spin_t1inp_2 == 0) cycle

     itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1) = czero
     call ccsd_den2_21_1(sh5,sp3,sp2,sh4,itm_hpph)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do h5 = 1,norb1
        i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact * itm_hpph(h5,p3,p2,h4) * t1inp(p1,h5,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hpph)
  end subroutine ccsd_den2_21_perm
  !--------------------------------------------
end subroutine ccsd_den2_21
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_21_1(sh5,sp3,sp2,sh4,i1)

!      i1 ( h5 p3 p2 h4 )_yt + = +1 * Sum ( h6 p7 ) * y ( h5 h6 p3 p7 )_y * t ( p2 p7 h4 h6 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh5,sp3,sp2,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)
  integer(c_int) :: h5,p3,p2,h4
  integer(c_int) :: h6,p7,sh6,sp7
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sp7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh6,sp3,sp7)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp7,sh4,sh6)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p2 = norb1+1,nact
     do h4 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i1(h5,p3,p2,h4) = i1(h5,p3,p2,h4) + fact * g2inp(h5,h6,p3,p7,spin_g2inp_1) * t2inp(p2,p7,h4,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_den2_21_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_22(sh1,sp2,sh3,sh4,i0)

!  i0 ( h1 p2 h3 h4 )_ytt + = -1 * P( h3 h4 ) * Sum ( p5 ) * i1 ( h1 p5 p2 h4 )_yt * t ( p5 h3 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_22_perm(sh1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  do h1 = 1,norb1
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact_p * i0_perm(h1,p2,h3,h4)
  end do
  end do
  end do
  end do

  if(.not. (sh1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccsd_den2_22_perm(sh1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h1 = 1,norb1
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact_p * i0_perm(h1,p2,h4,h3)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccsd_den2_22_perm(sh1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: h1,p2,h3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_hpph_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hpph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1))
  do sp5 = 1,2
     spin_itm_hpph_1 = tdcc_spin_dummy2(sh1,sp5,sp2,sh4)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh3)
     if(spin_itm_hpph_1 * spin_t1inp_2 == 0) cycle

     itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1) = czero
     call ccsd_den2_22_1(sh1,sp5,sp2,sh4,itm_hpph)

     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact * itm_hpph(h1,p5,p2,h4) * t1inp(p5,h3,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hpph)
  end subroutine ccsd_den2_22_perm
  !--------------------------------------------
end subroutine ccsd_den2_22
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_22_1(sh1,sp5,sp2,sh4,i1)

!      i1 ( h1 p5 p2 h4 )_yt + = +1 * Sum ( h6 p7 ) * y ( h1 h6 p5 p7 )_y * t ( p2 p7 h4 h6 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp5,sp2,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)
  integer(c_int) :: h1,p5,p2,h4
  integer(c_int) :: h6,p7,sh6,sp7
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sp7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh1,sh6,sp5,sp7)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp7,sh4,sh6)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h1 = 1,norb1
     do p5 = norb1+1,nact
     do p2 = norb1+1,nact
     do h4 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i1(h1,p5,p2,h4) = i1(h1,p5,p2,h4) + fact * g2inp(h1,h6,p5,p7,spin_g2inp_1) * t2inp(p2,p7,h4,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_den2_22_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_23(sp1,sp2,sp3,sh4,i0)

!  i0 ( p1 p2 p3 h4 )_ytt + = -1/2 * Sum ( p5 ) * i1 ( p3 p5 p1 p2 )_yt * t ( p5 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,p3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_pppp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_pppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_pppp_1 = tdcc_spin_dummy2(sp3,sp5,sp1,sp2)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh4)
     if(spin_itm_pppp_1 * spin_t1inp_2 == 0) cycle

     itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_23_1(sp3,sp5,sp1,sp2,itm_pppp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact * itm_pppp(p3,p5,p1,p2) * t1inp(p5,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pppp)
end subroutine ccsd_den2_23
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_23_1(sp3,sp5,sp1,sp2,i1)

!      i1 ( p3 p5 p1 p2 )_yt + = +1 * Sum ( h6 h7 ) * y ( h6 h7 p3 p5 )_y * t ( p1 p2 h6 h7 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp5,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p3,p5,p1,p2
  integer(c_int) :: h6,h7,sh6,sh7
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh6,sh7,sp3,sp5)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh6,sh7)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do p3 = norb1+1,nact
     do p5 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h6 = 1,norb1
     do h7 = 1,norb1
        i1(p3,p5,p1,p2) = i1(p3,p5,p1,p2) + fact * g2inp(h6,h7,p3,p5,spin_g2inp_1) * t2inp(p1,p2,h6,h7,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_den2_23_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_24(sh1,sp2,sh3,sh4,i0)

!  i0 ( h1 p2 h3 h4 )_ytt + = +1/2 * Sum ( h5 ) * i1 ( h1 h5 h3 h4 )_yt * t ( p2 h5 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,p2,h3,h4
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hhhh_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh5 = 1,2
     spin_itm_hhhh_1 = tdcc_spin_dummy2(sh1,sh5,sh3,sh4)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh5)
     if(spin_itm_hhhh_1 * spin_t1inp_2 == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsd_den2_24_1(sh1,sh5,sh3,sh4,itm_hhhh)

     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
        i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact * itm_hhhh(h1,h5,h3,h4) * t1inp(p2,h5,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsd_den2_24
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_24_1(sh1,sh5,sh3,sh4,i1)

!      i1 ( h1 h5 h3 h4 )_yt + = +1 * Sum ( p6 p7 ) * y ( h1 h5 p6 p7 )_y * t ( p6 p7 h3 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh5,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h1,h5,h3,h4
  integer(c_int) :: p6,p7,sp6,sp7
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp6 = 1,2
  do sp7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh1,sh5,sp6,sp7)
     spin_t2inp_2 = tdcc_spin_t2inp(sp6,sp7,sh3,sh4)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h1 = 1,norb1
     do h5 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h1,h5,h3,h4) = i1(h1,h5,h3,h4) + fact * g2inp(h1,h5,p6,p7,spin_g2inp_1) * t2inp(p6,p7,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_den2_24_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_25(sp1,sp2,sp3,sh4,i0)

!  i0 ( p1 p2 p3 h4 )_yttt + = -1 * Sum ( p5 ) * i1 ( p3 p5 p1 p2 )_ytt * t ( p5 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,p3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_pppp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_pppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_pppp_1 = tdcc_spin_dummy2(sp3,sp5,sp1,sp2)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh4)
     if(spin_itm_pppp_1 * spin_t1inp_2 == 0) cycle

     itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_25_1(sp3,sp5,sp1,sp2,itm_pppp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact * itm_pppp(p3,p5,p1,p2) * t1inp(p5,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_pppp)
end subroutine ccsd_den2_25
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_25_1(sp3,sp5,sp1,sp2,i1)

!      i1 ( p3 p5 p1 p2 )_ytt + = +1 * Sum ( h6 ) * i2 ( h6 p3 p5 p1 )_yt * t ( p2 h6 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp5,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p3,p5,p1,p2
  integer(c_int) :: h6,sh6
  integer(c_int) :: spin_itm_hppp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sh6 = 1,2
     spin_itm_hppp_1 = tdcc_spin_dummy2(sh6,sp3,sp5,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh6)
     if(spin_itm_hppp_1 * spin_t1inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_25_1_1(sh6,sp3,sp5,sp1,itm_hppp)

     do p3 = norb1+1,nact
     do p5 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h6 = 1,norb1
        i1(p3,p5,p1,p2) = i1(p3,p5,p1,p2) + fact * itm_hppp(h6,p3,p5,p1) * t1inp(p2,h6,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hppp)
end subroutine ccsd_den2_25_1
!##########################################################
!##########################################################
subroutine ccsd_den2_25_1_1(sh6,sp3,sp5,sp1,i2)

!          i2 ( h6 p3 p5 p1 )_yt + = +1 * Sum ( h7 ) * y ( h7 h6 p3 p5 )_y * t ( p1 h7 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh6,sp3,sp5,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h6,p3,p5,p1
  integer(c_int) :: h7,sh7
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh7,sh6,sp3,sp5)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh7)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do p5 = norb1+1,nact
     do p1 = norb1+1,nact
     do h7 = 1,norb1
        i2(h6,p3,p5,p1) = i2(h6,p3,p5,p1) + fact * g2inp(h7,h6,p3,p5,spin_g2inp_1) * t1inp(p1,h7,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_25_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_26(sh1,sp2,sh3,sh4,i0)

!  i0 ( h1 p2 h3 h4 )_yttt + = +1 * Sum ( h5 ) * i1 ( h1 h5 h3 h4 )_ytt * t ( p2 h5 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,p2,h3,h4
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hhhh_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh5 = 1,2
     spin_itm_hhhh_1 = tdcc_spin_dummy2(sh1,sh5,sh3,sh4)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh5)
     if(spin_itm_hhhh_1 * spin_t1inp_2 == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccsd_den2_26_1(sh1,sh5,sh3,sh4,itm_hhhh)

     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
        i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact * itm_hhhh(h1,h5,h3,h4) * t1inp(p2,h5,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhhh)
end subroutine ccsd_den2_26
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_26_1(sh1,sh5,sh3,sh4,i1)

!      i1 ( h1 h5 h3 h4 )_ytt + = +1 * Sum ( p6 ) * i2 ( h1 h5 p6 h3 )_yt * t ( p6 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh5,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h1,h5,h3,h4
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_itm_hhph_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sp6 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh1,sh5,sp6,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp6,sh4)
     if(spin_itm_hhph_1 * spin_t1inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call ccsd_den2_26_1_1(sh1,sh5,sp6,sh3,itm_hhph)

     do h1 = 1,norb1
     do h5 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p6 = norb1+1,nact
        i1(h1,h5,h3,h4) = i1(h1,h5,h3,h4) + fact * itm_hhph(h1,h5,p6,h3) * t1inp(p6,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhph)
end subroutine ccsd_den2_26_1
!##########################################################
!##########################################################
subroutine ccsd_den2_26_1_1(sh1,sh5,sp6,sh3,i2)

!          i2 ( h1 h5 p6 h3 )_yt + = +1 * Sum ( p7 ) * y ( h1 h5 p7 p6 )_y * t ( p7 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh5,sp6,sh3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer(c_int) :: h1,h5,p6,h3
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh1,sh5,sp7,sp6)
     spin_t1inp_2 = tdcc_spin_t1inp(sp7,sh3)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h1 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do h3 = 1,norb1
     do p7 = norb1+1,nact
        i2(h1,h5,p6,h3) = i2(h1,h5,p6,h3) + fact * g2inp(h1,h5,p7,p6,spin_g2inp_1) * t1inp(p7,h3,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_26_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_27(sp1,sh2,sp3,sp4,i0)

!  i0 ( p1 h2 p3 p4 )_yt + = +1 * Sum ( h5 ) * y ( h5 h2 p3 p4 )_y * t ( p1 h5 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sh2,sp3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,h2,p3,p4
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh5 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh2,sp3,sp4)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh5)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do p1 = norb1+1,nact
     do h2 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h5 = 1,norb1
        i0(p1,h2,p3,p4) = i0(p1,h2,p3,p4) + fact * g2inp(h5,h2,p3,p4,spin_g2inp_1) * t1inp(p1,h5,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_27
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_28(sh1,sh2,sh3,sp4,i0)

!  i0 ( h1 h2 h3 p4 )_yt + = -1 * Sum ( p5 ) * y ( h1 h2 p5 p4 )_y * t ( p5 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh2,sh3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,h2,h3,p4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh1,sh2,sp5,sp4)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh3)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i0(h1,h2,h3,p4) = i0(h1,h2,h3,p4) + fact * g2inp(h1,h2,p5,p4,spin_g2inp_1) * t1inp(p5,h3,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_28
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_29(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = +1 * P( h3 h4 ) * P( p1 p2 ) * i1 ( p1 h3 )_yt * t ( p2 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_29_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_29_perm(sp1,sp2,sh4,sh3,i0_perm)
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
     call ccsd_den2_29_perm(sp2,sp1,sh3,sh4,i0_perm)
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
     call ccsd_den2_29_perm(sp2,sp1,sh4,sh3,i0_perm)
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
  subroutine ccsd_den2_29_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_itm_ph_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_ph(:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_ph((norb1+1):nact,1:norb1))
  do sdum = 1,1
     spin_itm_ph_1 = tdcc_spin_dummy1(sp1,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_itm_ph_1 * spin_t1inp_2 == 0) cycle

     itm_ph((norb1+1):nact,1:norb1) = czero
     call ccsd_den2_29_1(sp1,sh3,itm_ph)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_ph(p1,h3) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_ph)
  end subroutine ccsd_den2_29_perm
  !--------------------------------------------
end subroutine ccsd_den2_29
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_29_1(sp1,sh3,i1)

!      i1 ( p1 h3 )_yt + = +1 * Sum ( h5 p6 ) * y ( h5 p6 )_y * t ( p6 p1 h5 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sh3
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer(c_int) :: p1,h3
  integer(c_int) :: h5,p6,sh5,sp6
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh5 = 1,2
  do sp6 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh5,sp6)
     spin_t2inp_2 = tdcc_spin_t2inp(sp6,sp1,sh5,sh3)
     if(spin_g1inp_1 * spin_t2inp_2 == 0) cycle

     do p1 = norb1+1,nact
     do h3 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
        i1(p1,h3) = i1(p1,h3) + fact * g1inp(h5,p6,spin_g1inp_1) * t2inp(p6,p1,h5,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_den2_29_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_30(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1 * P( h3 h4 ) * Sum ( h5 ) * i1 ( h5 h3 )_yt * t ( p1 p2 h5 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_30_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_30_perm(sp1,sp2,sh4,sh3,i0_perm)
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
  subroutine ccsd_den2_30_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hh_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh5 = 1,2
     spin_itm_hh_1 = tdcc_spin_dummy1(sh5,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh4)
     if(spin_itm_hh_1 * spin_t2inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_den2_30_1(sh5,sh3,itm_hh)

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
  end subroutine ccsd_den2_30_perm
  !--------------------------------------------
end subroutine ccsd_den2_30
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_30_1(sh5,sh3,i1)

!      i1 ( h5 h3 )_yt + = +1 * Sum ( p6 ) * y ( h5 p6 )_y * t ( p6 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh5,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h5,h3
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp6 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh5,sp6)
     spin_t1inp_2 = tdcc_spin_t1inp(sp6,sh3)
     if(spin_g1inp_1 * spin_t1inp_2 == 0) cycle

     do h5 = 1,norb1
     do h3 = 1,norb1
     do p6 = norb1+1,nact
        i1(h5,h3) = i1(h5,h3) + fact * g1inp(h5,p6,spin_g1inp_1) * t1inp(p6,h3,spin_t1inp_2)
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_30_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_31(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1 * P( p1 p2 ) * Sum ( p5 ) * i1 ( p5 p1 )_yt * t ( p5 p2 h3 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_31_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_31_perm(sp2,sp1,sh3,sh4,i0_perm)
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
  subroutine ccsd_den2_31_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_pp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_pp_1 = tdcc_spin_dummy1(sp5,sp1)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp2,sh3,sh4)
     if(spin_itm_pp_1 * spin_t2inp_2 == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_31_1(sp5,sp1,itm_pp)

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
  end subroutine ccsd_den2_31_perm
  !--------------------------------------------
end subroutine ccsd_den2_31
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_31_1(sp5,sp1,i1)

!      i1 ( p5 p1 )_yt + = +1 * Sum ( h6 ) * y ( h6 p5 )_y * t ( p1 h6 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sp5,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p5,p1
  integer(c_int) :: h6,sh6
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh6,sp5)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh6)
     if(spin_g1inp_1 * spin_t1inp_2 == 0) cycle

     do p5 = norb1+1,nact
     do p1 = norb1+1,nact
     do h6 = 1,norb1
        i1(p5,p1) = i1(p5,p1) + fact * g1inp(h6,p5,spin_g1inp_1) * t1inp(p1,h6,spin_t1inp_2)
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_31_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_32(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_yttt + = -1/2 * P( h3 h4 ) * P( p1 p2 ) * i1 ( h3 p1 )_ytt * t ( p2 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_32_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_32_perm(sp1,sp2,sh4,sh3,i0_perm)
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
     call ccsd_den2_32_perm(sp2,sp1,sh3,sh4,i0_perm)
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
     call ccsd_den2_32_perm(sp2,sp1,sh4,sh3,i0_perm)
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
  subroutine ccsd_den2_32_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_itm_hp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sdum = 1,1
     spin_itm_hp_1 = tdcc_spin_dummy1(sh3,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_itm_hp_1 * spin_t1inp_2 == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsd_den2_32_1(sh3,sp1,itm_hp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_hp(h3,p1) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hp)
  end subroutine ccsd_den2_32_perm
  !--------------------------------------------
end subroutine ccsd_den2_32
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_32_1(sh3,sp1,i1)

!      i1 ( h3 p1 )_ytt + = +1 * Sum ( h5 ) * i2 ( h5 h3 )_yt * t ( p1 h5 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh3,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,p1
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hh_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh5 = 1,2
     spin_itm_hh_1 = tdcc_spin_dummy1(sh5,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh5)
     if(spin_itm_hh_1 * spin_t1inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_den2_32_1_1(sh5,sh3,itm_hh)

     do h3 = 1,norb1
     do p1 = norb1+1,nact
     do h5 = 1,norb1
        i1(h3,p1) = i1(h3,p1) + fact * itm_hh(h5,h3) * t1inp(p1,h5,spin_t1inp_2)
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
end subroutine ccsd_den2_32_1
!##########################################################
!##########################################################
subroutine ccsd_den2_32_1_1(sh5,sh3,i2)

!          i2 ( h5 h3 )_yt + = +1 * Sum ( h6 p7 p8 ) * y ( h6 h5 p7 p8 )_y * t ( p7 p8 h6 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh5,sh3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1)
  integer(c_int) :: h5,h3
  integer(c_int) :: h6,p7,p8,sh6,sp7,sp8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
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
        i2(h5,h3) = i2(h5,h3) + fact * g2inp(h6,h5,p7,p8,spin_g2inp_1) * t2inp(p7,p8,h6,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_den2_32_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_33(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_yttt + = -1/2 * P( h3 h4 ) * P( p1 p2 ) * i1 ( p1 h3 )_ytt * t ( p2 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_33_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_33_perm(sp1,sp2,sh4,sh3,i0_perm)
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
     call ccsd_den2_33_perm(sp2,sp1,sh3,sh4,i0_perm)
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
     call ccsd_den2_33_perm(sp2,sp1,sh4,sh3,i0_perm)
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
  subroutine ccsd_den2_33_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: sdum
  integer(c_int) :: spin_itm_ph_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_ph(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_ph((norb1+1):nact,1:norb1))
  do sdum = 1,1
     spin_itm_ph_1 = tdcc_spin_dummy1(sp1,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh4)
     if(spin_itm_ph_1 * spin_t1inp_2 == 0) cycle

     itm_ph((norb1+1):nact,1:norb1) = czero
     call ccsd_den2_33_1(sp1,sh3,itm_ph)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_ph(p1,h3) * t1inp(p2,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_ph)
  end subroutine ccsd_den2_33_perm
  !--------------------------------------------
end subroutine ccsd_den2_33
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_33_1(sp1,sh3,i1)

!      i1 ( p1 h3 )_ytt + = +1 * Sum ( p5 ) * i2 ( p5 p1 )_yt * t ( p5 h3 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sh3
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer(c_int) :: p1,h3
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_pp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_pp_1 = tdcc_spin_dummy1(sp5,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh3)
     if(spin_itm_pp_1 * spin_t1inp_2 == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_33_1_1(sp5,sp1,itm_pp)

     do p1 = norb1+1,nact
     do h3 = 1,norb1
     do p5 = norb1+1,nact
        i1(p1,h3) = i1(p1,h3) + fact * itm_pp(p5,p1) * t1inp(p5,h3,spin_t1inp_2)
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
end subroutine ccsd_den2_33_1
!##########################################################
!##########################################################
subroutine ccsd_den2_33_1_1(sp5,sp1,i2)

!          i2 ( p5 p1 )_yt + = +1 * Sum ( h6 h7 p8 ) * y ( h6 h7 p8 p5 )_y * t ( p8 p1 h6 h7 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp5,sp1
  complex(kind(0d0)),intent(inout) :: i2((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p5,p1
  integer(c_int) :: h6,h7,p8,sh6,sh7,sp8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
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
        i2(p5,p1) = i2(p5,p1) + fact * g2inp(h6,h7,p8,p5,spin_g2inp_1) * t2inp(p8,p1,h6,h7,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_den2_33_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_34(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_yttt + = -1 * P( h3 h4 ) * P( p1 p2 ) * Sum ( h5 p6 ) * i1 ( h5 p6 p1 h3 )_ytt * t ( p2 p6 h4 h5 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccsd_den2_34_perm(sp1,sp2,sh3,sh4,i0_perm)
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
     call ccsd_den2_34_perm(sp1,sp2,sh4,sh3,i0_perm)
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
     call ccsd_den2_34_perm(sp2,sp1,sh3,sh4,i0_perm)
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
     call ccsd_den2_34_perm(sp2,sp1,sh4,sh3,i0_perm)
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
  subroutine ccsd_den2_34_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: h5,p6,sh5,sp6
  integer(c_int) :: spin_itm_hpph_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hpph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1))
  do sh5 = 1,2
  do sp6 = 1,2
     spin_itm_hpph_1 = tdcc_spin_dummy2(sh5,sp6,sp1,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp6,sh4,sh5)
     if(spin_itm_hpph_1 * spin_t2inp_2 == 0) cycle

     itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1) = czero
     call ccsd_den2_34_1(sh5,sp6,sp1,sh3,itm_hpph)

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
  end subroutine ccsd_den2_34_perm
  !--------------------------------------------
end subroutine ccsd_den2_34
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_34_1(sh5,sp6,sp1,sh3,i1)

!      i1 ( h5 p6 p1 h3 )_ytt + = +1 * Sum ( p7 ) * i2 ( h5 p7 p6 p1 )_yt * t ( p7 h3 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh5,sp6,sp1,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)
  integer(c_int) :: h5,p6,p1,h3
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_itm_hppp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
     spin_itm_hppp_1 = tdcc_spin_dummy2(sh5,sp7,sp6,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp7,sh3)
     if(spin_itm_hppp_1 * spin_t1inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_34_1_1(sh5,sp7,sp6,sp1,itm_hppp)

     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do h3 = 1,norb1
     do p7 = norb1+1,nact
        i1(h5,p6,p1,h3) = i1(h5,p6,p1,h3) + fact * itm_hppp(h5,p7,p6,p1) * t1inp(p7,h3,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hppp)
end subroutine ccsd_den2_34_1
!##########################################################
!##########################################################
subroutine ccsd_den2_34_1_1(sh5,sp7,sp6,sp1,i2)

!          i2 ( h5 p7 p6 p1 )_yt + = +1 * Sum ( h8 ) * y ( h8 h5 p7 p6 )_y * t ( p1 h8 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh5,sp7,sp6,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h5,p7,p6,p1
  integer(c_int) :: h8,sh8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh8,sh5,sp7,sp6)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh8)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h5 = 1,norb1
     do p7 = norb1+1,nact
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do h8 = 1,norb1
        i2(h5,p7,p6,p1) = i2(h5,p7,p6,p1) + fact * g2inp(h8,h5,p7,p6,spin_g2inp_1) * t1inp(p1,h8,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_34_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_35(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_yttt + = +1/2 * Sum ( p5 ) * i1 ( p5 h3 p1 p2 )_ytt * t ( p5 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_phpp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_phpp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  allocate(itm_phpp((norb1+1):nact,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_phpp_1 = tdcc_spin_dummy2(sp5,sh3,sp1,sp2)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh4)
     if(spin_itm_phpp_1 * spin_t1inp_2 == 0) cycle

     itm_phpp((norb1+1):nact,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_35_1(sp5,sh3,sp1,sp2,itm_phpp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_phpp(p5,h3,p1,p2) * t1inp(p5,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_phpp)
end subroutine ccsd_den2_35
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_35_1(sp5,sh3,sp1,sp2,i1)

!      i1 ( p5 h3 p1 p2 )_ytt + = +1 * Sum ( h6 h7 ) * i2 ( h6 h7 p5 h3 )_yt * t ( p1 p2 h6 h7 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp5,sh3,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p5,h3,p1,p2
  integer(c_int) :: h6,h7,sh6,sh7
  integer(c_int) :: spin_itm_hhph_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sh6 = 1,2
  do sh7 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh6,sh7,sp5,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh6,sh7)
     if(spin_itm_hhph_1 * spin_t2inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call ccsd_den2_35_1_1(sh6,sh7,sp5,sh3,itm_hhph)

     do p5 = norb1+1,nact
     do h3 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h6 = 1,norb1
     do h7 = 1,norb1
        i1(p5,h3,p1,p2) = i1(p5,h3,p1,p2) + fact * itm_hhph(h6,h7,p5,h3) * t2inp(p1,p2,h6,h7,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hhph)
end subroutine ccsd_den2_35_1
!##########################################################
!##########################################################
subroutine ccsd_den2_35_1_1(sh6,sh7,sp5,sh3,i2)

!          i2 ( h6 h7 p5 h3 )_yt + = +1 * Sum ( p8 ) * y ( h6 h7 p8 p5 )_y * t ( p8 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh6,sh7,sp5,sh3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer(c_int) :: h6,h7,p5,h3
  integer(c_int) :: p8,sp8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh6,sh7,sp8,sp5)
     spin_t1inp_2 = tdcc_spin_t1inp(sp8,sh3)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h6 = 1,norb1
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do h3 = 1,norb1
     do p8 = norb1+1,nact
        i2(h6,h7,p5,h3) = i2(h6,h7,p5,h3) + fact * g2inp(h6,h7,p8,p5,spin_g2inp_1) * t1inp(p8,h3,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_35_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_36(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_yttt + = +1/2 * Sum ( h5 ) * i1 ( h5 p1 h3 h4 )_ytt * t ( p2 h5 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hphh_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh5 = 1,2
     spin_itm_hphh_1 = tdcc_spin_dummy2(sh5,sp1,sh3,sh4)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh5)
     if(spin_itm_hphh_1 * spin_t1inp_2 == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccsd_den2_36_1(sh5,sp1,sh3,sh4,itm_hphh)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_hphh(h5,p1,h3,h4) * t1inp(p2,h5,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphh)
end subroutine ccsd_den2_36
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_36_1(sh5,sp1,sh3,sh4,i1)

!      i1 ( h5 p1 h3 h4 )_ytt + = +1 * Sum ( p6 p7 ) * i2 ( h5 p6 p7 p1 )_yt * t ( p6 p7 h3 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sh5,sp1,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h5,p1,h3,h4
  integer(c_int) :: p6,p7,sp6,sp7
  integer(c_int) :: spin_itm_hppp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp6 = 1,2
  do sp7 = 1,2
     spin_itm_hppp_1 = tdcc_spin_dummy2(sh5,sp6,sp7,sp1)
     spin_t2inp_2 = tdcc_spin_t2inp(sp6,sp7,sh3,sh4)
     if(spin_itm_hppp_1 * spin_t2inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_36_1_1(sh5,sp6,sp7,sp1,itm_hppp)

     do h5 = 1,norb1
     do p1 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h5,p1,h3,h4) = i1(h5,p1,h3,h4) + fact * itm_hppp(h5,p6,p7,p1) * t2inp(p6,p7,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  deallocate(itm_hppp)
end subroutine ccsd_den2_36_1
!##########################################################
!##########################################################
subroutine ccsd_den2_36_1_1(sh5,sp6,sp7,sp1,i2)

!          i2 ( h5 p6 p7 p1 )_yt + = +1 * Sum ( h8 ) * y ( h8 h5 p6 p7 )_y * t ( p1 h8 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh5,sp6,sp7,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h5,p6,p7,p1
  integer(c_int) :: h8,sh8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh8,sh5,sp6,sp7)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh8)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
     do p1 = norb1+1,nact
     do h8 = 1,norb1
        i2(h5,p6,p7,p1) = i2(h5,p6,p7,p1) + fact * g2inp(h8,h5,p6,p7,spin_g2inp_1) * t1inp(p1,h8,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_36_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_37(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytttt + = +1 * Sum ( p5 ) * i1 ( p5 h3 p1 p2 )_yttt * t ( p5 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: p1,p2,h3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_phpp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_phpp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_phpp((norb1+1):nact,1:norb1,(norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_phpp_1 = tdcc_spin_dummy2(sp5,sh3,sp1,sp2)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh4)
     if(spin_itm_phpp_1 * spin_t1inp_2 == 0) cycle

     itm_phpp((norb1+1):nact,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_37_1(sp5,sh3,sp1,sp2,itm_phpp)

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * itm_phpp(p5,h3,p1,p2) * t1inp(p5,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_phpp)
end subroutine ccsd_den2_37
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_37_1(sp5,sh3,sp1,sp2,i1)

!      i1 ( p5 h3 p1 p2 )_yttt + = +1 * Sum ( h6 ) * i2 ( h6 p5 h3 p1 )_ytt * t ( p2 h6 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp5,sh3,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p5,h3,p1,p2
  integer(c_int) :: h6,sh6
  integer(c_int) :: spin_itm_hphp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh6 = 1,2
     spin_itm_hphp_1 = tdcc_spin_dummy2(sh6,sp5,sh3,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh6)
     if(spin_itm_hphp_1 * spin_t1inp_2 == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccsd_den2_37_1_1(sh6,sp5,sh3,sp1,itm_hphp)

     do p5 = norb1+1,nact
     do h3 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h6 = 1,norb1
        i1(p5,h3,p1,p2) = i1(p5,h3,p1,p2) + fact * itm_hphp(h6,p5,h3,p1) * t1inp(p2,h6,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphp)
end subroutine ccsd_den2_37_1
!##########################################################
!##########################################################
subroutine ccsd_den2_37_1_1(sh6,sp5,sh3,sp1,i2)

!          i2 ( h6 p5 h3 p1 )_ytt + = +1 * Sum ( h7 ) * i3 ( h7 h6 p5 h3 )_yt * t ( p1 h7 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh6,sp5,sh3,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,p5,h3,p1
  integer(c_int) :: h7,sh7
  integer(c_int) :: spin_itm_hhph_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sh7 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh7,sh6,sp5,sh3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh7)
     if(spin_itm_hhph_1 * spin_t1inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call ccsd_den2_37_1_1_1(sh7,sh6,sp5,sh3,itm_hhph)

     do h6 = 1,norb1
     do p5 = norb1+1,nact
     do h3 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
        i2(h6,p5,h3,p1) = i2(h6,p5,h3,p1) + fact * itm_hhph(h7,h6,p5,h3) * t1inp(p1,h7,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hhph)
end subroutine ccsd_den2_37_1_1
!##########################################################
subroutine ccsd_den2_37_1_1_1(sh7,sh6,sp5,sh3,i3)

!              i3 ( h7 h6 p5 h3 )_yt + = +1 * Sum ( p8 ) * y ( h7 h6 p8 p5 )_y * t ( p8 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh7,sh6,sp5,sh3
  complex(kind(0d0)),intent(inout) :: i3(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer(c_int) :: h7,h6,p5,h3
  integer(c_int) :: p8,sp8
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh7,sh6,sp8,sp5)
     spin_t1inp_2 = tdcc_spin_t1inp(sp8,sh3)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h7 = 1,norb1
     do h6 = 1,norb1
     do p5 = norb1+1,nact
     do h3 = 1,norb1
     do p8 = norb1+1,nact
        i3(h7,h6,p5,h3) = i3(h7,h6,p5,h3) + fact * g2inp(h7,h6,p8,p5,spin_g2inp_1) * t1inp(p8,h3,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_37_1_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_38(sh1,sp2,sp3,sh4,i0)

!  i0 ( h1 p2 p3 h4 )_ytt + = -1 * Sum ( p5 ) * i1 ( h1 p3 p5 p2 )_yt * t ( p5 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: h1,p2,p3,h4
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_itm_hppp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_hppp_1 = tdcc_spin_dummy2(sh1,sp3,sp5,sp2)
     spin_t1inp_2 = tdcc_spin_t1inp(sp5,sh4)
     if(spin_itm_hppp_1 * spin_t1inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den2_38_1(sh1,sp3,sp5,sp2,itm_hppp)

     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(h1,p2,p3,h4) = i0(h1,p2,p3,h4) + fact * itm_hppp(h1,p3,p5,p2) * t1inp(p5,h4,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hppp)
end subroutine ccsd_den2_38
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den2_38_1(sh1,sp3,sp5,sp2,i1)

!      i1 ( h1 p3 p5 p2 )_yt + = +1 * Sum ( h6 ) * y ( h1 h6 p3 p5 )_y * t ( p2 h6 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp3,sp5,sp2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h1,p3,p5,p2
  integer(c_int) :: h6,sh6
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh1,sh6,sp3,sp5)
     spin_t1inp_2 = tdcc_spin_t1inp(sp2,sh6)
     if(spin_g2inp_1 * spin_t1inp_2 == 0) cycle

     do h1 = 1,norb1
     do p3 = norb1+1,nact
     do p5 = norb1+1,nact
     do p2 = norb1+1,nact
     do h6 = 1,norb1
        i1(h1,p3,p5,p2) = i1(h1,p3,p5,p2) + fact * g2inp(h1,h6,p3,p5,spin_g2inp_1) * t1inp(p2,h6,spin_t1inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccsd_den2_38_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccsd_den2_main(den)

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact

  implicit none
  complex(kind(0d0)),intent(inout) :: den(1:nact,1:nact,1:nact,1:nact,1:*)

  call ccsd_den2_1(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_2(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_3(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_4(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_5(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_6(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_7(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_8(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_9(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_10(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_11(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_12(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_13(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_14(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_15(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_16(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_17(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_18(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_19(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_20(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_21(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_22(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_23(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_24(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_25(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_26(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_27(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_28(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_29(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_30(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_31(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_32(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_33(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_34(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_35(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_36(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_37(1,1,1,1,den(1,1,1,1,1))
  call ccsd_den2_38(1,1,1,1,den(1,1,1,1,1))

  call ccsd_den2_1(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_2(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_3(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_4(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_5(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_6(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_7(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_8(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_9(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_10(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_11(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_12(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_13(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_14(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_15(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_16(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_17(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_18(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_19(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_20(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_21(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_22(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_23(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_24(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_25(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_26(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_27(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_28(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_29(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_30(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_31(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_32(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_33(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_34(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_35(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_36(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_37(1,2,1,2,den(1,1,1,1,3))
  call ccsd_den2_38(1,2,1,2,den(1,1,1,1,3))

  call ccsd_den2_1(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_2(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_3(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_4(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_5(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_6(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_7(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_8(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_9(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_10(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_11(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_12(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_13(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_14(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_15(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_16(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_17(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_18(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_19(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_20(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_21(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_22(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_23(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_24(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_25(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_26(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_27(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_28(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_29(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_30(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_31(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_32(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_33(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_34(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_35(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_36(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_37(1,2,2,1,den(1,1,1,1,5))
  call ccsd_den2_38(1,2,2,1,den(1,1,1,1,5))

end subroutine ccsd_den2_main
!**********************************************************
