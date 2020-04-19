!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_1(sp1,sp2,sp3,sp4,i0)

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

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
        i0(p1,p2,p3,p4) = i0(p1,p2,p3,p4) + fact * &
             g2inp(h5,h6,p3,p4,spin_g2inp_1) * t2inp(p1,p2,h5,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine ccdt_den2p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_2(sh1,sh2,sh3,sh4,i0)

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

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(h1,h2,h3,h4) = i0(h1,h2,h3,h4) + fact * &
             g2inp(h1,h2,p5,p6,spin_g2inp_1) * t2inp(p5,p6,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine ccdt_den2p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_3(sp1,sp2,sh3,sh4,i0)

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
  !$omp parallel default(shared) private(p1,p2,h3,h4)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
!     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * t2inp(p1,p2,h3,h4,spin_t2inp_1)
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + t2inp(p1,p2,h3,h4,spin_t2inp_1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_den2p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_4(sp1,sp2,sh3,sh4,i0)

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
  call ccdt_den2p_4_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_4_perm(sp1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h4,h3)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_4_perm(sp2,sp1,sh3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_4_perm(sp2,sp1,sh4,sh3,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h4,h3)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_4_perm(sp1,sp2,sh3,sh4,i0)

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
     call ccdt_den2p_4_1(sh5,sp6,sp1,sh3,itm_hpph)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_hpph(h5,p6,p1,h3) * t2inp(p2,p6,h4,h5,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hpph)
  end subroutine ccdt_den2p_4_perm
  !--------------------------------------------
end subroutine ccdt_den2p_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_4_1(sh5,sp6,sp1,sh3,i1)

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

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do h3 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
!        i1(h5,p6,p1,h3) = i1(h5,p6,p1,h3) + fact * &
!        g2inp(h7,h5,p8,p6,spin_g2inp_1) * t2inp(p8,p1,h7,h3,spin_t2inp_2)
        i1(h5,p6,p1,h3) = i1(h5,p6,p1,h3) + &
             g2inp(h7,h5,p8,p6,spin_g2inp_1) * t2inp(p8,p1,h7,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine ccdt_den2p_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_5(sp1,sp2,sh3,sh4,i0)

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
  call ccdt_den2p_5_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_5_perm(sp1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h4,h3)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_5_perm(sp1,sp2,sh3,sh4,i0)

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
     call ccdt_den2p_5_1(sh5,sh3,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_hh(h5,h3) * t2inp(p1,p2,h5,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hh)
  end subroutine ccdt_den2p_5_perm
  !--------------------------------------------
end subroutine ccdt_den2p_5
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_5_1(sh5,sh3,i1)

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

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h5 = 1,norb1
     do h3 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
!        i1(h5,h3) = i1(h5,h3) + fact * &
!        g2inp(h6,h5,p7,p8,spin_g2inp_1) * t2inp(p7,p8,h6,h3,spin_t2inp_2)
        i1(h5,h3) = i1(h5,h3) + &
             g2inp(h6,h5,p7,p8,spin_g2inp_1) * t2inp(p7,p8,h6,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_6(sp1,sp2,sh3,sh4,i0)

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
  call ccdt_den2p_6_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_6_perm(sp2,sp1,sh3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_6_perm(sp1,sp2,sh3,sh4,i0)

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
     call ccdt_den2p_6_1(sp5,sp1,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_pp(p5,p1) * t2inp(p5,p2,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_pp)
  end subroutine ccdt_den2p_6_perm
  !--------------------------------------------
end subroutine ccdt_den2p_6
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_6_1(sp5,sp1,i1)

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

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p5 = norb1+1,nact
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
!        i1(p5,p1) = i1(p5,p1) + fact * &
!        g2inp(h6,h7,p8,p5,spin_g2inp_1) * t2inp(p8,p1,h6,h7,spin_t2inp_2)
        i1(p5,p1) = i1(p5,p1) + &
             g2inp(h6,h7,p8,p5,spin_g2inp_1) * t2inp(p8,p1,h6,h7,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_7(sp1,sp2,sh3,sh4,i0)

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
     call ccdt_den2p_7_1(sh5,sh6,sh3,sh4,itm_hhhh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_hhhh(h5,h6,h3,h4) * t2inp(p1,p2,h5,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccdt_den2p_7
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_7_1(sh5,sh6,sh3,sh4,i1)

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

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
!        i1(h5,h6,h3,h4) = i1(h5,h6,h3,h4) + fact * &
!        g2inp(h5,h6,p7,p8,spin_g2inp_1) * t2inp(p7,p8,h3,h4,spin_t2inp_2)
        i1(h5,h6,h3,h4) = i1(h5,h6,h3,h4) + &
             g2inp(h5,h6,p7,p8,spin_g2inp_1) * t2inp(p7,p8,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine ccdt_den2p_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_8(sh1,sp2,sp3,sh4,i0)

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

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
!        i0(h1,p2,p3,h4) = i0(h1,p2,p3,h4) + fact * &
!        g2inp(h5,h1,p6,p3,spin_g2inp_1) * t2inp(p6,p2,h5,h4,spin_t2inp_2)
        i0(h1,p2,p3,h4) = i0(h1,p2,p3,h4) + &
             g2inp(h5,h1,p6,p3,spin_g2inp_1) * t2inp(p6,p2,h5,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine ccdt_den2p_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_9(sh1,sh2,sp3,sp4,i0)

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
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h1 = 1,norb1
  do h2 = 1,norb1
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
!     i0(h1,h2,p3,p4) = i0(h1,h2,p3,p4) + fact * &
!     g2inp(h1,h2,p3,p4,spin_g2inp_1)
     i0(h1,h2,p3,p4) = i0(h1,h2,p3,p4) + g2inp(h1,h2,p3,p4,spin_g2inp_1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_den2p_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_10(sp1,sp2,sp3,sp4,i0)

!  i0 ( p1 p2 p3 p4 )_yt + = +1/6 * Sum ( h5 h6 h7 p8 ) * y ( h5 h6 h7 p3 p4 p8 )_y * t ( p1 p2 p8 h5 h6 h7 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sp1,sp2,sp3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,p3,p4
  integer :: h5,h6,h7,p8,sh5,sh6,sh7,sp8
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 6.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sh7 = 1,2
  do sp8 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh5,sh6,sh7,sp3,sp4,sp8)
     spin_t3inp_2 = tdcc_spin_t3inp(sp1,sp2,sp8,sh5,sh6,sh7)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
        i0(p1,p2,p3,p4) = i0(p1,p2,p3,p4) + fact * &
             g3inp(h5,h6,h7,p3,p4,p8,spin_g3inp_1) * t3inp(p1,p2,p8,h5,h6,h7,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  end do
end subroutine ccdt_den2p_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_11(sh1,sh2,sh3,sh4,i0)

!  i0 ( h1 h2 h3 h4 )_yt + = +1/6 * Sum ( h5 p6 p7 p8 ) * y ( h1 h2 h5 p6 p7 p8 )_y * t ( p6 p7 p8 h3 h4 h5 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sh1,sh2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: h1,h2,h3,h4
  integer :: h5,p6,p7,p8,sh5,sp6,sp7,sp8
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 6.0d+0 * runit

  do sh5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh1,sh2,sh5,sp6,sp7,sp8)
     spin_t3inp_2 = tdcc_spin_t3inp(sp6,sp7,sp8,sh3,sh4,sh5)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(h1,h2,h3,h4) = i0(h1,h2,h3,h4) + fact * &
             g3inp(h1,h2,h5,p6,p7,p8,spin_g3inp_1) * t3inp(p6,p7,p8,h3,h4,h5,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  end do
end subroutine ccdt_den2p_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_12(sp1,sp2,sp3,sh4,i0)

!  i0 ( p1 p2 p3 h4 )_yt + = +1/2 * Sum ( h5 h6 p7 ) * y ( h5 h6 p7 p3 )_y * t ( p7 p1 p2 h5 h6 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,p3,h4
  integer :: h5,h6,p7,sh5,sh6,sp7
  integer :: spin_g2inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh6,sp7,sp3)
     spin_t3inp_2 = tdcc_spin_t3inp(sp7,sp1,sp2,sh5,sh6,sh4)
     if(spin_g2inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact * &
             g2inp(h5,h6,p7,p3,spin_g2inp_1) * t3inp(p7,p1,p2,h5,h6,h4,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_12
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_13(sh1,sp2,sh3,sh4,i0)

!  i0 ( h1 p2 h3 h4 )_yt + = -1/2 * Sum ( h5 p6 p7 ) * y ( h5 h1 p6 p7 )_y * t ( p6 p7 p2 h5 h3 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: h1,p2,h3,h4
  integer :: h5,p6,p7,sh5,sp6,sp7
  integer :: spin_g2inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh1,sp6,sp7)
     spin_t3inp_2 = tdcc_spin_t3inp(sp6,sp7,sp2,sh5,sh3,sh4)
     if(spin_g2inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact * &
             g2inp(h5,h1,p6,p7,spin_g2inp_1) * t3inp(p6,p7,p2,h5,h3,h4,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_13
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_14(sp1,sp2,sp3,sh4,i0)

!  i0 ( p1 p2 p3 h4 )_ytt + = -1/4 * Sum ( h5 h6 ) * i1 ( h5 h6 p3 h4 )_yt * t ( p1 p2 h5 h6 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,p3,h4
  integer :: h5,h6,sh5,sh6
  integer :: spin_itm_hhph_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sh5 = 1,2
  do sh6 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh5,sh6,sp3,sh4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh6)
     if(spin_itm_hhph_1 * spin_t2inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call ccdt_den2p_14_1(sh5,sh6,sp3,sh4,itm_hhph)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
        i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact * &
             itm_hhph(h5,h6,p3,h4) * t2inp(p1,p2,h5,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hhph)
end subroutine ccdt_den2p_14
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_14_1(sh5,sh6,sp3,sh4,i1)

!      i1 ( h5 h6 p3 h4 )_yt + = +1 * Sum ( h7 p8 p9 ) * y ( h5 h6 h7 p3 p8 p9 )_y * t ( p8 p9 h4 h7 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh5,sh6,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer :: h5,h6,p3,h4
  integer :: h7,p8,p9,sh7,sp8,sp9
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh7 = 1,2
  do sp8 = 1,2
  do sp9 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh5,sh6,sh7,sp3,sp8,sp9)
     spin_t2inp_2 = tdcc_spin_t2inp(sp8,sp9,sh4,sh7)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
!        i1(h5,h6,p3,h4) = i1(h5,h6,p3,h4) + fact * &
!        g3inp(h5,h6,h7,p3,p8,p9,spin_g3inp_1) * t2inp(p8,p9,h4,h7,spin_t2inp_2)
        i1(h5,h6,p3,h4) = i1(h5,h6,p3,h4) + &
             g3inp(h5,h6,h7,p3,p8,p9,spin_g3inp_1) * t2inp(p8,p9,h4,h7,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_14_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_15(sh1,sp2,sh3,sh4,i0)

!  i0 ( h1 p2 h3 h4 )_ytt + = +1/4 * Sum ( p5 p6 ) * i1 ( h1 p5 p6 p2 )_yt * t ( p5 p6 h3 h4 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: h1,p2,h3,h4
  integer :: p5,p6,sp5,sp6
  integer :: spin_itm_hppp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
  do sp6 = 1,2
     spin_itm_hppp_1 = tdcc_spin_dummy2(sh1,sp5,sp6,sp2)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp6,sh3,sh4)
     if(spin_itm_hppp_1 * spin_t2inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_den2p_15_1(sh1,sp5,sp6,sp2,itm_hppp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact * &
             itm_hppp(h1,p5,p6,p2) * t2inp(p5,p6,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hppp)
end subroutine ccdt_den2p_15
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_15_1(sh1,sp5,sp6,sp2,i1)

!      i1 ( h1 p5 p6 p2 )_yt + = +1 * Sum ( h7 h8 p9 ) * y ( h1 h7 h8 p5 p6 p9 )_y * t ( p2 p9 h7 h8 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh1,sp5,sp6,sp2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h1,p5,p6,p2
  integer :: h7,h8,p9,sh7,sh8,sp9
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp9 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh1,sh7,sh8,sp5,sp6,sp9)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp9,sh7,sh8)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do p2 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p9 = norb1+1,nact
!        i1(h1,p5,p6,p2) = i1(h1,p5,p6,p2) + fact * &
!        g3inp(h1,h7,h8,p5,p6,p9,spin_g3inp_1) * t2inp(p2,p9,h7,h8,spin_t2inp_2)
        i1(h1,p5,p6,p2) = i1(h1,p5,p6,p2) + &
             g3inp(h1,h7,h8,p5,p6,p9,spin_g3inp_1) * t2inp(p2,p9,h7,h8,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_15_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_16(sp1,sp2,sp3,sh4,i0)

!  i0 ( p1 p2 p3 h4 )_ytt + = +1/2 * P( p1 p2 ) * Sum ( h5 p6 ) * i1 ( h5 p3 p6 p1 )_yt * t ( p2 p6 h4 h5 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,p3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccdt_den2p_16_perm(sp1,sp2,sp3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
  do h4 = 1,norb1
     i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact_p * i0_perm(p1,p2,p3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sp3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_16_perm(sp2,sp1,sp3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
  do h4 = 1,norb1
     i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact_p * i0_perm(p2,p1,p3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_16_perm(sp1,sp2,sp3,sh4,i0)

  implicit none
  integer,intent(in) :: sp1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer :: p1,p2,p3,h4
  integer :: h5,p6,sh5,sp6
  integer :: spin_itm_hppp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sh5 = 1,2
  do sp6 = 1,2
     spin_itm_hppp_1 = tdcc_spin_dummy2(sh5,sp3,sp6,sp1)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp6,sh4,sh5)
     if(spin_itm_hppp_1 * spin_t2inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_den2p_16_1(sh5,sp3,sp6,sp1,itm_hppp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
        i0(p1,p2,p3,h4) = i0(p1,p2,p3,h4) + fact * &
             itm_hppp(h5,p3,p6,p1) * t2inp(p2,p6,h4,h5,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hppp)
  end subroutine ccdt_den2p_16_perm
  !--------------------------------------------
end subroutine ccdt_den2p_16
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_16_1(sh5,sp3,sp6,sp1,i1)

!      i1 ( h5 p3 p6 p1 )_yt + = +1 * Sum ( h7 h8 p9 ) * y ( h7 h8 h5 p9 p3 p6 )_y * t ( p9 p1 h7 h8 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh5,sp3,sp6,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h5,p3,p6,p1
  integer :: h7,h8,p9,sh7,sh8,sp9
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp9 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh7,sh8,sh5,sp9,sp3,sp6)
     spin_t2inp_2 = tdcc_spin_t2inp(sp9,sp1,sh7,sh8)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p9 = norb1+1,nact
!        i1(h5,p3,p6,p1) = i1(h5,p3,p6,p1) + fact * &
!        g3inp(h7,h8,h5,p9,p3,p6,spin_g3inp_1) * t2inp(p9,p1,h7,h8,spin_t2inp_2)
        i1(h5,p3,p6,p1) = i1(h5,p3,p6,p1) + &
             g3inp(h7,h8,h5,p9,p3,p6,spin_g3inp_1) * t2inp(p9,p1,h7,h8,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_16_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_17(sh1,sp2,sh3,sh4,i0)

!  i0 ( h1 p2 h3 h4 )_ytt + = -1/2 * P( h3 h4 ) * Sum ( h5 p6 ) * i1 ( h1 h5 p6 h3 )_yt * t ( p2 p6 h4 h5 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: h1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccdt_den2p_17_perm(sh1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h1 = 1,norb1
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact_p * i0_perm(h1,p2,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_17_perm(sh1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h1 = 1,norb1
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact_p * i0_perm(h1,p2,h4,h3)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_17_perm(sh1,sp2,sh3,sh4,i0)

  implicit none
  integer,intent(in) :: sh1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer :: h1,p2,h3,h4
  integer :: h5,p6,sh5,sp6
  integer :: spin_itm_hhph_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sh5 = 1,2
  do sp6 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh1,sh5,sp6,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp6,sh4,sh5)
     if(spin_itm_hhph_1 * spin_t2inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call ccdt_den2p_17_1(sh1,sh5,sp6,sh3,itm_hhph)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
        i0(h1,p2,h3,h4) = i0(h1,p2,h3,h4) + fact * &
             itm_hhph(h1,h5,p6,h3) * t2inp(p2,p6,h4,h5,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hhph)
  end subroutine ccdt_den2p_17_perm
  !--------------------------------------------
end subroutine ccdt_den2p_17
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_17_1(sh1,sh5,sp6,sh3,i1)

!      i1 ( h1 h5 p6 h3 )_yt + = +1 * Sum ( h7 p8 p9 ) * y ( h7 h1 h5 p8 p9 p6 )_y * t ( p8 p9 h7 h3 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh1,sh5,sp6,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer :: h1,h5,p6,h3
  integer :: h7,p8,p9,sh7,sp8,sp9
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh7 = 1,2
  do sp8 = 1,2
  do sp9 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh7,sh1,sh5,sp8,sp9,sp6)
     spin_t2inp_2 = tdcc_spin_t2inp(sp8,sp9,sh7,sh3)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do h3 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
!        i1(h1,h5,p6,h3) = i1(h1,h5,p6,h3) + fact * &
!        g3inp(h7,h1,h5,p8,p9,p6,spin_g3inp_1) * t2inp(p8,p9,h7,h3,spin_t2inp_2)
        i1(h1,h5,p6,h3) = i1(h1,h5,p6,h3) + &
             g3inp(h7,h1,h5,p8,p9,p6,spin_g3inp_1) * t2inp(p8,p9,h7,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_17_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_18(sp1,sh2,sp3,sp4,i0)

!  i0 ( p1 h2 p3 p4 )_yt + = +1/2 * Sum ( h5 h6 p7 ) * y ( h5 h6 h2 p7 p3 p4 )_y * t ( p7 p1 h5 h6 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sp1,sh2,sp3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,h2,p3,p4
  integer :: h5,h6,p7,sh5,sh6,sp7
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp7 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh5,sh6,sh2,sp7,sp3,sp4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp7,sp1,sh5,sh6)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do h2 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i0(p1,h2,p3,p4) = i0(p1,h2,p3,p4) + fact * &
             g3inp(h5,h6,h2,p7,p3,p4,spin_g3inp_1) * t2inp(p7,p1,h5,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_18
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_19(sh1,sh2,sh3,sp4,i0)

!  i0 ( h1 h2 h3 p4 )_yt + = -1/2 * Sum ( h5 p6 p7 ) * y ( h5 h1 h2 p6 p7 p4 )_y * t ( p6 p7 h5 h3 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh1,sh2,sh3,sp4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: h1,h2,h3,p4
  integer :: h5,p6,p7,sh5,sp6,sp7
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh5,sh1,sh2,sp6,sp7,sp4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp6,sp7,sh5,sh3)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p4 = norb1+1,nact
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h1,h2,h3,p4) = i0(h1,h2,h3,p4) + fact * &
             g3inp(h5,h1,h2,p6,p7,p4,spin_g3inp_1) * t2inp(p6,p7,h5,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_19
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_20(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = +1/4 * P( h3 h4 ) * P( p1 p2 ) * Sum ( h5 p6 ) * i1 ( h5 p6 p1 h3 )_yt * t ( p2 p6 h4 h5 )_t 1

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
  call ccdt_den2p_20_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_20_perm(sp1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h4,h3)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_20_perm(sp2,sp1,sh3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_20_perm(sp2,sp1,sh4,sh3,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h4,h3)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_20_perm(sp1,sp2,sh3,sh4,i0)

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
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  allocate(itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1))
  do sh5 = 1,2
  do sp6 = 1,2
     spin_itm_hpph_1 = tdcc_spin_dummy2(sh5,sp6,sp1,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp6,sh4,sh5)
     if(spin_itm_hpph_1 * spin_t2inp_2 == 0) cycle

     itm_hpph(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1) = czero
     call ccdt_den2p_20_1(sh5,sp6,sp1,sh3,itm_hpph)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_hpph(h5,p6,p1,h3) * t2inp(p2,p6,h4,h5,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hpph)
  end subroutine ccdt_den2p_20_perm
  !--------------------------------------------
end subroutine ccdt_den2p_20
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_20_1(sh5,sp6,sp1,sh3,i1)

!      i1 ( h5 p6 p1 h3 )_yt + = +1 * Sum ( h7 h8 p9 p10 ) * y ( h7 h8 h5 p9 p10 p6 )_y * t ( p9 p10 p1 h7 h8 h3 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sh5,sp6,sp1,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1)
  integer :: h5,p6,p1,h3
  integer :: h7,h8,p9,p10,sh7,sh8,sp9,sp10
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp9 = 1,2
  do sp10 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh7,sh8,sh5,sp9,sp10,sp6)
     spin_t3inp_2 = tdcc_spin_t3inp(sp9,sp10,sp1,sh7,sh8,sh3)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do h3 = 1,norb1
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
!        i1(h5,p6,p1,h3) = i1(h5,p6,p1,h3) + fact * &
!        g3inp(h7,h8,h5,p9,p10,p6,spin_g3inp_1) * t3inp(p9,p10,p1,h7,h8,h3,spin_t3inp_2)
        i1(h5,p6,p1,h3) = i1(h5,p6,p1,h3) + &
             g3inp(h7,h8,h5,p9,p10,p6,spin_g3inp_1) * t3inp(p9,p10,p1,h7,h8,h3,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  end do
end subroutine ccdt_den2p_20_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_21(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = +1/12 * Sum ( p5 p6 ) * i1 ( p5 p6 p1 p2 )_yt * t ( p5 p6 h3 h4 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,h3,h4
  integer :: p5,p6,sp5,sp6
  integer :: spin_itm_pppp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_pppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 12.0d+0 * runit

  allocate(itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
  do sp6 = 1,2
     spin_itm_pppp_1 = tdcc_spin_dummy2(sp5,sp6,sp1,sp2)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp6,sh3,sh4)
     if(spin_itm_pppp_1 * spin_t2inp_2 == 0) cycle

     itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_den2p_21_1(sp5,sp6,sp1,sp2,itm_pppp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_pppp(p5,p6,p1,p2) * t2inp(p5,p6,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_pppp)
end subroutine ccdt_den2p_21
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_21_1(sp5,sp6,sp1,sp2,i1)

!      i1 ( p5 p6 p1 p2 )_yt + = +1 * Sum ( h7 h8 h9 p10 ) * y ( h7 h8 h9 p10 p5 p6 )_y * t ( p10 p1 p2 h7 h8 h9 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sp5,sp6,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: p5,p6,p1,p2
  integer :: h7,h8,h9,p10,sh7,sh8,sh9,sp10
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sh9 = 1,2
  do sp10 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh7,sh8,sh9,sp10,sp5,sp6)
     spin_t3inp_2 = tdcc_spin_t3inp(sp10,sp1,sp2,sh7,sh8,sh9)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p10 = norb1+1,nact
!        i1(p5,p6,p1,p2) = i1(p5,p6,p1,p2) + fact * &
!        g3inp(h7,h8,h9,p10,p5,p6,spin_g3inp_1) * t3inp(p10,p1,p2,h7,h8,h9,spin_t3inp_2)
        i1(p5,p6,p1,p2) = i1(p5,p6,p1,p2) + &
             g3inp(h7,h8,h9,p10,p5,p6,spin_g3inp_1) * t3inp(p10,p1,p2,h7,h8,h9,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  end do
end subroutine ccdt_den2p_21_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_22(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = +1/12 * Sum ( h5 h6 ) * i1 ( h5 h6 h3 h4 )_yt * t ( p1 p2 h5 h6 )_t 1

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
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 12.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh5 = 1,2
  do sh6 = 1,2
     spin_itm_hhhh_1 = tdcc_spin_dummy2(sh5,sh6,sh3,sh4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh6)
     if(spin_itm_hhhh_1 * spin_t2inp_2 == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_den2p_22_1(sh5,sh6,sh3,sh4,itm_hhhh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_hhhh(h5,h6,h3,h4) * t2inp(p1,p2,h5,h6,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccdt_den2p_22
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_22_1(sh5,sh6,sh3,sh4,i1)

!      i1 ( h5 h6 h3 h4 )_yt + = +1 * Sum ( h7 p8 p9 p10 ) * y ( h7 h5 h6 p8 p9 p10 )_y * t ( p8 p9 p10 h7 h3 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sh5,sh6,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h5,h6,h3,h4
  integer :: h7,p8,p9,p10,sh7,sp8,sp9,sp10
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh7 = 1,2
  do sp8 = 1,2
  do sp9 = 1,2
  do sp10 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh7,sh5,sh6,sp8,sp9,sp10)
     spin_t3inp_2 = tdcc_spin_t3inp(sp8,sp9,sp10,sh7,sh3,sh4)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
!        i1(h5,h6,h3,h4) = i1(h5,h6,h3,h4) + fact * &
!        g3inp(h7,h5,h6,p8,p9,p10,spin_g3inp_1) * t3inp(p8,p9,p10,h7,h3,h4,spin_t3inp_2)
        i1(h5,h6,h3,h4) = i1(h5,h6,h3,h4) + &
             g3inp(h7,h5,h6,p8,p9,p10,spin_g3inp_1) * t3inp(p8,p9,p10,h7,h3,h4,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  end do
end subroutine ccdt_den2p_22_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_23(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/12 * P( h3 h4 ) * Sum ( h5 ) * i1 ( h5 h3 )_yt * t ( p1 p2 h5 h4 )_t 1

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
  call ccdt_den2p_23_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_23_perm(sp1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h4,h3)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_23_perm(sp1,sp2,sh3,sh4,i0)

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
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 12.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh5 = 1,2
     spin_itm_hh_1 = tdcc_spin_dummy1(sh5,sh3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh5,sh4)
     if(spin_itm_hh_1 * spin_t2inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccdt_den2p_23_1(sh5,sh3,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_hh(h5,h3) * t2inp(p1,p2,h5,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hh)
  end subroutine ccdt_den2p_23_perm
  !--------------------------------------------
end subroutine ccdt_den2p_23
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_23_1(sh5,sh3,i1)

!      i1 ( h5 h3 )_yt + = +1 * Sum ( h6 h7 p8 p9 p10 ) * y ( h6 h7 h5 p8 p9 p10 )_y * t ( p8 p9 p10 h6 h7 h3 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sh5,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h5,h3
  integer :: h6,h7,p8,p9,p10,sh6,sh7,sp8,sp9,sp10
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
  do sp8 = 1,2
  do sp9 = 1,2
  do sp10 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh6,sh7,sh5,sp8,sp9,sp10)
     spin_t3inp_2 = tdcc_spin_t3inp(sp8,sp9,sp10,sh6,sh7,sh3)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h5 = 1,norb1
     do h3 = 1,norb1
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
!        i1(h5,h3) = i1(h5,h3) + fact * &
!        g3inp(h6,h7,h5,p8,p9,p10,spin_g3inp_1) * t3inp(p8,p9,p10,h6,h7,h3,spin_t3inp_2)
        i1(h5,h3) = i1(h5,h3) + &
             g3inp(h6,h7,h5,p8,p9,p10,spin_g3inp_1) * t3inp(p8,p9,p10,h6,h7,h3,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  end do
  end do
end subroutine ccdt_den2p_23_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_24(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/12 * P( p1 p2 ) * Sum ( p5 ) * i1 ( p5 p1 )_yt * t ( p5 p2 h3 h4 )_t 1

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
  call ccdt_den2p_24_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_24_perm(sp2,sp1,sh3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_24_perm(sp1,sp2,sh3,sh4,i0)

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
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 12.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_itm_pp_1 = tdcc_spin_dummy1(sp5,sp1)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp2,sh3,sh4)
     if(spin_itm_pp_1 * spin_t2inp_2 == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_den2p_24_1(sp5,sp1,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_pp(p5,p1) * t2inp(p5,p2,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_pp)
  end subroutine ccdt_den2p_24_perm
  !--------------------------------------------
end subroutine ccdt_den2p_24
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_24_1(sp5,sp1,i1)

!      i1 ( p5 p1 )_yt + = +1 * Sum ( h6 h7 h8 p9 p10 ) * y ( h6 h7 h8 p9 p10 p5 )_y * t ( p9 p10 p1 h6 h7 h8 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sp5,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p5,p1
  integer :: h6,h7,h8,p9,p10,sh6,sh7,sh8,sp9,sp10
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
  do sh8 = 1,2
  do sp9 = 1,2
  do sp10 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh6,sh7,sh8,sp9,sp10,sp5)
     spin_t3inp_2 = tdcc_spin_t3inp(sp9,sp10,sp1,sh6,sh7,sh8)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p5 = norb1+1,nact
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
!        i1(p5,p1) = i1(p5,p1) + fact * &
!        g3inp(h6,h7,h8,p9,p10,p5,spin_g3inp_1) * t3inp(p9,p10,p1,h6,h7,h8,spin_t3inp_2)
        i1(p5,p1) = i1(p5,p1) + &
             g3inp(h6,h7,h8,p9,p10,p5,spin_g3inp_1) * t3inp(p9,p10,p1,h6,h7,h8,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  end do
  end do
end subroutine ccdt_den2p_24_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_25(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/4 * P( h3 h4 ) * Sum ( h5 h6 p7 ) * i1 ( h5 h6 p7 h3 )_yt * t ( p1 p2 p7 h5 h4 h6 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccdt_den2p_25_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_25_perm(sp1,sp2,sh4,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h4,h3)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_25_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer :: p1,p2,h3,h4
  integer :: h5,h6,p7,sh5,sh6,sp7
  integer :: spin_itm_hhph_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sh5 = 1,2
  do sh6 = 1,2
  do sp7 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh5,sh6,sp7,sh3)
     spin_t3inp_2 = tdcc_spin_t3inp(sp1,sp2,sp7,sh5,sh4,sh6)
     if(spin_itm_hhph_1 * spin_t3inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call ccdt_den2p_25_1(sh5,sh6,sp7,sh3,itm_hhph)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_hhph(h5,h6,p7,h3) * t3inp(p1,p2,p7,h5,h4,h6,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  deallocate(itm_hhph)
  end subroutine ccdt_den2p_25_perm
  !--------------------------------------------
end subroutine ccdt_den2p_25
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_25_1(sh5,sh6,sp7,sh3,i1)

!      i1 ( h5 h6 p7 h3 )_yt + = +1 * Sum ( h8 p9 p10 ) * y ( h8 h5 h6 p9 p10 p7 )_y * t ( p9 p10 h8 h3 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh5,sh6,sp7,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer :: h5,h6,p7,h3
  integer :: h8,p9,p10,sh8,sp9,sp10
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh8 = 1,2
  do sp9 = 1,2
  do sp10 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh8,sh5,sh6,sp9,sp10,sp7)
     spin_t2inp_2 = tdcc_spin_t2inp(sp9,sp10,sh8,sh3)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do h3 = 1,norb1
     do h8 = 1,norb1
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
!        i1(h5,h6,p7,h3) = i1(h5,h6,p7,h3) + fact * &
!        g3inp(h8,h5,h6,p9,p10,p7,spin_g3inp_1) * t2inp(p9,p10,h8,h3,spin_t2inp_2)
        i1(h5,h6,p7,h3) = i1(h5,h6,p7,h3) + &
             g3inp(h8,h5,h6,p9,p10,p7,spin_g3inp_1) * t2inp(p9,p10,h8,h3,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_25_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_26(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/4 * P( p1 p2 ) * Sum ( h5 p6 p7 ) * i1 ( h5 p6 p7 p1 )_yt * t ( p6 p2 p7 h3 h4 h5 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: p1,p2,h3,h4
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:nact,1:nact,1:nact,1:nact))

  i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
  call ccdt_den2p_26_perm(sp1,sp2,sh3,sh4,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p1,p2,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp1 * sp2 * sh3 * sh4 == 1)) then
     i0_perm(1:nact,1:nact,1:nact,1:nact) = czero
     call ccdt_den2p_26_perm(sp2,sp1,sh3,sh4,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do h3 = 1,norb1
  do h4 = 1,norb1
     i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact_p * i0_perm(p2,p1,h3,h4)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_den2p_26_perm(sp1,sp2,sh3,sh4,i0)

  implicit none
  integer,intent(in) :: sp1,sp2,sh3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)

  integer :: p1,p2,h3,h4
  integer :: h5,p6,p7,sh5,sp6,sp7
  integer :: spin_itm_hppp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sh5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_itm_hppp_1 = tdcc_spin_dummy2(sh5,sp6,sp7,sp1)
     spin_t3inp_2 = tdcc_spin_t3inp(sp6,sp2,sp7,sh3,sh4,sh5)
     if(spin_itm_hppp_1 * spin_t3inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_den2p_26_1(sh5,sp6,sp7,sp1,itm_hppp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_hppp(h5,p6,p7,p1) * t3inp(p6,p2,p7,h3,h4,h5,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  deallocate(itm_hppp)
  end subroutine ccdt_den2p_26_perm
  !--------------------------------------------
end subroutine ccdt_den2p_26
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_26_1(sh5,sp6,sp7,sp1,i1)

!      i1 ( h5 p6 p7 p1 )_yt + = +1 * Sum ( h8 h9 p10 ) * y ( h8 h9 h5 p10 p6 p7 )_y * t ( p10 p1 h8 h9 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh5,sp6,sp7,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h5,p6,p7,p1
  integer :: h8,h9,p10,sh8,sh9,sp10
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh8 = 1,2
  do sh9 = 1,2
  do sp10 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh8,sh9,sh5,sp10,sp6,sp7)
     spin_t2inp_2 = tdcc_spin_t2inp(sp10,sp1,sh8,sh9)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p10 = norb1+1,nact
!        i1(h5,p6,p7,p1) = i1(h5,p6,p7,p1) + fact * &
!        g3inp(h8,h9,h5,p10,p6,p7,spin_g3inp_1) * t2inp(p10,p1,h8,h9,spin_t2inp_2)
        i1(h5,p6,p7,p1) = i1(h5,p6,p7,p1) + &
             g3inp(h8,h9,h5,p10,p6,p7,spin_g3inp_1) * t2inp(p10,p1,h8,h9,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_den2p_26_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2p_27(sh1,sp2,sp3,sh4,i0)

!  i0 ( h1 p2 p3 h4 )_yt + = +1/4 * Sum ( h5 h6 p7 p8 ) * y ( h5 h6 h1 p7 p8 p3 )_y * t ( p7 p8 p2 h5 h6 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sh1,sp2,sp3,sh4
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact,1:nact,1:nact)
  integer :: h1,p2,p3,h4
  integer :: h5,h6,p7,p8,sh5,sh6,sp7,sp8
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh5,sh6,sh1,sp7,sp8,sp3)
     spin_t3inp_2 = tdcc_spin_t3inp(sp7,sp8,sp2,sh5,sh6,sh4)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(h1,p2,p3,h4) = i0(h1,p2,p3,h4) + fact * &
             g3inp(h5,h6,h1,p7,p8,p3,spin_g3inp_1) * t3inp(p7,p8,p2,h5,h6,h4,spin_t3inp_2)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  end do
end subroutine ccdt_den2p_27
!##########################################################
!##########################################################
!##########################################################
