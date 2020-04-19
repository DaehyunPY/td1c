!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2xp_1(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/2 * Sum ( p5 ) * i1 ( p5 p2 )_yt * t ( p1 p5 h3 h4 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

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
     spin_itm_pp_1 = tdcc_spin_dummy1(sp5,sp2)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp5,sh3,sh4)
     if(spin_itm_pp_1 * spin_t2inp_2 == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_den2xp_1_1(sp5,sp2,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_pp(p5,p2) * t2inp(p1,p5,h3,h4,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_pp)
end subroutine ccdt_den2xp_1
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2xp_1_1(sp5,sp2,i1)

!      i1 ( p5 p2 )_yt + = +1 * Sum ( h6 h7 p8 ) * y ( h6 h7 p5 p8 )_y * t ( p2 p8 h6 h7 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sp5,sp2
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p5,p2
  integer :: h6,h7,p8,sh6,sh7,sp8
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
  do sp8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh6,sh7,sp5,sp8)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp8,sh6,sh7)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p5 = norb1+1,nact
     do p2 = norb1+1,nact
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
!        i1(p5,p2) = i1(p5,p2) + fact * &
!        g2inp(h6,h7,p5,p8,spin_g2inp_1) * t2inp(p2,p8,h6,h7,spin_t2inp_2)
        i1(p5,p2) = i1(p5,p2) + &
             g2inp(h6,h7,p5,p8,spin_g2inp_1) * t2inp(p2,p8,h6,h7,spin_t2inp_2)
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
end subroutine ccdt_den2xp_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2xp_2(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/2 * Sum ( h5 ) * i1 ( h5 h4 )_yt * t ( p1 p2 h3 h5 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

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
     spin_itm_hh_1 = tdcc_spin_dummy1(sh5,sh4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp2,sh3,sh5)
     if(spin_itm_hh_1 * spin_t2inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccdt_den2xp_2_1(sh5,sh4,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
        i0(p1,p2,h3,h4) = i0(p1,p2,h3,h4) + fact * &
             itm_hh(h5,h4) * t2inp(p1,p2,h3,h5,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hh)
end subroutine ccdt_den2xp_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2xp_2_1(sh5,sh4,i1)

!      i1 ( h5 h4 )_yt + = +1 * Sum ( h6 p7 p8 ) * y ( h5 h6 p7 p8 )_y * t ( p7 p8 h4 h6 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh5,sh4
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h5,h4
  integer :: h6,p7,p8,sh6,sp7,sp8
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh5,sh6,sp7,sp8)
     spin_t2inp_2 = tdcc_spin_t2inp(sp7,sp8,sh4,sh6)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h5 = 1,norb1
     do h4 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
!        i1(h5,h4) = i1(h5,h4) + fact * &
!        g2inp(h5,h6,p7,p8,spin_g2inp_1) * t2inp(p7,p8,h4,h6,spin_t2inp_2)
        i1(h5,h4) = i1(h5,h4) + &
             g2inp(h5,h6,p7,p8,spin_g2inp_1) * t2inp(p7,p8,h4,h6,spin_t2inp_2)
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
end subroutine ccdt_den2xp_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2xp_3(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/4 * Sum ( h5 p6 p7 ) * i1 ( h5 p6 p7 p1 )_yt * t ( p2 p6 p7 h4 h5 h3 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

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
     spin_t3inp_2 = tdcc_spin_t3inp(sp2,sp6,sp7,sh4,sh5,sh3)
     if(spin_itm_hppp_1 * spin_t3inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_den2xp_3_1(sh5,sp6,sp7,sp1,itm_hppp)

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
             itm_hppp(h5,p6,p7,p1) * t3inp(p2,p6,p7,h4,h5,h3,spin_t3inp_2)
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
end subroutine ccdt_den2xp_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2xp_3_1(sh5,sp6,sp7,sp1,i1)

!      i1 ( h5 p6 p7 p1 )_yt + = +1 * Sum ( h8 h9 p10 ) * y ( h5 h8 h9 p6 p7 p10 )_y * t ( p1 p10 h8 h9 )_t 0

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
     spin_g3inp_1 = tdcc_spin_g3inp(sh5,sh8,sh9,sp6,sp7,sp10)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp10,sh8,sh9)
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
!        g3inp(h5,h8,h9,p6,p7,p10,spin_g3inp_1) * t2inp(p1,p10,h8,h9,spin_t2inp_2)
        i1(h5,p6,p7,p1) = i1(h5,p6,p7,p1) + &
             g3inp(h5,h8,h9,p6,p7,p10,spin_g3inp_1) * t2inp(p1,p10,h8,h9,spin_t2inp_2)
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
end subroutine ccdt_den2xp_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2xp_4(sp1,sp2,sh3,sh4,i0)

!  i0 ( p1 p2 h3 h4 )_ytt + = -1/4 * Sum ( h5 h6 p7 ) * i1 ( h5 h6 p7 h3 )_yt * t ( p2 p7 p1 h4 h5 h6 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp

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
     spin_t3inp_2 = tdcc_spin_t3inp(sp2,sp7,sp1,sh4,sh5,sh6)
     if(spin_itm_hhph_1 * spin_t3inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call ccdt_den2xp_4_1(sh5,sh6,sp7,sh3,itm_hhph)

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
             itm_hhph(h5,h6,p7,h3) * t3inp(p2,p7,p1,h4,h5,h6,spin_t3inp_2)
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
end subroutine ccdt_den2xp_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den2xp_4_1(sh5,sh6,sp7,sh3,i1)

!      i1 ( h5 h6 p7 h3 )_yt + = +1 * Sum ( h8 p9 p10 ) * y ( h5 h6 h8 p7 p9 p10 )_y * t ( p9 p10 h3 h8 )_t 0

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
     spin_g3inp_1 = tdcc_spin_g3inp(sh5,sh6,sh8,sp7,sp9,sp10)
     spin_t2inp_2 = tdcc_spin_t2inp(sp9,sp10,sh3,sh8)
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
!        g3inp(h5,h6,h8,p7,p9,p10,spin_g3inp_1) * t2inp(p9,p10,h3,h8,spin_t2inp_2)
        i1(h5,h6,p7,h3) = i1(h5,h6,p7,h3) + &
             g3inp(h5,h6,h8,p7,p9,p10,spin_g3inp_1) * t2inp(p9,p10,h3,h8,spin_t2inp_2)
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
end subroutine ccdt_den2xp_4_1
!##########################################################
!##########################################################
!##########################################################
