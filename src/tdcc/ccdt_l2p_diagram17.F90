!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_diagram17_1(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_vty + = 1/12 * Sum ( p5 h8 h9 h10 ) * y ( h8 h9 h10 p1 p2 p5 )_y * i1 ( h3 h4 p5 h8 h9 h10 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : g3inp,norb1

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: p5,h8,h9,h10,sp5,sh8,sh9,sh10
  integer :: spin_g3inp
  integer :: spin_itm_hhphhh
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhphhh(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 12.0d+0 * runit
!
! Demanding: ORDER-8
!
  allocate(itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1))
  do sp5 = 1,2
  do sh8 = 1,2
  do sh9 = 1,2
  do sh10 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh8,sh9,sh10,sp1,sp2,sp5)
     spin_itm_hhphhh = tdcc_spin_dummy3(sh3,sh4,sp5,sh8,sh9,sh10)
     if(spin_g3inp * spin_itm_hhphhh == 0) cycle

     itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_l2p_diagram17_1_1(sh3,sh4,sp5,sh8,sh9,sh10,itm_hhphhh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
  end do
  deallocate(itm_hhphhh)
end subroutine ccdt_l2p_diagram17_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_diagram17_1_1(sh3,sh4,sp5,sh8,sh9,sh10,i1)

!     i1 ( h3 h4 p5 h8 h9 h10 )_vt + = 1 * Sum ( p6 p7 ) * t ( p5 p6 p7 h8 h9 h10 )_t * v ( h3 h4 p6 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : t3inp,norb1
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
!
! Demanding: ORDER-8
!
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp7,sh8,sh9,sh10)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp6,sp7)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(6)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine ccdt_l2p_diagram17_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_diagram17_2(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_vty + = 1/12 * Sum ( p6 p7 ) * i1 ( p6 p7 p1 p2 )_yt * v ( h3 h4 p6 p7 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : int2x,norb1

  implicit none
  integer,intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,h4,p1,p2
  integer :: p6,p7,sp6,sp7
  integer :: spin_itm_pppp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 12.0d+0 * runit

  allocate(itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp6 = 1,2
  do sp7 = 1,2
     spin_itm_pppp = tdcc_spin_dummy2(sp6,sp7,sp1,sp2)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp6,sp7)
     if(spin_itm_pppp * spin_int2x == 0) cycle

     itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l2p_diagram17_2_1(sp6,sp7,sp1,sp2,itm_pppp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_pppp(p6,p7,p1,p2) * int2x(h3,h4,p6,p7,spin_int2x)
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
end subroutine ccdt_l2p_diagram17_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_diagram17_2_1(sp6,sp7,sp1,sp2,i1)

!     i1 ( p6 p7 p1 p2 )_yt + = 1 * Sum ( p5 h8 h9 h10 ) * y ( h8 h9 h10 p1 p2 p5 )_y * t ( p5 p6 p7 h8 h9 h10 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : g3inp,norb1
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sp6,sp7,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: p6,p7,p1,p2
  integer :: p5,h8,h9,h10,sp5,sh8,sh9,sh10
  integer :: spin_g3inp
  integer :: spin_t3inp
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit
!
! Demanding: ORDER-8
!
  do sp5 = 1,2
  do sh8 = 1,2
  do sh9 = 1,2
  do sh10 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh8,sh9,sh10,sp1,sp2,sp5)
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp7,sh8,sh9,sh10)
     if(spin_g3inp * spin_t3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p5 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
        i1(p6,p7,p1,p2) = i1(p6,p7,p1,p2) + fact * &
             g3inp(h8,h9,h10,p1,p2,p5,spin_g3inp) * t3inp(p5,p6,p7,h8,h9,h10,spin_t3inp)
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
end subroutine ccdt_l2p_diagram17_2_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccdt_l2p_diagram17_main()

  use mod_ormas,only : nact
  use mod_cc,only : g2out,norb1

  implicit none

  call ccdt_l2p_diagram17_1(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
  call ccdt_l2p_diagram17_2(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))

  call ccdt_l2p_diagram17_1(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  call ccdt_l2p_diagram17_2(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))

end subroutine ccdt_l2p_diagram17_main
!**********************************************************
