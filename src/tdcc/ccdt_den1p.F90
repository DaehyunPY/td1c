!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den1p_1(sp1,sp2,i0)

!  i0 ( p1 p2 )_yt + = +1/2 * Sum ( h3 h4 p5 ) * y ( h3 h4 p5 p2 )_y * t ( p5 p1 h3 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: p1,p2
  integer :: h3,h4,p5,sh3,sh4,sp5
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh3,sh4,sp5,sp2)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp1,sh3,sh4)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2) = i0(p1,p2) + fact * g2inp(h3,h4,p5,p2,spin_g2inp_1) * t2inp(p5,p1,h3,h4,spin_t2inp_2)
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
end subroutine ccdt_den1p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den1p_2(sh1,sh2,i0)

!  i0 ( h1 h2 )_yt + = -1/2 * Sum ( h3 p4 p5 ) * y ( h3 h1 p4 p5 )_y * t ( p4 p5 h3 h2 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: h1,h2
  integer :: h3,p4,p5,sh3,sp4,sp5
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh3 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh3,sh1,sp4,sp5)
     spin_t2inp_2 = tdcc_spin_t2inp(sp4,sp5,sh3,sh2)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i0(h1,h2) = i0(h1,h2) + fact * g2inp(h3,h1,p4,p5,spin_g2inp_1) * t2inp(p4,p5,h3,h2,spin_t2inp_2)
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
end subroutine ccdt_den1p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den1p_3(sp1,sh2,i0)

!  i0 ( p1 h2 )_yt + = +1/4 * Sum ( h3 h4 p5 p6 ) * y ( h3 h4 p5 p6 )_y * t ( p5 p6 p1 h3 h4 h2 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sp1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: p1,h2
  integer :: h3,h4,p5,p6,sh3,sh4,sp5,sp6
  integer :: spin_g2inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     spin_t3inp_2 = tdcc_spin_t3inp(sp5,sp6,sp1,sh3,sh4,sh2)
     if(spin_g2inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p1 = norb1+1,nact
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(p1,h2) = i0(p1,h2) + fact * g2inp(h3,h4,p5,p6,spin_g2inp_1) * t3inp(p5,p6,p1,h3,h4,h2,spin_t3inp_2)
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
end subroutine ccdt_den1p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den1p_4(sp1,sh2,i0)

!  i0 ( p1 h2 )_ytt + = -1/4 * Sum ( h3 p4 p5 ) * i1 ( h3 p4 p5 p1 )_yt * t ( p4 p5 h3 h2 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: p1,h2
  integer :: h3,p4,p5,sh3,sp4,sp5
  integer :: spin_itm_hppp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sh3 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
     spin_itm_hppp_1 = tdcc_spin_dummy2(sh3,sp4,sp5,sp1)
     spin_t2inp_2 = tdcc_spin_t2inp(sp4,sp5,sh3,sh2)
     if(spin_itm_hppp_1 * spin_t2inp_2 == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_den1p_4_1(sh3,sp4,sp5,sp1,itm_hppp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p1 = norb1+1,nact
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i0(p1,h2) = i0(p1,h2) + fact * itm_hppp(h3,p4,p5,p1) * t2inp(p4,p5,h3,h2,spin_t2inp_2)
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
end subroutine ccdt_den1p_4
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den1p_4_1(sh3,sp4,sp5,sp1,i1)

!      i1 ( h3 p4 p5 p1 )_yt + = +1 * Sum ( h6 h7 p8 ) * y ( h3 h6 h7 p4 p5 p8 )_y * t ( p1 p8 h6 h7 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh3,sp4,sp5,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h3,p4,p5,p1
  integer :: h6,h7,p8,sh6,sh7,sp8
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
  do sp8 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh3,sh6,sh7,sp4,sp5,sp8)
     spin_t2inp_2 = tdcc_spin_t2inp(sp1,sp8,sh6,sh7)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p8 = norb1+1,nact
!        i1(h3,p4,p5,p1) = i1(h3,p4,p5,p1) + fact * g3inp(h3,h6,h7,p4,p5,p8,spin_g3inp_1) * t2inp(p1,p8,h6,h7,spin_t2inp_2)
        i1(h3,p4,p5,p1) = i1(h3,p4,p5,p1) + g3inp(h3,h6,h7,p4,p5,p8,spin_g3inp_1) * t2inp(p1,p8,h6,h7,spin_t2inp_2)
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
end subroutine ccdt_den1p_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den1p_5(sp1,sp2,i0)

!  i0 ( p1 p2 )_yt + = +1/12 * Sum ( h3 h4 h5 p6 p7 ) * y ( h3 h4 h5 p6 p7 p2 )_y * t ( p6 p7 p1 h3 h4 h5 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: p1,p2
  integer :: h3,h4,h5,p6,p7,sh3,sh4,sh5,sp6,sp7
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 12.0d+0 * runit

  do sh3 = 1,2
  do sh4 = 1,2
  do sh5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh3,sh4,sh5,sp6,sp7,sp2)
     spin_t3inp_2 = tdcc_spin_t3inp(sp6,sp7,sp1,sh3,sh4,sh5)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(p1,p2) = i0(p1,p2) + fact * g3inp(h3,h4,h5,p6,p7,p2,spin_g3inp_1) * t3inp(p6,p7,p1,h3,h4,h5,spin_t3inp_2)
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
end subroutine ccdt_den1p_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_den1p_6(sh1,sh2,i0)

!  i0 ( h1 h2 )_yt + = -1/12 * Sum ( h3 h4 p5 p6 p7 ) * y ( h3 h4 h1 p5 p6 p7 )_y * t ( p5 p6 p7 h3 h4 h2 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t3inp

  implicit none
  integer,intent(in) :: sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: h1,h2
  integer :: h3,h4,p5,p6,p7,sh3,sh4,sp5,sp6,sp7
  integer :: spin_g3inp_1
  integer :: spin_t3inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 12.0d+0 * runit

  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh3,sh4,sh1,sp5,sp6,sp7)
     spin_t3inp_2 = tdcc_spin_t3inp(sp5,sp6,sp7,sh3,sh4,sh2)
     if(spin_g3inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h1,h2) = i0(h1,h2) + fact * g3inp(h3,h4,h1,p5,p6,p7,spin_g3inp_1) * t3inp(p5,p6,p7,h3,h4,h2,spin_t3inp_2)
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
end subroutine ccdt_den1p_6
!##########################################################
!##########################################################
!##########################################################
