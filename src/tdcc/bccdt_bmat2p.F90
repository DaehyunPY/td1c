!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_bmat2p_1(sh1,sp2,i0)

!  i0 ( h1 p2 )_ydt + = +1 * Sum ( h3 p4 ) * y ( h3 p4 )_y * dt ( p2 p4 h1 h3 )_dt 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : dt2inp

  implicit none
  integer,intent(in) :: sh1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: h1,p2
  integer :: h3,p4,sh3,sp4
  integer :: spin_g1inp_1
  integer :: spin_dt2inp_2
  integer,external :: tdcc_spin_g1inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh3 = 1,2
  do sp4 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh3,sp4)
     spin_dt2inp_2 = tdcc_spin_t2inp(sp2,sp4,sh1,sh3)
     if(spin_g1inp_1 * spin_dt2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do p4 = norb1+1,nact
        i0(h1,p2) = i0(h1,p2) + fact * &
             g1inp(h3,p4,spin_g1inp_1) * dt2inp(p2,p4,h1,h3,spin_dt2inp_2)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccdt_bmat2p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_bmat2p_2(sh1,sp2,i0)

!  i0 ( h1 p2 )_ydt + = +1/4 * Sum ( h3 h4 p5 p6 ) * y ( h3 h4 p5 p6 )_y * dt ( p2 p5 p6 h1 h3 h4 )_dt 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : dt3inp

  implicit none
  integer,intent(in) :: sh1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: h1,p2
  integer :: h3,h4,p5,p6,sh3,sh4,sp5,sp6
  integer :: spin_g2inp_1
  integer :: spin_dt3inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     spin_dt3inp_2 = tdcc_spin_t3inp(sp2,sp5,sp6,sh1,sh3,sh4)
     if(spin_g2inp_1 * spin_dt3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(h1,p2) = i0(h1,p2) + fact * &
             g2inp(h3,h4,p5,p6,spin_g2inp_1) * dt3inp(p2,p5,p6,h1,h3,h4,spin_dt3inp_2)
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
end subroutine bccdt_bmat2p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_bmat2p_3(sh1,sp2,i0)

!  i0 ( h1 p2 )_ydtt + = -1/4 * Sum ( h3 h4 p5 ) * i1 ( h3 h4 p5 h1 )_ydt * t ( p2 p5 h3 h4 )_t 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: h1,p2
  integer :: h3,h4,p5,sh3,sh4,sp5
  integer :: spin_itm_hhph_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh3,sh4,sp5,sh1)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp5,sh3,sh4)
     if(spin_itm_hhph_1 * spin_t2inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call bccdt_bmat2p_3_1(sh3,sh4,sp5,sh1,itm_hhph)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(h1,p2) = i0(h1,p2) + fact * &
             itm_hhph(h3,h4,p5,h1) * t2inp(p2,p5,h3,h4,spin_t2inp_2)
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
end subroutine bccdt_bmat2p_3
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_bmat2p_3_1(sh3,sh4,sp5,sh1,i1)

!      i1 ( h3 h4 p5 h1 )_ydt + = +1 * Sum ( h6 p7 p8 ) * y ( h6 h3 h4 p7 p8 p5 )_y * dt ( p8 p7 h1 h6 )_dt 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : dt2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp5,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer :: h3,h4,p5,h1
  integer :: h6,p7,p8,sh6,sp7,sp8
  integer :: spin_g3inp_1
  integer :: spin_dt2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh6,sh3,sh4,sp7,sp8,sp5)
     spin_dt2inp_2 = tdcc_spin_t2inp(sp8,sp7,sh1,sh6)
     if(spin_g3inp_1 * spin_dt2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h3,h4,p5,h1) = i1(h3,h4,p5,h1) + fact * &
             g3inp(h6,h3,h4,p7,p8,p5,spin_g3inp_1) * dt2inp(p8,p7,h1,h6,spin_dt2inp_2)
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
end subroutine bccdt_bmat2p_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_bmat2p_4(sh1,sp2,i0)

!  i0 ( h1 p2 )_ytdt + = -1/4 * Sum ( h3 h4 p5 ) * i1 ( h3 h4 p5 h1 )_yt * dt ( p5 p2 h3 h4 )_dt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, dt2inp

  implicit none
  integer,intent(in) :: sh1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: h1,p2
  integer :: h3,h4,p5,sh3,sh4,sp5
  integer :: spin_itm_hhph_1
  integer :: spin_dt2inp_2
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh3,sh4,sp5,sh1)
     spin_dt2inp_2 = tdcc_spin_t2inp(sp5,sp2,sh3,sh4)
     if(spin_itm_hhph_1 * spin_dt2inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call bccdt_bmat2p_4_1(sh3,sh4,sp5,sh1,itm_hhph)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(h1,p2) = i0(h1,p2) + fact * &
             itm_hhph(h3,h4,p5,h1) * dt2inp(p5,p2,h3,h4,spin_dt2inp_2)
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
end subroutine bccdt_bmat2p_4
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_bmat2p_4_1(sh3,sh4,sp5,sh1,i1)

!      i1 ( h3 h4 p5 h1 )_yt + = +1 * Sum ( h6 p7 p8 ) * y ( h3 h4 h6 p5 p7 p8 )_y * t ( p7 p8 h1 h6 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh3,sh4,sp5,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer :: h3,h4,p5,h1
  integer :: h6,p7,p8,sh6,sp7,sp8
  integer :: spin_g3inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_g3inp_1 = tdcc_spin_g3inp(sh3,sh4,sh6,sp5,sp7,sp8)
     spin_t2inp_2 = tdcc_spin_t2inp(sp7,sp8,sh1,sh6)
     if(spin_g3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h3,h4,p5,h1) = i1(h3,h4,p5,h1) + fact * &
             g3inp(h3,h4,h6,p5,p7,p8,spin_g3inp_1) * t2inp(p7,p8,h1,h6,spin_t2inp_2)
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
end subroutine bccdt_bmat2p_4_1
!##########################################################
!##########################################################
!##########################################################
