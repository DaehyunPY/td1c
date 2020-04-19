!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_1(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_v + = 1 * v ( h3 h4 p1 p2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : int2x,norb1

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh4,sp1,sp2)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * int2x(h3,h4,p1,p2,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccd_l2p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_2(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = 1 * P( 4 ) * y ( h3 p1 )_y * f ( h4 p2 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : fock

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call bccd_l2p_2_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_2_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_2_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_2_perm(sh4,sh3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p2,p1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine bccd_l2p_2_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: sdum
  integer(c_int) :: spin_g1inp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sdum = 1,1
     spin_g1inp = tdcc_spin_g1inp(sh3,sp1)
     spin_fock = tdcc_spin_fock(sh4,sp2)
     if(spin_g1inp * spin_fock == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g1inp(h3,p1,spin_g1inp) * fock(h4,p2,spin_fock)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end subroutine bccd_l2p_2_perm
  !--------------------------------------------
end subroutine bccd_l2p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_3(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1 * P( 2 ) * Sum ( h7 ) * y ( h7 p1 )_y * v ( h3 h4 h7 p2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call bccd_l2p_3_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_3_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine bccd_l2p_3_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h7,sh7
  integer(c_int) :: spin_g1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh7,sp1)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sh7,sp2)
     if(spin_g1inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h7 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g1inp(h7,p1,spin_g1inp) * int2x(h3,h4,h7,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end subroutine bccd_l2p_3_perm
  !--------------------------------------------
end subroutine bccd_l2p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_4(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1 * P( 2 ) * Sum ( p5 ) * y ( h3 p5 )_y * v ( h4 p5 p1 p2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call bccd_l2p_4_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_4_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine bccd_l2p_4_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_g1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh3,sp5)
     spin_int2x = tdcc_spin_int2x(sh4,sp5,sp1,sp2)
     if(spin_g1inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end subroutine bccd_l2p_4_perm
  !--------------------------------------------
end subroutine bccd_l2p_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_5(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = -1 * P( 2 ) * Sum ( h9 ) * y ( h3 h9 p1 p2 )_y * i1 ( h4 h9 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call bccd_l2p_5_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_5_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine bccd_l2p_5_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h9,sh9
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_hh
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh9 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh9,sp1,sp2)
     spin_itm_hh = tdcc_spin_fock(sh4,sh9)
     if(spin_g2inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call bccd_l2p_5_1(sh4,sh9,itm_hh)
     call bccd_l2p_5_2(sh4,sh9,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hh)
  end subroutine bccd_l2p_5_perm
  !--------------------------------------------
end subroutine bccd_l2p_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_5_1(sh3,sh9,i1)

!     i1 ( h3 h9 )_f + = 1 * f ( h3 h9 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh3,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h3,h9
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh3,sh9)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do h3 = 1,norb1
  do h9 = 1,norb1
     i1(h3,h9) = i1(h3,h9) + fact * fock(h3,h9,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccd_l2p_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_5_2(sh3,sh9,i1)

!     i1 ( h3 h9 )_vt + = -1/2 * Sum ( h8 p5 p6 ) * t ( p5 p6 h8 h9 )_t * v ( h3 h8 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h3,h9
  integer(c_int) :: h8,p5,p6,sh8,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh3,sh8,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine bccd_l2p_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_6(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = 1 * P( 2 ) * Sum ( p10 ) * y ( h3 h4 p1 p10 )_y * i1 ( p10 p2 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call bccd_l2p_6_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_6_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine bccd_l2p_6_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: p10,sp10
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_pp
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp10 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp1,sp10)
     spin_itm_pp = tdcc_spin_fock(sp10,sp2)
     if(spin_g2inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_6_1(sp10,sp2,itm_pp)
     call bccd_l2p_6_2(sp10,sp2,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_pp)
  end subroutine bccd_l2p_6_perm
  !--------------------------------------------
end subroutine bccd_l2p_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_6_1(sp10,sp1,i1)

!     i1 ( p10 p1 )_f + = 1 * f ( p10 p1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sp10,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p10,p1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp10,sp1)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do p10 = norb1+1,nact
  do p1 = norb1+1,nact
     i1(p10,p1) = i1(p10,p1) + fact * fock(p10,p1,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccd_l2p_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_6_2(sp10,sp1,i1)

!     i1 ( p10 p1 )_vt + = 1/2 * Sum ( h7 h8 p6 ) * t ( p6 p10 h7 h8 )_t * v ( h7 h8 p1 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp10,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p10,p1
  integer(c_int) :: h7,h8,p6,sh7,sh8,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp10,sh7,sh8)
     spin_int2x = tdcc_spin_int2x(sh7,sh8,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine bccd_l2p_6_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_7(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = 1/2 * Sum ( h9 h10 ) * y ( h9 h10 p1 p2 )_y * i1 ( h3 h4 h9 h10 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h9,h10,sh9,sh10
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_hhhh
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh10 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh9,sh10,sp1,sp2)
     spin_itm_hhhh = tdcc_spin_int2x(sh3,sh4,sh9,sh10)
     if(spin_g2inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call bccd_l2p_7_1(sh3,sh4,sh9,sh10,itm_hhhh)
     call bccd_l2p_7_2(sh3,sh4,sh9,sh10,itm_hhhh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hhhh)
end subroutine bccd_l2p_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_7_1(sh3,sh4,sh9,sh10,i1)

!     i1 ( h3 h4 h9 h10 )_v + = 1 * v ( h3 h4 h9 h10 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sh9,sh10
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h3,h4,h9,h10
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh4,sh9,sh10)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do h9 = 1,norb1
  do h10 = 1,norb1
     i1(h3,h4,h9,h10) = i1(h3,h4,h9,h10) + fact * int2x(h3,h4,h9,h10,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccd_l2p_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_7_2(sh3,sh4,sh9,sh10,i1)

!     i1 ( h3 h4 h9 h10 )_vt + = 1/2 * Sum ( p5 p6 ) * t ( p5 p6 h9 h10 )_t * v ( h3 h4 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sh9,sh10
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h3,h4,h9,h10
  integer(c_int) :: p5,p6,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh9,sh10)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccd_l2p_7_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_8(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1 * P( 4 ) * Sum ( h9 p7 ) * y ( h3 h9 p1 p7 )_y * i1 ( h4 p7 h9 p2 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call bccd_l2p_8_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_8_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_8_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_8_perm(sh4,sh3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p2,p1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine bccd_l2p_8_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h9,p7,sh9,sp7
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_hphp
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sp7 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh9,sp1,sp7)
     spin_itm_hphp = tdcc_spin_int2x(sh4,sp7,sh9,sp2)
     if(spin_g2inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call bccd_l2p_8_1(sh4,sp7,sh9,sp2,itm_hphp)
     call bccd_l2p_8_2(sh4,sp7,sh9,sp2,itm_hphp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hphp)
  end subroutine bccd_l2p_8_perm
  !--------------------------------------------
end subroutine bccd_l2p_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_8_1(sh3,sp7,sh9,sp1,i1)

!     i1 ( h3 p7 h9 p1 )_v + = 1 * v ( h3 p7 h9 p1 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sp7,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,p7,h9,p1
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sp7,sh9,sp1)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do p7 = norb1+1,nact
  do h9 = 1,norb1
  do p1 = norb1+1,nact
     i1(h3,p7,h9,p1) = i1(h3,p7,h9,p1) + fact * int2x(h3,p7,h9,p1,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccd_l2p_8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_8_2(sh3,sp7,sh9,sp1,i1)

!     i1 ( h3 p7 h9 p1 )_vt + = -1 * Sum ( h8 p6 ) * t ( p6 p7 h8 h9 )_t * v ( h3 h8 p1 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sp7,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,p7,h9,p1
  integer(c_int) :: h8,p6,sh8,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp7,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh3,sh8,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccd_l2p_8_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_9(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = 1/2 * Sum ( p5 p6 ) * y ( h3 h4 p5 p6 )_y * v ( p5 p6 p1 p2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: p5,p6,sp5,sp6
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     spin_int2x = tdcc_spin_int2x(sp5,sp6,sp1,sp2)
     if(spin_g2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccd_l2p_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_10(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1 * P( 2 ) * Sum ( h9 ) * i1 ( h3 h9 )_yt * v ( h4 h9 p1 p2 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call bccd_l2p_10_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_10_perm(sh4,sh3,sp1,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h4,h3,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine bccd_l2p_10_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h9,sh9
  integer(c_int) :: spin_itm_hh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh9 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh3,sh9)
     spin_int2x = tdcc_spin_int2x(sh4,sh9,sp1,sp2)
     if(spin_itm_hh * spin_int2x == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call bccd_l2p_10_1(sh3,sh9,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hh)
  end subroutine bccd_l2p_10_perm
  !--------------------------------------------
end subroutine bccd_l2p_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_10_1(sh3,sh9,i1)

!     i1 ( h3 h9 )_yt + = -1/2 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h9 )_t * y ( h3 h7 p5 p6 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h3,h9
  integer(c_int) :: h7,p5,p6,sh7,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh9)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh7,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine bccd_l2p_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_11(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/2 * P( 2 ) * Sum ( p6 ) * i1 ( p6 p1 )_yt * v ( h3 h4 p2 p6 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call bccd_l2p_11_perm(sh3,sh4,sp1,sp2,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p1,p2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh3 * sh4 * sp1 * sp2 == 1)) then
     i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_11_perm(sh3,sh4,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
     i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact_p * i0_perm(h3,h4,p2,p1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine bccd_l2p_11_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: p6,sp6
  integer(c_int) :: spin_itm_pp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp6 = 1,2
     spin_itm_pp = tdcc_spin_dummy1(sp6,sp1)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp2,sp6)
     if(spin_itm_pp * spin_int2x == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l2p_11_1(sp6,sp1,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_pp)
  end subroutine bccd_l2p_11_perm
  !--------------------------------------------
end subroutine bccd_l2p_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_11_1(sp6,sp1,i1)

!     i1 ( p6 p1 )_yt + = -1 * Sum ( h7 h8 p5 ) * t ( p5 p6 h7 h8 )_t * y ( h7 h8 p1 p5 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sp6,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p6,p1
  integer(c_int) :: h7,h8,p5,sh7,sh8,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh7,sh8,sp1,sp5)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine bccd_l2p_11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_12(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1/4 * Sum ( h9 h8 ) * i1 ( h3 h4 h8 h9 )_yt * v ( h8 h9 p1 p2 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h9,h8,sh9,sh8
  integer(c_int) :: spin_itm_hhhh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh8 = 1,2
     spin_itm_hhhh = tdcc_spin_dummy2(sh3,sh4,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh8,sh9,sp1,sp2)
     if(spin_itm_hhhh * spin_int2x == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call bccd_l2p_12_1(sh3,sh4,sh8,sh9,itm_hhhh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hhhh)
end subroutine bccd_l2p_12
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l2p_12_1(sh3,sh4,sh8,sh9,i1)

!     i1 ( h3 h4 h8 h9 )_yt + = -1 * Sum ( p5 p6 ) * t ( p5 p6 h8 h9 )_t * y ( h3 h4 p5 p6 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sh8,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h3,h4,h8,h9
  integer(c_int) :: p5,p6,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh8,sh9)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccd_l2p_12_1
!##########################################################
!##########################################################
!##########################################################
