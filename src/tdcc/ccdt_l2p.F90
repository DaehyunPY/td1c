!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_1(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_v + = 1 * v ( h3 h4 p1 p2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
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
end subroutine ccdt_l2p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_2(sh3,sh4,sp1,sp2,i0)

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
  call ccdt_l2p_2_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_2_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccdt_l2p_2_perm(sh3,sh4,sp1,sp2,i0)

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
     call ccdt_l2p_2_1(sh4,sh9,itm_hh)
     call ccdt_l2p_2_2(sh4,sh9,itm_hh)

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
  end subroutine ccdt_l2p_2_perm
  !--------------------------------------------
end subroutine ccdt_l2p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_2_1(sh3,sh9,i1)

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
end subroutine ccdt_l2p_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_2_2(sh3,sh9,i1)

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
end subroutine ccdt_l2p_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_3(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yf + = 1 * P( 2 ) * Sum ( p5 ) * y ( h3 h4 p1 p5 )_y * i1 ( p5 p2 )_f 2

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
  call ccdt_l2p_3_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_3_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2p_3_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_pp
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp5 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp1,sp5)
     spin_itm_pp = tdcc_spin_fock(sp5,sp2)
     if(spin_g2inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l2p_3_1(sp5,sp2,itm_pp)
     call ccdt_l2p_3_2(sp5,sp2,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_pp)
  end subroutine ccdt_l2p_3_perm
  !--------------------------------------------
end subroutine ccdt_l2p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_3_1(sp5,sp1,i1)

!     i1 ( p5 p1 )_f + = 1 * f ( p5 p1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sp5,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p5,p1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp5,sp1)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do p5 = norb1+1,nact
  do p1 = norb1+1,nact
     i1(p5,p1) = i1(p5,p1) + fact * fock(p5,p1,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l2p_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_3_2(sp5,sp1,i1)

!     i1 ( p5 p1 )_vt + = -1/2 * Sum ( h7 h8 p6 ) * t ( p5 p6 h7 h8 )_t * v ( h7 h8 p1 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp5,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p5,p1
  integer(c_int) :: h7,h8,p6,sh7,sh8,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_int2x = tdcc_spin_int2x(sh7,sh8,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_l2p_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_4(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1/2 * Sum ( h9 h7 ) * y ( h7 h9 p1 p2 )_y * i1 ( h3 h4 h7 h9 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h9,h7,sh9,sh7
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_hhhh
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh7 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh7,sh9,sp1,sp2)
     spin_itm_hhhh = tdcc_spin_int2x(sh3,sh4,sh7,sh9)
     if(spin_g2inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_l2p_4_1(sh3,sh4,sh7,sh9,itm_hhhh)
     call ccdt_l2p_4_2(sh3,sh4,sh7,sh9,itm_hhhh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hhhh)
end subroutine ccdt_l2p_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_4_1(sh3,sh4,sh7,sh9,i1)

!     i1 ( h3 h4 h7 h9 )_v + = -1 * v ( h3 h4 h7 h9 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sh7,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h3,h4,h7,h9
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sh4,sh7,sh9)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do h4 = 1,norb1
  do h7 = 1,norb1
  do h9 = 1,norb1
     i1(h3,h4,h7,h9) = i1(h3,h4,h7,h9) + fact * int2x(h3,h4,h7,h9,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l2p_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_4_2(sh3,sh4,sh7,sh9,i1)

!     i1 ( h3 h4 h7 h9 )_vt + = -1/2 * Sum ( p5 p6 ) * t ( p5 p6 h7 h9 )_t * v ( h3 h4 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sh7,sh9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h3,h4,h7,h9
  integer(c_int) :: p5,p6,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh9)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine ccdt_l2p_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_5(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1 * P( 4 ) * Sum ( h9 p5 ) * y ( h3 h9 p1 p5 )_y * i1 ( h4 p5 h9 p2 )_v 2

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
  call ccdt_l2p_5_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_5_perm(sh4,sh3,sp1,sp2,i0_perm)
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
     call ccdt_l2p_5_perm(sh3,sh4,sp2,sp1,i0_perm)
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
     call ccdt_l2p_5_perm(sh4,sh3,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2p_5_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h9,p5,sh9,sp5
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_hphp
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sp5 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh3,sh9,sp1,sp5)
     spin_itm_hphp = tdcc_spin_int2x(sh4,sp5,sh9,sp2)
     if(spin_g2inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccdt_l2p_5_1(sh4,sp5,sh9,sp2,itm_hphp)
     call ccdt_l2p_5_2(sh4,sp5,sh9,sp2,itm_hphp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hphp)
  end subroutine ccdt_l2p_5_perm
  !--------------------------------------------
end subroutine ccdt_l2p_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_5_1(sh3,sp5,sh9,sp1,i1)

!     i1 ( h3 p5 h9 p1 )_v + = 1 * v ( h3 p5 h9 p1 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sp5,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,p5,h9,p1
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sp5,sh9,sp1)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do p5 = norb1+1,nact
  do h9 = 1,norb1
  do p1 = norb1+1,nact
     i1(h3,p5,h9,p1) = i1(h3,p5,h9,p1) + fact * int2x(h3,p5,h9,p1,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l2p_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_5_2(sh3,sp5,sh9,sp1,i1)

!     i1 ( h3 p5 h9 p1 )_vt + = 1 * Sum ( h8 p6 ) * t ( p5 p6 h8 h9 )_t * v ( h3 h8 p1 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sp5,sh9,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,p5,h9,p1
  integer(c_int) :: h8,p6,sh8,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh3,sh8,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine ccdt_l2p_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_6(sh3,sh4,sp1,sp2,i0)

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
end subroutine ccdt_l2p_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_7(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1/2 * P( 2 ) * Sum ( h11 h12 p13 ) * y ( h3 h11 h12 p1 p2 p13 )_y * i1 ( h4 p13 h11 h12 )_v 5

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l2p_7_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_7_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccdt_l2p_7_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h11,h12,p13,sh11,sh12,sp13
  integer(c_int) :: spin_g3inp
  integer(c_int) :: spin_itm_hphh
  integer(c_int),external :: tdcc_spin_g3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh11 = 1,2
  do sh12 = 1,2
  do sp13 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh3,sh11,sh12,sp1,sp2,sp13)
     spin_itm_hphh = tdcc_spin_int2x(sh4,sp13,sh11,sh12)
     if(spin_g3inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccdt_l2p_7_1(sh4,sp13,sh11,sh12,itm_hphh)
     call ccdt_l2p_7_2(sh4,sp13,sh11,sh12,itm_hphh)

     call ccdt_l2p_7_3(sh4,sp13,sh11,sh12,itm_hphh)
!test     if (sh4*sp13*sh11*sh12==1) then
!test        call ccdt_l2p_7_3(sh4,sp13,sh11,sh12,itm_hphh)
!test     end if

     call ccdt_l2p_7_4(sh4,sp13,sh11,sh12,itm_hphh)
     call ccdt_l2p_7_5(sh4,sp13,sh11,sh12,itm_hphh)

!
! Demanding: ORDER-7
!

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p13 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * g3inp(h3,h11,h12,p1,p2,p13,spin_g3inp) * itm_hphh(h4,p13,h11,h12)
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
  deallocate(itm_hphh)
  end subroutine ccdt_l2p_7_perm
  !--------------------------------------------
end subroutine ccdt_l2p_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_7_1(sh3,sp13,sh11,sh12,i1)

!     i1 ( h3 p13 h11 h12 )_v + = 1 * v ( h3 p13 h11 h12 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sp13,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h3,p13,h11,h12
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh3,sp13,sh11,sh12)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h3 = 1,norb1
  do p13 = norb1+1,nact
  do h11 = 1,norb1
  do h12 = 1,norb1
     i1(h3,p13,h11,h12) = i1(h3,p13,h11,h12) + fact * int2x(h3,p13,h11,h12,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l2p_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_7_2(sh3,sp13,sh11,sh12,i1)

!     i1 ( h3 p13 h11 h12 )_ft + = 1 * Sum ( p5 ) * t ( p5 p13 h11 h12 )_t * f ( h3 p5 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : fock

  implicit none
  integer(c_int),intent(in) :: sh3,sp13,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h3,p13,h11,h12
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp13,sh11,sh12)
     spin_fock = tdcc_spin_fock(sh3,sp5)
     if(spin_t2inp * spin_fock == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,p13,h11,h12) = i1(h3,p13,h11,h12) + fact * t2inp(p5,p13,h11,h12,spin_t2inp) * fock(h3,p5,spin_fock)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
end subroutine ccdt_l2p_7_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_7_3(sh3,sp13,sh11,sh12,i1)

!     i1 ( h3 p13 h11 h12 )_vt + = -2 * Sum ( h8 p6 ) * t ( p6 p13 h8 h11 )_t * v ( h3 h8 h12 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sp13,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h3,p13,h11,h12
  integer(c_int) :: h8,p6,sh8,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp13,sh8,sh11)
     spin_int2x = tdcc_spin_int2x(sh3,sh8,sh12,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i1(h3,p13,h11,h12) = i1(h3,p13,h11,h12) + fact * &
             t2inp(p6,p13,h8,h11,spin_t2inp) * int2x(h3,h8,h12,p6,spin_int2x)
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
end subroutine ccdt_l2p_7_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_7_4(sh3,sp13,sh11,sh12,i1)

!     i1 ( h3 p13 h11 h12 )_vt + = 1/2 * Sum ( p5 p6 ) * t ( p5 p6 h11 h12 )_t * v ( h3 p13 p5 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sp13,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h3,p13,h11,h12
  integer(c_int) :: p5,p6,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh11,sh12)
     spin_int2x = tdcc_spin_int2x(sh3,sp13,sp5,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,p13,h11,h12) = i1(h3,p13,h11,h12) + fact * &
             t2inp(p5,p6,h11,h12,spin_t2inp) * int2x(h3,p13,p5,p6,spin_int2x)
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
end subroutine ccdt_l2p_7_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_7_5(sh3,sp13,sh11,sh12,i1)

!     i1 ( h3 p13 h11 h12 )_vt + = -1/2 * Sum ( h10 p6 p7 ) * t ( p6 p7 p13 h10 h11 h12 )_t * v ( h3 h10 p6 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sp13,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h3,p13,h11,h12
  integer(c_int) :: h10,p6,p7,sh10,sp6,sp7
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit
!
! Demanding: ORDER-7
!
  do sh10 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp6,sp7,sp13,sh10,sh11,sh12)
     spin_int2x = tdcc_spin_int2x(sh3,sh10,sp6,sp7)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do h10 = 1,norb1
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h3,p13,h11,h12) = i1(h3,p13,h11,h12) + fact * &
             t3inp(p6,p7,p13,h10,h11,h12,spin_t3inp) * int2x(h3,h10,p6,p7,spin_int2x)
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
end subroutine ccdt_l2p_7_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_8(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_yv + = -1/2 * P( 2 ) * Sum ( h12 p5 p11 ) 
!  * y ( h3 h4 h12 p1 p5 p11 )_y * i1 ( p5 p11 h12 p2 )_v 4

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l2p_8_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_8_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2p_8_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h12,p5,p11,sh12,sp5,sp11
  integer(c_int) :: spin_g3inp
  integer(c_int) :: spin_itm_pphp
  integer(c_int),external :: tdcc_spin_g3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh12 = 1,2
  do sp5 = 1,2
  do sp11 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh12,sp1,sp5,sp11)
     spin_itm_pphp = tdcc_spin_int2x(sp5,sp11,sh12,sp2)
     if(spin_g3inp * spin_itm_pphp == 0) cycle

     itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccdt_l2p_8_1(sp5,sp11,sh12,sp2,itm_pphp)
     call ccdt_l2p_8_2(sp5,sp11,sh12,sp2,itm_pphp)
     call ccdt_l2p_8_3(sp5,sp11,sh12,sp2,itm_pphp)
     call ccdt_l2p_8_4(sp5,sp11,sh12,sp2,itm_pphp)
!
! Demanding: ORDER-7
!
     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h12 = 1,norb1
     do p5 = norb1+1,nact
     do p11 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             g3inp(h3,h4,h12,p1,p5,p11,spin_g3inp) * itm_pphp(p5,p11,h12,p2)
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
  deallocate(itm_pphp)
  end subroutine ccdt_l2p_8_perm
  !--------------------------------------------
end subroutine ccdt_l2p_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_8_1(sp5,sp11,sh12,sp1,i1)

!     i1 ( p5 p11 h12 p1 )_v + = 1 * v ( p5 p11 h12 p1 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sp5,sp11,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: p5,p11,h12,p1
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sp5,sp11,sh12,sp1)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p5 = norb1+1,nact
  do p11 = norb1+1,nact
  do h12 = 1,norb1
  do p1 = norb1+1,nact
     i1(p5,p11,h12,p1) = i1(p5,p11,h12,p1) + fact * int2x(p5,p11,h12,p1,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l2p_8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_8_2(sp5,sp11,sh12,sp1,i1)

!     i1 ( p5 p11 h12 p1 )_vt + = 1/2 * Sum ( h7 h8 ) * t ( p5 p11 h7 h8 )_t * v ( h7 h8 h12 p1 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp5,sp11,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: p5,p11,h12,p1
  integer(c_int) :: h7,h8,sh7,sh8
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp11,sh7,sh8)
     spin_int2x = tdcc_spin_int2x(sh7,sh8,sh12,sp1)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p5 = norb1+1,nact
     do p11 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
        i1(p5,p11,h12,p1) = i1(p5,p11,h12,p1) + fact * &
             t2inp(p5,p11,h7,h8,spin_t2inp) * int2x(h7,h8,h12,p1,spin_int2x)
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
end subroutine ccdt_l2p_8_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_8_3(sp5,sp11,sh12,sp1,i1)

!     i1 ( p5 p11 h12 p1 )_vt + = 2 * Sum ( h8 p6 ) * t ( p5 p6 h8 h12 )_t * v ( h8 p11 p1 p6 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp5,sp11,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: p5,p11,h12,p1
  integer(c_int) :: h8,p6,sh8,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 2.0d+0 * runit

  do sh8 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh8,sh12)
     spin_int2x = tdcc_spin_int2x(sh8,sp11,sp1,sp6)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p5 = norb1+1,nact
     do p11 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do p6 = norb1+1,nact
        i1(p5,p11,h12,p1) = i1(p5,p11,h12,p1) + fact * &
             t2inp(p5,p6,h8,h12,spin_t2inp) * int2x(h8,p11,p1,p6,spin_int2x)
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
end subroutine ccdt_l2p_8_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_8_4(sp5,sp11,sh12,sp1,i1)

!     i1 ( p5 p11 h12 p1 )_vt + = 1/2 * Sum ( h9 h10 p7 ) * t ( p5 p7 p11 h9 h10 h12 )_t * v ( h9 h10 p1 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp5,sp11,sh12,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: p5,p11,h12,p1
  integer(c_int) :: h9,h10,p7,sh9,sh10,sp7
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit
!
! Demanding: ORDER-7
!
  do sh9 = 1,2
  do sh10 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp7,sp11,sh9,sh10,sh12)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp1,sp7)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p5 = norb1+1,nact
     do p11 = norb1+1,nact
     do h12 = 1,norb1
     do p1 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p7 = norb1+1,nact
        i1(p5,p11,h12,p1) = i1(p5,p11,h12,p1) + fact * &
             t3inp(p5,p7,p11,h9,h10,h12,spin_t3inp) * int2x(h9,h10,p1,p7,spin_int2x)
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
end subroutine ccdt_l2p_8_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_9(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/2 * P( 2 ) * Sum ( p11 ) * i1 ( p11 p1 )_yt * v ( h3 h4 p2 p11 )_v 2

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
  call ccdt_l2p_9_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_9_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2p_9_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: p11,sp11
  integer(c_int) :: spin_itm_pp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp11 = 1,2
     spin_itm_pp = tdcc_spin_dummy1(sp11,sp1)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp2,sp11)
     if(spin_itm_pp * spin_int2x == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l2p_9_1(sp11,sp1,itm_pp)
     call ccdt_l2p_9_2(sp11,sp1,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p11 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_pp(p11,p1) * int2x(h3,h4,p2,p11,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_pp)
  end subroutine ccdt_l2p_9_perm
  !--------------------------------------------
end subroutine ccdt_l2p_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_9_1(sp11,sp1,i1)

!     i1 ( p11 p1 )_yt + = -1 * Sum ( h7 h8 p5 ) * t ( p5 p11 h7 h8 )_t * y ( h7 h8 p1 p5 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sp11,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p11,p1
  integer(c_int) :: h7,h8,p5,sh7,sh8,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp11,sh7,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh7,sh8,sp1,sp5)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p11 = norb1+1,nact
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(p11,p1) = i1(p11,p1) + fact * &
             t2inp(p5,p11,h7,h8,spin_t2inp) * g2inp(h7,h8,p1,p5,spin_g2inp)
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
end subroutine ccdt_l2p_9_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_9_2(sp11,sp1,i1)

!     i1 ( p11 p1 )_yt + = 1/6 * Sum ( h8 h9 h10 p5 p6 ) * t ( p5 p6 p11 h8 h9 h10 )_t * y ( h8 h9 h10 p1 p5 p6 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sp11,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p11,p1
  integer(c_int) :: h8,h9,h10,p5,p6,sh8,sh9,sh10,sp5,sp6
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 6.0d+0 * runit
!
! Demanding: ORDER-7
!
  do sh8 = 1,2
  do sh9 = 1,2
  do sh10 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp11,sh8,sh9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh8,sh9,sh10,sp1,sp5,sp6)
     if(spin_t3inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p11 = norb1+1,nact
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(p11,p1) = i1(p11,p1) + fact * &
             t3inp(p5,p6,p11,h8,h9,h10,spin_t3inp) * g3inp(h8,h9,h10,p1,p5,p6,spin_g3inp)
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
end subroutine ccdt_l2p_9_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_10(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/2 * P( 2 ) * Sum ( h11 ) * i1 ( h3 h11 )_yt * v ( h4 h11 p1 p2 )_v 2

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
  call ccdt_l2p_10_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_10_perm(sh4,sh3,sp1,sp2,i0_perm)
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
  subroutine ccdt_l2p_10_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h11,sh11
  integer(c_int) :: spin_itm_hh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh11 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh3,sh11)
     spin_int2x = tdcc_spin_int2x(sh4,sh11,sp1,sp2)
     if(spin_itm_hh * spin_int2x == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccdt_l2p_10_1(sh3,sh11,itm_hh)
     call ccdt_l2p_10_2(sh3,sh11,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h11 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hh(h3,h11) * int2x(h4,h11,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hh)
  end subroutine ccdt_l2p_10_perm
  !--------------------------------------------
end subroutine ccdt_l2p_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_10_1(sh3,sh11,i1)

!     i1 ( h3 h11 )_yt + = -1 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h11 )_t * y ( h3 h7 p5 p6 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h3,h11
  integer(c_int) :: h7,p5,p6,sh7,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh11)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh7,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h3 = 1,norb1
     do h11 = 1,norb1
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h11) = i1(h3,h11) + fact * &
             t2inp(p5,p6,h7,h11,spin_t2inp) * g2inp(h3,h7,p5,p6,spin_g2inp)
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
end subroutine ccdt_l2p_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_10_2(sh3,sh11,i1)

!     i1 ( h3 h11 )_yt + = 1/6 * Sum ( h8 h9 p5 p6 p7 ) * t ( p5 p6 p7 h8 h9 h11 )_t * y ( h3 h8 h9 p5 p6 p7 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h3,h11
  integer(c_int) :: h8,h9,p5,p6,p7,sh8,sh9,sp5,sp6,sp7
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 6.0d+0 * runit
!
! Demanding: ORDER-7
!
  do sh8 = 1,2
  do sh9 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp7,sh8,sh9,sh11)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh8,sh9,sp5,sp6,sp7)
     if(spin_t3inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h3 = 1,norb1
     do h11 = 1,norb1
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h3,h11) = i1(h3,h11) + fact * &
             t3inp(p5,p6,p7,h8,h9,h11,spin_t3inp) * g3inp(h3,h8,h9,p5,p6,p7,spin_g3inp)
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
end subroutine ccdt_l2p_10_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_11(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1/4 * Sum ( h11 h9 ) * i1 ( h3 h4 h9 h11 )_yt * v ( h9 h11 p1 p2 )_v 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h11,h9,sh11,sh9
  integer(c_int) :: spin_itm_hhhh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh11 = 1,2
  do sh9 = 1,2
     spin_itm_hhhh = tdcc_spin_dummy2(sh3,sh4,sh9,sh11)
     spin_int2x = tdcc_spin_int2x(sh9,sh11,sp1,sp2)
     if(spin_itm_hhhh * spin_int2x == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_l2p_11_1(sh3,sh4,sh9,sh11,itm_hhhh)
     call ccdt_l2p_11_2(sh3,sh4,sh9,sh11,itm_hhhh)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h11 = 1,norb1
     do h9 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hhhh(h3,h4,h9,h11) * int2x(h9,h11,p1,p2,spin_int2x)
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
end subroutine ccdt_l2p_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_11_1(sh3,sh4,sh9,sh11,i1)

!     i1 ( h3 h4 h9 h11 )_yt + = -1 * Sum ( p5 p6 ) * t ( p5 p6 h9 h11 )_t * y ( h3 h4 p5 p6 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sh9,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h3,h4,h9,h11
  integer(c_int) :: p5,p6,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh9,sh11)
     spin_g2inp = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h9 = 1,norb1
     do h11 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h4,h9,h11) = i1(h3,h4,h9,h11) + fact * &
             t2inp(p5,p6,h9,h11,spin_t2inp) * g2inp(h3,h4,p5,p6,spin_g2inp)
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
end subroutine ccdt_l2p_11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_11_2(sh3,sh4,sh9,sh11,i1)

! i1 ( h3 h4 h9 h11 )_yt + = -1/3 * Sum ( h8 p5 p6 p7 ) 
!  * t ( p5 p6 p7 h8 h9 h11 )_t * y ( h3 h4 h8 p5 p6 p7 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sh9,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h3,h4,h9,h11
  integer(c_int) :: h8,p5,p6,p7,sh8,sp5,sp6,sp7
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 3.0d+0 * runit
!
! Demanding: ORDER-8
!
  do sh8 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp7,sh8,sh9,sh11)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh8,sp5,sp6,sp7)
     if(spin_t3inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h9 = 1,norb1
     do h11 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do p7 = norb1+1,nact
        i1(h3,h4,h9,h11) = i1(h3,h4,h9,h11) + fact * &
             t3inp(p5,p6,p7,h8,h9,h11,spin_t3inp) * g3inp(h3,h4,h8,p5,p6,p7,spin_g3inp)
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
end subroutine ccdt_l2p_11_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_12(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytf + = 1/2 * P( 2 ) * Sum ( h5 ) * i1 ( h3 h4 h5 p1 )_yt * f ( h5 p2 )_f 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l2p_12_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_12_perm(sh3,sh4,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2p_12_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h5,sh5
  integer(c_int) :: spin_itm_hhhp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh5 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh5,sp1)
     spin_fock = tdcc_spin_fock(sh5,sp2)
     if(spin_itm_hhhp * spin_fock == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccdt_l2p_12_1(sh3,sh4,sh5,sp1,itm_hhhp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hhhp)
  end subroutine ccdt_l2p_12_perm
  !--------------------------------------------
end subroutine ccdt_l2p_12
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_12_1(sh3,sh4,sh5,sp1,i1)

!     i1 ( h3 h4 h5 p1 )_yt + = -1 * Sum ( h8 p6 p7 ) * t ( p6 p7 h5 h8 )_t * y ( h3 h4 h8 p1 p6 p7 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sh5,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,h4,h5,p1
  integer(c_int) :: h8,p6,p7,sh8,sp6,sp7
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit
!
! Demanding: ORDER-7
!
  do sh8 = 1,2
  do sp6 = 1,2
  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp6,sp7,sh5,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh8,sp1,sp6,sp7)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
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
     !$omp end do
     !$omp end parallel
  end do
  end do
  end do
end subroutine ccdt_l2p_12_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_13(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1/2 * Sum ( h9 p6 ) * i1 ( h9 p6 p1 p2 )_yt * v ( h3 h4 h9 p6 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h9,p6,sh9,sp6
  integer(c_int) :: spin_itm_hppp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sh9 = 1,2
  do sp6 = 1,2
     spin_itm_hppp = tdcc_spin_dummy2(sh9,sp6,sp1,sp2)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sh9,sp6)
     if(spin_itm_hppp * spin_int2x == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l2p_13_1(sh9,sp6,sp1,sp2,itm_hppp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do p6 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hppp(h9,p6,p1,p2) * int2x(h3,h4,h9,p6,spin_int2x)
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
end subroutine ccdt_l2p_13
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_13_1(sh9,sp6,sp1,sp2,i1)

!     i1 ( h9 p6 p1 p2 )_yt + = 1 * Sum ( h7 h8 p5 ) * t ( p5 p6 h7 h8 )_t * y ( h7 h8 h9 p1 p2 p5 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sh9,sp6,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h9,p6,p1,p2
  integer(c_int) :: h7,h8,p5,sh7,sh8,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit
!
! Demanding: ORDER-7
!
  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh7,sh8,sh9,sp1,sp2,sp5)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h9 = 1,norb1
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(h9,p6,p1,p2) = i1(h9,p6,p1,p2) + fact * &
             t2inp(p5,p6,h7,h8,spin_t2inp) * g3inp(h7,h8,h9,p1,p2,p5,spin_g3inp)
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
end subroutine ccdt_l2p_13_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_14(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1/2 * P( 4 ) * Sum ( h9 h8 ) * i1 ( h3 h9 h8 p1 )_yt * v ( h4 h8 h9 p2 )_v 1

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
  call ccdt_l2p_14_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_14_perm(sh4,sh3,sp1,sp2,i0_perm)
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
     call ccdt_l2p_14_perm(sh3,sh4,sp2,sp1,i0_perm)
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
     call ccdt_l2p_14_perm(sh4,sh3,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2p_14_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h9,h8,sh9,sh8
  integer(c_int) :: spin_itm_hhhp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sh8 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh9,sh8,sp1)
     spin_int2x = tdcc_spin_int2x(sh4,sh8,sh9,sp2)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccdt_l2p_14_1(sh3,sh9,sh8,sp1,itm_hhhp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do h8 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hhhp(h3,h9,h8,p1) * int2x(h4,h8,h9,p2,spin_int2x)
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
  deallocate(itm_hhhp)
  end subroutine ccdt_l2p_14_perm
  !--------------------------------------------
end subroutine ccdt_l2p_14
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_14_1(sh3,sh9,sh8,sp1,i1)

!     i1 ( h3 h9 h8 p1 )_yt + = -1 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h8 )_t * y ( h3 h7 h9 p1 p5 p6 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh9,sh8,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,h9,h8,p1
  integer(c_int) :: h7,p5,p6,sh7,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit
!
! Demanding: ORDER-7
!
  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh7,sh9,sp1,sp5,sp6)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h9 = 1,norb1
     do h8 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h9,h8,p1) = i1(h3,h9,h8,p1) + fact * &
             t2inp(p5,p6,h7,h8,spin_t2inp) * g3inp(h3,h7,h9,p1,p5,p6,spin_g3inp)
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
end subroutine ccdt_l2p_14_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_15(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1/2 * P( 4 ) * Sum ( p6 p9 ) * i1 ( h3 p6 p1 p9 )_yt * v ( h4 p9 p2 p6 )_v 1

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
  call ccdt_l2p_15_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_15_perm(sh4,sh3,sp1,sp2,i0_perm)
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
     call ccdt_l2p_15_perm(sh3,sh4,sp2,sp1,i0_perm)
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
     call ccdt_l2p_15_perm(sh4,sh3,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2p_15_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: p6,p9,sp6,sp9
  integer(c_int) :: spin_itm_hppp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp6 = 1,2
  do sp9 = 1,2
     spin_itm_hppp = tdcc_spin_dummy2(sh3,sp6,sp1,sp9)
     spin_int2x = tdcc_spin_int2x(sh4,sp9,sp2,sp6)
     if(spin_itm_hppp * spin_int2x == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l2p_15_1(sh3,sp6,sp1,sp9,itm_hppp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p6 = norb1+1,nact
     do p9 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hppp(h3,p6,p1,p9) * int2x(h4,p9,p2,p6,spin_int2x)
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
  end subroutine ccdt_l2p_15_perm
  !--------------------------------------------
end subroutine ccdt_l2p_15
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_15_1(sh3,sp6,sp1,sp9,i1)

!     i1 ( h3 p6 p1 p9 )_yt + = -1 * Sum ( h7 h8 p5 ) * t ( p5 p6 h7 h8 )_t * y ( h3 h7 h8 p1 p5 p9 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sp6,sp1,sp9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,p6,p1,p9
  integer(c_int) :: h7,h8,p5,sh7,sh8,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit
!
! Demanding: ORDER-7
!
  do sh7 = 1,2
  do sh8 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh7,sh8,sp1,sp5,sp9)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do p6 = norb1+1,nact
     do p1 = norb1+1,nact
     do p9 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(h3,p6,p1,p9) = i1(h3,p6,p1,p9) + fact * &
             t2inp(p5,p6,h7,h8,spin_t2inp) * g3inp(h3,h7,h8,p1,p5,p9,spin_g3inp)
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
end subroutine ccdt_l2p_15_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_16(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = -1/2 * Sum ( h8 p9 ) * i1 ( h3 h4 h8 p9 )_yt * v ( h8 p9 p1 p2 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1,int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: h8,p9,sh8,sp9
  integer(c_int) :: spin_itm_hhhp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh8 = 1,2
  do sp9 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh3,sh4,sh8,sp9)
     spin_int2x = tdcc_spin_int2x(sh8,sp9,sp1,sp2)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccdt_l2p_16_1(sh3,sh4,sh8,sp9,itm_hhhp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h8 = 1,norb1
     do p9 = norb1+1,nact
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hhhp(h3,h4,h8,p9) * int2x(h8,p9,p1,p2,spin_int2x)
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
  deallocate(itm_hhhp)
end subroutine ccdt_l2p_16
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_16_1(sh3,sh4,sh8,sp9,i1)

!     i1 ( h3 h4 h8 p9 )_yt + = 1 * Sum ( h7 p5 p6 ) * t ( p5 p6 h7 h8 )_t * y ( h3 h4 h7 p5 p6 p9 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sh8,sp9
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,h4,h8,p9
  integer(c_int) :: h7,p5,p6,sh7,sp5,sp6
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit
!
! Demanding: ORDER-7
!
  do sh7 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh4,sh7,sp5,sp6,sp9)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do h8 = 1,norb1
     do p9 = norb1+1,nact
     do h7 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h4,h8,p9) = i1(h3,h4,h8,p9) + fact * &
             t2inp(p5,p6,h7,h8,spin_t2inp) * g3inp(h3,h4,h7,p5,p6,p9,spin_g3inp)
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
end subroutine ccdt_l2p_16_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_17(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_vty + = 1/12 * Sum ( p5 h8 h9 h10 ) * y ( h8 h9 h10 p1 p2 p5 )_y * i1 ( h3 h4 p5 h8 h9 h10 )_vt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: p5,h8,h9,h10,sp5,sh8,sh9,sh10
  integer(c_int) :: spin_g3inp
  integer(c_int) :: spin_itm_hhphhh
  integer(c_int),external :: tdcc_spin_g3inp
  integer(c_int),external :: tdcc_spin_dummy3
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
     call ccdt_l2p_17_1(sh3,sh4,sp5,sh8,sh9,sh10,itm_hhphhh)

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
end subroutine ccdt_l2p_17
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_17_1(sh3,sh4,sp5,sh8,sh9,sh10,i1)

!     i1 ( h3 h4 p5 h8 h9 h10 )_vt + = 1 * Sum ( p6 p7 ) * t ( p5 p6 p7 h8 h9 h10 )_t * v ( h3 h4 p6 p7 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp5,sh8,sh9,sh10
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h3,h4,p5,h8,h9,h10
  integer(c_int) :: p6,p7,sp6,sp7
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
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
end subroutine ccdt_l2p_17_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_18(sh3,sh4,sp1,sp2,i0)

! i0 ( h3 h4 p1 p2 )_ytv + = 1/4 * P( 4 ) * Sum ( p7 h10 ) * i1 ( h3 p7 h10 p1 )_yt * v ( h4 h10 p2 p7 )_v 1

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
  call ccdt_l2p_18_perm(sh3,sh4,sp1,sp2,i0_perm)
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
     call ccdt_l2p_18_perm(sh4,sh3,sp1,sp2,i0_perm)
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
     call ccdt_l2p_18_perm(sh3,sh4,sp2,sp1,i0_perm)
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
     call ccdt_l2p_18_perm(sh4,sh3,sp2,sp1,i0_perm)
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
  subroutine ccdt_l2p_18_perm(sh3,sh4,sp1,sp2,i0)

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h3,h4,p1,p2
  integer(c_int) :: p7,h10,sp7,sh10
  integer(c_int) :: spin_itm_hphp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp7 = 1,2
  do sh10 = 1,2
     spin_itm_hphp = tdcc_spin_dummy2(sh3,sp7,sh10,sp1)
     spin_int2x = tdcc_spin_int2x(sh4,sh10,sp2,sp7)
     if(spin_itm_hphp * spin_int2x == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccdt_l2p_18_1(sh3,sp7,sh10,sp1,itm_hphp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p7 = norb1+1,nact
     do h10 = 1,norb1
        i0(h3,h4,p1,p2) = i0(h3,h4,p1,p2) + fact * &
             itm_hphp(h3,p7,h10,p1) * int2x(h4,h10,p2,p7,spin_int2x)
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
  end subroutine ccdt_l2p_18_perm
  !--------------------------------------------
end subroutine ccdt_l2p_18
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l2p_18_1(sh3,sp7,sh10,sp1,i1)

!     i1 ( h3 p7 h10 p1 )_yt + = 1 * Sum ( h8 h9 p5 p6 ) * t ( p5 p6 p7 h8 h9 h10 )_t * y ( h3 h8 h9 p1 p5 p6 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sh3,sp7,sh10,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h3,p7,h10,p1
  integer(c_int) :: h8,h9,p5,p6,sh8,sh9,sp5,sp6
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit
!
! Demanding: ORDER-8
!
  do sh8 = 1,2
  do sh9 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp6,sp7,sh8,sh9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh3,sh8,sh9,sp1,sp5,sp6)
     if(spin_t3inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do p7 = norb1+1,nact
     do h10 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,p7,h10,p1) = i1(h3,p7,h10,p1) + fact * &
             t3inp(p5,p6,p7,h8,h9,h10,spin_t3inp) * g3inp(h3,h8,h9,p1,p5,p6,spin_g3inp)
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
end subroutine ccdt_l2p_18_1
!##########################################################
!##########################################################
!##########################################################
