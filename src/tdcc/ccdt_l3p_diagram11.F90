!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_diagram11_1(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = 1/4 * P( 3 ) * Sum ( h9 h10 ) * i1 ( h4 h5 h6 h9 h10 p1 )_yt * v ( h9 h10 p2 p3 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : int2x,norb1

  implicit none
  integer(c_int),intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_diagram11_1_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_diagram11_1_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_diagram11_1_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3p_diagram11_1_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer(c_int),intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h4,h5,h6,p1,p2,p3
  integer(c_int) :: h9,h10,sh9,sh10
  integer(c_int) :: spin_itm_hhhhhp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy3
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhhhp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sh10 = 1,2
     spin_itm_hhhhhp = tdcc_spin_dummy3(sh4,sh5,sh6,sh9,sh10,sp1)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp2,sp3)
     if(spin_itm_hhhhhp * spin_int2x == 0) cycle

     itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccdt_l3p_diagram11_1_1(sh4,sh5,sh6,sh9,sh10,sp1,itm_hhhhhp)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             itm_hhhhhp(h4,h5,h6,h9,h10,p1) * int2x(h9,h10,p2,p3,spin_int2x)
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
  deallocate(itm_hhhhhp)
  end subroutine ccdt_l3p_diagram11_1_perm
  !--------------------------------------------
end subroutine ccdt_l3p_diagram11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_diagram11_1_1(sh4,sh5,sh6,sh9,sh10,sp1,i1)

!     i1 ( h4 h5 h6 h9 h10 p1 )_yt + = 1 * Sum ( p7 p8 ) * t ( p7 p8 h9 h10 )_t * y ( h4 h5 h6 p1 p7 p8 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : t2inp,norb1
  use mod_cc,only : g3inp

  implicit none
  integer(c_int),intent(in) :: sh4,sh5,sh6,sh9,sh10,sp1
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer(c_int) :: h4,h5,h6,h9,h10,p1
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g3inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp7,sp8)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h4,h5,h6,h9,h10,p1) = i1(h4,h5,h6,h9,h10,p1) + fact * &
             t2inp(p7,p8,h9,h10,spin_t2inp) * g3inp(h4,h5,h6,p1,p7,p8,spin_g3inp)
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
end subroutine ccdt_l3p_diagram11_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_diagram11_2(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = 1/4 * P( 3 ) * Sum ( p7 p8 ) * y ( h4 h5 h6 p1 p7 p8 )_y * i1 ( p7 p8 p2 p3 )_vt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : g3inp,norb1

  implicit none
  integer(c_int),intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_diagram11_2_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_diagram11_2_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p2,p1,p3)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_diagram11_2_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3p_diagram11_2_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer(c_int),intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer(c_int) :: h4,h5,h6,p1,p2,p3
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_g3inp
  integer(c_int) :: spin_itm_pppp
  integer(c_int),external :: tdcc_spin_g3inp
  integer(c_int),external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_pppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
  do sp8 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp7,sp8)
     spin_itm_pppp = tdcc_spin_dummy2(sp7,sp8,sp2,sp3)
     if(spin_g3inp * spin_itm_pppp == 0) cycle

     itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_diagram11_2_1(sp7,sp8,sp2,sp3,itm_pppp)

!!!!!DEMANDING!!!!!

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h6,p1,p7,p8,spin_g3inp) * itm_pppp(p7,p8,p2,p3)
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
  deallocate(itm_pppp)
  end subroutine ccdt_l3p_diagram11_2_perm
  !--------------------------------------------
end subroutine ccdt_l3p_diagram11_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_diagram11_2_1(sp7,sp8,sp2,sp3,i1)

!     i1 ( p7 p8 p2 p3 )_vt + = 1 * Sum ( h9 h10 ) * t ( p7 p8 h9 h10 )_t * v ( h9 h10 p2 p3 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : t2inp,norb1
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp7,sp8,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p7,p8,p2,p3
  integer(c_int) :: h9,h10,sh9,sh10
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh9 = 1,2
  do sh10 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh10)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp2,sp3)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
        i1(p7,p8,p2,p3) = i1(p7,p8,p2,p3) + fact * &
             t2inp(p7,p8,h9,h10,spin_t2inp) * int2x(h9,h10,p2,p3,spin_int2x)
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
end subroutine ccdt_l3p_diagram11_2_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccdt_l3p_diagram11_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : g3out,norb1

  implicit none

  call ccdt_l3p_diagram11_1(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
  call ccdt_l3p_diagram11_2(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))

  call ccdt_l3p_diagram11_1(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  call ccdt_l3p_diagram11_2(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))

end subroutine ccdt_l3p_diagram11_main
!**********************************************************
