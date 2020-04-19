!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_diagram8_1(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vtt + = 1/4 * P( 3 ) 
!  * Sum ( h9 h10 ) * t ( p4 p5 h9 h10 )_t * i1 ( h9 h10 p6 h1 h2 h3 )_vt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : t2inp,norb1

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_diagram8_1_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_diagram8_1_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_diagram8_1_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h1,h2,h3)
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
  subroutine ccdt_t3p_diagram8_1_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: p4,p5,p6,h1,h2,h3
  integer(c_int) :: h9,h10,sh9,sh10
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hhphhh
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhphhh(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh10 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh9,sh10)
     spin_itm_hhphhh = tdcc_spin_dummy3(sh9,sh10,sp6,sh1,sh2,sh3)
     if(spin_t2inp * spin_itm_hhphhh == 0) cycle

     itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_diagram8_1_1(sh9,sh10,sp6,sh1,sh2,sh3,itm_hhphhh)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h9 = 1,norb1
     do h10 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t2inp(p4,p5,h9,h10,spin_t2inp) * itm_hhphhh(h9,h10,p6,h1,h2,h3)
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
  deallocate(itm_hhphhh)
  end subroutine ccdt_t3p_diagram8_1_perm
  !--------------------------------------------
end subroutine ccdt_t3p_diagram8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_diagram8_1_1(sh9,sh10,sp4,sh1,sh2,sh3,i1)

!     i1 ( h9 h10 p4 h1 h2 h3 )_vt + = 1 * Sum ( p7 p8 ) * t ( p4 p7 p8 h1 h2 h3 )_t * v ( h9 h10 p7 p8 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : t3inp,norb1
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh9,sh10,sp4,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h9,h10,p4,h1,h2,h3
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp7,sp8,sh1,sh2,sh3)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp7,sp8)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h9,h10,p4,h1,h2,h3) = i1(h9,h10,p4,h1,h2,h3) + fact * &
             t3inp(p4,p7,p8,h1,h2,h3,spin_t3inp) * int2x(h9,h10,p7,p8,spin_int2x)
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
end subroutine ccdt_t3p_diagram8_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_diagram8_2(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vtt + = 1/4 * P( 3 ) * Sum ( p7 p8 ) * t ( p6 p7 p8 h1 h2 h3 )_t * i1 ( p4 p5 p7 p8 )_vt 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : t3inp,norb1

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_diagram8_2_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_diagram8_2_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h2,h3)
  end do
  end do
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_diagram8_2_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h1,h2,h3)
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
  subroutine ccdt_t3p_diagram8_2_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: p4,p5,p6,h1,h2,h3
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_itm_pppp
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_pppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
  do sp8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp6,sp7,sp8,sh1,sh2,sh3)
     spin_itm_pppp = tdcc_spin_dummy2(sp4,sp5,sp7,sp8)
     if(spin_t3inp * spin_itm_pppp == 0) cycle

     itm_pppp((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_t3p_diagram8_2_1(sp4,sp5,sp7,sp8,itm_pppp)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p6,p7,p8,h1,h2,h3,spin_t3inp) * itm_pppp(p4,p5,p7,p8)
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
  end subroutine ccdt_t3p_diagram8_2_perm
  !--------------------------------------------
end subroutine ccdt_t3p_diagram8_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_diagram8_2_1(sp4,sp5,sp7,sp8,i1)

!     i1 ( p4 p5 p7 p8 )_vt + = 1 * Sum ( h9 h10 ) * v ( h9 h10 p7 p8 )_v * t ( p4 p5 h9 h10 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : int2x,norb1
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp7,sp8
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p4,p5,p7,p8
  integer(c_int) :: h9,h10,sh9,sh10
  integer(c_int) :: spin_int2x
  integer(c_int) :: spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh9 = 1,2
  do sh10 = 1,2
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp7,sp8)
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh9,sh10)
     if(spin_int2x * spin_t2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
        i1(p4,p5,p7,p8) = i1(p4,p5,p7,p8) + fact * &
             int2x(h9,h10,p7,p8,spin_int2x) * t2inp(p4,p5,h9,h10,spin_t2inp)
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
end subroutine ccdt_t3p_diagram8_2_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccdt_t3p_diagram8_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : t3out,norb1

  implicit none

  call ccdt_t3p_diagram8_1(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
  call ccdt_t3p_diagram8_2(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))

  call ccdt_t3p_diagram8_1(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
  call ccdt_t3p_diagram8_2(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))

end subroutine ccdt_t3p_diagram8_main
!**********************************************************
