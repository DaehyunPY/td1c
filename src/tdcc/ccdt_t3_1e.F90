!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3_1e_1(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_tf + = -1 * P( 3 ) * Sum ( h7 ) * t ( p4 p5 p6 h1 h2 h7 )_t * f ( h7 h3 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : fock

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3_1e_1_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
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

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_1_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_1_perm(sp4,sp5,sp6,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h3,h2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_t3_1e_1_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: p4,p5,p6,h1,h2,h3
  integer(c_int) :: h7,sh7
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp6,sh1,sh2,sh7)
     spin_fock = tdcc_spin_fock(sh7,sh3)
     if(spin_t3inp * spin_fock == 0) cycle

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h7 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p6,h1,h2,h7,spin_t3inp) * fock(h7,h3,spin_fock)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccdt_t3_1e_1_perm
  !--------------------------------------------
end subroutine ccdt_t3_1e_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3_1e_2(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_tf + = 1 * P( 3 ) * Sum ( p7 ) * t ( p4 p5 p7 h1 h2 h3 )_t * f ( p6 p7 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : fock

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3_1e_2_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
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

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_2_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
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

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_2_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
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

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_t3_1e_2_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: p4,p5,p6,h1,h2,h3
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp7,sh1,sh2,sh3)
     spin_fock = tdcc_spin_fock(sp6,sp7)
     if(spin_t3inp * spin_fock == 0) cycle

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p7 = norb1+1,nact
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p7,h1,h2,h3,spin_t3inp) * fock(p6,p7,spin_fock)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccdt_t3_1e_2_perm
  !--------------------------------------------
end subroutine ccdt_t3_1e_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3_1e_3(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_ftt + = -1 * P( 9 ) * Sum ( h8 ) * t ( p4 p5 h1 h8 )_t * i1 ( h8 p6 h2 h3 )_ft 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3_1e_3_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
  fact_p = +1.0d+0
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

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_3_perm(sp4,sp5,sp6,sh2,sh1,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h2,h1,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_3_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_3_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
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

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_3_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
  end if
  fact_p = -1.0d+0
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

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_3_perm(sp6,sp5,sp4,sh2,sh1,sh3,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h2,h1,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_3_perm(sp4,sp6,sp5,sh2,sh1,sh3,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h2,h1,h3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_3_perm(sp6,sp5,sp4,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_3_perm(sp4,sp6,sp5,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h3,h2,h1)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_t3_1e_3_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer(c_int),intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer(c_int) :: p4,p5,p6,h1,h2,h3
  integer(c_int) :: h8,sh8
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hphh
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh1,sh8)
     spin_itm_hphh = tdcc_spin_dummy2(sh8,sp6,sh2,sh3)
     if(spin_t2inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccdt_t3_1e_3_1(sh8,sp6,sh2,sh3,itm_hphh)

     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h8 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t2inp(p4,p5,h1,h8,spin_t2inp) * itm_hphh(h8,p6,h2,h3)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  deallocate(itm_hphh)
  end subroutine ccdt_t3_1e_3_perm
  !--------------------------------------------
end subroutine ccdt_t3_1e_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3_1e_3_1(sh8,sp4,sh1,sh2,i1)

!     i1 ( h8 p4 h1 h2 )_ft + = -1 * Sum ( p7 ) * t ( p4 p7 h1 h2 )_t * f ( h8 p7 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : fock

  implicit none
  integer(c_int),intent(in) :: sh8,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h8,p4,h1,h2
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp7,sh1,sh2)
     spin_fock = tdcc_spin_fock(sh8,sp7)
     if(spin_t2inp * spin_fock == 0) cycle

     do h8 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p7 = norb1+1,nact
        i1(h8,p4,h1,h2) = i1(h8,p4,h1,h2) + fact * &
             t2inp(p4,p7,h1,h2,spin_t2inp) * fock(h8,p7,spin_fock)
     end do
     end do
     end do
     end do
     end do
  end do
end subroutine ccdt_t3_1e_3_1
!##########################################################
!##########################################################
!##########################################################
