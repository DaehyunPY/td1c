!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t2_1e_1(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_tf + = -1 * P( 2 ) * Sum ( h5 ) * t ( p3 p4 h1 h5 )_t * f ( h5 h2 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer :: p3,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccdt_t2_1e_1_perm(sp3,sp4,sh1,sh2,i0_perm)
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccdt_t2_1e_1_perm(sp3,sp4,sh2,sh1,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h2,h1)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_t2_1e_1_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer,intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer :: p3,p4,h1,h2
  integer :: h5,sh5
  integer :: spin_t2inp
  integer :: spin_fock
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh1,sh5)
     spin_fock = tdcc_spin_fock(sh5,sh2)
     if(spin_t2inp * spin_fock == 0) cycle

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h5 = 1,norb1
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * t2inp(p3,p4,h1,h5,spin_t2inp) * fock(h5,h2,spin_fock)
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccdt_t2_1e_1_perm
  !--------------------------------------------
end subroutine ccdt_t2_1e_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t2_1e_2(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_tf + = 1 * P( 2 ) * Sum ( p5 ) * t ( p3 p5 h1 h2 )_t * f ( p4 p5 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer :: p3,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccdt_t2_1e_2_perm(sp3,sp4,sh1,sh2,i0_perm)
  fact_p = +1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p3,p4,h1,h2)
  end do
  end do
  end do
  end do

  if(.not. (sp3 * sp4 * sh1 * sh2 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccdt_t2_1e_2_perm(sp4,sp3,sh1,sh2,i0_perm)
  end if
  fact_p = -1.0d+0
  do p3 = norb1+1,nact
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact_p * i0_perm(p4,p3,h1,h2)
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_t2_1e_2_perm(sp3,sp4,sh1,sh2,i0)

  implicit none
  integer,intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)

  integer :: p3,p4,h1,h2
  integer :: p5,sp5
  integer :: spin_t2inp
  integer :: spin_fock
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp5,sh1,sh2)
     spin_fock = tdcc_spin_fock(sp4,sp5)
     if(spin_t2inp * spin_fock == 0) cycle

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p5 = norb1+1,nact
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * t2inp(p3,p5,h1,h2,spin_t2inp) * fock(p4,p5,spin_fock)
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccdt_t2_1e_2_perm
  !--------------------------------------------
end subroutine ccdt_t2_1e_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t2_1e_3(sp3,sp4,sh1,sh2,i0)

! i0 ( p3 p4 h1 h2 )_tf + = 1 * Sum ( p6 h5 ) * t ( p3 p4 p6 h1 h2 h5 )_t * f ( h5 p6 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sp3,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1)
  integer :: p3,p4,h1,h2
  integer :: p6,h5,sp6,sh5
  integer :: spin_t3inp
  integer :: spin_fock
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp6 = 1,2
  do sh5 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp4,sp6,sh1,sh2,sh5)
     spin_fock = tdcc_spin_fock(sh5,sp6)
     if(spin_t3inp * spin_fock == 0) cycle

     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p6 = norb1+1,nact
     do h5 = 1,norb1
        i0(p3,p4,h1,h2) = i0(p3,p4,h1,h2) + fact * t3inp(p3,p4,p6,h1,h2,h5,spin_t3inp) * fock(h5,p6,spin_fock)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccdt_t2_1e_3
!##########################################################
!##########################################################
!##########################################################
