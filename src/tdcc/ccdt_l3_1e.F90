!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3_1e_1(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yf + = 1 * P( 9 ) * y ( h4 h5 p1 p2 )_y * f ( h6 p3 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3_1e_1_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_1_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_1_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_1_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_1_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_1_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_1_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_1_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p3,p2,p1)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_1_perm(sh4,sh6,sh5,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3_1e_1_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: sdum
  integer :: spin_g2inp
  integer :: spin_fock
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sdum = 1,1
     spin_g2inp = tdcc_spin_g2inp(sh4,sh5,sp1,sp2)
     spin_fock = tdcc_spin_fock(sh6,sp3)
     if(spin_g2inp * spin_fock == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g2inp(h4,h5,p1,p2,spin_g2inp) * fock(h6,p3,spin_fock)
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccdt_l3_1e_1_perm
  !--------------------------------------------
end subroutine ccdt_l3_1e_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3_1e_2(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yf + = -1 * P( 3 ) * Sum ( h7 ) * y ( h4 h5 h7 p1 p2 p3 )_y * f ( h6 h7 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3_1e_2_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_2_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_2_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3_1e_2_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h7,sh7
  integer :: spin_g3inp
  integer :: spin_fock
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh7,sp1,sp2,sp3)
     spin_fock = tdcc_spin_fock(sh6,sh7)
     if(spin_g3inp * spin_fock == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h7 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h7,p1,p2,p3,spin_g3inp) * fock(h6,h7,spin_fock)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccdt_l3_1e_2_perm
  !--------------------------------------------
end subroutine ccdt_l3_1e_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3_1e_3(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yf + = 1 * P( 3 ) * Sum ( p7 ) * y ( h4 h5 h6 p1 p2 p7 )_y * f ( p7 p3 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3_1e_3_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
  fact_p = +1.0d+0
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_3_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = -1.0d+0
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3_1e_3_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = -1.0d+0
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
  end do
  end do
  end do
  end do
  end do
  end do

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3_1e_3_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: p7,sp7
  integer :: spin_g3inp
  integer :: spin_fock
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp2,sp7)
     spin_fock = tdcc_spin_fock(sp7,sp3)
     if(spin_g3inp * spin_fock == 0) cycle

     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h6,p1,p2,p7,spin_g3inp) * fock(p7,p3,spin_fock)
     end do
     end do
     end do
     end do
     end do
     end do
     end do
  end do
  end subroutine ccdt_l3_1e_3_perm
  !--------------------------------------------
end subroutine ccdt_l3_1e_3
!##########################################################
!##########################################################
!##########################################################
