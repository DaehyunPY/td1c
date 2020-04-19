!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_1(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( h12 ) * t ( p4 p5 h1 h12 )_t * i1 ( h12 p6 h2 h3 )_v 5

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_1_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_1_perm(sp4,sp5,sp6,sh2,sh1,sh3,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h2,h1,h3)
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
     call ccdt_t3p_1_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
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
     call ccdt_t3p_1_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_1_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
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

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_1_perm(sp6,sp5,sp4,sh2,sh1,sh3,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_1_perm(sp4,sp6,sp5,sh2,sh1,sh3,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_1_perm(sp6,sp5,sp4,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_1_perm(sp4,sp6,sp5,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_t3p_1_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h12,sh12
  integer :: spin_t2inp
  integer :: spin_itm_hphh
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh12 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh1,sh12)
     spin_itm_hphh = tdcc_spin_int2x(sh12,sp6,sh2,sh3)
     if(spin_t2inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccdt_t3p_1_1(sh12,sp6,sh2,sh3,itm_hphh)
     call ccdt_t3p_1_2(sh12,sp6,sh2,sh3,itm_hphh)
     call ccdt_t3p_1_3(sh12,sp6,sh2,sh3,itm_hphh)
     call ccdt_t3p_1_4(sh12,sp6,sh2,sh3,itm_hphh)
     call ccdt_t3p_1_5(sh12,sp6,sh2,sh3,itm_hphh)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h12 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t2inp(p4,p5,h1,h12,spin_t2inp) * itm_hphh(h12,p6,h2,h3)
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
  deallocate(itm_hphh)

  end subroutine ccdt_t3p_1_perm
  !--------------------------------------------
end subroutine ccdt_t3p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_1_1(sh12,sp4,sh1,sh2,i1)

!     i1 ( h12 p4 h1 h2 )_v + = 1 * v ( h12 p4 h1 h2 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h12,p4,h1,h2
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh12,sp4,sh1,sh2)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h12,p4,h1,h2) = i1(h12,p4,h1,h2) + fact * int2x(h12,p4,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_t3p_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_1_2(sh12,sp4,sh1,sh2,i1)

!     i1 ( h12 p4 h1 h2 )_ft + = -1 * Sum ( p7 ) * t ( p4 p7 h1 h2 )_t * f ( h12 p7 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h12,p4,h1,h2
  integer :: p7,sp7
  integer :: spin_t2inp
  integer :: spin_fock
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp7,sh1,sh2)
     spin_fock = tdcc_spin_fock(sh12,sp7)
     if(spin_t2inp * spin_fock == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p7 = norb1+1,nact
        i1(h12,p4,h1,h2) = i1(h12,p4,h1,h2) + fact * &
             t2inp(p4,p7,h1,h2,spin_t2inp) * fock(h12,p7,spin_fock)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
end subroutine ccdt_t3p_1_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_1_3(sh12,sp4,sh1,sh2,i1)

!     i1 ( h12 p4 h1 h2 )_vt + = 1 * P( 2 ) * Sum ( h9 p8 ) * t ( p4 p8 h1 h9 )_t * v ( h9 h12 h2 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h12,p4,h1,h2
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:)

  allocate(i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1))

  i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
  call ccdt_t3p_1_3_perm(sh12,sp4,sh1,sh2,i1_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h12,p4,h1,h2) = i1(h12,p4,h1,h2) + fact_p * i1_perm(h12,p4,h1,h2)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sh12 * sp4 * sh1 * sh2 == 1)) then
     i1_perm(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call ccdt_t3p_1_3_perm(sh12,sp4,sh2,sh1,i1_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h12 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h12,p4,h1,h2) = i1(h12,p4,h1,h2) + fact_p * i1_perm(h12,p4,h2,h1)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_t3p_1_3_perm(sh12,sp4,sh1,sh2,i1)

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)

  integer :: h12,p4,h1,h2
  integer :: h9,p8,sh9,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh9 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp8,sh1,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sh12,sh2,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h9 = 1,norb1
     do p8 = norb1+1,nact
        i1(h12,p4,h1,h2) = i1(h12,p4,h1,h2) + fact * &
             t2inp(p4,p8,h1,h9,spin_t2inp) * int2x(h9,h12,h2,p8,spin_int2x)
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
  end subroutine ccdt_t3p_1_3_perm
  !--------------------------------------------
end subroutine ccdt_t3p_1_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_1_4(sh12,sp4,sh1,sh2,i1)

!     i1 ( h12 p4 h1 h2 )_vt + = 1/2 * Sum ( p8 p9 ) * t ( p8 p9 h1 h2 )_t * v ( h12 p4 p8 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h12,p4,h1,h2
  integer :: p8,p9,sp8,sp9
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp8 = 1,2
  do sp9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp8,sp9,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sh12,sp4,sp8,sp9)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
        i1(h12,p4,h1,h2) = i1(h12,p4,h1,h2) + fact * &
             t2inp(p8,p9,h1,h2,spin_t2inp) * int2x(h12,p4,p8,p9,spin_int2x)
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
end subroutine ccdt_t3p_1_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_1_5(sh12,sp4,sh1,sh2,i1)

!     i1 ( h12 p4 h1 h2 )_vt + = 1/2 * Sum ( h9 p7 p8 ) * t ( p4 p7 p8 h1 h2 h9 )_t * v ( h9 h12 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh12,sp4,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h12,p4,h1,h2
  integer :: h9,p7,p8,sh9,sp7,sp8
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh9 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp7,sp8,sh1,sh2,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sh12,sp7,sp8)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h12 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h9 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h12,p4,h1,h2) = i1(h12,p4,h1,h2) + fact * &
             t3inp(p4,p7,p8,h1,h2,h9,spin_t3inp) * int2x(h9,h12,p7,p8,spin_int2x)
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
end subroutine ccdt_t3p_1_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_2(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( p11 ) * t ( p4 p11 h1 h2 )_t * i1 ( p5 p6 h3 p11 )_v 4

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_2_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_2_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
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
     call ccdt_t3p_2_perm(sp4,sp5,sp6,sh1,sh3,sh2,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h3,h2)
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
     call ccdt_t3p_2_perm(sp5,sp4,sp6,sh1,sh2,sh3,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p5,p4,p6,h1,h2,h3)
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
     call ccdt_t3p_2_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_2_perm(sp5,sp4,sp6,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p5,p4,p6,h3,h2,h1)
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
     call ccdt_t3p_2_perm(sp6,sp5,sp4,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_2_perm(sp5,sp4,sp6,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p5,p4,p6,h1,h3,h2)
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
     call ccdt_t3p_2_perm(sp6,sp5,sp4,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h3,h2)
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
  subroutine ccdt_t3p_2_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: p11,sp11
  integer :: spin_t2inp
  integer :: spin_itm_pphp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp11 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp11,sh1,sh2)
     spin_itm_pphp = tdcc_spin_int2x(sp5,sp6,sh3,sp11)
     if(spin_t2inp * spin_itm_pphp == 0) cycle

     itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccdt_t3p_2_1(sp5,sp6,sh3,sp11,itm_pphp)
     call ccdt_t3p_2_2(sp5,sp6,sh3,sp11,itm_pphp)
     call ccdt_t3p_2_3(sp5,sp6,sh3,sp11,itm_pphp)
     call ccdt_t3p_2_4(sp5,sp6,sh3,sp11,itm_pphp)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p11 = norb1+1,nact
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t2inp(p4,p11,h1,h2,spin_t2inp) * itm_pphp(p5,p6,h3,p11)
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
  deallocate(itm_pphp)
  end subroutine ccdt_t3p_2_perm
  !--------------------------------------------
end subroutine ccdt_t3p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_2_1(sp4,sp5,sh1,sp11,i1)

!     i1 ( p4 p5 h1 p11 )_v + = 1 * v ( p4 p5 h1 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p4,p5,h1,p11
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sp4,sp5,sh1,sp11)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do p11 = norb1+1,nact
     i1(p4,p5,h1,p11) = i1(p4,p5,h1,p11) + fact * int2x(p4,p5,h1,p11,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_t3p_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_2_2(sp4,sp5,sh1,sp11,i1)

!     i1 ( p4 p5 h1 p11 )_vt + = 1/2 * Sum ( h8 h9 ) * t ( p4 p5 h8 h9 )_t * v ( h8 h9 h1 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p4,p5,h1,p11
  integer :: h8,h9,sh8,sh9
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sh9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh8,sh9,sh1,sp11)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do p11 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
        i1(p4,p5,h1,p11) = i1(p4,p5,h1,p11) + fact * &
             t2inp(p4,p5,h8,h9,spin_t2inp) * int2x(h8,h9,h1,p11,spin_int2x)
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
end subroutine ccdt_t3p_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_2_3(sp4,sp5,sh1,sp11,i1)

!     i1 ( p4 p5 h1 p11 )_vt + = 1 * P( 2 ) * Sum ( h9 p8 ) * t ( p4 p8 h1 h9 )_t * v ( h9 p5 p8 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p4,p5,h1,p11
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i1_perm(:,:,:,:)

  allocate(i1_perm((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact))

  i1_perm((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
  call ccdt_t3p_2_3_perm(sp4,sp5,sh1,sp11,i1_perm)
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do p11 = norb1+1,nact
     i1(p4,p5,h1,p11) = i1(p4,p5,h1,p11) + fact_p * i1_perm(p4,p5,h1,p11)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sh1 * sp11 == 1)) then
     i1_perm((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccdt_t3p_2_3_perm(sp5,sp4,sh1,sp11,i1_perm)
  end if
  fact_p = -1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do h1 = 1,norb1
  do p11 = norb1+1,nact
     i1(p4,p5,h1,p11) = i1(p4,p5,h1,p11) + fact_p * i1_perm(p5,p4,h1,p11)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(i1_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_t3p_2_3_perm(sp4,sp5,sh1,sp11,i1)

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)

  integer :: p4,p5,h1,p11
  integer :: h9,p8,sh9,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh9 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp8,sh1,sh9)
     spin_int2x = tdcc_spin_int2x(sh9,sp5,sp8,sp11)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do p11 = norb1+1,nact
     do h9 = 1,norb1
     do p8 = norb1+1,nact
        i1(p4,p5,h1,p11) = i1(p4,p5,h1,p11) + fact * &
             t2inp(p4,p8,h1,h9,spin_t2inp) * int2x(h9,p5,p8,p11,spin_int2x)
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
  end subroutine ccdt_t3p_2_3_perm
  !--------------------------------------------
end subroutine ccdt_t3p_2_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_2_4(sp4,sp5,sh1,sp11,i1)

!     i1 ( p4 p5 h1 p11 )_vt + = 1/2 * Sum ( h8 h9 p7 ) * t ( p4 p5 p7 h1 h8 h9 )_t * v ( h8 h9 p7 p11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sh1,sp11
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p4,p5,h1,p11
  integer :: h8,h9,p7,sh8,sh9,sp7
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sh9 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp7,sh1,sh8,sh9)
     spin_int2x = tdcc_spin_int2x(sh8,sh9,sp7,sp11)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do p11 = norb1+1,nact
     do h8 = 1,norb1
     do h9 = 1,norb1
     do p7 = norb1+1,nact
        i1(p4,p5,h1,p11) = i1(p4,p5,h1,p11) + fact * &
             t3inp(p4,p5,p7,h1,h8,h9,spin_t3inp) * int2x(h8,h9,p7,p11,spin_int2x)
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
end subroutine ccdt_t3p_2_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_3(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_tf + = -1 * P( 3 ) * Sum ( h7 ) * t ( p4 p5 p6 h1 h2 h7 )_t * i1 ( h7 h3 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_3_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_3_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
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
     call ccdt_t3p_3_perm(sp4,sp5,sp6,sh1,sh3,sh2,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h3,h2)
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
  subroutine ccdt_t3p_3_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h7,sh7
  integer :: spin_t3inp
  integer :: spin_itm_hh
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp6,sh1,sh2,sh7)
     spin_itm_hh = tdcc_spin_fock(sh7,sh3)
     if(spin_t3inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccdt_t3p_3_1(sh7,sh3,itm_hh)
     call ccdt_t3p_3_2(sh7,sh3,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h7 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p6,h1,h2,h7,spin_t3inp) * itm_hh(h7,h3)
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
  deallocate(itm_hh)
  end subroutine ccdt_t3p_3_perm
  !--------------------------------------------
end subroutine ccdt_t3p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_3_1(sh7,sh1,i1)

!     i1 ( h7 h1 )_f + = 1 * f ( h7 h1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh7,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h7,h1
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh7,sh1)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do h7 = 1,norb1
  do h1 = 1,norb1
     i1(h7,h1) = i1(h7,h1) + fact * fock(h7,h1,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_t3p_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_3_2(sh7,sh1,i1)

!     i1 ( h7 h1 )_vt + = 1/2 * Sum ( h10 p8 p9 ) * t ( p8 p9 h1 h10 )_t * v ( h7 h10 p8 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh7,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h7,h1
  integer :: h10,p8,p9,sh10,sp8,sp9
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh10 = 1,2
  do sp8 = 1,2
  do sp9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp8,sp9,sh1,sh10)
     spin_int2x = tdcc_spin_int2x(sh7,sh10,sp8,sp9)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h7 = 1,norb1
     do h1 = 1,norb1
     do h10 = 1,norb1
     do p8 = norb1+1,nact
     do p9 = norb1+1,nact
        i1(h7,h1) = i1(h7,h1) + fact * t2inp(p8,p9,h1,h10,spin_t2inp) * int2x(h7,h10,p8,p9,spin_int2x)
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
end subroutine ccdt_t3p_3_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_4(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_tf + = 1 * P( 3 ) * Sum ( p7 ) * t ( p4 p5 p7 h1 h2 h3 )_t * i1 ( p6 p7 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_4_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_4_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_4_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
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
  subroutine ccdt_t3p_4_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: p7,sp7
  integer :: spin_t3inp
  integer :: spin_itm_pp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp7,sh1,sh2,sh3)
     spin_itm_pp = tdcc_spin_fock(sp6,sp7)
     if(spin_t3inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_t3p_4_1(sp6,sp7,itm_pp)
     call ccdt_t3p_4_2(sp6,sp7,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p7 = norb1+1,nact
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p7,h1,h2,h3,spin_t3inp) * itm_pp(p6,p7)
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
  deallocate(itm_pp)
  end subroutine ccdt_t3p_4_perm
  !--------------------------------------------
end subroutine ccdt_t3p_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_4_1(sp4,sp7,i1)

!     i1 ( p4 p7 )_f + = 1 * f ( p4 p7 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sp4,sp7
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p4,p7
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp4,sp7)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do p4 = norb1+1,nact
  do p7 = norb1+1,nact
     i1(p4,p7) = i1(p4,p7) + fact * fock(p4,p7,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_t3p_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_4_2(sp4,sp7,i1)

!     i1 ( p4 p7 )_vt + = -1/2 * Sum ( h9 h10 p8 ) * t ( p4 p8 h9 h10 )_t * v ( h9 h10 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp7
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p4,p7
  integer :: h9,h10,p8,sh9,sh10,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh9 = 1,2
  do sh10 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp8,sh9,sh10)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp7,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p4 = norb1+1,nact
     do p7 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p8 = norb1+1,nact
        i1(p4,p7) = i1(p4,p7) + fact * t2inp(p4,p8,h9,h10,spin_t2inp) * int2x(h9,h10,p7,p8,spin_int2x)
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
end subroutine ccdt_t3p_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_5(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = 1/2 * P( 3 ) * Sum ( h7 h8 ) * t ( p4 p5 p6 h1 h7 h8 )_t * i1 ( h7 h8 h2 h3 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_5_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_5_perm(sp4,sp5,sp6,sh2,sh1,sh3,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h2,h1,h3)
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
     call ccdt_t3p_5_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
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
  subroutine ccdt_t3p_5_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h7,h8,sh7,sh8
  integer :: spin_t3inp
  integer :: spin_itm_hhhh
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh7 = 1,2
  do sh8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp6,sh1,sh7,sh8)
     spin_itm_hhhh = tdcc_spin_int2x(sh7,sh8,sh2,sh3)
     if(spin_t3inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_5_1(sh7,sh8,sh2,sh3,itm_hhhh)
     call ccdt_t3p_5_2(sh7,sh8,sh2,sh3,itm_hhhh)

!!!!!DEMANDING!!!!!

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h7 = 1,norb1
     do h8 = 1,norb1
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p6,h1,h7,h8,spin_t3inp) * itm_hhhh(h7,h8,h2,h3)
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
  deallocate(itm_hhhh)
  end subroutine ccdt_t3p_5_perm
  !--------------------------------------------
end subroutine ccdt_t3p_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_5_1(sh7,sh8,sh1,sh2,i1)

!     i1 ( h7 h8 h1 h2 )_v + = 1 * v ( h7 h8 h1 h2 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh7,sh8,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h7,h8,h1,h2
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh7,sh8,sh1,sh2)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h7 = 1,norb1
  do h8 = 1,norb1
  do h1 = 1,norb1
  do h2 = 1,norb1
     i1(h7,h8,h1,h2) = i1(h7,h8,h1,h2) + fact * int2x(h7,h8,h1,h2,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_t3p_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_5_2(sh7,sh8,sh1,sh2,i1)

!     i1 ( h7 h8 h1 h2 )_vt + = 1/2 * Sum ( p9 p10 ) * t ( p9 p10 h1 h2 )_t * v ( h7 h8 p9 p10 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh7,sh8,sh1,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h7,h8,h1,h2
  integer :: p9,p10,sp9,sp10
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp9 = 1,2
  do sp10 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp9,sp10,sh1,sh2)
     spin_int2x = tdcc_spin_int2x(sh7,sh8,sp9,sp10)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h7 = 1,norb1
     do h8 = 1,norb1
     do h1 = 1,norb1
     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
        i1(h7,h8,h1,h2) = i1(h7,h8,h1,h2) + fact * &
             t2inp(p9,p10,h1,h2,spin_t2inp) * int2x(h7,h8,p9,p10,spin_int2x)
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
end subroutine ccdt_t3p_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_6(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = -1 * P( 9 ) * Sum ( h8 p7 ) * t ( p4 p5 p7 h1 h2 h8 )_t * i1 ( h8 p6 h3 p7 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t3inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_6_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_6_perm(sp4,sp5,sp6,sh3,sh2,sh1,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h3,h2,h1)
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
     call ccdt_t3p_6_perm(sp4,sp5,sp6,sh1,sh3,sh2,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p5,p6,h1,h3,h2)
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
     call ccdt_t3p_6_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_6_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
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

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_6_perm(sp6,sp5,sp4,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_6_perm(sp4,sp6,sp5,sh3,sh2,sh1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sp4 * sp5 * sp6 * sh1 * sh2 * sh3 == 1)) then
     i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_6_perm(sp6,sp5,sp4,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p6,p5,p4,h1,h3,h2)
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
     call ccdt_t3p_6_perm(sp4,sp6,sp5,sh1,sh3,sh2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do p4 = norb1+1,nact
  do p5 = norb1+1,nact
  do p6 = norb1+1,nact
  do h1 = 1,norb1
  do h2 = 1,norb1
  do h3 = 1,norb1
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p4,p6,p5,h1,h3,h2)
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
  subroutine ccdt_t3p_6_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h8,p7,sh8,sp7
  integer :: spin_t3inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh8 = 1,2
  do sp7 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp7,sh1,sh2,sh8)
     spin_itm_hphp = tdcc_spin_int2x(sh8,sp6,sh3,sp7)
     if(spin_t3inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccdt_t3p_6_1(sh8,sp6,sh3,sp7,itm_hphp)
     call ccdt_t3p_6_2(sh8,sp6,sh3,sp7,itm_hphp)

!!!!!DEMANDING!!!!!

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h8 = 1,norb1
     do p7 = norb1+1,nact
        i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact * &
             t3inp(p4,p5,p7,h1,h2,h8,spin_t3inp) * itm_hphp(h8,p6,h3,p7)
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
  deallocate(itm_hphp)
  end subroutine ccdt_t3p_6_perm
  !--------------------------------------------
end subroutine ccdt_t3p_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_6_1(sh8,sp4,sh1,sp7,i1)

!     i1 ( h8 p4 h1 p7 )_v + = 1 * v ( h8 p4 h1 p7 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh8,sp4,sh1,sp7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h8,p4,h1,p7
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh8,sp4,sh1,sp7)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h8 = 1,norb1
  do p4 = norb1+1,nact
  do h1 = 1,norb1
  do p7 = norb1+1,nact
     i1(h8,p4,h1,p7) = i1(h8,p4,h1,p7) + fact * int2x(h8,p4,h1,p7,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_t3p_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_6_2(sh8,sp4,sh1,sp7,i1)

!     i1 ( h8 p4 h1 p7 )_vt + = -1 * Sum ( h10 p9 ) * t ( p4 p9 h1 h10 )_t * v ( h8 h10 p7 p9 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : t2inp
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh8,sp4,sh1,sp7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h8,p4,h1,p7
  integer :: h10,p9,sh10,sp9
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh10 = 1,2
  do sp9 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp9,sh1,sh10)
     spin_int2x = tdcc_spin_int2x(sh8,sh10,sp7,sp9)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h8 = 1,norb1
     do p4 = norb1+1,nact
     do h1 = 1,norb1
     do p7 = norb1+1,nact
     do h10 = 1,norb1
     do p9 = norb1+1,nact
        i1(h8,p4,h1,p7) = i1(h8,p4,h1,p7) + fact * &
             t2inp(p4,p9,h1,h10,spin_t2inp) * int2x(h8,h10,p7,p9,spin_int2x)
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
end subroutine ccdt_t3p_6_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_7(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vt + = 1/2 * P( 3 ) * Sum ( p7 p8 ) * t ( p4 p7 p8 h1 h2 h3 )_t * v ( p5 p6 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_7_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_7_perm(sp5,sp4,sp6,sh1,sh2,sh3,i0_perm)
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
     i0(p4,p5,p6,h1,h2,h3) = i0(p4,p5,p6,h1,h2,h3) + fact_p * i0_perm(p5,p4,p6,h1,h2,h3)
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
     call ccdt_t3p_7_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
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

  deallocate(i0_perm)
!DEBUG
!  do p4 = norb1+1,nact
!  do p5 = norb1+1,p4-1
!  do p6 = norb1+1,nact
!  do h1 = 1,norb1
!  do h2 = 1,h1-1
!  do h3 = 1,norb1
!     write(6,"('ccdt_t3p_7: ',6i5,f20.10)") p4,p5,p6,h1,h2,h3,dble(i0(p4,p5,p6,h1,h2,h3))
!  end do
!  end do
!  end do
!  end do
!  end do
!  end do
!DEBUG

  contains
  !--------------------------------------------
  subroutine ccdt_t3p_7_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: p7,p8,sp7,sp8
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp7,sp8,sh1,sh2,sh3)
     spin_int2x = tdcc_spin_int2x(sp5,sp6,sp7,sp8)
     if(spin_t3inp * spin_int2x == 0) cycle

!!!!!DEMANDING!!!!!

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
             t3inp(p4,p7,p8,h1,h2,h3,spin_t3inp) * int2x(p5,p6,p7,p8,spin_int2x)
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
  end subroutine ccdt_t3p_7_perm
  !--------------------------------------------
end subroutine ccdt_t3p_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_8(sp4,sp5,sp6,sh1,sh2,sh3,i0)

! i0 ( p4 p5 p6 h1 h2 h3 )_vtt + = 1/4 * P( 3 ) * Sum ( h9 h10 ) * t ( p4 p5 h9 h10 )_t * i1 ( h9 h10 p6 h1 h2 h3 )_vt 1

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: p4,p5,p6,h1,h2,h3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))

  i0_perm((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
  call ccdt_t3p_8_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_8_perm(sp6,sp5,sp4,sh1,sh2,sh3,i0_perm)
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
     call ccdt_t3p_8_perm(sp4,sp6,sp5,sh1,sh2,sh3,i0_perm)
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
  subroutine ccdt_t3p_8_perm(sp4,sp5,sp6,sh1,sh2,sh3,i0)

  implicit none
  integer,intent(in) :: sp4,sp5,sp6,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i0((norb1+1):nact,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)

  integer :: p4,p5,p6,h1,h2,h3
  integer :: h9,h10,sh9,sh10
  integer :: spin_t2inp
  integer :: spin_itm_hhphhh
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hhphhh(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1))
  do sh9 = 1,2
  do sh10 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh9,sh10)
     spin_itm_hhphhh = tdcc_spin_dummy3(sh9,sh10,sp6,sh1,sh2,sh3)
     if(spin_t2inp * spin_itm_hhphhh == 0) cycle

     itm_hhphhh(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_t3p_8_1(sh9,sh10,sp6,sh1,sh2,sh3,itm_hhphhh)

!!!!!DEMANDING!!!!!

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
  end subroutine ccdt_t3p_8_perm
  !--------------------------------------------
end subroutine ccdt_t3p_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_t3p_8_1(sh9,sh10,sp4,sh1,sh2,sh3,i1)

!     i1 ( h9 h10 p4 h1 h2 h3 )_vt + = 1 * Sum ( p7 p8 ) * t ( p4 p7 p8 h1 h2 h3 )_t * v ( h9 h10 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only :nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh9,sh10,sp4,sh1,sh2,sh3
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h9,h10,p4,h1,h2,h3
  integer :: p7,p8,sp7,sp8
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp7,sp8,sh1,sh2,sh3)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp7,sp8)
     if(spin_t3inp * spin_int2x == 0) cycle

!!!!!DEMANDING!!!!!

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
end subroutine ccdt_t3p_8_1
!##########################################################
!##########################################################
!##########################################################
