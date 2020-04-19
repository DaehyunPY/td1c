!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_1(sh4,sh5,sh6,sp1,sp2,sp3,i0)

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
  call ccdt_l3p_1_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_1_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
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
     call ccdt_l3p_1_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
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
     call ccdt_l3p_1_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_1_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
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
     call ccdt_l3p_1_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_1_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_1_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_1_perm(sh4,sh6,sh5,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3p_1_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

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

     !$omp parallel default(shared)
     !$omp do collapse(6)
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
     !$omp end do
     !$omp end parallel
  end do
  end subroutine ccdt_l3p_1_perm
  !--------------------------------------------
end subroutine ccdt_l3p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_2(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = -1 * P( 9 ) * Sum ( h7 ) * y ( h4 h7 p1 p2 )_y * v ( h5 h6 h7 p3 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_2_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_2_perm(sh5,sh4,sh6,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p2,p3)
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
     call ccdt_l3p_2_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
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
     call ccdt_l3p_2_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_2_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
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
     call ccdt_l3p_2_perm(sh5,sh4,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p3,p2,p1)
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
     call ccdt_l3p_2_perm(sh5,sh4,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p3,p2)
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
     call ccdt_l3p_2_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_2_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3p_2_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h7,sh7
  integer :: spin_g2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh7 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh4,sh7,sp1,sp2)
     spin_int2x = tdcc_spin_int2x(sh5,sh6,sh7,sp3)
     if(spin_g2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h7 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g2inp(h4,h7,p1,p2,spin_g2inp) * int2x(h5,h6,h7,p3,spin_int2x)
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
  end subroutine ccdt_l3p_2_perm
  !--------------------------------------------
end subroutine ccdt_l3p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_3(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = -1 * P( 9 ) * Sum ( p7 ) * y ( h4 h5 p1 p7 )_y * v ( h6 p7 p2 p3 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_3_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_3_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
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
     call ccdt_l3p_3_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
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
     call ccdt_l3p_3_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
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
     call ccdt_l3p_3_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_3_perm(sh6,sh5,sh4,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p2,p1,p3)
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
     call ccdt_l3p_3_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_3_perm(sh4,sh6,sh5,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p2,p1,p3)
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
     call ccdt_l3p_3_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3p_3_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: p7,sp7
  integer :: spin_g2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp7 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh4,sh5,sp1,sp7)
     spin_int2x = tdcc_spin_int2x(sh6,sp7,sp2,sp3)
     if(spin_g2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g2inp(h4,h5,p1,p7,spin_g2inp) * int2x(h6,p7,p2,p3,spin_int2x)
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
  end subroutine ccdt_l3p_3_perm
  !--------------------------------------------
end subroutine ccdt_l3p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_4(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yf + = -1 * P( 3 ) * Sum ( h11 ) * y ( h4 h5 h11 p1 p2 p3 )_y * i1 ( h6 h11 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_4_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_4_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
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
     call ccdt_l3p_4_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
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
  subroutine ccdt_l3p_4_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h11,sh11
  integer :: spin_g3inp
  integer :: spin_itm_hh
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh11 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh11,sp1,sp2,sp3)
     spin_itm_hh = tdcc_spin_fock(sh6,sh11)
     if(spin_g3inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccdt_l3p_4_1(sh6,sh11,itm_hh)
     call ccdt_l3p_4_2(sh6,sh11,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h11 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h11,p1,p2,p3,spin_g3inp) * itm_hh(h6,h11)
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
  end subroutine ccdt_l3p_4_perm
  !--------------------------------------------
end subroutine ccdt_l3p_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_4_1(sh4,sh11,i1)

!     i1 ( h4 h11 )_f + = 1 * f ( h4 h11 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh4,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h4,h11
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh4,sh11)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do h4 = 1,norb1
  do h11 = 1,norb1
     i1(h4,h11) = i1(h4,h11) + fact * fock(h4,h11,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_4_2(sh4,sh11,i1)

!     i1 ( h4 h11 )_vt + = -1/2 * Sum ( h10 p7 p8 ) * t ( p7 p8 h10 h11 )_t * v ( h4 h10 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h4,h11
  integer :: h10,p7,p8,sh10,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh10 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh10,sh11)
     spin_int2x = tdcc_spin_int2x(sh4,sh10,sp7,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h4 = 1,norb1
     do h11 = 1,norb1
     do h10 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h4,h11) = i1(h4,h11) + fact * &
             t2inp(p7,p8,h10,h11,spin_t2inp) * int2x(h4,h10,p7,p8,spin_int2x)
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
end subroutine ccdt_l3p_4_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_5(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yf + = 1 * P( 3 ) * Sum ( p7 ) * y ( h4 h5 h6 p1 p2 p7 )_y * i1 ( p7 p3 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_5_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_5_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_5_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
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
  subroutine ccdt_l3p_5_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: p7,sp7
  integer :: spin_g3inp
  integer :: spin_itm_pp
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp2,sp7)
     spin_itm_pp = tdcc_spin_fock(sp7,sp3)
     if(spin_g3inp * spin_itm_pp == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_5_1(sp7,sp3,itm_pp)
     call ccdt_l3p_5_2(sp7,sp3,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h6,p1,p2,p7,spin_g3inp) * itm_pp(p7,p3)
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
  end subroutine ccdt_l3p_5_perm
  !--------------------------------------------
end subroutine ccdt_l3p_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_5_1(sp7,sp1,i1)

!     i1 ( p7 p1 )_f + = 1 * f ( p7 p1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sp7,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p7,p1
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp7,sp1)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do p7 = norb1+1,nact
  do p1 = norb1+1,nact
     i1(p7,p1) = i1(p7,p1) + fact * fock(p7,p1,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_5_2(sp7,sp1,i1)

!     i1 ( p7 p1 )_vt + = -1/2 * Sum ( h9 h10 p8 ) * t ( p7 p8 h9 h10 )_t * v ( h9 h10 p1 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sp7,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p7,p1
  integer :: h9,h10,p8,sh9,sh10,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh9 = 1,2
  do sh10 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh10)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp1,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p7 = norb1+1,nact
     do p1 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p8 = norb1+1,nact
        i1(p7,p1) = i1(p7,p1) + fact * &
             t2inp(p7,p8,h9,h10,spin_t2inp) * int2x(h9,h10,p1,p8,spin_int2x)
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
end subroutine ccdt_l3p_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_6(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = -1/2 * P( 3 ) * Sum ( h11 h9 ) * y ( h4 h9 h11 p1 p2 p3 )_y * i1 ( h5 h6 h9 h11 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_6_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_6_perm(sh5,sh4,sh6,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p2,p3)
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
     call ccdt_l3p_6_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
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
  subroutine ccdt_l3p_6_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h11,h9,sh11,sh9
  integer :: spin_g3inp
  integer :: spin_itm_hhhh
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh11 = 1,2
  do sh9 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh9,sh11,sp1,sp2,sp3)
     spin_itm_hhhh = tdcc_spin_int2x(sh5,sh6,sh9,sh11)
     if(spin_g3inp * spin_itm_hhhh == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call ccdt_l3p_6_1(sh5,sh6,sh9,sh11,itm_hhhh)
     call ccdt_l3p_6_2(sh5,sh6,sh9,sh11,itm_hhhh)

!!!!!DEMANDING!!!!!

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h11 = 1,norb1
     do h9 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h9,h11,p1,p2,p3,spin_g3inp) * itm_hhhh(h5,h6,h9,h11)
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
  end subroutine ccdt_l3p_6_perm
  !--------------------------------------------
end subroutine ccdt_l3p_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_6_1(sh4,sh5,sh9,sh11,i1)

!     i1 ( h4 h5 h9 h11 )_v + = -1 * v ( h4 h5 h9 h11 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh9,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h4,h5,h9,h11
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh4,sh5,sh9,sh11)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h9 = 1,norb1
  do h11 = 1,norb1
     i1(h4,h5,h9,h11) = i1(h4,h5,h9,h11) + fact * int2x(h4,h5,h9,h11,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_6_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_6_2(sh4,sh5,sh9,sh11,i1)

!     i1 ( h4 h5 h9 h11 )_vt + = -1/2 * Sum ( p7 p8 ) * t ( p7 p8 h9 h11 )_t * v ( h4 h5 p7 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh9,sh11
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h4,h5,h9,h11
  integer :: p7,p8,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh11)
     spin_int2x = tdcc_spin_int2x(sh4,sh5,sp7,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h9 = 1,norb1
     do h11 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h4,h5,h9,h11) = i1(h4,h5,h9,h11) + fact * &
             t2inp(p7,p8,h9,h11,spin_t2inp) * int2x(h4,h5,p7,p8,spin_int2x)
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
end subroutine ccdt_l3p_6_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_7(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = -1 * P( 9 ) * Sum ( h11 p7 ) * y ( h4 h5 h11 p1 p2 p7 )_y * i1 ( h6 p7 h11 p3 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_7_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_7_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
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
     call ccdt_l3p_7_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
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
     call ccdt_l3p_7_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_7_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
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
     call ccdt_l3p_7_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_7_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_7_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_7_perm(sh4,sh6,sh5,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3p_7_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h11,p7,sh11,sp7
  integer :: spin_g3inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sh11 = 1,2
  do sp7 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh11,sp1,sp2,sp7)
     spin_itm_hphp = tdcc_spin_int2x(sh6,sp7,sh11,sp3)
     if(spin_g3inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call ccdt_l3p_7_1(sh6,sp7,sh11,sp3,itm_hphp)
     call ccdt_l3p_7_2(sh6,sp7,sh11,sp3,itm_hphp)

!!!!!DEMANDING!!!!!

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h11 = 1,norb1
     do p7 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             g3inp(h4,h5,h11,p1,p2,p7,spin_g3inp) * itm_hphp(h6,p7,h11,p3)
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
  end subroutine ccdt_l3p_7_perm
  !--------------------------------------------
end subroutine ccdt_l3p_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_7_1(sh4,sp7,sh11,sp1,i1)

!     i1 ( h4 p7 h11 p1 )_v + = 1 * v ( h4 p7 h11 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sp7,sh11,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h4,p7,h11,p1
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh4,sp7,sh11,sp1)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h4 = 1,norb1
  do p7 = norb1+1,nact
  do h11 = 1,norb1
  do p1 = norb1+1,nact
     i1(h4,p7,h11,p1) = i1(h4,p7,h11,p1) + fact * int2x(h4,p7,h11,p1,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine ccdt_l3p_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_7_2(sh4,sp7,sh11,sp1,i1)

!     i1 ( h4 p7 h11 p1 )_vt + = 1 * Sum ( h10 p8 ) * t ( p7 p8 h10 h11 )_t * v ( h4 h10 p1 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sp7,sh11,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h4,p7,h11,p1
  integer :: h10,p8,sh10,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh10 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh10,sh11)
     spin_int2x = tdcc_spin_int2x(sh4,sh10,sp1,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h4 = 1,norb1
     do p7 = norb1+1,nact
     do h11 = 1,norb1
     do p1 = norb1+1,nact
     do h10 = 1,norb1
     do p8 = norb1+1,nact
        i1(h4,p7,h11,p1) = i1(h4,p7,h11,p1) + fact * &
             t2inp(p7,p8,h10,h11,spin_t2inp) * int2x(h4,h10,p1,p8,spin_int2x)
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
end subroutine ccdt_l3p_7_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_8(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_yv + = 1/2 * P( 3 ) * Sum ( p7 p8 ) * y ( h4 h5 h6 p1 p7 p8 )_y * v ( p7 p8 p2 p3 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_8_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_8_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
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
     call ccdt_l3p_8_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
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
  subroutine ccdt_l3p_8_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: p7,p8,sp7,sp8
  integer :: spin_g3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp7,sp8)
     spin_int2x = tdcc_spin_int2x(sp7,sp8,sp2,sp3)
     if(spin_g3inp * spin_int2x == 0) cycle

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
             g3inp(h4,h5,h6,p1,p7,p8,spin_g3inp) * int2x(p7,p8,p2,p3,spin_int2x)
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
  end subroutine ccdt_l3p_8_perm
  !--------------------------------------------
end subroutine ccdt_l3p_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_9(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = -1/2 * P( 9 ) * Sum ( p8 ) * i1 ( h4 p8 p1 p2 )_yt * v ( h5 h6 p3 p8 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_9_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_9_perm(sh5,sh4,sh6,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p2,p3)
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
     call ccdt_l3p_9_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
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
     call ccdt_l3p_9_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_9_perm(sh4,sh5,sh6,sp1,sp3,sp2,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h5,h6,p1,p3,p2)
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
     call ccdt_l3p_9_perm(sh5,sh4,sh6,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p3,p2,p1)
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
     call ccdt_l3p_9_perm(sh5,sh4,sh6,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h5,h4,h6,p1,p3,p2)
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
     call ccdt_l3p_9_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_9_perm(sh6,sh5,sh4,sp1,sp3,sp2,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3p_9_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: p8,sp8
  integer :: spin_itm_hppp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp8 = 1,2
     spin_itm_hppp = tdcc_spin_dummy2(sh4,sp8,sp1,sp2)
     spin_int2x = tdcc_spin_int2x(sh5,sh6,sp3,sp8)
     if(spin_itm_hppp * spin_int2x == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_9_1(sh4,sp8,sp1,sp2,itm_hppp)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             itm_hppp(h4,p8,p1,p2) * int2x(h5,h6,p3,p8,spin_int2x)
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
  deallocate(itm_hppp)
  end subroutine ccdt_l3p_9_perm
  !--------------------------------------------
end subroutine ccdt_l3p_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_9_1(sh4,sp8,sp1,sp2,i1)

!     i1 ( h4 p8 p1 p2 )_yt + = 1 * Sum ( h9 h10 p7 ) * t ( p7 p8 h9 h10 )_t * y ( h4 h9 h10 p1 p2 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sp8,sp1,sp2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,p8,p1,p2
  integer :: h9,h10,p7,sh9,sh10,sp7
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh9 = 1,2
  do sh10 = 1,2
  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh9,sh10,sp1,sp2,sp7)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h4 = 1,norb1
     do p8 = norb1+1,nact
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
     do p7 = norb1+1,nact
        i1(h4,p8,p1,p2) = i1(h4,p8,p1,p2) + fact * &
             t2inp(p7,p8,h9,h10,spin_t2inp) * g3inp(h4,h9,h10,p1,p2,p7,spin_g3inp)
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
end subroutine ccdt_l3p_9_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_10(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = -1/2 * P( 9 ) * Sum ( h10 ) * i1 ( h4 h5 h10 p1 )_yt * v ( h6 h10 p2 p3 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_10_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_10_perm(sh6,sh5,sh4,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p1,p2,p3)
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
     call ccdt_l3p_10_perm(sh4,sh6,sh5,sp1,sp2,sp3,i0_perm)
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
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p1,p2,p3)
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
     call ccdt_l3p_10_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
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
     call ccdt_l3p_10_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
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

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_10_perm(sh6,sh5,sh4,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h6,h5,h4,p2,p1,p3)
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
     call ccdt_l3p_10_perm(sh6,sh5,sh4,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  if(.not. (sh4 * sh5 * sh6 * sp1 * sp2 * sp3 == 1)) then
     i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call ccdt_l3p_10_perm(sh4,sh6,sh5,sp2,sp1,sp3,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
  do h4 = 1,norb1
  do h5 = 1,norb1
  do h6 = 1,norb1
  do p1 = norb1+1,nact
  do p2 = norb1+1,nact
  do p3 = norb1+1,nact
     i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact_p * i0_perm(h4,h6,h5,p2,p1,p3)
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
     call ccdt_l3p_10_perm(sh4,sh6,sh5,sp3,sp2,sp1,i0_perm)
  end if
  fact_p = +1.0d+0
  !$omp parallel default(shared)
  !$omp do collapse(6)
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
  !$omp end do
  !$omp end parallel

  deallocate(i0_perm)

  contains
  !--------------------------------------------
  subroutine ccdt_l3p_10_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h10,sh10
  integer :: spin_itm_hhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh10 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh4,sh5,sh10,sp1)
     spin_int2x = tdcc_spin_int2x(sh6,sh10,sp2,sp3)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccdt_l3p_10_1(sh4,sh5,sh10,sp1,itm_hhhp)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do p3 = norb1+1,nact
     do h10 = 1,norb1
        i0(h4,h5,h6,p1,p2,p3) = i0(h4,h5,h6,p1,p2,p3) + fact * &
             itm_hhhp(h4,h5,h10,p1) * int2x(h6,h10,p2,p3,spin_int2x)
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
  deallocate(itm_hhhp)
  end subroutine ccdt_l3p_10_perm
  !--------------------------------------------
end subroutine ccdt_l3p_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_10_1(sh4,sh5,sh10,sp1,i1)

!     i1 ( h4 h5 h10 p1 )_yt + = 1 * Sum ( h9 p7 p8 ) * t ( p7 p8 h9 h10 )_t * y ( h4 h5 h9 p1 p7 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh10,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h10,p1
  integer :: h9,p7,p8,sh9,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh9 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh9,sp1,sp7,sp8)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h4 = 1,norb1
     do h5 = 1,norb1
     do h10 = 1,norb1
     do p1 = norb1+1,nact
     do h9 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h4,h5,h10,p1) = i1(h4,h5,h10,p1) + fact * &
             t2inp(p7,p8,h9,h10,spin_t2inp) * g3inp(h4,h5,h9,p1,p7,p8,spin_g3inp)
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
end subroutine ccdt_l3p_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_11(sh4,sh5,sh6,sp1,sp2,sp3,i0)

! i0 ( h4 h5 h6 p1 p2 p3 )_ytv + = 1/4 * P( 3 ) * Sum ( h9 h10 ) * i1 ( h4 h5 h6 h9 h10 p1 )_yt * v ( h9 h10 p2 p3 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h4,h5,h6,p1,p2,p3
  complex(kind(0d0)) :: fact_p
  complex(kind(0d0)),allocatable :: i0_perm(:,:,:,:,:,:)

  allocate(i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))

  i0_perm(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
  call ccdt_l3p_11_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0_perm)
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
     call ccdt_l3p_11_perm(sh4,sh5,sh6,sp2,sp1,sp3,i0_perm)
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
     call ccdt_l3p_11_perm(sh4,sh5,sh6,sp3,sp2,sp1,i0_perm)
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
  subroutine ccdt_l3p_11_perm(sh4,sh5,sh6,sp1,sp2,sp3,i0)

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sp1,sp2,sp3
  complex(kind(0d0)),intent(inout) :: &
       i0(1:norb1,1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)

  integer :: h4,h5,h6,p1,p2,p3
  integer :: h9,h10,sh9,sh10
  integer :: spin_itm_hhhhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy3
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhhhp(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sh10 = 1,2
     spin_itm_hhhhhp = tdcc_spin_dummy3(sh4,sh5,sh6,sh9,sh10,sp1)
     spin_int2x = tdcc_spin_int2x(sh9,sh10,sp2,sp3)
     if(spin_itm_hhhhhp * spin_int2x == 0) cycle

     itm_hhhhhp(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call ccdt_l3p_11_1(sh4,sh5,sh6,sh9,sh10,sp1,itm_hhhhhp)

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
  end subroutine ccdt_l3p_11_perm
  !--------------------------------------------
end subroutine ccdt_l3p_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_l3p_11_1(sh4,sh5,sh6,sh9,sh10,sp1,i1)

!     i1 ( h4 h5 h6 h9 h10 p1 )_yt + = 1 * Sum ( p7 p8 ) * t ( p7 p8 h9 h10 )_t * y ( h4 h5 h6 p1 p7 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh4,sh5,sh6,sh9,sh10,sp1
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,1:norb1,1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h4,h5,h6,h9,h10,p1
  integer :: p7,p8,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh4,sh5,sh6,sp1,sp7,sp8)
     if(spin_t2inp * spin_g3inp == 0) cycle

!!!!!DEMANDING!!!!!

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
end subroutine ccdt_l3p_11_1
!##########################################################
!##########################################################
!##########################################################
