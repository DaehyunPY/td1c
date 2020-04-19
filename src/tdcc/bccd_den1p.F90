!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_den1p_1(sh1,sp2,i0)

!  i0 ( h1 p2 )_y + = +1 * y ( h1 p2 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp

  implicit none
  integer,intent(in) :: sh1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: h1,p2
  integer :: sdum
  integer :: spin_g1inp_1
  integer,external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  spin_g1inp_1 = tdcc_spin_g1inp(sh1,sp2)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do h1 = 1,norb1
  do p2 = norb1+1,nact
     i0(h1,p2) = i0(h1,p2) + fact * g1inp(h1,p2,spin_g1inp_1)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccd_den1p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_den1p_2(sp1,sp2,i0)

!  i0 ( p1 p2 )_yt + = +1/2 * Sum ( h3 h4 p5 ) * y ( h3 h4 p5 p2 )_y * t ( p5 p1 h3 h4 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: p1,p2
  integer :: h3,h4,p5,sh3,sh4,sp5
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh3,sh4,sp5,sp2)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp1,sh3,sh4)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(p1,p2) = i0(p1,p2) + fact * g2inp(h3,h4,p5,p2,spin_g2inp_1) * t2inp(p5,p1,h3,h4,spin_t2inp_2)
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
end subroutine bccd_den1p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_den1p_3(sh1,sh2,i0)

!  i0 ( h1 h2 )_yt + = -1/2 * Sum ( h3 p4 p5 ) * y ( h3 h1 p4 p5 )_y * t ( p4 p5 h3 h2 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: h1,h2
  integer :: h3,p4,p5,sh3,sp4,sp5
  integer :: spin_g2inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh3 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh3,sh1,sp4,sp5)
     spin_t2inp_2 = tdcc_spin_t2inp(sp4,sp5,sh3,sh2)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1,norb1
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i0(h1,h2) = i0(h1,h2) + fact * g2inp(h3,h1,p4,p5,spin_g2inp_1) * t2inp(p4,p5,h3,h2,spin_t2inp_2)
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
end subroutine bccd_den1p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_den1p_4(sp1,sh2,i0)

!  i0 ( p1 h2 )_yt + = +1 * Sum ( h3 p4 ) * y ( h3 p4 )_y * t ( p4 p1 h3 h2 )_t 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t2inp

  implicit none
  integer,intent(in) :: sp1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer :: p1,h2
  integer :: h3,p4,sh3,sp4
  integer :: spin_g1inp_1
  integer :: spin_t2inp_2
  integer,external :: tdcc_spin_g1inp
  integer,external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh3 = 1,2
  do sp4 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh3,sp4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp4,sp1,sh3,sh2)
     if(spin_g1inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p1 = norb1+1,nact
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p4 = norb1+1,nact
        i0(p1,h2) = i0(p1,h2) + fact * g1inp(h3,p4,spin_g1inp_1) * t2inp(p4,p1,h3,h2,spin_t2inp_2)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccd_den1p_4
!##########################################################
!##########################################################
!##########################################################
