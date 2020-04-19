!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_1(sh2,sp1,i0)

! i0 ( h2 p1 )_f + = 1 * f ( h2 p1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh2,sp1)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do h2 = 1,norb1
  do p1 = norb1+1,nact
     i0(h2,p1) = i0(h2,p1) + fact * fock(h2,p1,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccdt_l1p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_2(sh2,sp1,i0)

! i0 ( h2 p1 )_yf + = -1 * Sum ( h7 ) * y ( h7 p1 )_y * i1 ( h2 h7 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h7,sh7
  integer :: spin_g1inp
  integer :: spin_itm_hh
  integer,external :: tdcc_spin_g1inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh7 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh7,sp1)
     spin_itm_hh = tdcc_spin_fock(sh2,sh7)
     if(spin_g1inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call bccdt_l1p_2_1(sh2,sh7,itm_hh)
     call bccdt_l1p_2_2(sh2,sh7,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * g1inp(h7,p1,spin_g1inp) * itm_hh(h2,h7)
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hh)
end subroutine bccdt_l1p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_2_1(sh2,sh7,i1)

!     i1 ( h2 h7 )_f + = 1 * f ( h2 h7 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh2,sh7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h2,h7
  integer :: sdum
  integer :: spin_fock
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh2,sh7)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do h2 = 1,norb1
  do h7 = 1,norb1
     i1(h2,h7) = i1(h2,h7) + fact * fock(h2,h7,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccdt_l1p_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_2_2(sh2,sh7,i1)

!     i1 ( h2 h7 )_vt + = -1/2 * Sum ( h6 p3 p4 ) * t ( p3 p4 h6 h7 )_t * v ( h2 h6 p3 p4 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh2,sh7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h2,h7
  integer :: h6,p3,p4,sh6,sp3,sp4
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh6 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh6,sh7)
     spin_int2x = tdcc_spin_int2x(sh2,sh6,sp3,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do h7 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,h7) = i1(h2,h7) + fact * &
             t2inp(p3,p4,h6,h7,spin_t2inp) * int2x(h2,h6,p3,p4,spin_int2x)
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
end subroutine bccdt_l1p_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_3(sh2,sp1,i0)

! i0 ( h2 p1 )_yf + = 1 * Sum ( p7 ) * y ( h2 p7 )_y * f ( p7 p1 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: p7,sp7
  integer :: spin_g1inp
  integer :: spin_fock
  integer,external :: tdcc_spin_g1inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh2,sp7)
     spin_fock = tdcc_spin_fock(sp7,sp1)
     if(spin_g1inp * spin_fock == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * &
             g1inp(h2,p7,spin_g1inp) * fock(p7,p1,spin_fock)
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
end subroutine bccdt_l1p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_4(sh2,sp1,i0)

! i0 ( h2 p1 )_yv + = -1 * Sum ( h4 p3 ) * y ( h4 p3 )_y * v ( h2 p3 h4 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h4,p3,sh4,sp3
  integer :: spin_g1inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh4 = 1,2
  do sp3 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh4,sp3)
     spin_int2x = tdcc_spin_int2x(sh2,sp3,sh4,sp1)
     if(spin_g1inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h4 = 1,norb1
     do p3 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * &
             g1inp(h4,p3,spin_g1inp) * int2x(h2,p3,h4,p1,spin_int2x)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccdt_l1p_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_5(sh2,sp1,i0)

! i0 ( h2 p1 )_yv + = 1/2 * Sum ( h12 h9 p11 ) * y ( h9 h12 p1 p11 )_y * i1 ( h2 p11 h9 h12 )_v 5

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h12,h9,p11,sh12,sh9,sp11
  integer :: spin_g2inp
  integer :: spin_itm_hphh
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh12 = 1,2
  do sh9 = 1,2
  do sp11 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh9,sh12,sp1,sp11)
     spin_itm_hphh = tdcc_spin_int2x(sh2,sp11,sh9,sh12)
     if(spin_g2inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call bccdt_l1p_5_1(sh2,sp11,sh9,sh12,itm_hphh)
     call bccdt_l1p_5_2(sh2,sp11,sh9,sh12,itm_hphh)
     call bccdt_l1p_5_3(sh2,sp11,sh9,sh12,itm_hphh)
     call bccdt_l1p_5_4(sh2,sp11,sh9,sh12,itm_hphh)
     call bccdt_l1p_5_5(sh2,sp11,sh9,sh12,itm_hphh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h12 = 1,norb1
     do h9 = 1,norb1
     do p11 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * &
             g2inp(h9,h12,p1,p11,spin_g2inp) * itm_hphh(h2,p11,h9,h12)
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
end subroutine bccdt_l1p_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_5_1(sh2,sp11,sh9,sh12,i1)

!     i1 ( h2 p11 h9 h12 )_v + = -1 * v ( h2 p11 h9 h12 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp11,sh9,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h2,p11,h9,h12
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh2,sp11,sh9,sh12)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h2 = 1,norb1
  do p11 = norb1+1,nact
  do h9 = 1,norb1
  do h12 = 1,norb1
     i1(h2,p11,h9,h12) = i1(h2,p11,h9,h12) + fact * &
          int2x(h2,p11,h9,h12,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccdt_l1p_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_5_2(sh2,sp11,sh9,sh12,i1)

!     i1 ( h2 p11 h9 h12 )_ft + = -1 * Sum ( p5 ) * t ( p5 p11 h9 h12 )_t * f ( h2 p5 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh2,sp11,sh9,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h2,p11,h9,h12
  integer :: p5,sp5
  integer :: spin_t2inp
  integer :: spin_fock
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp11,sh9,sh12)
     spin_fock = tdcc_spin_fock(sh2,sp5)
     if(spin_t2inp * spin_fock == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do h9 = 1,norb1
     do h12 = 1,norb1
     do p5 = norb1+1,nact
        i1(h2,p11,h9,h12) = i1(h2,p11,h9,h12) + fact * &
             t2inp(p5,p11,h9,h12,spin_t2inp) * fock(h2,p5,spin_fock)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
end subroutine bccdt_l1p_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_5_3(sh2,sp11,sh9,sh12,i1)

!     i1 ( h2 p11 h9 h12 )_vt + = -2 * Sum ( h6 p4 ) * t ( p4 p11 h6 h12 )_t * v ( h2 h6 h9 p4 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh2,sp11,sh9,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h2,p11,h9,h12
  integer :: h6,p4,sh6,sp4
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  do sh6 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp11,sh6,sh12)
     spin_int2x = tdcc_spin_int2x(sh2,sh6,sh9,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do h9 = 1,norb1
     do h12 = 1,norb1
     do h6 = 1,norb1
     do p4 = norb1+1,nact
        i1(h2,p11,h9,h12) = i1(h2,p11,h9,h12) + fact * &
             t2inp(p4,p11,h6,h12,spin_t2inp) * int2x(h2,h6,h9,p4,spin_int2x)
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
end subroutine bccdt_l1p_5_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_5_4(sh2,sp11,sh9,sh12,i1)

!     i1 ( h2 p11 h9 h12 )_vt + = -1/2 * Sum ( p3 p4 ) * t ( p3 p4 h9 h12 )_t * v ( h2 p11 p3 p4 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh2,sp11,sh9,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h2,p11,h9,h12
  integer :: p3,p4,sp3,sp4
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh9,sh12)
     spin_int2x = tdcc_spin_int2x(sh2,sp11,sp3,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do h9 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,p11,h9,h12) = i1(h2,p11,h9,h12) + fact * &
             t2inp(p3,p4,h9,h12,spin_t2inp) * int2x(h2,p11,p3,p4,spin_int2x)
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
end subroutine bccdt_l1p_5_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_5_5(sh2,sp11,sh9,sh12,i1)

!     i1 ( h2 p11 h9 h12 )_vt + = 1/2 * Sum ( h8 p4 p5 ) * t ( p4 p5 p11 h8 h9 h12 )_t * v ( h2 h8 p4 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh2,sp11,sh9,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h2,p11,h9,h12
  integer :: h8,p4,p5,sh8,sp4,sp5
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp11,sh8,sh9,sh12)
     spin_int2x = tdcc_spin_int2x(sh2,sh8,sp4,sp5)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do h9 = 1,norb1
     do h12 = 1,norb1
     do h8 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h2,p11,h9,h12) = i1(h2,p11,h9,h12) + fact * &
             t3inp(p4,p5,p11,h8,h9,h12,spin_t3inp) * int2x(h2,h8,p4,p5,spin_int2x)
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
end subroutine bccdt_l1p_5_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_6(sh2,sp1,i0)

! i0 ( h2 p1 )_yv + = -1/2 * Sum ( h7 p8 p5 ) * y ( h2 h7 p5 p8 )_y * v ( p5 p8 h7 p1 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h7,p8,p5,sh7,sp8,sp5
  integer :: spin_g2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_g2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh7 = 1,2
  do sp8 = 1,2
  do sp5 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh2,sh7,sp5,sp8)
     spin_int2x = tdcc_spin_int2x(sp5,sp8,sh7,sp1)
     if(spin_g2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do p8 = norb1+1,nact
     do p5 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * &
             g2inp(h2,h7,p5,p8,spin_g2inp) * int2x(p5,p8,h7,p1,spin_int2x)
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
end subroutine bccdt_l1p_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_7(sh2,sp1,i0)

! i0 ( h2 p1 )_vt + = 1 * Sum ( h13 p11 ) * i1 ( p11 h13 )_t * v ( h2 h13 p1 p11 )_v 3

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h13,p11,sh13,sp11
  integer :: spin_itm_ph
  integer :: spin_int2x
  integer,external :: tdcc_spin_t1inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_ph(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_ph((norb1+1):nact,1:norb1))
  do sh13 = 1,2
  do sp11 = 1,2
     spin_itm_ph = tdcc_spin_t1inp(sp11,sh13)
     spin_int2x = tdcc_spin_int2x(sh2,sh13,sp1,sp11)
     if(spin_itm_ph * spin_int2x == 0) cycle

     itm_ph((norb1+1):nact,1:norb1) = czero
     call bccdt_l1p_7_1(sp11,sh13,itm_ph)
     call bccdt_l1p_7_2(sp11,sh13,itm_ph)
     call bccdt_l1p_7_3(sp11,sh13,itm_ph)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h13 = 1,norb1
     do p11 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_ph(p11,h13) * int2x(h2,h13,p1,p11,spin_int2x)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_ph)
end subroutine bccdt_l1p_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_7_1(sp11,sh13,i1)

!     i1 ( p11 h13 )_yt + = 1 * Sum ( h5 p3 ) * t ( p3 p11 h5 h13 )_t * y ( h5 p3 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g1inp

  implicit none
  integer,intent(in) :: sp11,sh13
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer :: p11,h13
  integer :: h5,p3,sh5,sp3
  integer :: spin_t2inp
  integer :: spin_g1inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp11,sh5,sh13)
     spin_g1inp = tdcc_spin_g1inp(sh5,sp3)
     if(spin_t2inp * spin_g1inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p11 = norb1+1,nact
     do h13 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
        i1(p11,h13) = i1(p11,h13) + fact * &
             t2inp(p3,p11,h5,h13,spin_t2inp) * g1inp(h5,p3,spin_g1inp)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccdt_l1p_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_7_2(sp11,sh13,i1)

!     i1 ( p11 h13 )_yt + = 1/4 * Sum ( h6 h7 p3 p4 ) * t ( p3 p4 p11 h6 h7 h13 )_t * y ( h6 h7 p3 p4 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sp11,sh13
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer :: p11,h13
  integer :: h6,h7,p3,p4,sh6,sh7,sp3,sp4
  integer :: spin_t3inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp4,sp11,sh6,sh7,sh13)
     spin_g2inp = tdcc_spin_g2inp(sh6,sh7,sp3,sp4)
     if(spin_t3inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p11 = norb1+1,nact
     do h13 = 1,norb1
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(p11,h13) = i1(p11,h13) + fact * &
             t3inp(p3,p4,p11,h6,h7,h13,spin_t3inp) * g2inp(h6,h7,p3,p4,spin_g2inp)
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
end subroutine bccdt_l1p_7_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_7_3(sp11,sh13,i1)

!     i1 ( p11 h13 )_ytt + = 1/2 * Sum ( h5 h6 p3 ) * t ( p3 p11 h5 h6 )_t * i2 ( h5 h6 h13 p3 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp11,sh13
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer :: p11,h13
  integer :: h5,h6,p3,sh5,sh6,sp3
  integer :: spin_t2inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh5 = 1,2
  do sh6 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp11,sh5,sh6)
     spin_itm_hhhp = tdcc_spin_dummy2(sh5,sh6,sh13,sp3)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call bccdt_l1p_7_3_1(sh5,sh6,sh13,sp3,itm_hhhp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p11 = norb1+1,nact
     do h13 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
        i1(p11,h13) = i1(p11,h13) + fact * &
             t2inp(p3,p11,h5,h6,spin_t2inp) * itm_hhhp(h5,h6,h13,p3)
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
  deallocate(itm_hhhp)
end subroutine bccdt_l1p_7_3
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_7_3_1(sh5,sh6,sh13,sp3,i2)

!         i2 ( h5 h6 h13 p3 )_yt + = 1/2 * Sum ( h9 p7 p8 ) * t ( p7 p8 h9 h13 )_t * y ( h5 h6 h9 p3 p7 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh5,sh6,sh13,sp3
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h5,h6,h13,p3
  integer :: h9,p7,p8,sh9,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh9 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh9,sh13)
     spin_g3inp = tdcc_spin_g3inp(sh5,sh6,sh9,sp3,sp7,sp8)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h5 = 1,norb1
     do h6 = 1,norb1
     do h13 = 1,norb1
     do p3 = norb1+1,nact
     do h9 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i2(h5,h6,h13,p3) = i2(h5,h6,h13,p3) + fact * &
             t2inp(p7,p8,h9,h13,spin_t2inp) * g3inp(h5,h6,h9,p3,p7,p8,spin_g3inp)
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
end subroutine bccdt_l1p_7_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_8(sh2,sp1,i0)

! i0 ( h2 p1 )_ytf + = -1 * Sum ( h3 ) * i1 ( h2 h3 )_yt * f ( h3 p1 )_f 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h3,sh3
  integer :: spin_itm_hh
  integer :: spin_fock
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh3 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh2,sh3)
     spin_fock = tdcc_spin_fock(sh3,sp1)
     if(spin_itm_hh * spin_fock == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call bccdt_l1p_8_1(sh2,sh3,itm_hh)
     call bccdt_l1p_8_2(sh2,sh3,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h3 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_hh(h2,h3) * fock(h3,p1,spin_fock)
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hh)
end subroutine bccdt_l1p_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_8_1(sh2,sh3,i1)

!     i1 ( h2 h3 )_yt + = 1/2 * Sum ( h6 p4 p5 ) * t ( p4 p5 h3 h6 )_t * y ( h2 h6 p4 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh2,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h2,h3
  integer :: h6,p4,p5,sh6,sp4,sp5
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh6 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp5,sh3,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh6,sp4,sp5)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h6 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h2,h3) = i1(h2,h3) + fact * &
             t2inp(p4,p5,h3,h6,spin_t2inp) * g2inp(h2,h6,p4,p5,spin_g2inp)
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
end subroutine bccdt_l1p_8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_8_2(sh2,sh3,i1)

!     i1 ( h2 h3 )_yt + = 1/12 * Sum ( h7 h8 p4 p5 p6 ) * t ( p4 p5 p6 h3 h7 h8 )_t * y ( h2 h7 h8 p4 p5 p6 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh2,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h2,h3
  integer :: h7,h8,p4,p5,p6,sh7,sh8,sp4,sp5,sp6
  integer :: spin_t3inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 12.0d+0 * runit

  do sh7 = 1,2
  do sh8 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp6,sh3,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh2,sh7,sh8,sp4,sp5,sp6)
     if(spin_t3inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do h3 = 1,norb1
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h2,h3) = i1(h2,h3) + fact * &
             t3inp(p4,p5,p6,h3,h7,h8,spin_t3inp) * g3inp(h2,h7,h8,p4,p5,p6,spin_g3inp)
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
end subroutine bccdt_l1p_8_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_9(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( h10 h8 ) * i1 ( h10 h8 )_yt * v ( h2 h8 h10 p1 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h10,h8,sh10,sh8
  integer :: spin_itm_hh
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh10 = 1,2
  do sh8 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh10,sh8)
     spin_int2x = tdcc_spin_int2x(sh2,sh8,sh10,sp1)
     if(spin_itm_hh * spin_int2x == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call bccdt_l1p_9_1(sh10,sh8,itm_hh)
     call bccdt_l1p_9_2(sh10,sh8,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h10 = 1,norb1
     do h8 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_hh(h10,h8) * int2x(h2,h8,h10,p1,spin_int2x)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hh)
end subroutine bccdt_l1p_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_9_1(sh10,sh8,i1)

!     i1 ( h10 h8 )_yt + = 1/2 * Sum ( h5 p3 p4 ) * t ( p3 p4 h5 h8 )_t * y ( h5 h10 p3 p4 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh10,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h10,h8
  integer :: h5,p3,p4,sh5,sp3,sp4
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh5,sh10,sp3,sp4)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h10 = 1,norb1
     do h8 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h10,h8) = i1(h10,h8) + fact * &
             t2inp(p3,p4,h5,h8,spin_t2inp) * g2inp(h5,h10,p3,p4,spin_g2inp)
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
end subroutine bccdt_l1p_9_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_9_2(sh10,sh8,i1)

!     i1 ( h10 h8 )_yt + = 1/12 * Sum ( h6 h7 p3 p4 p5 ) * t ( p3 p4 p5 h6 h7 h8 )_t * y ( h6 h7 h10 p3 p4 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh10,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer :: h10,h8
  integer :: h6,h7,p3,p4,p5,sh6,sh7,sp3,sp4,sp5
  integer :: spin_t3inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 12.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp4,sp5,sh6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh6,sh7,sh10,sp3,sp4,sp5)
     if(spin_t3inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h10 = 1,norb1
     do h8 = 1,norb1
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h10,h8) = i1(h10,h8) + fact * &
             t3inp(p3,p4,p5,h6,h7,h8,spin_t3inp) * g3inp(h6,h7,h10,p3,p4,p5,spin_g3inp)
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
end subroutine bccdt_l1p_9_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_10(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( p9 p10 ) * i1 ( p9 p10 )_yt * v ( h2 p10 p1 p9 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: p9,p10,sp9,sp10
  integer :: spin_itm_pp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy1
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp9 = 1,2
  do sp10 = 1,2
     spin_itm_pp = tdcc_spin_dummy1(sp9,sp10)
     spin_int2x = tdcc_spin_int2x(sh2,sp10,sp1,sp9)
     if(spin_itm_pp * spin_int2x == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call bccdt_l1p_10_1(sp9,sp10,itm_pp)
     call bccdt_l1p_10_2(sp9,sp10,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_pp(p9,p10) * int2x(h2,p10,p1,p9,spin_int2x)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_pp)
end subroutine bccdt_l1p_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_10_1(sp9,sp10,i1)

!     i1 ( p9 p10 )_yt + = 1/2 * Sum ( h5 h6 p3 ) * t ( p3 p9 h5 h6 )_t * y ( h5 h6 p3 p10 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sp9,sp10
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p9,p10
  integer :: h5,h6,p3,sh5,sh6,sp3
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp9,sh5,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh5,sh6,sp3,sp10)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
        i1(p9,p10) = i1(p9,p10) + fact * &
             t2inp(p3,p9,h5,h6,spin_t2inp) * g2inp(h5,h6,p3,p10,spin_g2inp)
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
end subroutine bccdt_l1p_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_10_2(sp9,sp10,i1)

!     i1 ( p9 p10 )_yt + = 1/12 * Sum ( h6 h7 h8 p3 p4 ) * t ( p3 p4 p9 h6 h7 h8 )_t * y ( h6 h7 h8 p3 p4 p10 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sp9,sp10
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer :: p9,p10
  integer :: h6,h7,h8,p3,p4,sh6,sh7,sh8,sp3,sp4
  integer :: spin_t3inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 12.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
  do sh8 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp4,sp9,sh6,sh7,sh8)
     spin_g3inp = tdcc_spin_g3inp(sh6,sh7,sh8,sp3,sp4,sp10)
     if(spin_t3inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p9 = norb1+1,nact
     do p10 = norb1+1,nact
     do h6 = 1,norb1
     do h7 = 1,norb1
     do h8 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(p9,p10) = i1(p9,p10) + fact * &
             t3inp(p3,p4,p9,h6,h7,h8,spin_t3inp) * g3inp(h6,h7,h8,p3,p4,p10,spin_g3inp)
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
end subroutine bccdt_l1p_10_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_11(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( h7 h9 p10 ) * i1 ( h2 h7 h9 p10 )_yt * v ( h9 p10 h7 p1 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h7,h9,p10,sh7,sh9,sp10
  integer :: spin_itm_hhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh7 = 1,2
  do sh9 = 1,2
  do sp10 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh2,sh7,sh9,sp10)
     spin_int2x = tdcc_spin_int2x(sh9,sp10,sh7,sp1)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call bccdt_l1p_11_1(sh2,sh7,sh9,sp10,itm_hhhp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h9 = 1,norb1
     do p10 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_hhhp(h2,h7,h9,p10) * int2x(h9,p10,h7,p1,spin_int2x)
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
  deallocate(itm_hhhp)
end subroutine bccdt_l1p_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_11_1(sh2,sh7,sh9,sp10,i1)

!     i1 ( h2 h7 h9 p10 )_yt + = -1/2 * Sum ( h5 p3 p4 ) * t ( p3 p4 h5 h9 )_t * y ( h2 h5 h7 p3 p4 p10 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh2,sh7,sh9,sp10
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h2,h7,h9,p10
  integer :: h5,p3,p4,sh5,sp3,sp4
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh9)
     spin_g3inp = tdcc_spin_g3inp(sh2,sh5,sh7,sp3,sp4,sp10)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do h7 = 1,norb1
     do h9 = 1,norb1
     do p10 = norb1+1,nact
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,h7,h9,p10) = i1(h2,h7,h9,p10) + fact * &
             t2inp(p3,p4,h5,h9,spin_t2inp) * g3inp(h2,h5,h7,p3,p4,p10,spin_g3inp)
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
end subroutine bccdt_l1p_11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_12(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1/2 * Sum ( p14 h13 h11 ) * i1 ( h2 p14 h11 h13 )_yt * v ( h11 h13 p1 p14 )_v 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: p14,h13,h11,sp14,sh13,sh11
  integer :: spin_itm_hphh
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sp14 = 1,2
  do sh13 = 1,2
  do sh11 = 1,2
     spin_itm_hphh = tdcc_spin_dummy2(sh2,sp14,sh11,sh13)
     spin_int2x = tdcc_spin_int2x(sh11,sh13,sp1,sp14)
     if(spin_itm_hphh * spin_int2x == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call bccdt_l1p_12_1(sh2,sp14,sh11,sh13,itm_hphh)
     call bccdt_l1p_12_2(sh2,sp14,sh11,sh13,itm_hphh)
     call bccdt_l1p_12_3(sh2,sp14,sh11,sh13,itm_hphh)
     call bccdt_l1p_12_4(sh2,sp14,sh11,sh13,itm_hphh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p14 = norb1+1,nact
     do h13 = 1,norb1
     do h11 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_hphh(h2,p14,h11,h13) * int2x(h11,h13,p1,p14,spin_int2x)
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
end subroutine bccdt_l1p_12
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_12_1(sh2,sp14,sh11,sh13,i1)

!     i1 ( h2 p14 h11 h13 )_yt + = -1 * Sum ( p3 ) * t ( p3 p14 h11 h13 )_t * y ( h2 p3 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g1inp

  implicit none
  integer,intent(in) :: sh2,sp14,sh11,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h2,p14,h11,h13
  integer :: p3,sp3
  integer :: spin_t2inp
  integer :: spin_g1inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp14,sh11,sh13)
     spin_g1inp = tdcc_spin_g1inp(sh2,sp3)
     if(spin_t2inp * spin_g1inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p14 = norb1+1,nact
     do h11 = 1,norb1
     do h13 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,p14,h11,h13) = i1(h2,p14,h11,h13) + fact * &
             t2inp(p3,p14,h11,h13,spin_t2inp) * g1inp(h2,p3,spin_g1inp)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
end subroutine bccdt_l1p_12_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_12_2(sh2,sp14,sh11,sh13,i1)

!     i1 ( h2 p14 h11 h13 )_yt + = 1/2 * Sum ( h6 p3 p4 ) * t ( p3 p4 p14 h6 h11 h13 )_t * y ( h2 h6 p3 p4 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh2,sp14,sh11,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h2,p14,h11,h13
  integer :: h6,p3,p4,sh6,sp3,sp4
  integer :: spin_t3inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh6 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp4,sp14,sh6,sh11,sh13)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh6,sp3,sp4)
     if(spin_t3inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p14 = norb1+1,nact
     do h11 = 1,norb1
     do h13 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,p14,h11,h13) = i1(h2,p14,h11,h13) + fact * &
             t3inp(p3,p4,p14,h6,h11,h13,spin_t3inp) * g2inp(h2,h6,p3,p4,spin_g2inp)
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
end subroutine bccdt_l1p_12_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_12_3(sh2,sp14,sh11,sh13,i1)

!     i1 ( h2 p14 h11 h13 )_ytt + = -1/4 * Sum ( p7 p8 ) * t ( p7 p8 h11 h13 )_t * i2 ( h2 p14 p7 p8 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh2,sp14,sh11,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h2,p14,h11,h13
  integer :: p7,p8,sp7,sp8
  integer :: spin_t2inp
  integer :: spin_itm_hppp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp8,sh11,sh13)
     spin_itm_hppp = tdcc_spin_dummy2(sh2,sp14,sp7,sp8)
     if(spin_t2inp * spin_itm_hppp == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call bccdt_l1p_12_3_1(sh2,sp14,sp7,sp8,itm_hppp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p14 = norb1+1,nact
     do h11 = 1,norb1
     do h13 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i1(h2,p14,h11,h13) = i1(h2,p14,h11,h13) + fact * &
             t2inp(p7,p8,h11,h13,spin_t2inp) * itm_hppp(h2,p14,p7,p8)
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
end subroutine bccdt_l1p_12_3
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_12_3_1(sh2,sp14,sp7,sp8,i2)

!         i2 ( h2 p14 p7 p8 )_yt + = 1 * Sum ( h5 h6 p3 ) * t ( p3 p14 h5 h6 )_t * y ( h2 h5 h6 p3 p7 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh2,sp14,sp7,sp8
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h2,p14,p7,p8
  integer :: h5,h6,p3,sh5,sh6,sp3
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp14,sh5,sh6)
     spin_g3inp = tdcc_spin_g3inp(sh2,sh5,sh6,sp3,sp7,sp8)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p14 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
        i2(h2,p14,p7,p8) = i2(h2,p14,p7,p8) + fact * &
             t2inp(p3,p14,h5,h6,spin_t2inp) * g3inp(h2,h5,h6,p3,p7,p8,spin_g3inp)
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
end subroutine bccdt_l1p_12_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_12_4(sh2,sp14,sh11,sh13,i1)

!     i1 ( h2 p14 h11 h13 )_ytt + = 1 * Sum ( h9 p7 ) * t ( p7 p14 h9 h11 )_t * i2 ( h2 h9 h13 p7 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh2,sp14,sh11,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer :: h2,p14,h11,h13
  integer :: h9,p7,sh9,sp7
  integer :: spin_t2inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sp7 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp7,sp14,sh9,sh11)
     spin_itm_hhhp = tdcc_spin_dummy2(sh2,sh9,sh13,sp7)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call bccdt_l1p_12_4_1(sh2,sh9,sh13,sp7,itm_hhhp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p14 = norb1+1,nact
     do h11 = 1,norb1
     do h13 = 1,norb1
     do h9 = 1,norb1
     do p7 = norb1+1,nact
        i1(h2,p14,h11,h13) = i1(h2,p14,h11,h13) + fact * &
             t2inp(p7,p14,h9,h11,spin_t2inp) * itm_hhhp(h2,h9,h13,p7)
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
end subroutine bccdt_l1p_12_4
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_12_4_1(sh2,sh9,sh13,sp7,i2)

!         i2 ( h2 h9 h13 p7 )_yt + = -1 * Sum ( h5 p3 p4 ) * t ( p3 p4 h5 h13 )_t * y ( h2 h5 h9 p3 p4 p7 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh2,sh9,sh13,sp7
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h2,h9,h13,p7
  integer :: h5,p3,p4,sh5,sp3,sp4
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh13)
     spin_g3inp = tdcc_spin_g3inp(sh2,sh5,sh9,sp3,sp4,sp7)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do h9 = 1,norb1
     do h13 = 1,norb1
     do p7 = norb1+1,nact
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i2(h2,h9,h13,p7) = i2(h2,h9,h13,p7) + fact * &
             t2inp(p3,p4,h5,h13,spin_t2inp) * g3inp(h2,h5,h9,p3,p4,p7,spin_g3inp)
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
end subroutine bccdt_l1p_12_4_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_13(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1/4 * Sum ( h8 h12 h13 ) * i1 ( h2 h8 h12 h13 )_yt * v ( h12 h13 h8 p1 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h8,h12,h13,sh8,sh12,sh13
  integer :: spin_itm_hhhh
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh8 = 1,2
  do sh12 = 1,2
  do sh13 = 1,2
     spin_itm_hhhh = tdcc_spin_dummy2(sh2,sh8,sh12,sh13)
     spin_int2x = tdcc_spin_int2x(sh12,sh13,sh8,sp1)
     if(spin_itm_hhhh * spin_int2x == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call bccdt_l1p_13_1(sh2,sh8,sh12,sh13,itm_hhhh)
     call bccdt_l1p_13_2(sh2,sh8,sh12,sh13,itm_hhhh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h8 = 1,norb1
     do h12 = 1,norb1
     do h13 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_hhhh(h2,h8,h12,h13) * int2x(h12,h13,h8,p1,spin_int2x)
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
  deallocate(itm_hhhh)
end subroutine bccdt_l1p_13
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_13_1(sh2,sh8,sh12,sh13,i1)

!     i1 ( h2 h8 h12 h13 )_yt + = -1 * Sum ( p3 p4 ) * t ( p3 p4 h12 h13 )_t * y ( h2 h8 p3 p4 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh2,sh8,sh12,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h2,h8,h12,h13
  integer :: p3,p4,sp3,sp4
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh12,sh13)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh8,sp3,sp4)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do h8 = 1,norb1
     do h12 = 1,norb1
     do h13 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,h8,h12,h13) = i1(h2,h8,h12,h13) + fact * &
             t2inp(p3,p4,h12,h13,spin_t2inp) * g2inp(h2,h8,p3,p4,spin_g2inp)
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
end subroutine bccdt_l1p_13_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_13_2(sh2,sh8,sh12,sh13,i1)

!     i1 ( h2 h8 h12 h13 )_yt + = 1/3 * Sum ( h6 p3 p4 p5 ) * t ( p3 p4 p5 h6 h12 h13 )_t * y ( h2 h6 h8 p3 p4 p5 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh2,sh8,sh12,sh13
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer :: h2,h8,h12,h13
  integer :: h6,p3,p4,p5,sh6,sp3,sp4,sp5
  integer :: spin_t3inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 3.0d+0 * runit

  do sh6 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp4,sp5,sh6,sh12,sh13)
     spin_g3inp = tdcc_spin_g3inp(sh2,sh6,sh8,sp3,sp4,sp5)
     if(spin_t3inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do h8 = 1,norb1
     do h12 = 1,norb1
     do h13 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h2,h8,h12,h13) = i1(h2,h8,h12,h13) + fact * &
             t3inp(p3,p4,p5,h6,h12,h13,spin_t3inp) * g3inp(h2,h6,h8,p3,p4,p5,spin_g3inp)
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
end subroutine bccdt_l1p_13_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_14(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( p11 h10 p12 ) * i1 ( h2 p11 h10 p12 )_yt * v ( h10 p12 p1 p11 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: p11,h10,p12,sp11,sh10,sp12
  integer :: spin_itm_hphp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp11 = 1,2
  do sh10 = 1,2
  do sp12 = 1,2
     spin_itm_hphp = tdcc_spin_dummy2(sh2,sp11,sh10,sp12)
     spin_int2x = tdcc_spin_int2x(sh10,sp12,sp1,sp11)
     if(spin_itm_hphp * spin_int2x == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call bccdt_l1p_14_1(sh2,sp11,sh10,sp12,itm_hphp)
     call bccdt_l1p_14_2(sh2,sp11,sh10,sp12,itm_hphp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p11 = norb1+1,nact
     do h10 = 1,norb1
     do p12 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_hphp(h2,p11,h10,p12) * int2x(h10,p12,p1,p11,spin_int2x)
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
  deallocate(itm_hphp)
end subroutine bccdt_l1p_14
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_14_1(sh2,sp11,sh10,sp12,i1)

!     i1 ( h2 p11 h10 p12 )_yt + = 1 * Sum ( h5 p3 ) * t ( p3 p11 h5 h10 )_t * y ( h2 h5 p3 p12 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer,intent(in) :: sh2,sp11,sh10,sp12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h2,p11,h10,p12
  integer :: h5,p3,sh5,sp3
  integer :: spin_t2inp
  integer :: spin_g2inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp11,sh5,sh10)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh5,sp3,sp12)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do h10 = 1,norb1
     do p12 = norb1+1,nact
     do h5 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,p11,h10,p12) = i1(h2,p11,h10,p12) + fact * &
             t2inp(p3,p11,h5,h10,spin_t2inp) * g2inp(h2,h5,p3,p12,spin_g2inp)
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
end subroutine bccdt_l1p_14_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_14_2(sh2,sp11,sh10,sp12,i1)

!     i1 ( h2 p11 h10 p12 )_yt + = -1/4 * Sum ( h6 h7 p3 p4 ) * t ( p3 p4 p11 h6 h7 h10 )_t * y ( h2 h6 h7 p3 p4 p12 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh2,sp11,sh10,sp12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h2,p11,h10,p12
  integer :: h6,h7,p3,p4,sh6,sh7,sp3,sp4
  integer :: spin_t3inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  do sh6 = 1,2
  do sh7 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp3,sp4,sp11,sh6,sh7,sh10)
     spin_g3inp = tdcc_spin_g3inp(sh2,sh6,sh7,sp3,sp4,sp12)
     if(spin_t3inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do h10 = 1,norb1
     do p12 = norb1+1,nact
     do h6 = 1,norb1
     do h7 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,p11,h10,p12) = i1(h2,p11,h10,p12) + fact * &
             t3inp(p3,p4,p11,h6,h7,h10,spin_t3inp) * g3inp(h2,h6,h7,p3,p4,p12,spin_g3inp)
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
end subroutine bccdt_l1p_14_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_15(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = -1/4 * Sum ( h7 h8 h6 ) * i1 ( h7 h8 h6 p1 )_yt * v ( h2 h6 h7 h8 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: h7,h8,h6,sh7,sh8,sh6
  integer :: spin_itm_hhhp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh7 = 1,2
  do sh8 = 1,2
  do sh6 = 1,2
     spin_itm_hhhp = tdcc_spin_dummy2(sh7,sh8,sh6,sp1)
     spin_int2x = tdcc_spin_int2x(sh2,sh6,sh7,sh8)
     if(spin_itm_hhhp * spin_int2x == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call bccdt_l1p_15_1(sh7,sh8,sh6,sp1,itm_hhhp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do h6 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_hhhp(h7,h8,h6,p1) * int2x(h2,h6,h7,h8,spin_int2x)
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
  deallocate(itm_hhhp)
end subroutine bccdt_l1p_15
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_15_1(sh7,sh8,sh6,sp1,i1)

!     i1 ( h7 h8 h6 p1 )_yt + = 1 * Sum ( h5 p3 p4 ) * t ( p3 p4 h5 h6 )_t * y ( h5 h7 h8 p1 p3 p4 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh7,sh8,sh6,sp1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h7,h8,h6,p1
  integer :: h5,p3,p4,sh5,sp3,sp4
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh6)
     spin_g3inp = tdcc_spin_g3inp(sh5,sh7,sh8,sp1,sp3,sp4)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h7 = 1,norb1
     do h8 = 1,norb1
     do h6 = 1,norb1
     do p1 = norb1+1,nact
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h7,h8,h6,p1) = i1(h7,h8,h6,p1) + fact * &
             t2inp(p3,p4,h5,h6,spin_t2inp) * g3inp(h5,h7,h8,p1,p3,p4,spin_g3inp)
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
end subroutine bccdt_l1p_15_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_16(sh2,sp1,i0)

! i0 ( h2 p1 )_vty + = -1/2 * Sum ( p15 p11 h14 h13 h16 ) * y ( h13 h14 h16 p1 p11 p15 )_y * i1 ( h2 p11 p15 h13 h14 h16 )_vt 4

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: p15,p11,h14,h13,h16,sp15,sp11,sh14,sh13,sh16
  integer :: spin_g3inp
  integer :: spin_itm_hpphhh
  integer,external :: tdcc_spin_g3inp
  integer,external :: tdcc_spin_dummy3
  complex(kind(0d0)),allocatable :: itm_hpphhh(:,:,:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hpphhh(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1))
  do sp15 = 1,2
  do sp11 = 1,2
  do sh14 = 1,2
  do sh13 = 1,2
  do sh16 = 1,2
     spin_g3inp = tdcc_spin_g3inp(sh13,sh14,sh16,sp1,sp11,sp15)
     spin_itm_hpphhh = tdcc_spin_dummy3(sh2,sp11,sp15,sh13,sh14,sh16)
     if(spin_g3inp * spin_itm_hpphhh == 0) cycle

     itm_hpphhh(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1) = czero
     call bccdt_l1p_16_1(sh2,sp11,sp15,sh13,sh14,sh16,itm_hpphhh)
     call bccdt_l1p_16_2(sh2,sp11,sp15,sh13,sh14,sh16,itm_hpphhh)
     call bccdt_l1p_16_3(sh2,sp11,sp15,sh13,sh14,sh16,itm_hpphhh)
     call bccdt_l1p_16_4(sh2,sp11,sp15,sh13,sh14,sh16,itm_hpphhh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p15 = norb1+1,nact
     do p11 = norb1+1,nact
     do h14 = 1,norb1
     do h13 = 1,norb1
     do h16 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * &
             g3inp(h13,h14,h16,p1,p11,p15,spin_g3inp) * itm_hpphhh(h2,p11,p15,h13,h14,h16)
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
  deallocate(itm_hpphhh)
end subroutine bccdt_l1p_16
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_16_1(sh2,sp11,sp15,sh13,sh14,sh16,i1)

!     i1 ( h2 p11 p15 h13 h14 h16 )_vt + = 1 * Sum ( p4 ) * t ( p4 p15 h13 h14 )_t * i2 ( h2 p11 h16 p4 )_v 2

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sh2,sp11,sp15,sh13,sh14,sh16
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h2,p11,p15,h13,h14,h16
  integer :: p4,sp4
  integer :: spin_t2inp
  integer :: spin_itm_hphp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp15,sh13,sh14)
     spin_itm_hphp = tdcc_spin_int2x(sh2,sp11,sh16,sp4)
     if(spin_t2inp * spin_itm_hphp == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call bccdt_l1p_16_1_1(sh2,sp11,sh16,sp4,itm_hphp)
     call bccdt_l1p_16_1_2(sh2,sp11,sh16,sp4,itm_hphp)

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do p15 = norb1+1,nact
     do h13 = 1,norb1
     do h14 = 1,norb1
     do h16 = 1,norb1
     do p4 = norb1+1,nact
        i1(h2,p11,p15,h13,h14,h16) = i1(h2,p11,p15,h13,h14,h16) + fact * &
             t2inp(p4,p15,h13,h14,spin_t2inp) * itm_hphp(h2,p11,h16,p4)
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
  deallocate(itm_hphp)
end subroutine bccdt_l1p_16_1
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_16_1_1(sh2,sp11,sh16,sp4,i2)

!         i2 ( h2 p11 h16 p4 )_v + = 1 * v ( h2 p11 h16 p4 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp11,sh16,sp4
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h2,p11,h16,p4
  integer :: sdum
  integer :: spin_int2x
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh2,sp11,sh16,sp4)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h2 = 1,norb1
  do p11 = norb1+1,nact
  do h16 = 1,norb1
  do p4 = norb1+1,nact
     i2(h2,p11,h16,p4) = i2(h2,p11,h16,p4) + fact * &
          int2x(h2,p11,h16,p4,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccdt_l1p_16_1_1
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_16_1_2(sh2,sp11,sh16,sp4,i2)

!         i2 ( h2 p11 h16 p4 )_vt + = -1 * Sum ( h10 p8 ) * t ( p8 p11 h10 h16 )_t * v ( h2 h10 p4 p8 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh2,sp11,sh16,sp4
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: h2,p11,h16,p4
  integer :: h10,p8,sh10,sp8
  integer :: spin_t2inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sh10 = 1,2
  do sp8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp8,sp11,sh10,sh16)
     spin_int2x = tdcc_spin_int2x(sh2,sh10,sp4,sp8)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do h16 = 1,norb1
     do p4 = norb1+1,nact
     do h10 = 1,norb1
     do p8 = norb1+1,nact
        i2(h2,p11,h16,p4) = i2(h2,p11,h16,p4) + fact * &
             t2inp(p8,p11,h10,h16,spin_t2inp) * int2x(h2,h10,p4,p8,spin_int2x)
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
end subroutine bccdt_l1p_16_1_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_16_2(sh2,sp11,sp15,sh13,sh14,sh16,i1)

!     i1 ( h2 p11 p15 h13 h14 h16 )_ft + = 1/6 * Sum ( p4 ) * t ( p4 p11 p15 h13 h14 h16 )_t * f ( h2 p4 )_f 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : fock

  implicit none
  integer,intent(in) :: sh2,sp11,sp15,sh13,sh14,sh16
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h2,p11,p15,h13,h14,h16
  integer :: p4,sp4
  integer :: spin_t3inp
  integer :: spin_fock
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 6.0d+0 * runit

  do sp4 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp11,sp15,sh13,sh14,sh16)
     spin_fock = tdcc_spin_fock(sh2,sp4)
     if(spin_t3inp * spin_fock == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do p15 = norb1+1,nact
     do h13 = 1,norb1
     do h14 = 1,norb1
     do h16 = 1,norb1
     do p4 = norb1+1,nact
        i1(h2,p11,p15,h13,h14,h16) = i1(h2,p11,p15,h13,h14,h16) + fact * &
             t3inp(p4,p11,p15,h13,h14,h16,spin_t3inp) * fock(h2,p4,spin_fock)
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
end subroutine bccdt_l1p_16_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_16_3(sh2,sp11,sp15,sh13,sh14,sh16,i1)

!     i1 ( h2 p11 p15 h13 h14 h16 )_vt + = 1/2 * Sum ( h8 p5 ) * t ( p5 p11 p15 h8 h13 h14 )_t * v ( h2 h8 h16 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh2,sp11,sp15,sh13,sh14,sh16
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h2,p11,p15,h13,h14,h16
  integer :: h8,p5,sh8,sp5
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh8 = 1,2
  do sp5 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp5,sp11,sp15,sh8,sh13,sh14)
     spin_int2x = tdcc_spin_int2x(sh2,sh8,sh16,sp5)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do p15 = norb1+1,nact
     do h13 = 1,norb1
     do h14 = 1,norb1
     do h16 = 1,norb1
     do h8 = 1,norb1
     do p5 = norb1+1,nact
        i1(h2,p11,p15,h13,h14,h16) = i1(h2,p11,p15,h13,h14,h16) + fact * &
             t3inp(p5,p11,p15,h8,h13,h14,spin_t3inp) * int2x(h2,h8,h16,p5,spin_int2x)
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
end subroutine bccdt_l1p_16_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_16_4(sh2,sp11,sp15,sh13,sh14,sh16,i1)

!     i1 ( h2 p11 p15 h13 h14 h16 )_vt + = 1/6 * Sum ( p4 p5 ) * t ( p4 p5 p15 h13 h14 h16 )_t * v ( h2 p11 p4 p5 )_v 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer,intent(in) :: sh2,sp11,sp15,sh13,sh14,sh16
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer :: h2,p11,p15,h13,h14,h16
  integer :: p4,p5,sp4,sp5
  integer :: spin_t3inp
  integer :: spin_int2x
  integer,external :: tdcc_spin_t3inp
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 6.0d+0 * runit

  do sp4 = 1,2
  do sp5 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp4,sp5,sp15,sh13,sh14,sh16)
     spin_int2x = tdcc_spin_int2x(sh2,sp11,sp4,sp5)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(6)
     do h2 = 1,norb1
     do p11 = norb1+1,nact
     do p15 = norb1+1,nact
     do h13 = 1,norb1
     do h14 = 1,norb1
     do h16 = 1,norb1
     do p4 = norb1+1,nact
     do p5 = norb1+1,nact
        i1(h2,p11,p15,h13,h14,h16) = i1(h2,p11,p15,h13,h14,h16) + fact * &
             t3inp(p4,p5,p15,h13,h14,h16,spin_t3inp) * int2x(h2,p11,p4,p5,spin_int2x)
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
end subroutine bccdt_l1p_16_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_17(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = -1/4 * Sum ( p4 p7 p8 ) * i1 ( h2 p4 p7 p8 )_yt * v ( p7 p8 p1 p4 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: p4,p7,p8,sp4,sp7,sp8
  integer :: spin_itm_hppp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hppp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact))
  do sp4 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_itm_hppp = tdcc_spin_dummy2(sh2,sp4,sp7,sp8)
     spin_int2x = tdcc_spin_int2x(sp7,sp8,sp1,sp4)
     if(spin_itm_hppp * spin_int2x == 0) cycle

     itm_hppp(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact) = czero
     call bccdt_l1p_17_1(sh2,sp4,sp7,sp8,itm_hppp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p4 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_hppp(h2,p4,p7,p8) * int2x(p7,p8,p1,p4,spin_int2x)
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
  deallocate(itm_hppp)
end subroutine bccdt_l1p_17
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_17_1(sh2,sp4,sp7,sp8,i1)

!     i1 ( h2 p4 p7 p8 )_yt + = 1 * Sum ( h5 h6 p3 ) * t ( p3 p4 h5 h6 )_t * y ( h2 h5 h6 p3 p7 p8 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh2,sp4,sp7,sp8
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,(norb1+1):nact)
  integer :: h2,p4,p7,p8
  integer :: h5,h6,p3,sh5,sh6,sp3
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh6)
     spin_g3inp = tdcc_spin_g3inp(sh2,sh5,sh6,sp3,sp7,sp8)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p4 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,p4,p7,p8) = i1(h2,p4,p7,p8) + fact * &
             t2inp(p3,p4,h5,h6,spin_t2inp) * g3inp(h2,h5,h6,p3,p7,p8,spin_g3inp)
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
end subroutine bccdt_l1p_17_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_18(sh2,sp1,i0)

! i0 ( h2 p1 )_yttv + = -1/8 * Sum ( p12 p13 h11 ) * i1 ( p12 p13 h11 p1 )_ytt * v ( h2 h11 p12 p13 )_v 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer,intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer :: h2,p1
  integer :: p12,p13,h11,sp12,sp13,sh11
  integer :: spin_itm_pphp
  integer :: spin_int2x
  integer,external :: tdcc_spin_dummy2
  integer,external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 8.0d+0 * runit

  allocate(itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp12 = 1,2
  do sp13 = 1,2
  do sh11 = 1,2
     spin_itm_pphp = tdcc_spin_dummy2(sp12,sp13,sh11,sp1)
     spin_int2x = tdcc_spin_int2x(sh2,sh11,sp12,sp13)
     if(spin_itm_pphp * spin_int2x == 0) cycle

     itm_pphp((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call bccdt_l1p_18_1(sp12,sp13,sh11,sp1,itm_pphp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p12 = norb1+1,nact
     do p13 = norb1+1,nact
     do h11 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * &
             itm_pphp(p12,p13,h11,p1) * int2x(h2,h11,p12,p13,spin_int2x)
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
end subroutine bccdt_l1p_18
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_18_1(sp12,sp13,sh11,sp1,i1)

!     i1 ( p12 p13 h11 p1 )_ytt + = 1 * Sum ( h9 h10 ) * t ( p12 p13 h9 h10 )_t * i2 ( h9 h10 h11 p1 )_yt 1

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer,intent(in) :: sp12,sp13,sh11,sp1
  complex(kind(0d0)),intent(inout) :: &
       i1((norb1+1):nact,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer :: p12,p13,h11,p1
  integer :: h9,h10,sh9,sh10
  integer :: spin_t2inp
  integer :: spin_itm_hhhp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_dummy2
  complex(kind(0d0)),allocatable :: itm_hhhp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact))
  do sh9 = 1,2
  do sh10 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp12,sp13,sh9,sh10)
     spin_itm_hhhp = tdcc_spin_dummy2(sh9,sh10,sh11,sp1)
     if(spin_t2inp * spin_itm_hhhp == 0) cycle

     itm_hhhp(1:norb1,1:norb1,1:norb1,(norb1+1):nact) = czero
     call bccdt_l1p_18_1_1(sh9,sh10,sh11,sp1,itm_hhhp)

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do p12 = norb1+1,nact
     do p13 = norb1+1,nact
     do h11 = 1,norb1
     do p1 = norb1+1,nact
     do h9 = 1,norb1
     do h10 = 1,norb1
        i1(p12,p13,h11,p1) = i1(p12,p13,h11,p1) + fact * &
             t2inp(p12,p13,h9,h10,spin_t2inp) * itm_hhhp(h9,h10,h11,p1)
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
end subroutine bccdt_l1p_18_1
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_l1p_18_1_1(sh9,sh10,sh11,sp1,i2)

!         i2 ( h9 h10 h11 p1 )_yt + = 1 * Sum ( h5 p3 p4 ) * t ( p3 p4 h5 h11 )_t * y ( h5 h9 h10 p1 p3 p4 )_y 0

  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g3inp

  implicit none
  integer,intent(in) :: sh9,sh10,sh11,sp1
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,1:norb1,1:norb1,(norb1+1):nact)
  integer :: h9,h10,h11,p1
  integer :: h5,p3,p4,sh5,sp3,sp4
  integer :: spin_t2inp
  integer :: spin_g3inp
  integer,external :: tdcc_spin_t2inp
  integer,external :: tdcc_spin_g3inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh11)
     spin_g3inp = tdcc_spin_g3inp(sh5,sh9,sh10,sp1,sp3,sp4)
     if(spin_t2inp * spin_g3inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h9 = 1,norb1
     do h10 = 1,norb1
     do h11 = 1,norb1
     do p1 = norb1+1,nact
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i2(h9,h10,h11,p1) = i2(h9,h10,h11,p1) + fact * &
             t2inp(p3,p4,h5,h11,spin_t2inp) * g3inp(h5,h9,h10,p1,p3,p4,spin_g3inp)
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
end subroutine bccdt_l1p_18_1_1
!##########################################################
!##########################################################
!##########################################################
