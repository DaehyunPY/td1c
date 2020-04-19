!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_1(sh2,sp1,i0)

! i0 ( h2 p1 )_f + = 1 * f ( h2 p1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
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
end subroutine bccd_l1p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_2(sh2,sp1,i0)

! i0 ( h2 p1 )_yf + = -1 * Sum ( h7 ) * y ( h7 p1 )_y * i1 ( h2 h7 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h7,sh7
  integer(c_int) :: spin_g1inp
  integer(c_int) :: spin_itm_hh
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh7 = 1,2
     spin_g1inp = tdcc_spin_g1inp(sh7,sp1)
     spin_itm_hh = tdcc_spin_fock(sh2,sh7)
     if(spin_g1inp * spin_itm_hh == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call bccd_l1p_2_1(sh2,sh7,itm_hh)
     call bccd_l1p_2_2(sh2,sh7,itm_hh)

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
end subroutine bccd_l1p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_2_1(sh2,sh7,i1)

!     i1 ( h2 h7 )_f + = 1 * f ( h2 h7 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh2,sh7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h2,h7
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
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
end subroutine bccd_l1p_2_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_2_2(sh2,sh7,i1)

!     i1 ( h2 h7 )_vt + = -1/2 * Sum ( h6 p3 p4 ) * t ( p3 p4 h6 h7 )_t * v ( h2 h6 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sh7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h2,h7
  integer(c_int) :: h6,p3,p4,sh6,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
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
        i1(h2,h7) = i1(h2,h7) + fact * t2inp(p3,p4,h6,h7,spin_t2inp) * int2x(h2,h6,p3,p4,spin_int2x)
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
end subroutine bccd_l1p_2_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_3(sh2,sp1,i0)

! i0 ( h2 p1 )_yf + = 1 * Sum ( p7 ) * y ( h2 p7 )_y * f ( p7 p1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : fock

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: p7,sp7
  integer(c_int) :: spin_g1inp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_fock
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
        i0(h2,p1) = i0(h2,p1) + fact * g1inp(h2,p7,spin_g1inp) * fock(p7,p1,spin_fock)
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
end subroutine bccd_l1p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_4(sh2,sp1,i0)

! i0 ( h2 p1 )_yv + = -1 * Sum ( h4 p3 ) * y ( h4 p3 )_y * v ( h2 p3 h4 p1 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h4,p3,sh4,sp3
  integer(c_int) :: spin_g1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_int2x
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
        i0(h2,p1) = i0(h2,p1) + fact * g1inp(h4,p3,spin_g1inp) * int2x(h2,p3,h4,p1,spin_int2x)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccd_l1p_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_5(sh2,sp1,i0)

! i0 ( h2 p1 )_yv + = -1/2 * Sum ( h11 h12 p9 ) * y ( h11 h12 p1 p9 )_y * i1 ( h2 p9 h11 h12 )_v 4

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h11,h12,p9,sh11,sh12,sp9
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_itm_hphh
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sh11 = 1,2
  do sh12 = 1,2
  do sp9 = 1,2
     spin_g2inp = tdcc_spin_g2inp(sh11,sh12,sp1,sp9)
     spin_itm_hphh = tdcc_spin_int2x(sh2,sp9,sh11,sh12)
     if(spin_g2inp * spin_itm_hphh == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call bccd_l1p_5_1(sh2,sp9,sh11,sh12,itm_hphh)
     call bccd_l1p_5_2(sh2,sp9,sh11,sh12,itm_hphh)
     call bccd_l1p_5_3(sh2,sp9,sh11,sh12,itm_hphh)
     call bccd_l1p_5_4(sh2,sp9,sh11,sh12,itm_hphh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p9 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * g2inp(h11,h12,p1,p9,spin_g2inp) * itm_hphh(h2,p9,h11,h12)
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
end subroutine bccd_l1p_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_5_1(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_v + = 1 * v ( h2 p9 h11 h12 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_int2x = tdcc_spin_int2x(sh2,sp9,sh11,sh12)
  !$omp parallel default(shared)
  !$omp do collapse(4)
  do h2 = 1,norb1
  do p9 = norb1+1,nact
  do h11 = 1,norb1
  do h12 = 1,norb1
     i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * int2x(h2,p9,h11,h12,spin_int2x)
  end do
  end do
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccd_l1p_5_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_5_2(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_ft + = 1 * Sum ( p5 ) * t ( p5 p9 h11 h12 )_t * f ( h2 p5 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : fock

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: p5,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp5 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp5,sp9,sh11,sh12)
     spin_fock = tdcc_spin_fock(sh2,sp5)
     if(spin_t2inp * spin_fock == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p5 = norb1+1,nact
        i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * t2inp(p5,p9,h11,h12,spin_t2inp) * fock(h2,p5,spin_fock)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
end subroutine bccd_l1p_5_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_5_3(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_vt + = -2 * Sum ( h6 p4 ) * t ( p4 p9 h6 h11 )_t * v ( h2 h6 h12 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: h6,p4,sh6,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -2.0d+0 * runit

  do sh6 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp4,sp9,sh6,sh11)
     spin_int2x = tdcc_spin_int2x(sh2,sh6,sh12,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do h6 = 1,norb1
     do p4 = norb1+1,nact
        i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * t2inp(p4,p9,h6,h11,spin_t2inp) * int2x(h2,h6,h12,p4,spin_int2x)
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
end subroutine bccd_l1p_5_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_5_4(sh2,sp9,sh11,sh12,i1)

!     i1 ( h2 p9 h11 h12 )_vt + = 1/2 * Sum ( p3 p4 ) * t ( p3 p4 h11 h12 )_t * v ( h2 p9 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh11,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h11,h12
  integer(c_int) :: p3,p4,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh11,sh12)
     spin_int2x = tdcc_spin_int2x(sh2,sp9,sp3,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h11 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,p9,h11,h12) = i1(h2,p9,h11,h12) + fact * t2inp(p3,p4,h11,h12,spin_t2inp) * int2x(h2,p9,p3,p4,spin_int2x)
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
end subroutine bccd_l1p_5_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_6(sh2,sp1,i0)

! i0 ( h2 p1 )_yv + = -1/2 * Sum ( h7 p8 p5 ) * y ( h2 h7 p5 p8 )_y * v ( p5 p8 h7 p1 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h7,p8,p5,sh7,sp8,sp5
  integer(c_int) :: spin_g2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_int2x
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
        i0(h2,p1) = i0(h2,p1) + fact * g2inp(h2,h7,p5,p8,spin_g2inp) * int2x(p5,p8,h7,p1,spin_int2x)
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
end subroutine bccd_l1p_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_7(sh2,sp1,i0)

! i0 ( h2 p1 )_vt + = 1 * Sum ( h10 p9 ) * i1 ( p9 h10 )_t * v ( h2 h10 p1 p9 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h10,p9,sh10,sp9
  integer(c_int) :: spin_itm_ph
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_ph(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_ph((norb1+1):nact,1:norb1))
  do sh10 = 1,2
  do sp9 = 1,2
     spin_itm_ph = tdcc_spin_t1inp(sp9,sh10)
     spin_int2x = tdcc_spin_int2x(sh2,sh10,sp1,sp9)
     if(spin_itm_ph * spin_int2x == 0) cycle

     itm_ph((norb1+1):nact,1:norb1) = czero
     call bccd_l1p_7_1(sp9,sh10,itm_ph)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h10 = 1,norb1
     do p9 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * itm_ph(p9,h10) * int2x(h2,h10,p1,p9,spin_int2x)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_ph)
end subroutine bccd_l1p_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_7_1(sp9,sh10,i1)

!     i1 ( p9 h10 )_yt + = 1 * Sum ( h5 p3 ) * t ( p3 p9 h5 h10 )_t * y ( h5 p3 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g1inp

  implicit none
  integer(c_int),intent(in) :: sp9,sh10
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,1:norb1)
  integer(c_int) :: p9,h10
  integer(c_int) :: h5,p3,sh5,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g1inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp9,sh5,sh10)
     spin_g1inp = tdcc_spin_g1inp(sh5,sp3)
     if(spin_t2inp * spin_g1inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p9 = norb1+1,nact
     do h10 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
        i1(p9,h10) = i1(p9,h10) + fact * t2inp(p3,p9,h5,h10,spin_t2inp) * g1inp(h5,p3,spin_g1inp)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccd_l1p_7_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_8(sh2,sp1,i0)

! i0 ( h2 p1 )_ytf + = -1 * Sum ( h3 ) * i1 ( h2 h3 )_yt * f ( h3 p1 )_f 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h3,sh3
  integer(c_int) :: spin_itm_hh
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh3 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh2,sh3)
     spin_fock = tdcc_spin_fock(sh3,sp1)
     if(spin_itm_hh * spin_fock == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call bccd_l1p_8_1(sh2,sh3,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h3 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * itm_hh(h2,h3) * fock(h3,p1,spin_fock)
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  deallocate(itm_hh)
end subroutine bccd_l1p_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_8_1(sh2,sh3,i1)

!     i1 ( h2 h3 )_yt + = 1/2 * Sum ( h6 p4 p5 ) * t ( p4 p5 h3 h6 )_t * y ( h2 h6 p4 p5 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh3
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h2,h3
  integer(c_int) :: h6,p4,p5,sh6,sp4,sp5
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
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
        i1(h2,h3) = i1(h2,h3) + fact * t2inp(p4,p5,h3,h6,spin_t2inp) * g2inp(h2,h6,p4,p5,spin_g2inp)
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
end subroutine bccd_l1p_8_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_9(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( h6 h8 ) * i1 ( h6 h8 )_yt * v ( h2 h8 h6 p1 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h6,h8,sh6,sh8
  integer(c_int) :: spin_itm_hh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh6 = 1,2
  do sh8 = 1,2
     spin_itm_hh = tdcc_spin_dummy1(sh6,sh8)
     spin_int2x = tdcc_spin_int2x(sh2,sh8,sh6,sp1)
     if(spin_itm_hh * spin_int2x == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call bccd_l1p_9_1(sh6,sh8,itm_hh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h6 = 1,norb1
     do h8 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * itm_hh(h6,h8) * int2x(h2,h8,h6,p1,spin_int2x)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_hh)
end subroutine bccd_l1p_9
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_9_1(sh6,sh8,i1)

!     i1 ( h6 h8 )_yt + = 1/2 * Sum ( h5 p3 p4 ) * t ( p3 p4 h5 h8 )_t * y ( h5 h6 p3 p4 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh6,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h6,h8
  integer(c_int) :: h5,p3,p4,sh5,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh5,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh5,sh6,sp3,sp4)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h6 = 1,norb1
     do h8 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h6,h8) = i1(h6,h8) + fact * t2inp(p3,p4,h5,h8,spin_t2inp) * g2inp(h5,h6,p3,p4,spin_g2inp)
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
end subroutine bccd_l1p_9_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_10(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( p7 p8 ) * i1 ( p7 p8 )_yt * v ( h2 p8 p1 p7 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: p7,p8,sp7,sp8
  integer(c_int) :: spin_itm_pp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp7 = 1,2
  do sp8 = 1,2
     spin_itm_pp = tdcc_spin_dummy1(sp7,sp8)
     spin_int2x = tdcc_spin_int2x(sh2,sp8,sp1,sp7)
     if(spin_itm_pp * spin_int2x == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call bccd_l1p_10_1(sp7,sp8,itm_pp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * itm_pp(p7,p8) * int2x(h2,p8,p1,p7,spin_int2x)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
  deallocate(itm_pp)
end subroutine bccd_l1p_10
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_10_1(sp7,sp8,i1)

!     i1 ( p7 p8 )_yt + = 1/2 * Sum ( h5 h6 p3 ) * t ( p3 p7 h5 h6 )_t * y ( h5 h6 p3 p8 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sp7,sp8
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p7,p8
  integer(c_int) :: h5,h6,p3,sh5,sh6,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp7,sh5,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh5,sh6,sp3,sp8)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
        i1(p7,p8) = i1(p7,p8) + fact * t2inp(p3,p7,h5,h6,spin_t2inp) * g2inp(h5,h6,p3,p8,spin_g2inp)
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
end subroutine bccd_l1p_10_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_11(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1/2 * Sum ( p9 h12 h6 ) * i1 ( h2 p9 h6 h12 )_yt * v ( h6 h12 p1 p9 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: p9,h12,h6,sp9,sh12,sh6
  integer(c_int) :: spin_itm_hphh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 2.0d+0 * runit

  allocate(itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1))
  do sp9 = 1,2
  do sh12 = 1,2
  do sh6 = 1,2
     spin_itm_hphh = tdcc_spin_dummy2(sh2,sp9,sh6,sh12)
     spin_int2x = tdcc_spin_int2x(sh6,sh12,sp1,sp9)
     if(spin_itm_hphh * spin_int2x == 0) cycle

     itm_hphh(1:norb1,(norb1+1):nact,1:norb1,1:norb1) = czero
     call bccd_l1p_11_1(sh2,sp9,sh6,sh12,itm_hphh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p9 = norb1+1,nact
     do h12 = 1,norb1
     do h6 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * itm_hphh(h2,p9,h6,h12) * int2x(h6,h12,p1,p9,spin_int2x)
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
end subroutine bccd_l1p_11
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_11_1(sh2,sp9,sh6,sh12,i1)

!     i1 ( h2 p9 h6 h12 )_yt + = -1 * Sum ( p3 ) * t ( p3 p9 h6 h12 )_t * y ( h2 p3 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g1inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp9,sh6,sh12
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,1:norb1)
  integer(c_int) :: h2,p9,h6,h12
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g1inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp9,sh6,sh12)
     spin_g1inp = tdcc_spin_g1inp(sh2,sp3)
     if(spin_t2inp * spin_g1inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p9 = norb1+1,nact
     do h6 = 1,norb1
     do h12 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,p9,h6,h12) = i1(h2,p9,h6,h12) + fact * t2inp(p3,p9,h6,h12,spin_t2inp) * g1inp(h2,p3,spin_g1inp)
     end do
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
end subroutine bccd_l1p_11_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_12(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = -1/4 * Sum ( h7 h8 h6 ) * i1 ( h2 h7 h6 h8 )_yt * v ( h6 h8 h7 p1 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: h7,h8,h6,sh7,sh8,sh6
  integer(c_int) :: spin_itm_hhhh
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hhhh(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1))
  do sh7 = 1,2
  do sh8 = 1,2
  do sh6 = 1,2
     spin_itm_hhhh = tdcc_spin_dummy2(sh2,sh7,sh6,sh8)
     spin_int2x = tdcc_spin_int2x(sh6,sh8,sh7,sp1)
     if(spin_itm_hhhh * spin_int2x == 0) cycle

     itm_hhhh(1:norb1,1:norb1,1:norb1,1:norb1) = czero
     call bccd_l1p_12_1(sh2,sh7,sh6,sh8,itm_hhhh)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do h7 = 1,norb1
     do h8 = 1,norb1
     do h6 = 1,norb1
        i0(h2,p1) = i0(h2,p1) + fact * itm_hhhh(h2,h7,h6,h8) * int2x(h6,h8,h7,p1,spin_int2x)
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
end subroutine bccd_l1p_12
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_12_1(sh2,sh7,sh6,sh8,i1)

!     i1 ( h2 h7 h6 h8 )_yt + = 1 * Sum ( p3 p4 ) * t ( p3 p4 h6 h8 )_t * y ( h2 h7 p3 p4 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sh7,sh6,sh8
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,h7,h6,h8
  integer(c_int) :: p3,p4,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh6,sh8)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh7,sp3,sp4)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do h7 = 1,norb1
     do h6 = 1,norb1
     do h8 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i1(h2,h7,h6,h8) = i1(h2,h7,h6,h8) + fact * t2inp(p3,p4,h6,h8,spin_t2inp) * g2inp(h2,h7,p3,p4,spin_g2inp)
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
end subroutine bccd_l1p_12_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_13(sh2,sp1,i0)

! i0 ( h2 p1 )_ytv + = 1 * Sum ( p8 h6 p7 ) * i1 ( h2 p8 h6 p7 )_yt * v ( h6 p7 p1 p8 )_v 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: p8,h6,p7,sp8,sh6,sp7
  integer(c_int) :: spin_itm_hphp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),allocatable :: itm_hphp(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  allocate(itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact))
  do sp8 = 1,2
  do sh6 = 1,2
  do sp7 = 1,2
     spin_itm_hphp = tdcc_spin_dummy2(sh2,sp8,sh6,sp7)
     spin_int2x = tdcc_spin_int2x(sh6,sp7,sp1,sp8)
     if(spin_itm_hphp * spin_int2x == 0) cycle

     itm_hphp(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact) = czero
     call bccd_l1p_13_1(sh2,sp8,sh6,sp7,itm_hphp)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h2 = 1,norb1
     do p1 = norb1+1,nact
     do p8 = norb1+1,nact
     do h6 = 1,norb1
     do p7 = norb1+1,nact
        i0(h2,p1) = i0(h2,p1) + fact * itm_hphp(h2,p8,h6,p7) * int2x(h6,p7,p1,p8,spin_int2x)
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
end subroutine bccd_l1p_13
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccd_l1p_13_1(sh2,sp8,sh6,sp7,i1)

!     i1 ( h2 p8 h6 p7 )_yt + = 1 * Sum ( h5 p3 ) * t ( p3 p8 h5 h6 )_t * y ( h2 h5 p3 p7 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : g2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp8,sh6,sp7
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p8,h6,p7
  integer(c_int) :: h5,p3,sh5,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_g2inp
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp8,sh5,sh6)
     spin_g2inp = tdcc_spin_g2inp(sh2,sh5,sp3,sp7)
     if(spin_t2inp * spin_g2inp == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h2 = 1,norb1
     do p8 = norb1+1,nact
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do h5 = 1,norb1
     do p3 = norb1+1,nact
        i1(h2,p8,h6,p7) = i1(h2,p8,h6,p7) + fact * t2inp(p3,p8,h5,h6,spin_t2inp) * g2inp(h2,h5,p3,p7,spin_g2inp)
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
end subroutine bccd_l1p_13_1
!##########################################################
!##########################################################
!##########################################################
