!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_t1p_1(sp2,sh1,i0)

! i0 ( p2 h1 )_f + = 1 * f ( p2 h1 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sp2,sh1)
  !$omp parallel default(shared)
  !$omp do collapse(2)
  do p2 = norb1+1,nact
  do h1 = 1,norb1
     i0(p2,h1) = i0(p2,h1) + fact * fock(p2,h1,spin_fock)
  end do
  end do
  !$omp end do
  !$omp end parallel
end subroutine bccdt_t1p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_t1p_2(sp2,sh1,i0)

! i0 ( p2 h1 )_tf + = 1 * Sum ( p7 h8 ) * t ( p2 p7 h1 h8 )_t * f ( h8 p7 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : fock

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: p7,h8,sp7,sh8
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = 1.0d+0 * runit

  do sp7 = 1,2
  do sh8 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp2,sp7,sh1,sh8)
     spin_fock = tdcc_spin_fock(sh8,sp7)
     if(spin_t2inp * spin_fock == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do p7 = norb1+1,nact
     do h8 = 1,norb1
        i0(p2,h1) = i0(p2,h1) + fact * t2inp(p2,p7,h1,h8,spin_t2inp) * fock(h8,p7,spin_fock)
     end do
     end do
     end do
     end do
     !$omp end do
     !$omp end parallel
  end do
  end do
end subroutine bccdt_t1p_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_t1p_3(sp2,sh1,i0)

! i0 ( p2 h1 )_vt + = -1/2 * Sum ( h4 h5 p3 ) * t ( p2 p3 h4 h5 )_t * v ( h4 h5 h1 p3 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: h4,h5,p3,sh4,sh5,sp3
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh4 = 1,2
  do sh5 = 1,2
  do sp3 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp2,sp3,sh4,sh5)
     spin_int2x = tdcc_spin_int2x(sh4,sh5,sh1,sp3)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p3 = norb1+1,nact
        i0(p2,h1) = i0(p2,h1) + fact * t2inp(p2,p3,h4,h5,spin_t2inp) * int2x(h4,h5,h1,p3,spin_int2x)
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
end subroutine bccdt_t1p_3
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_t1p_4(sp2,sh1,i0)

! i0 ( p2 h1 )_vt + = -1/2 * Sum ( h5 p3 p4 ) * t ( p3 p4 h1 h5 )_t * v ( h5 p2 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: h5,p3,p4,sh5,sp3,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh5 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp3,sp4,sh1,sh5)
     spin_int2x = tdcc_spin_int2x(sh5,sp2,sp3,sp4)
     if(spin_t2inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p2 = norb1+1,nact
     do h1 = 1,norb1
        do h5 = 1,norb1
        do p3 = norb1+1,nact
        do p4 = norb1+1,nact
           i0(p2,h1) = i0(p2,h1) + fact * &
                t2inp(p3,p4,h1,h5,spin_t2inp) * int2x(h5,p2,p3,p4,spin_int2x)
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
end subroutine bccdt_t1p_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine bccdt_t1p_5(sp2,sh1,i0)

! i0 ( p2 h1 )_vt + = 1/4 * Sum ( h5 h6 p3 p4 ) * t ( p2 p3 p4 h1 h5 h6 )_t * v ( h5 h6 p3 p4 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sp2,sh1
  complex(kind(0d0)),intent(inout) :: i0((norb1+1):nact,1:norb1)
  integer(c_int) :: p2,h1
  integer(c_int) :: h5,h6,p3,p4,sh5,sh6,sp3,sp4
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = 1.0d+0 / 4.0d+0 * runit

  do sh5 = 1,2
  do sh6 = 1,2
  do sp3 = 1,2
  do sp4 = 1,2
     spin_t3inp = tdcc_spin_t3inp(sp2,sp3,sp4,sh1,sh5,sh6)
     spin_int2x = tdcc_spin_int2x(sh5,sh6,sp3,sp4)
     if(spin_t3inp * spin_int2x == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do p2 = norb1+1,nact
     do h1 = 1,norb1
     do h5 = 1,norb1
     do h6 = 1,norb1
     do p3 = norb1+1,nact
     do p4 = norb1+1,nact
        i0(p2,h1) = i0(p2,h1) &
             + fact*t3inp(p2,p3,p4,h1,h5,h6,spin_t3inp)*int2x(h5,h6,p3,p4,spin_int2x)
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
end subroutine bccdt_t1p_5
!##########################################################
!##########################################################
!##########################################################
