!##########################################################
subroutine bccdt_l1p_16(sh2,sp1,i0)

!16: i0 ( i a )_vty + = -1/2 * Sum ( b c j k l ) * y ( k j l a c b )_y * i1 ( i c b k j l )_vt 4

!16_1: i1 ( i a b j k l )_vt + = 1 * Sum ( c ) * t ( c b j k )_t * i2 ( i a l c )_v 2
!16_2: i1 ( i a b j k l )_ft + = 1/6 * Sum ( c ) * t ( c a b j k l )_t * f ( i c )_f 0
!16_3: i1 ( i a b j k l )_vt + = 1/2 * Sum ( m c ) * t ( c a b m j k )_t * v ( i m l c )_v 0
!16_4: i1 ( i a b j k l )_vt + = 1/6 * Sum ( c d ) * t ( c d b j k l )_t * v ( i a c d )_v 0

!16_1_1: i2 ( i a j b )_v + = 1 * v ( i a j b )_v 0
!16_1_2: i2 ( i a j b )_vt + = -1 * Sum ( k c ) * t ( c a k j )_t * v ( i k b c )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g3inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp1
  complex(kind(0d0)),intent(inout) :: i0(1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p1
  integer(c_int) :: p15,p11,h14,h13,h16,sp15,sp11,sh14,sh13,sh16
  integer(c_int) :: spin_g3inp
  integer(c_int) :: spin_itm_hpphhh
  integer(c_int),external :: tdcc_spin_g3inp
  integer(c_int),external :: tdcc_spin_dummy3
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

!16_1: i1 ( i a b j k l )_vt + = 1 * Sum ( c ) * t ( c b j k )_t * i2 ( i a l c )_v 2
!16_2: i1 ( i a b j k l )_ft + = 1/6 * Sum ( c ) * t ( c a b j k l )_t * f ( i c )_f 0
!16_3: i1 ( i a b j k l )_vt + = 1/2 * Sum ( m c ) * t ( c a b m j k )_t * v ( i m l c )_v 0
!16_4: i1 ( i a b j k l )_vt + = 1/6 * Sum ( c d ) * t ( c d b j k l )_t * v ( i a c d )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sh2,sp11,sp15,sh13,sh14,sh16
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,p11,p15,h13,h14,h16
  integer(c_int) :: p4,sp4
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_itm_hphp
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
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

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp11,sh16,sp4
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p11,h16,p4
  integer(c_int) :: sdum
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_int2x
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

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp11,sh16,sp4
  complex(kind(0d0)),intent(inout) :: i2(1:norb1,(norb1+1):nact,1:norb1,(norb1+1):nact)
  integer(c_int) :: h2,p11,h16,p4
  integer(c_int) :: h10,p8,sh10,sp8
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
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

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : fock

  implicit none
  integer(c_int),intent(in) :: sh2,sp11,sp15,sh13,sh14,sh16
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,p11,p15,h13,h14,h16
  integer(c_int) :: p4,sp4
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_fock
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

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp11,sp15,sh13,sh14,sh16
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,p11,p15,h13,h14,h16
  integer(c_int) :: h8,p5,sh8,sp5
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
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

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t3inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh2,sp11,sp15,sh13,sh14,sh16
  complex(kind(0d0)),intent(inout) :: &
       i1(1:norb1,(norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:norb1)
  integer(c_int) :: h2,p11,p15,h13,h14,h16
  integer(c_int) :: p4,p5,sp4,sp5
  integer(c_int) :: spin_t3inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t3inp
  integer(c_int),external :: tdcc_spin_int2x
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
