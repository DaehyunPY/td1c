!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_ene_1(i0)

! i0 ( )_tf + = +1 * Sum ( p5 h6 ) * t ( p5 h6 )_t * i1 ( h6 p5 )_f 2

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  complex(kind(0d0)),intent(inout) :: i0
  integer(c_int) :: p5,h6,sp5,sh6
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_itm_hp
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),allocatable :: itm_hp(:,:)
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  allocate(itm_hp(1:norb1,(norb1+1):nact))
  do sp5 = 1,2
  do sh6 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp5,sh6)
     spin_itm_hp = tdcc_spin_fock(sh6,sp5)
     if(spin_t1inp * spin_itm_hp == 0) cycle

     itm_hp(1:norb1,(norb1+1):nact) = czero
     call ccsd_ene_1_1(sh6,sp5,itm_hp)
     call ccsd_ene_1_2(sh6,sp5,itm_hp)

     do p5 = norb1+1,nact
     do h6 = 1,norb1
        i0 = i0 + fact * t1inp(p5,h6,spin_t1inp) * itm_hp(h6,p5)
     end do
     end do
  end do
  end do
  deallocate(itm_hp)
end subroutine ccsd_ene_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_ene_1_1(sh6,sp5,i1)

!     i1 ( h6 p5 )_f + = +1 * f ( h6 p5 )_f 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, fock

  implicit none
  integer(c_int),intent(in) :: sh6,sp5
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,p5
  integer(c_int) :: sdum
  integer(c_int) :: spin_fock
  integer(c_int),external :: tdcc_spin_fock
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  spin_fock = tdcc_spin_fock(sh6,sp5)
  do h6 = 1,norb1
  do p5 = norb1+1,nact
     i1(h6,p5) = i1(h6,p5) + fact * fock(h6,p5,spin_fock)
  end do
  end do
end subroutine ccsd_ene_1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_ene_1_2(sh6,sp5,i1)

!     i1 ( h6 p5 )_vt + = +1/2 * Sum ( h4 p3 ) * t ( p3 h4 )_t * v ( h4 h6 p3 p5 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp
  use mod_cc,only : int2x

  implicit none
  integer(c_int),intent(in) :: sh6,sp5
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,(norb1+1):nact)
  integer(c_int) :: h6,p5
  integer(c_int) :: h4,p3,sh4,sp3
  integer(c_int) :: spin_t1inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t1inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  do sh4 = 1,2
  do sp3 = 1,2
     spin_t1inp = tdcc_spin_t1inp(sp3,sh4)
     spin_int2x = tdcc_spin_int2x(sh4,sh6,sp3,sp5)
     if(spin_t1inp * spin_int2x == 0) cycle

     do h6 = 1,norb1
     do p5 = norb1+1,nact
     do h4 = 1,norb1
     do p3 = norb1+1,nact
        i1(h6,p5) = i1(h6,p5) + fact * t1inp(p3,h4,spin_t1inp) * int2x(h4,h6,p3,p5,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_ene_1_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_ene_2(i0)

! i0 ( )_vt + = +1/4 * Sum ( h3 h4 p1 p2 ) * t ( p1 p2 h3 h4 )_t * v ( h3 h4 p1 p2 )_v 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp
  use mod_cc,only : int2x

  implicit none
  complex(kind(0d0)),intent(inout) :: i0
  integer(c_int) :: h3,h4,p1,p2,sh3,sh4,sp1,sp2
  integer(c_int) :: spin_t2inp
  integer(c_int) :: spin_int2x
  integer(c_int),external :: tdcc_spin_t2inp
  integer(c_int),external :: tdcc_spin_int2x
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  do sh3 = 1,2
  do sh4 = 1,2
  do sp1 = 1,2
  do sp2 = 1,2
     spin_t2inp = tdcc_spin_t2inp(sp1,sp2,sh3,sh4)
     spin_int2x = tdcc_spin_int2x(sh3,sh4,sp1,sp2)
     if(spin_t2inp * spin_int2x == 0) cycle

     do h3 = 1,norb1
     do h4 = 1,norb1
     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
        i0 = i0 + fact * t2inp(p1,p2,h3,h4,spin_t2inp) * int2x(h3,h4,p1,p2,spin_int2x)
     end do
     end do
     end do
     end do
  end do
  end do
  end do
  end do
end subroutine ccsd_ene_2
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
complex(kind(0d0)) function ccsd_ene_main()

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero

  implicit none
  complex(kind(0d0)) :: cene

  cene = czero
  call ccsd_ene_1(cene)
  call ccsd_ene_2(cene)
  ccsd_ene_main = cene

end function ccsd_ene_main
!**********************************************************
