!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den1_1(sp1,sp2,i0)

!  i0 ( p1 p2 )_yt + = +1/2 * Sum ( h3 h4 p5 ) * y ( h3 h4 p5 p2 )_y * t ( p5 p1 h3 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: p1,p2
  integer(c_int) :: h3,h4,p5,sh3,sh4,sp5
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 2.0d+0 * runit

  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh3,sh4,sp5,sp2)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp1,sh3,sh4)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

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
  end do
  end do
  end do
end subroutine ccd_den1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccd_den1_2(sh1,sh2,i0)

!  i0 ( h1 h2 )_yt + = -1/2 * Sum ( h3 p4 p5 ) * y ( h3 h1 p4 p5 )_y * t ( p4 p5 h3 h2 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: h1,h2
  integer(c_int) :: h3,p4,p5,sh3,sp4,sp5
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  do sh3 = 1,2
  do sp4 = 1,2
  do sp5 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh3,sh1,sp4,sp5)
     spin_t2inp_2 = tdcc_spin_t2inp(sp4,sp5,sh3,sh2)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

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
  end do
  end do
  end do
end subroutine ccd_den1_2
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccd_den1_main(den)

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact

  implicit none
  complex(kind(0d0)),intent(inout) :: den(1:nact,1:nact,1:*)

  call ccd_den1_1(1,1,den(1,1,1))
  call ccd_den1_2(1,1,den(1,1,1))

end subroutine ccd_den1_main
!**********************************************************
