!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_1(sh1,sp2,i0)

!  i0 ( h1 p2 )_y + = +1 * y ( h1 p2 )_y 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: h1,p2
  integer(c_int) :: sdum
  integer(c_int) :: spin_g1inp_1
  integer(c_int),external :: tdcc_spin_g1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  spin_g1inp_1 = tdcc_spin_g1inp(sh1,sp2)
  do h1 = 1,norb1
  do p2 = norb1+1,nact
     i0(h1,p2) = i0(h1,p2) + fact * g1inp(h1,p2,spin_g1inp_1)
  end do
  end do
end subroutine ccsd_den1_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_2(sp1,sh2,i0)

!  i0 ( p1 h2 )_t + = +1 * t ( p1 h2 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: p1,h2
  integer(c_int) :: sdum
  integer(c_int) :: spin_t1inp_1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  spin_t1inp_1 = tdcc_spin_t1inp(sp1,sh2)
  do p1 = norb1+1,nact
  do h2 = 1,norb1
     i0(p1,h2) = i0(p1,h2) + fact * t1inp(p1,h2,spin_t1inp_1)
  end do
  end do
end subroutine ccsd_den1_2
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_3(sp1,sh2,i0)

!  i0 ( p1 h2 )_ytt + = -1 * Sum ( h3 ) * i1 ( h3 h2 )_yt * t ( p1 h3 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: p1,h2
  integer(c_int) :: h3,sh3
  integer(c_int) :: spin_itm_hh_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh3 = 1,2
     spin_itm_hh_1 = tdcc_spin_dummy1(sh3,sh2)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh3)
     if(spin_itm_hh_1 * spin_t1inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_den1_3_1(sh3,sh2,itm_hh)

     do p1 = norb1+1,nact
     do h2 = 1,norb1
     do h3 = 1,norb1
        i0(p1,h2) = i0(p1,h2) + fact * itm_hh(h3,h2) * t1inp(p1,h3,spin_t1inp_2)
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
end subroutine ccsd_den1_3
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_3_1(sh3,sh2,i1)

!      i1 ( h3 h2 )_yt + = +1 * Sum ( p4 ) * y ( h3 p4 )_y * t ( p4 h2 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h3,h2
  integer(c_int) :: p4,sp4
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sp4 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh3,sp4)
     spin_t1inp_2 = tdcc_spin_t1inp(sp4,sh2)
     if(spin_g1inp_1 * spin_t1inp_2 == 0) cycle

     do h3 = 1,norb1
     do h2 = 1,norb1
     do p4 = norb1+1,nact
        i1(h3,h2) = i1(h3,h2) + fact * g1inp(h3,p4,spin_g1inp_1) * t1inp(p4,h2,spin_t1inp_2)
     end do
     end do
     end do
  end do
end subroutine ccsd_den1_3_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_4(sp1,sp2,i0)

!  i0 ( p1 p2 )_yt + = +1 * Sum ( h3 ) * y ( h3 p2 )_y * t ( p1 h3 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: p1,p2
  integer(c_int) :: h3,sh3
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh3 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh3,sp2)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh3)
     if(spin_g1inp_1 * spin_t1inp_2 == 0) cycle

     do p1 = norb1+1,nact
     do p2 = norb1+1,nact
     do h3 = 1,norb1
        i0(p1,p2) = i0(p1,p2) + fact * g1inp(h3,p2,spin_g1inp_1) * t1inp(p1,h3,spin_t1inp_2)
     end do
     end do
     end do
  end do
end subroutine ccsd_den1_4
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_5(sh1,sh2,i0)

!  i0 ( h1 h2 )_yt + = -1 * Sum ( p3 ) * y ( h1 p3 )_y * t ( p3 h2 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t1inp

  implicit none
  integer(c_int),intent(in) :: sh1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: h1,h2
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),parameter :: fact = -1.0d+0 * runit

  do sp3 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh1,sp3)
     spin_t1inp_2 = tdcc_spin_t1inp(sp3,sh2)
     if(spin_g1inp_1 * spin_t1inp_2 == 0) cycle

     do h1 = 1,norb1
     do h2 = 1,norb1
     do p3 = norb1+1,nact
        i0(h1,h2) = i0(h1,h2) + fact * g1inp(h1,p3,spin_g1inp_1) * t1inp(p3,h2,spin_t1inp_2)
     end do
     end do
     end do
  end do
end subroutine ccsd_den1_5
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_6(sp1,sp2,i0)

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
end subroutine ccsd_den1_6
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_7(sh1,sh2,i0)

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
end subroutine ccsd_den1_7
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_8(sp1,sh2,i0)

!  i0 ( p1 h2 )_yt + = +1 * Sum ( h3 p4 ) * y ( h3 p4 )_y * t ( p4 p1 h3 h2 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g1inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: p1,h2
  integer(c_int) :: h3,p4,sh3,sp4
  integer(c_int) :: spin_g1inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g1inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh3 = 1,2
  do sp4 = 1,2
     spin_g1inp_1 = tdcc_spin_g1inp(sh3,sp4)
     spin_t2inp_2 = tdcc_spin_t2inp(sp4,sp1,sh3,sh2)
     if(spin_g1inp_1 * spin_t2inp_2 == 0) cycle

     do p1 = norb1+1,nact
     do h2 = 1,norb1
     do h3 = 1,norb1
     do p4 = norb1+1,nact
        i0(p1,h2) = i0(p1,h2) + fact * g1inp(h3,p4,spin_g1inp_1) * t2inp(p4,p1,h3,h2,spin_t2inp_2)
     end do
     end do
     end do
     end do
  end do
  end do
end subroutine ccsd_den1_8
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_9(sp1,sh2,i0)

!  i0 ( p1 h2 )_ytt + = -1/2 * Sum ( h3 ) * i1 ( h3 h2 )_yt * t ( p1 h3 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: p1,h2
  integer(c_int) :: h3,sh3
  integer(c_int) :: spin_itm_hh_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_hh(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_hh(1:norb1,1:norb1))
  do sh3 = 1,2
     spin_itm_hh_1 = tdcc_spin_dummy1(sh3,sh2)
     spin_t1inp_2 = tdcc_spin_t1inp(sp1,sh3)
     if(spin_itm_hh_1 * spin_t1inp_2 == 0) cycle

     itm_hh(1:norb1,1:norb1) = czero
     call ccsd_den1_9_1(sh3,sh2,itm_hh)

     do p1 = norb1+1,nact
     do h2 = 1,norb1
     do h3 = 1,norb1
        i0(p1,h2) = i0(p1,h2) + fact * itm_hh(h3,h2) * t1inp(p1,h3,spin_t1inp_2)
     end do
     end do
     end do
  end do
  deallocate(itm_hh)
end subroutine ccsd_den1_9
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_9_1(sh3,sh2,i1)

!      i1 ( h3 h2 )_yt + = +1 * Sum ( h4 p5 p6 ) * y ( h4 h3 p5 p6 )_y * t ( p5 p6 h4 h2 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh2
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1)
  integer(c_int) :: h3,h2
  integer(c_int) :: h4,p5,p6,sh4,sp5,sp6
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh4 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh4,sh3,sp5,sp6)
     spin_t2inp_2 = tdcc_spin_t2inp(sp5,sp6,sh4,sh2)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do h3 = 1,norb1
     do h2 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i1(h3,h2) = i1(h3,h2) + fact * g2inp(h4,h3,p5,p6,spin_g2inp_1) * t2inp(p5,p6,h4,h2,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_den1_9_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_10(sp1,sh2,i0)

!  i0 ( p1 h2 )_ytt + = -1/2 * Sum ( p3 ) * i1 ( p3 p1 )_yt * t ( p3 h2 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t1inp

  implicit none
  integer(c_int),intent(in) :: sp1,sh2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: p1,h2
  integer(c_int) :: p3,sp3
  integer(c_int) :: spin_itm_pp_1
  integer(c_int) :: spin_t1inp_2
  integer(c_int),external :: tdcc_spin_dummy1
  integer(c_int),external :: tdcc_spin_t1inp
  complex(kind(0d0)),allocatable :: itm_pp(:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 2.0d+0 * runit

  allocate(itm_pp((norb1+1):nact,(norb1+1):nact))
  do sp3 = 1,2
     spin_itm_pp_1 = tdcc_spin_dummy1(sp3,sp1)
     spin_t1inp_2 = tdcc_spin_t1inp(sp3,sh2)
     if(spin_itm_pp_1 * spin_t1inp_2 == 0) cycle

     itm_pp((norb1+1):nact,(norb1+1):nact) = czero
     call ccsd_den1_10_1(sp3,sp1,itm_pp)

     do p1 = norb1+1,nact
     do h2 = 1,norb1
     do p3 = norb1+1,nact
        i0(p1,h2) = i0(p1,h2) + fact * itm_pp(p3,p1) * t1inp(p3,h2,spin_t1inp_2)
     end do
     end do
     end do
  end do
  deallocate(itm_pp)
end subroutine ccsd_den1_10
!##########################################################
!##########################################################
!##########################################################
subroutine ccsd_den1_10_1(sp3,sp1,i1)

!      i1 ( p3 p1 )_yt + = +1 * Sum ( h4 h5 p6 ) * y ( h4 h5 p6 p3 )_y * t ( p6 p1 h4 h5 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, g2inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sp3,sp1
  complex(kind(0d0)),intent(inout) :: i1((norb1+1):nact,(norb1+1):nact)
  integer(c_int) :: p3,p1
  integer(c_int) :: h4,h5,p6,sh4,sh5,sp6
  integer(c_int) :: spin_g2inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh4 = 1,2
  do sh5 = 1,2
  do sp6 = 1,2
     spin_g2inp_1 = tdcc_spin_g2inp(sh4,sh5,sp6,sp3)
     spin_t2inp_2 = tdcc_spin_t2inp(sp6,sp1,sh4,sh5)
     if(spin_g2inp_1 * spin_t2inp_2 == 0) cycle

     do p3 = norb1+1,nact
     do p1 = norb1+1,nact
     do h4 = 1,norb1
     do h5 = 1,norb1
     do p6 = norb1+1,nact
        i1(p3,p1) = i1(p3,p1) + fact * g2inp(h4,h5,p6,p3,spin_g2inp_1) * t2inp(p6,p1,h4,h5,spin_t2inp_2)
     end do
     end do
     end do
     end do
     end do
  end do
  end do
  end do
end subroutine ccsd_den1_10_1
!##########################################################
!##########################################################
!##########################################################
!**********************************************************
subroutine ccsd_den1_main(den)

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact

  implicit none
  complex(kind(0d0)),intent(inout) :: den(1:nact,1:nact,1:*)

  call ccsd_den1_1(1,1,den(1,1,1))
  call ccsd_den1_2(1,1,den(1,1,1))
  call ccsd_den1_3(1,1,den(1,1,1))
  call ccsd_den1_4(1,1,den(1,1,1))
  call ccsd_den1_5(1,1,den(1,1,1))
  call ccsd_den1_6(1,1,den(1,1,1))
  call ccsd_den1_7(1,1,den(1,1,1))
  call ccsd_den1_8(1,1,den(1,1,1))
  call ccsd_den1_9(1,1,den(1,1,1))
  call ccsd_den1_10(1,1,den(1,1,1))

end subroutine ccsd_den1_main
!**********************************************************
