!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_bmat3p_1(sh1,sp2,i0)

!  i0 ( h1 p2 )_dyt + = +1/4 * Sum ( h3 h4 p5 p6 ) * dy ( h3 h4 p5 p6 )_dy * t ( p2 p5 p6 h1 h3 h4 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, dg2inp
  use mod_cc,only : t3inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: h1,p2
  integer(c_int) :: h3,h4,p5,p6,sh3,sh4,sp5,sp6
  integer(c_int) :: spin_dg2inp_1
  integer(c_int) :: spin_t3inp_2
  integer(c_int),external :: tdcc_spin_g2inp
  integer(c_int),external :: tdcc_spin_t3inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 / 4.0d+0 * runit

  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
  do sp6 = 1,2
     spin_dg2inp_1 = tdcc_spin_g2inp(sh3,sh4,sp5,sp6)
     spin_t3inp_2 = tdcc_spin_t3inp(sp2,sp5,sp6,sh1,sh3,sh4)
     if(spin_dg2inp_1 * spin_t3inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do p6 = norb1+1,nact
        i0(h1,p2) = i0(h1,p2) + fact * &
             dg2inp(h3,h4,p5,p6,spin_dg2inp_1) * t3inp(p2,p5,p6,h1,h3,h4,spin_t3inp_2)
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
end subroutine ccdt_bmat3p_1
!##########################################################
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_bmat3p_2(sh1,sp2,i0)

!  i0 ( h1 p2 )_dytt + = -1/4 * Sum ( h3 h4 p5 ) * i1 ( h3 h4 p5 h1 )_dyt * t ( p2 p5 h3 h4 )_t 1

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, t2inp

  implicit none
  integer(c_int),intent(in) :: sh1,sp2
  complex(kind(0d0)),intent(inout) :: i0(1:nact,1:nact)
  integer(c_int) :: h1,p2
  integer(c_int) :: h3,h4,p5,sh3,sh4,sp5
  integer(c_int) :: spin_itm_hhph_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_dummy2
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),allocatable :: itm_hhph(:,:,:,:)
  complex(kind(0d0)),parameter :: fact = -1.0d+0 / 4.0d+0 * runit

  allocate(itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1))
  do sh3 = 1,2
  do sh4 = 1,2
  do sp5 = 1,2
     spin_itm_hhph_1 = tdcc_spin_dummy2(sh3,sh4,sp5,sh1)
     spin_t2inp_2 = tdcc_spin_t2inp(sp2,sp5,sh3,sh4)
     if(spin_itm_hhph_1 * spin_t2inp_2 == 0) cycle

     itm_hhph(1:norb1,1:norb1,(norb1+1):nact,1:norb1) = czero
     call ccdt_bmat3p_2_1(sh3,sh4,sp5,sh1,itm_hhph)

     !$omp parallel default(shared)
     !$omp do collapse(2)
     do h1 = 1,norb1
     do p2 = norb1+1,nact
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
        i0(h1,p2) = i0(h1,p2) + fact * &
             itm_hhph(h3,h4,p5,h1) * t2inp(p2,p5,h3,h4,spin_t2inp_2)
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
  deallocate(itm_hhph)
end subroutine ccdt_bmat3p_2
!##########################################################
!##########################################################
!##########################################################
subroutine ccdt_bmat3p_2_1(sh3,sh4,sp5,sh1,i1)

!      i1 ( h3 h4 p5 h1 )_dyt + = +1 * Sum ( h6 p7 p8 ) * dy ( h6 h3 h4 p7 p8 p5 )_dy * t ( p8 p7 h1 h6 )_t 0

  use, intrinsic :: iso_c_binding
  use mod_const,only : czero,runit
  use mod_ormas,only : nact
  use mod_cc,only : norb1, dg3inp
  use mod_cc,only : t2inp

  implicit none
  integer(c_int),intent(in) :: sh3,sh4,sp5,sh1
  complex(kind(0d0)),intent(inout) :: i1(1:norb1,1:norb1,(norb1+1):nact,1:norb1)
  integer(c_int) :: h3,h4,p5,h1
  integer(c_int) :: h6,p7,p8,sh6,sp7,sp8
  integer(c_int) :: spin_dg3inp_1
  integer(c_int) :: spin_t2inp_2
  integer(c_int),external :: tdcc_spin_g3inp
  integer(c_int),external :: tdcc_spin_t2inp
  complex(kind(0d0)),parameter :: fact = +1.0d+0 * runit

  do sh6 = 1,2
  do sp7 = 1,2
  do sp8 = 1,2
     spin_dg3inp_1 = tdcc_spin_g3inp(sh6,sh3,sh4,sp7,sp8,sp5)
     spin_t2inp_2 = tdcc_spin_t2inp(sp8,sp7,sh1,sh6)
     if(spin_dg3inp_1 * spin_t2inp_2 == 0) cycle

     !$omp parallel default(shared)
     !$omp do collapse(4)
     do h3 = 1,norb1
     do h4 = 1,norb1
     do p5 = norb1+1,nact
     do h1 = 1,norb1
     do h6 = 1,norb1
     do p7 = norb1+1,nact
     do p8 = norb1+1,nact
!        i1(h3,h4,p5,h1) = i1(h3,h4,p5,h1) + fact * &
!        dg3inp(h6,h3,h4,p7,p8,p5,spin_dg3inp_1) * t2inp(p8,p7,h1,h6,spin_t2inp_2)
        i1(h3,h4,p5,h1) = i1(h3,h4,p5,h1) + &
             dg3inp(h6,h3,h4,p7,p8,p5,spin_dg3inp_1) * t2inp(p8,p7,h1,h6,spin_t2inp_2)
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
end subroutine ccdt_bmat3p_2_1
!##########################################################
!##########################################################
!##########################################################
