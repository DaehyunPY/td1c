!################################################################################
subroutine tdcc_print(cic)

  use, intrinsic :: iso_c_binding
  use mod_cc, only : cc_rank,t0inp,g0inp,t1inp,g1inp,t2inp,g2inp,t3inp,g3inp

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)

  call tdcc_getcc(cic)
  write(6, "('CC0 amplitude: ',4e14.5)") t0inp,g0inp
  if (cc_rank >= 1) call tdcc_print_cc1(t1inp,g1inp)
  if (cc_rank >= 2) call tdcc_print_cc2(t2inp,g2inp)
  if (cc_rank >= 3) call tdcc_print_cc3(t3inp,g3inp)

end subroutine tdcc_print
!######################################################################
subroutine tdcc_print_cc1(tamp, gamp)

  use mod_ormas, only : nact
  use mod_cc, only : norb1, fock

  implicit none
  complex(kind(0d0)), intent(in) :: tamp((norb1+1):nact,1:norb1,1:2)
  complex(kind(0d0)), intent(in) :: gamp(1:norb1,(norb1+1):nact,1:2)

  integer(c_int) :: iact, aact

  write(6, "('CC1 amplitude: aa')")
  do iact = 1, norb1
  do aact = norb1+1, nact
     write(6, "(2i5,4e14.5)") aact,iact,&
          tamp(aact,iact,1), &
          gamp(iact,aact,1)
  end do
  end do

end subroutine tdcc_print_cc1
!######################################################################
subroutine tdcc_print_cc2(tamp, gamp)

  use mod_ormas, only : nact
  use mod_cc, only : norb1, int2x

  implicit none
  complex(kind(0d0)), intent(in) :: tamp((norb1+1):nact,(norb1+1):nact,1:norb1,1:norb1,1:6)
  complex(kind(0d0)), intent(in) :: gamp(1:norb1,1:norb1,(norb1+1):nact,(norb1+1):nact,1:6)

  integer(c_int) :: iact, jact, aact, bact

  write(6, "('CC2 amplitude: aaaa')")
  do iact = 1, norb1
  do jact = 1, norb1
  do aact = norb1+1, nact
  do bact = norb1+1, nact
     write(6, "(4i5,4e14.5)") aact,bact,iact,jact,&
          tamp(aact,bact,iact,jact,1), &
          gamp(iact,jact,aact,bact,1)
  end do
  end do
  end do
  end do

  write(6, "('CC2 amplitude: abab')")
  do iact = 1, norb1
  do jact = 1, norb1
  do aact = norb1+1, nact
  do bact = norb1+1, nact
     write(6, "(4i5,4e14.5)") aact,bact,iact,jact,&
          tamp(aact,bact,iact,jact,3), &
          gamp(iact,jact,aact,bact,3)
  end do
  end do
  end do
  end do

  !redundant
  !write(6, "('CC2 amplitude: abba')")
  !do iact = 1, norb1
  !do jact = 1, norb1
  !do aact = norb1+1, nact
  !do bact = norb1+1, nact
  !   write(6, "(4i5,4e14.5)") aact,bact,iact,jact,&
  !        tamp(aact,bact,iact,jact,5), &
  !        gamp(iact,jact,aact,bact,5)
  !end do
  !end do
  !end do
  !end do
  !redundant

end subroutine tdcc_print_cc2
!######################################################################
subroutine tdcc_print_cc3(tamp, gamp)

  use mod_ormas, only : nact,act1_ll,act1_ul
  use mod_cc, only : norb1

  implicit none
  complex(kind(0d0)), intent(in) :: &
       tamp((norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,&
       act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,1:20)
  complex(kind(0d0)), intent(in) :: &
       gamp(act1_ll:norb1,act1_ll:norb1,act1_ll:norb1,&
       (norb1+1):act1_ul,(norb1+1):act1_ul,(norb1+1):act1_ul,1:20)

  integer(c_int) :: iact, jact, kact, aact, bact, cact

  write(6, "('CC3 amplitude: aaaaaa')")
  do iact = act1_ll, norb1
  do jact = act1_ll, norb1
  do kact = act1_ll, norb1
  do aact = norb1+1, act1_ul
  do bact = norb1+1, act1_ul
  do cact = norb1+1, act1_ul
     write(6, "(6i5,4e14.5)") aact,bact,cact,iact,jact,kact,&
          tamp(aact,bact,cact,iact,jact,kact,1), &
          gamp(iact,jact,kact,aact,bact,cact,1)
  end do
  end do
  end do
  end do
  end do
  end do

  write(6, "('CC3 amplitude: aabaab')")
  do iact = act1_ll, norb1
  do jact = act1_ll, norb1
  do kact = act1_ll, norb1
  do aact = norb1+1, act1_ul
  do bact = norb1+1, act1_ul
  do cact = norb1+1, act1_ul
     write(6, "(6i5,4e14.5)") aact,bact,cact,iact,jact,kact,&
          tamp(aact,bact,cact,iact,jact,kact,3), &
          gamp(iact,jact,kact,aact,bact,cact,3)
  end do
  end do
  end do
  end do
  end do
  end do

end subroutine tdcc_print_cc3
!################################################################################
