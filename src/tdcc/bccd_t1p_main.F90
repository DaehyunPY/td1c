!**********************************************************
subroutine bccd_t1p_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,cc_code,t1out
  use mod_cc2

  implicit none
!debug
!write(6,"('bccdt_t1p_main: skip t1p.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
!1
     call bccdt_t1p_1(1,1,t1out)
     call bccdt_t1p_2(1,1,t1out)
     call bccdt_t1p_3(1,1,t1out)
     call bccdt_t1p_4(1,1,t1out)
!2
! same as 1     call bccd_t1p_1(1,1,t1out)
! same as 1     call bccd_t1p_2(1,1,t1out)
! same as 1     call bccd_t1p_3(1,1,t1out)
! same as 1     call bccd_t1p_4(1,1,t1out)
!3
! same as 1     call ccsd_t1_1(1,1,t1out)
! same as 1!     call ccsd_t1_2(1,1,t1out) ! t1 term
! same as 1!     call ccsd_t1_3(1,1,t1out) ! t1 term
! same as 1!     call ccsd_t1_4(1,1,t1out) ! t1 term
! same as 1     call ccsd_t1_5(1,1,t1out) ! t1 terms have to be partially disabled
! same as 1     call ccsd_t1_6(1,1,t1out) ! t1 terms have to be partially disabled
! same as 1     call ccsd_t1_7(1,1,t1out)
  else
     call bccdt_t1p_man01(t1out)
  end if

end subroutine bccd_t1p_main
!**********************************************************
