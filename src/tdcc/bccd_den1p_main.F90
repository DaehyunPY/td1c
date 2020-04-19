!**********************************************************
subroutine bccd_den1p_main(den)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : cc_code
  use mod_cc2

  implicit none
  complex(kind(0d0)), intent(inout) :: den(*)
!debug
!write(6,"('ccd_den1p_main: skip den1.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     call bccd_den1p_1(1,1,den)
     call bccd_den1p_2(1,1,den)
     call bccd_den1p_3(1,1,den)
     call bccd_den1p_4(1,1,den)
  else
     call ccdt_den1p_man01(den,cc_work1,cc_work2)
     call ccdt_den1p_man02(den,cc_work1,cc_work2)
     ! L1 contribution
     call bccdt_den1p_man01(den,cc_work1,cc_work2)
  end if

end subroutine bccd_den1p_main
!**********************************************************
