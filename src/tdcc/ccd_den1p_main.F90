!**********************************************************
subroutine ccd_den1p_main(den)

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
     call ccd_den1_1(1,1,den)
     call ccd_den1_2(1,1,den)
  else
     call ccdt_den1p_man01(den,cc_work1,cc_work2,cc_work3)
     call ccdt_den1p_man02(den,cc_work1,cc_work2,cc_work3)
    !call ccdt_den1p_man03(den,cc_work1,cc_work2,cc_work3) ! T3
    !call ccdt_den1p_man04(den,cc_work1,cc_work2,cc_work3) ! L3
    !call ccdt_den1p_man05(den,cc_work1,cc_work2,cc_work3) ! T3,L3
    !call ccdt_den1p_man06(den,cc_work1,cc_work2,cc_work3) ! T3,L3
  end if

end subroutine ccd_den1p_main
!**********************************************************
