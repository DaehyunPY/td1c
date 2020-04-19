!**********************************************************
subroutine ccd_den2p_main(den)

  use mod_bas, only : smul
  use mod_ormas,only : nact,den2_abonly,nelact
  use mod_cc,only : cc_code
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: den(1:nact,1:nact,1:nact,1:nact,1:*)
!debug
!write(6,"('ccdt_den2p_main: skip den2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (.not.(smul==1.and.den2_abonly).and.nelact(1)>=2) then
        call ccd_den2_1(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_2(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_3(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_4(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_5(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_6(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_7(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_8(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_9(1,1,1,1,den(1,1,1,1,1))
     end if
     if (nelact(2)>=1) then
        call ccd_den2_1(1,2,1,2,den(1,1,1,1,3));call ccd_den2_1(1,2,2,1,den(1,1,1,1,5))
        call ccd_den2_2(1,2,1,2,den(1,1,1,1,3));call ccd_den2_2(1,2,2,1,den(1,1,1,1,5))
        call ccd_den2_3(1,2,1,2,den(1,1,1,1,3));call ccd_den2_3(1,2,2,1,den(1,1,1,1,5))
        call ccd_den2_4(1,2,1,2,den(1,1,1,1,3));call ccd_den2_4(1,2,2,1,den(1,1,1,1,5))
        call ccd_den2_5(1,2,1,2,den(1,1,1,1,3));call ccd_den2_5(1,2,2,1,den(1,1,1,1,5))
        call ccd_den2_6(1,2,1,2,den(1,1,1,1,3));call ccd_den2_6(1,2,2,1,den(1,1,1,1,5))
        call ccd_den2_7(1,2,1,2,den(1,1,1,1,3));call ccd_den2_7(1,2,2,1,den(1,1,1,1,5))
        call ccd_den2_8(1,2,1,2,den(1,1,1,1,3));call ccd_den2_8(1,2,2,1,den(1,1,1,1,5))
        call ccd_den2_9(1,2,1,2,den(1,1,1,1,3));call ccd_den2_9(1,2,2,1,den(1,1,1,1,5))
     end if
  else
     if (.not.(smul==1.and.den2_abonly).and.nelact(1)>=2) then
        call ccd_den2_1(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_2(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_3(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_4(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_5(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_6(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_7(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_8(1,1,1,1,den(1,1,1,1,1))
        call ccd_den2_9(1,1,1,1,den(1,1,1,1,1))
     else
        !write(6,"('ccdt_den2p_main: skip den2aa.')")
     end if
     if (nelact(2)>=1) then
        call ccdt_den2p_man01(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        call ccdt_den2p_man02(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        call ccdt_den2p_man03(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        call ccdt_den2p_man04(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        call ccdt_den2p_man05(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        call ccdt_den2p_man06(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        call ccdt_den2p_man07(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        call ccdt_den2p_man08(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        call ccdt_den2p_man09(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man10(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man11(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man12(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man13(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man14(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man15(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man16(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man17(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man18(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man19(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man20(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man21(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man22(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man23(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man24(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man25(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man26(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
        !call ccdt_den2p_man27(den(1,1,1,1,3),den(1,1,1,1,5),cc_work1,cc_work2,cc_work3)
     end if
  end if

end subroutine ccd_den2p_main
!**********************************************************
