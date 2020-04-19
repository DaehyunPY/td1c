!**********************************************************
subroutine bccdt_l2p_main()

  use mod_ormas,only : nact,nelact
  use mod_cc,only : norb1,cc_code,g2out
  use mod_cc2

  implicit none
!debug
!write(6,"('ccdt_l2p_main: skip l2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (nelact(1)>=2.and.norb1 >= 2) then
        call bccdt_l2p_1(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_2(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_3(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_4(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_5(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_6(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_7(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_8(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_9(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_10(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_11(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_12(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_13(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_14(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_15(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_16(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_17(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_18(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_19(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_20(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccdt_l2p_21(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
     end if
     if (nelact(2)>=1) then
        call bccdt_l2p_1(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_2(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_3(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_4(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_5(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_6(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_7(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_8(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_9(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_10(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_11(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_12(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_13(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_14(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_15(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_16(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_17(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_18(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_19(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_20(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccdt_l2p_21(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
     end if
  else
    !call ccdt_l2p_man01(g2out,cc_work1,cc_work2,cc_work3) ! merged with 6
    !call ccdt_l2p_man02(g2out,cc_work1,cc_work2,cc_work3) ! merged with 3
     call ccdt_l2p_man03(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man04(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man05(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man06(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man07(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man08(g2out,cc_work1,cc_work2,cc_work3)
    !call ccdt_l2p_man09(g2out,cc_work1,cc_work2,cc_work3) ! merged with 10
     call ccdt_l2p_man10(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man11(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man12(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man13(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man14(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man15(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man16(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man17(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man18(g2out,cc_work1,cc_work2,cc_work3)
     ! L1 contributions
     call bccdt_l2p_man01(g2out,cc_work1,cc_work2,cc_work3)
  end if

end subroutine bccdt_l2p_main
!**********************************************************
