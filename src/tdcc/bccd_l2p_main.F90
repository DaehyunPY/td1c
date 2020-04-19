!**********************************************************
subroutine bccd_l2p_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,nelact
  use mod_cc,only : norb1,cc_code,g2out
  use mod_cc2

  implicit none
!debug
!write(6,"('ccdt_l2p_main: skip l2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (nelact(1)>=2) then
        call bccd_l2p_1(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_2(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_3(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_4(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_5(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_6(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_7(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_8(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_9(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_10(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_11(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call bccd_l2p_12(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
     end if
     if (nelact(2)>=1) then
        call bccd_l2p_1(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_2(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_3(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_4(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_5(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_6(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_7(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_8(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_9(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_10(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_11(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
        call bccd_l2p_12(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
     end if
  else
     call ccdt_l2p_man03(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man04(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man05(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man06(g2out,cc_work1,cc_work2,cc_work3) ! 1 & 9
     call ccdt_l2p_man10(g2out,cc_work1,cc_work2,cc_work3)
     call ccdt_l2p_man11(g2out,cc_work1,cc_work2,cc_work3)
     ! L1 contribution
     call bccdt_l2p_man01(g2out,cc_work1,cc_work2,cc_work3)
  end if

end subroutine bccd_l2p_main
!**********************************************************
