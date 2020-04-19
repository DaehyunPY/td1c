!**********************************************************
subroutine ccdt_l2p_1e_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,cc_code,g2out
  use mod_cc2

  implicit none
!debug
!write(6,"('ccdt_l2p_1e_main: skip t2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (norb1 >= 2) then
        call ccdt_l2_1e_1(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call ccdt_l2_1e_2(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call ccdt_l2_1e_3(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
        call ccdt_l2_1e_4(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1))
     end if
     call ccdt_l2_1e_1(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
     call ccdt_l2_1e_2(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
     call ccdt_l2_1e_3(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
     call ccdt_l2_1e_4(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  else
     call ccdt_l2p_1e_man01(g2out,cc_work1,cc_work2)
     call ccdt_l2p_1e_man02(g2out,cc_work1,cc_work2)
     call ccdt_l2p_1e_man03(g2out,cc_work1,cc_work2)
     call ccdt_l2p_1e_man04(g2out,cc_work1,cc_work2)
!     call ccdt_l2_1e_1(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1)); 
!     call ccdt_l2_1e_1(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
!     call ccdt_l2_1e_2(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1)); 
!     call ccdt_l2_1e_2(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
!     call ccdt_l2_1e_3(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1)); 
!     call ccdt_l2_1e_3(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
!     call ccdt_l2_1e_4(1,1,1,1,g2out(1,1,norb1+1,norb1+1,1)); 
!     call ccdt_l2_1e_4(1,2,1,2,g2out(1,1,norb1+1,norb1+1,2))
  end if

end subroutine ccdt_l2p_1e_main
!**********************************************************
