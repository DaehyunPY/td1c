!**********************************************************
subroutine ccd_t2p_main()

  use mod_ormas, only : nact,nelact
  use mod_cc, only : norb1,cc_code,t2out
  use mod_cc2

  implicit none
!debug
!write(6,"('ccd_t2p_main: skip t2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (nelact(1)>=2) then
        call ccd_t2_1(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
        call ccd_t2_2(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
        call ccd_t2_3(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
        call ccd_t2_4(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
        call ccd_t2_5(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
        call ccd_t2_6(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
     end if
     if (nelact(2)>=1) then
        call ccd_t2_1(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
        call ccd_t2_2(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
        call ccd_t2_3(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
        call ccd_t2_4(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
        call ccd_t2_5(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
        call ccd_t2_6(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
     end if
  else
     call ccdt_t2p_man01(t2out,cc_work1,cc_work2,cc_work3)
     call ccdt_t2p_man02(t2out,cc_work1,cc_work2,cc_work3)
     call ccdt_t2p_man03(t2out,cc_work1,cc_work2,cc_work3)
     call ccdt_t2p_man04(t2out,cc_work1,cc_work2,cc_work3)
!     call ccd_t2_1(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
!     call ccd_t2_2(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
!     call ccd_t2_3(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
!     call ccd_t2_4(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
!     call ccd_t2_5(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
!     call ccd_t2_6(1,1,1,1,t2out(norb1+1,norb1+1,1,1,1))
!   
!     call ccd_t2_1(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
!     call ccd_t2_2(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
!     call ccd_t2_3(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
!     call ccd_t2_4(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
!     call ccd_t2_5(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
!     call ccd_t2_6(1,2,1,2,t2out(norb1+1,norb1+1,1,1,2))
  end if

end subroutine ccd_t2p_main
!**********************************************************
