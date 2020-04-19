!**********************************************************
subroutine ccdt_t2p_main()

  use mod_ormas, only : nact,nelact
  use mod_cc, only : norb1,cc_code,t2out
  use mod_cc2

  implicit none
!debug
!write(6,"('ccdt_t2p_main: skip t2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (nelact(1)>=2.and.norb1 >= 2) then
        call ccdt_t2p_1(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2p_2(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2p_3(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2p_4(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2p_5(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2p_6(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2p_7(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2p_8(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2p_9(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
     end if
     if (nelact(2)>=1) then
        call ccdt_t2p_1(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
        call ccdt_t2p_2(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
        call ccdt_t2p_3(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
        call ccdt_t2p_4(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
        call ccdt_t2p_5(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
        call ccdt_t2p_6(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
        call ccdt_t2p_7(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
        call ccdt_t2p_8(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
        call ccdt_t2p_9(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
     end if
  else
     call ccdt_t2p_man01(t2out,cc_work1,cc_work2,cc_work3)
     call ccdt_t2p_man02(t2out,cc_work1,cc_work2,cc_work3)
     call ccdt_t2p_man03(t2out,cc_work1,cc_work2,cc_work3)
     call ccdt_t2p_man04(t2out,cc_work1,cc_work2,cc_work3)
  end if

end subroutine ccdt_t2p_main
!**********************************************************
