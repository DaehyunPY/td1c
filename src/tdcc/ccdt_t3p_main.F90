!**********************************************************
subroutine ccdt_t3p_main()

  use mod_ormas,only : nact,nelact
  use mod_cc,only : norb1,cc_code,t3out
  use mod_cc2

  implicit none
!debug
!write(6,"('ccdt_t3p_main: skip t3.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (nelact(1)>=3.and.norb1>=3) then
        call ccdt_t3p_1(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
        call ccdt_t3p_2(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
        call ccdt_t3p_3(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
        call ccdt_t3p_4(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
        call ccdt_t3p_5(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
        call ccdt_t3p_6(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
        call ccdt_t3p_7(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
        call ccdt_t3p_diagram8_2(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
     end if
     if (nelact(1)>=2.and.nelact(2)>=1.and.norb1>=2) then
        call ccdt_t3p_1(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
        call ccdt_t3p_2(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
        call ccdt_t3p_3(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
        call ccdt_t3p_4(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
        call ccdt_t3p_5(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
        call ccdt_t3p_6(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
        call ccdt_t3p_7(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
        call ccdt_t3p_diagram8_2(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
     end if
  else
     call ccdt_t3p_man01(t3out,cc_work1,cc_work2,cc_work3)
     call ccdt_t3p_man02(t3out,cc_work1,cc_work2,cc_work3)
     call ccdt_t3p_man03(t3out,cc_work1,cc_work2,cc_work3)
     call ccdt_t3p_man04(t3out,cc_work1,cc_work2,cc_work3)
     call ccdt_t3p_man05(t3out,cc_work1,cc_work2,cc_work3)
     call ccdt_t3p_man06(t3out,cc_work1,cc_work2,cc_work3)
     call ccdt_t3p_man07(t3out,cc_work1,cc_work2,cc_work3)
  end if

end subroutine ccdt_t3p_main
!**********************************************************
