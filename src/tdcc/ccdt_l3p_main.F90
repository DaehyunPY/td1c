!**********************************************************
subroutine ccdt_l3p_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact,nelact
  use mod_cc,only : norb1,cc_code,g3out
  use mod_cc2

  implicit none
!debug
!write(6,"('ccdt_l3p_main: skip l3.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (nelact(1)>=3.and.norb1>=3) then
        call ccdt_l3p_1(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_2(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_3(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_4(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_5(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_6(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_7(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_8(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_9(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_10(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3p_diagram11_2(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
     end if
     if (nelact(1)>=2.and.nelact(2)>=1.and.norb1>=2) then
        call ccdt_l3p_1(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_2(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_3(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_4(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_5(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_6(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_7(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_8(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_9(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_10(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3p_diagram11_2(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
     end if
  else
    !call ccdt_l3p_man01(g3out,cc_work1,cc_work2,cc_work3) ! merged with 8
    !call ccdt_l3p_man02(g3out,cc_work1,cc_work2,cc_work3) ! merged with 8
    !call ccdt_l3p_man03(g3out,cc_work1,cc_work2,cc_work3) ! merged with 8
     call ccdt_l3p_man04(g3out,cc_work1,cc_work2,cc_work3)
     call ccdt_l3p_man05(g3out,cc_work1,cc_work2,cc_work3)
     call ccdt_l3p_man06(g3out,cc_work1,cc_work2,cc_work3)
     call ccdt_l3p_man07(g3out,cc_work1,cc_work2,cc_work3)
     call ccdt_l3p_man08(g3out,cc_work1,cc_work2,cc_work3)
     call ccdt_l3p_man09(g3out,cc_work1,cc_work2,cc_work3)
     call ccdt_l3p_man10(g3out,cc_work1,cc_work2,cc_work3)
     call ccdt_l3p_man11(g3out,cc_work1,cc_work2,cc_work3)
  end if
  ! following diagrams involve T3 or L3 contributions
  !  ccdt_l3p_man04.F90
  !  ccdt_l3p_man05.F90
  !  ccdt_l3p_man06.F90
  !  ccdt_l3p_man07.F90
  !  ccdt_l3p_man08.F90
  !  ccdt_l3p_man09.F90
  !  ccdt_l3p_man10.F90
  !  ccdt_l3p_man11.F90

end subroutine ccdt_l3p_main
!**********************************************************
