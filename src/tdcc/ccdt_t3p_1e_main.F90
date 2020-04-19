!**********************************************************
subroutine ccdt_t3p_1e_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,cc_code,t3out
  use mod_cc2

  implicit none
!debug
!write(6,"('ccdt_t3p_1e_main: skip t2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (norb1 > 2) then
        call ccdt_t3_1e_1(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
        call ccdt_t3_1e_2(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
        call ccdt_t3_1e_3(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
     end if
     if (norb1 >= 2) then
        call ccdt_t3_1e_1(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
        call ccdt_t3_1e_2(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
        call ccdt_t3_1e_3(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
     end if
  else
     call ccdt_t3p_1e_man01(t3out,cc_work1,cc_work2)
     call ccdt_t3p_1e_man02(t3out,cc_work1,cc_work2)
     call ccdt_t3p_1e_man03(t3out,cc_work1,cc_work2)
!     call ccdt_t3_1e_1(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1)); 
!     call ccdt_t3_1e_1(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
!     call ccdt_t3_1e_2(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
!     call ccdt_t3_1e_2(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))
!     call ccdt_t3_1e_3(1,1,1,1,1,1,t3out(norb1+1,norb1+1,norb1+1,1,1,1,1))
!     call ccdt_t3_1e_3(1,1,2,1,1,2,t3out(norb1+1,norb1+1,norb1+1,1,1,1,2))     
  end if

end subroutine ccdt_t3p_1e_main
!**********************************************************
