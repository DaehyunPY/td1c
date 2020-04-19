!**********************************************************
subroutine ccdt_l3p_1e_main()

  use, intrinsic :: iso_c_binding
  use mod_ormas,only : nact
  use mod_cc,only : norb1,cc_code,g3out

  implicit none
!debug
!write(6,"('ccdt_l3p_1e_main: skip t2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (norb1 > 2) then
        call ccdt_l3_1e_1(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3_1e_2(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
        call ccdt_l3_1e_3(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
     end if
     if (norb1 >= 2) then
        call ccdt_l3_1e_1(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3_1e_2(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
        call ccdt_l3_1e_3(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
     end if
  else
     call ccdt_l3p_1e_man01(g3out)
     call ccdt_l3p_1e_man02(g3out)
     call ccdt_l3p_1e_man03(g3out)
!     call ccdt_l3_1e_1(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1)) 
!     call ccdt_l3_1e_1(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
!     call ccdt_l3_1e_2(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1))
!     call ccdt_l3_1e_2(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
!     call ccdt_l3_1e_3(1,1,1,1,1,1,g3out(1,1,1,norb1+1,norb1+1,norb1+1,1)) 
!     call ccdt_l3_1e_3(1,1,2,1,1,2,g3out(1,1,1,norb1+1,norb1+1,norb1+1,2))
  end if

end subroutine ccdt_l3p_1e_main
!**********************************************************
