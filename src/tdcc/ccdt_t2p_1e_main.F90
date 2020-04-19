!**********************************************************
subroutine ccdt_t2p_1e_main()

  use mod_ormas, only : nact
  use mod_cc, only : norb1,cc_code,t2out

  implicit none
!debug
!write(6,"('ccdt_t2p_1e_main: skip t2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     if (norb1 >= 2) then
        call ccdt_t2_1e_1(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2_1e_2(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
        call ccdt_t2_1e_3(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
     end if
     call ccdt_t2_1e_1(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
     call ccdt_t2_1e_2(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
     call ccdt_t2_1e_3(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
  else
     call ccdt_t2p_1e_man01(t2out)
     call ccdt_t2p_1e_man02(t2out)
     call ccdt_t2p_1e_man03(t2out)
!     call ccdt_t2_1e_1(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
!     call ccdt_t2_1e_2(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
!     call ccdt_t2_1e_3(1, 1, 1, 1, t2out(norb1+1,norb1+1,1,1,1))
!     call ccdt_t2_1e_1(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
!     call ccdt_t2_1e_2(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
!     call ccdt_t2_1e_3(1, 2, 1, 2, t2out(norb1+1,norb1+1,1,1,2))
  end if

end subroutine ccdt_t2p_1e_main
!**********************************************************
