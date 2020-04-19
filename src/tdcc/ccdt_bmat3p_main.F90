!**********************************************************
subroutine ccdt_bmat3p_main(den)

  use mod_ormas,only : nact
  use mod_cc,only : cc_code
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: den(1:nact,1:nact,1:*)
!debug
!write(6,"('ccdt_bmat3p_main: skip b2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     call ccdt_bmat3p_1(1,1,den)
     call ccdt_bmat3p_2(1,1,den)
  else
     call ccdt_bmat3p_man01(den,cc_work1,cc_work2)
     call ccdt_bmat3p_man02(den,cc_work1,cc_work2)
!     call ccdt_bmat3p_1(1,1,den)
!     call ccdt_bmat3p_2(1,1,den)
  end if

end subroutine ccdt_bmat3p_main
!**********************************************************
