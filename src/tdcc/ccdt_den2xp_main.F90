!**********************************************************
subroutine ccdt_den2xp_main(den)

  use mod_ormas,only : nact
  use mod_cc,only : cc_code
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: den(1:nact,1:nact,1:nact,1:nact,1:*)
!debug
!write(6,"('ccdt_den2xp_main: skip d2x.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     call ccdt_den2xp_1(1,1,1,1,den(1,1,1,1,1))
     call ccdt_den2xp_2(1,1,1,1,den(1,1,1,1,1))
     call ccdt_den2xp_3(1,1,1,1,den(1,1,1,1,1))
     call ccdt_den2xp_4(1,1,1,1,den(1,1,1,1,1))
   
     call ccdt_den2xp_1(1,2,1,2,den(1,1,1,1,2))
     call ccdt_den2xp_2(1,2,1,2,den(1,1,1,1,2))
     call ccdt_den2xp_3(1,2,1,2,den(1,1,1,1,2))
     call ccdt_den2xp_4(1,2,1,2,den(1,1,1,1,2))
  else
     call ccdt_den2xp_man01(den,cc_work1,cc_work2)
     call ccdt_den2xp_man02(den,cc_work1,cc_work2)
     call ccdt_den2xp_man03(den,cc_work1,cc_work2)
     call ccdt_den2xp_man04(den,cc_work1,cc_work2)
!     call ccdt_den2xp_1(1,1,1,1,den(1,1,1,1,1));call ccdt_den2xp_1(1,2,1,2,den(1,1,1,1,2))
!     call ccdt_den2xp_2(1,1,1,1,den(1,1,1,1,1));call ccdt_den2xp_2(1,2,1,2,den(1,1,1,1,2))
!     call ccdt_den2xp_3(1,1,1,1,den(1,1,1,1,1));call ccdt_den2xp_3(1,2,1,2,den(1,1,1,1,2))
!     call ccdt_den2xp_4(1,1,1,1,den(1,1,1,1,1));call ccdt_den2xp_4(1,2,1,2,den(1,1,1,1,2))
  end if

end subroutine ccdt_den2xp_main
!**********************************************************
