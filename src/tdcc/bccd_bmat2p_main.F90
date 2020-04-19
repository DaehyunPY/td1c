!**********************************************************
subroutine bccd_bmat2p_main(den)

  use mod_ormas,only : nact
  use mod_cc,only : cc_code
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: den(1:nact,1:nact,1:*)
!debug
!write(6,"('bccdt_bmat2p_main: skip b2.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     call bccd_bmat2p_1(1,1,den)
  else
     ! L1 contribution
     call bccdt_bmat2p_man01(den,cc_work1,cc_work2)
  end if

end subroutine bccd_bmat2p_main
!**********************************************************
