!**********************************************************
subroutine bccdt_t1p_main()

  use mod_cc,only : cc_code,t1out

  implicit none
!debug
!write(6,"('bccdt_t1p_main: skip t1p.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     call bccdt_t1p_1(1,1,t1out)
     call bccdt_t1p_2(1,1,t1out)
     call bccdt_t1p_3(1,1,t1out)
     call bccdt_t1p_4(1,1,t1out)
     call bccdt_t1p_5(1,1,t1out)
  else
     call bccdt_t1p_man01(t1out)
  end if

end subroutine bccdt_t1p_main
!**********************************************************
