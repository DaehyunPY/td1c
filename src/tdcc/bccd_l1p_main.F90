!**********************************************************
subroutine bccd_l1p_main()

  use mod_cc,only : cc_code,g1out

  implicit none
!debug
!write(6,"('bccdt_l1p_main: skip l1p.')")
!return
!debug

  if (cc_code(1:4) == 'auto') then
     call bccd_l1p_1(1,1,g1out)
     call bccd_l1p_2(1,1,g1out)
     call bccd_l1p_3(1,1,g1out)
     call bccd_l1p_4(1,1,g1out)
     call bccd_l1p_5(1,1,g1out)
     call bccd_l1p_6(1,1,g1out)
     call bccd_l1p_7(1,1,g1out)
     call bccd_l1p_8(1,1,g1out)
     call bccd_l1p_9(1,1,g1out)
     call bccd_l1p_10(1,1,g1out)
     call bccd_l1p_11(1,1,g1out)
     call bccd_l1p_12(1,1,g1out)
     call bccd_l1p_13(1,1,g1out)
  else
     write(6,"('bccdt_l1p_main: manual code nyi.')")
     call bccd_l1p_1(1,1,g1out)
     call bccd_l1p_2(1,1,g1out)
     call bccd_l1p_3(1,1,g1out)
     call bccd_l1p_4(1,1,g1out)
     call bccd_l1p_5(1,1,g1out)
     call bccd_l1p_6(1,1,g1out)
     call bccd_l1p_7(1,1,g1out)
     call bccd_l1p_8(1,1,g1out)
     call bccd_l1p_9(1,1,g1out)
     call bccd_l1p_10(1,1,g1out)
     call bccd_l1p_11(1,1,g1out)
     call bccd_l1p_12(1,1,g1out)
     call bccd_l1p_13(1,1,g1out)
  end if

end subroutine bccd_l1p_main
!**********************************************************
