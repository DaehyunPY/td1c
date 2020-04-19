!**********************************************************
subroutine ccdt_amat_main(amat_spin)

  use, intrinsic :: iso_c_binding
  use mod_control, only : xact2_type
  use mod_ormas,only : nact,nact1
  use mod_cc,only : cc_code,cc_rank
  use mod_cc2

  implicit none
  complex(kind(0d0)),intent(inout) :: amat_spin(1:nact**4,1:*)

  if (cc_code(1:4) == 'auto') then
     stop 'ccdt_amat_main: cc_code=auto not supported.'
  else
     call ccdt_amat_man01(amat_spin,cc_work1,cc_work2)
     call ccdt_amat_man02(amat_spin,cc_work1,cc_work2)
     call ccdt_amat_man03(amat_spin,cc_work1,cc_work2)
     call ccdt_amat_man04(amat_spin,cc_work1,cc_work2)
     !##### 2019/6/7 ##### 05 & 06 are symmetric, thus vanishing
     !call ccdt_amat_man05(amat_spin,cc_work1,cc_work2,cc_work3)
     !call ccdt_amat_man06(amat_spin,cc_work1,cc_work2,cc_work3)
     !##### 2019/6/7 #####
  end if

  if (cc_rank>=3 .and. nact.ne.nact1 .and. xact2_type==3) then
     call ccdt_amat_IIJJ(amat_spin,cc_work1,cc_work2,cc_work3)
     call ccdt_amat_AABB(amat_spin,cc_work1,cc_work2,cc_work3)
     call ccdt_amat_IAJJ(amat_spin,cc_work1,cc_work2,cc_work3)
     call ccdt_amat_IABB(amat_spin,cc_work1,cc_work2,cc_work3)
     call ccdt_amat_IIJB(amat_spin,cc_work1,cc_work2,cc_work3)
     call ccdt_amat_AAJB(amat_spin,cc_work1,cc_work2,cc_work3)
  end if

!debug
!stop 'for debug @ ccdt_amat_main.'
!debug

end subroutine ccdt_amat_main
!**********************************************************
