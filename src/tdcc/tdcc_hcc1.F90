!######################################################################
subroutine tdcc_hcc1(int1e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only :nact,lcic
  use mod_cc, only : cc_rank,optcc,fock

  implicit none
  complex(kind(0d0)), intent(in) :: int1e(1:nact,1:nact)
  complex(kind(0d0)), intent(in) :: cic(1:lcic)
  complex(kind(0d0)), intent(out) :: hcic(1:lcic)
!debug
!hcic(1:lcic) = 0d0
!write(6,"('tdcc_hcc1: skip hcc1.')")
!return
!debug

  if (cc_rank < 3 .or. .not.optcc) stop 'tdcc_hcc1 supports only OCCDT.'

  fock(:,:,1) = int1e(:,:)
  fock(:,:,2) = int1e(:,:)
  call tdcc_getcc(cic)
  call tdcc_zeroout

  call ccdt_t2p_1e_main
  call ccdt_t3p_1e_main
  call ccdt_l2p_1e_main
  call ccdt_l3p_1e_main

!!  call ccdt_t2_1ep_main(work(ind_tcc2))
!!  call ccdt_t3_1ep_main(work(ind_tcc3))
!!  call ccdt_l2_1ep_main(work(ind_gcc2))
!!  call ccdt_l3_1ep_main(work(ind_gcc3))
!  call ccnew_tcc2_clean(work(ind_tcc2))
!  call ccnew_gcc2_clean(work(ind_gcc2))
!  call ccnew_tcc3_clean(work(ind_tcc3))
!  call ccnew_gcc3_clean(work(ind_gcc3))

  hcic(1:lcic) = 0d0
  call tdcc_putcc(hcic)

end subroutine tdcc_hcc1
!######################################################################
