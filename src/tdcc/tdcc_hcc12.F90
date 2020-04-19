!######################################################################
subroutine tdcc_hcc12(int1e, int2e, cic, hcic)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only :nact,lcic
  use mod_cc, only : norb1,cc_rank,tonly,dot1,optcc,bcc,optbcc,fock,int2x
  use mod_cc, only : t0inp,g0inp,t1inp,g1inp,t2inp,g2inp,t3inp,g3inp
  use mod_cc, only : t0out,g0out,t1out,g1out,t2out,g2out,t3out,g3out

  implicit none
  complex(kind(0d0)), intent(in) :: int1e(1:*)
  complex(kind(0d0)), intent(in) :: int2e(1:*)
  complex(kind(0d0)), intent(in) :: cic(1:lcic)
  complex(kind(0d0)), intent(out) :: hcic(1:lcic)

  call tdcc_mkint1x(int1e, int2e, fock)
  call tdcc_mkint2x(int2e, int2x)
  call tdcc_getcc(cic)
  call tdcc_zeroout

  if (cc_rank == 2 .and. optcc) then
     call ccd_t2p_main
     if (.not. tonly) then
        call ccd_l2p_main
     end if
  else if (cc_rank == 2 .and. bcc) then
     call ccd_t2p_main
     if (.not. tonly) then
        call bccd_l1p_main ! L1: Not necessary for OBCC
        call bccd_l2p_main
     end if
  else if (cc_rank == 2 .and. optbcc) then
     call ccd_t2p_main
     if (.not. tonly) then
        call bccd_l2p_main
     end if
  else if (cc_rank == 2) then
     call ccsd_t1_main
     call ccsd_t2_main
     if (.not. tonly) then
        call ccsd_l1_main
        call ccsd_l2_main
     end if
  else if (cc_rank == 3 .and. optcc) then
     call ccdt_t2p_main
     call ccdt_t3p_main
!     call tdcc_tcc2_clean(hcca(ind_tcc2))
!     call tdcc_tcc3_clean(hcca(ind_tcc3))
     if (.not. tonly) then
        call ccdt_l2p_main
        call ccdt_l3p_main
!        call tdcc_gcc2_clean(hcca(ind_gcc2))
!        call tdcc_gcc3_clean(hcca(ind_gcc3))
     end if
  else if (cc_rank == 3 .and. bcc) then
     call ccdt_t2p_main
     call ccdt_t3p_main
     if (.not. tonly) then
        call bccdt_l1p_main ! L1: Not necessary for OBCC
        call bccdt_l2p_main
        call bccdt_l3p_main
     end if
  else if (cc_rank == 3 .and. optbcc) then
     call ccdt_t2p_main
     call ccdt_t3p_main
!     call tdcc_tcc2_clean(hcca(ind_tcc2))
!     call tdcc_tcc3_clean(hcca(ind_tcc3))
     if (.not. tonly) then
        call bccdt_l2p_main
        call bccdt_l3p_main
!        call tdcc_gcc2_clean(hcca(ind_gcc2))
!        call tdcc_gcc3_clean(hcca(ind_gcc3))
     end if
  else if (cc_rank == 3) then
     call ccsdt_t1_main
     call ccsdt_t2_main
     call ccsdt_t3_main
     if (.not. tonly) then
        call ccsdt_l1_main
        call ccsdt_l2_main
        call ccsdt_l3_main
     end if
  end if

  hcic(1:lcic) = 0d0
  call tdcc_putcc(hcic)

  !debug
  !write(6,"('tdcc_hcc12: output-1')")
  !call ormas_cic_printx(hcic, 6)
  !write(6,"('tdcc_hcc12: output-2')")
  !call tdcc_print(hcic)
  !stop
  !debug

end subroutine tdcc_hcc12
!######################################################################
