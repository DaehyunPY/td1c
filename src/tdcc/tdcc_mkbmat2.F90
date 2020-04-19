!################################################################################
subroutine tdcc_mkbmat2(fac,cic,dcic,bmat)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit
  use mod_ormas, only : nact, iprint
  use mod_cc, only : cc_rank, norb1, optcc, bcc, optbcc
  use mod_cc, only : t1inp,g1inp,t2inp,g2inp,t3inp,g3inp
  use mod_cc, only : dt2inp,dg2inp,dt3inp,dg3inp

  implicit none
  complex(kind(0d0)),intent(in) :: fac
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(in) :: dcic(1:*)
  complex(kind(0d0)), intent(inout) :: bmat(1:nact,1:nact)

  if (cc_rank == 2 .and. optbcc) then
     stop 'tdcc_mkbmat2: nyi for BCCD'
     t1inp = 0d0
     call tdcc_getgcc1(cic, g1inp)
     call tdcc_gettcc2(cic, t2inp)
     call tdcc_getgcc2(cic, g2inp)
     call tdcc_gettcc2(dcic,dt2inp)
     call tdcc_getgcc2(dcic, dg2inp)
!nyi     call bccd_bmat2p_main(bmat)
  else if (cc_rank == 3 .and. optcc) then
     call tdcc_gettcc2(cic,t2inp)
     call tdcc_getgcc2(cic,g2inp)
     call tdcc_gettcc2(dcic,dt2inp)
     call tdcc_getgcc2(dcic,dg2inp)
     call tdcc_gettcc3(cic, t3inp)
     call tdcc_getgcc3(cic, g3inp)
     call tdcc_gettcc3(dcic,dt3inp)
     call tdcc_getgcc3(dcic, dg3inp)
     call ccdt_bmat2_main(fac, bmat)
  else if (cc_rank == 3 .and. optbcc) then
     stop 'tdcc_mkbmat2: nyi for BCCDT'
     t1inp = 0d0
     call tdcc_getgcc1(cic, g1inp)
     call tdcc_gettcc2( cic, t2inp)
     call tdcc_getgcc2( cic, g2inp)
     call tdcc_gettcc2(dcic,dt2inp)
     call tdcc_getgcc2(dcic,dg2inp)
     call tdcc_gettcc3( cic, t3inp)
     call tdcc_getgcc3( cic, g3inp)
     call tdcc_gettcc3(dcic,dt3inp)
     call tdcc_getgcc3(dcic,dg3inp)
!nyi     call bccdt_bmat2p_main(bmat)
!nyi     call ccdt_bmat3p_main(bmat)
  end if

end subroutine tdcc_mkbmat2
!######################################################################
