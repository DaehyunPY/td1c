!################################################################################
subroutine tdcc_mkbmat(int1e,int2e,den1,den2,cic,dcic,bmat)

  use mod_const, only : czero
  use mod_ormas, only : nfun,nact,iprint
  use mod_cc, only : cc_rank,optcc,bcc,optbcc

  implicit none
  complex(kind(0d0)), intent(in) :: int1e(1:nact,1:nact)
  complex(kind(0d0)), intent(in) :: int2e(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)), intent(in) :: den1(1:nact,1:nact)
  complex(kind(0d0)), intent(in) :: den2(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(in) :: dcic(1:*)
  complex(kind(0d0)), intent(out) :: bmat(1:nact,1:nact)
  !--------------------------------------------------------------------
  complex(kind(0d0)), allocatable :: bmat1(:,:)
  complex(kind(0d0)), allocatable :: bmat2(:,:)
  complex(kind(0d0)), allocatable :: bmat3(:,:)

  allocate(bmat1(1:nact,1:nact))
  allocate(bmat2(1:nact,1:nact))
  allocate(bmat3(1:nact,1:nact))

  bmat(1:nact,1:nact) = czero
  bmat1(1:nact,1:nact) = czero
  bmat2(1:nact,1:nact) = czero
  bmat3(1:nact,1:nact) = czero

  if (cc_rank==2 .and. optcc) then
     call tdcc_mkbmat_1(int1e,int2e,den1,den2,bmat1)
  else if (cc_rank==2 .and. optbcc) then
     call tdcc_mkbmat_1(int1e,int2e,den1,den2,bmat1)
     call tdcc_mkbmat_2(cic,dcic,bmat2)
  else if (cc_rank==2) then
     stop 'tdcc_mkbmat: CCSD nyi.'
  else if (cc_rank==3 .and. optcc) then
     call tdcc_mkbmat_1(int1e,int2e,den1,den2,bmat1)
     call tdcc_mkbmat_2(cic,dcic,bmat2)
     call tdcc_mkbmat_3(cic,dcic,bmat3)
  else if (cc_rank==3 .and. optbcc) then
     call tdcc_mkbmat_1(int1e,int2e,den1,den2,bmat1)
     call tdcc_mkbmat_2(cic,dcic,bmat2)
     call tdcc_mkbmat_3(cic,dcic,bmat3)
  else if (cc_rank==3) then
     stop 'tdcc_mkbmat: CCSDT nyi.'
  else
     stop 'tdcc_mkbmat supports CCD/CCDT/OBCCD/OBCCDT only'
  end if
  bmat = bmat1 - bmat2 + bmat3

  if (iprint > 4) then
     write(6, "('# tdcc_mkbmat: bmat-1')")
     call util_matoutc(6,nact,bmat1)
     write(6, "('# tdcc_mkbmat: bmat-2')")
     call util_matoutc(6,nact,bmat2)
     write(6, "('# tdcc_mkbmat: bmat-3')")
     call util_matoutc(6,nact,bmat3)
     write(6, "('# tdcc_mkbmat: bmat-total')")
     call util_matoutc(6,nact,bmat)
  end if
  !stop 'for debug @ tdcc_mkbmat.'

  deallocate(bmat1)
  deallocate(bmat2)
  deallocate(bmat3)

end subroutine tdcc_mkbmat
!######################################################################
subroutine tdcc_mkbmat_1(int1e, int2e, den1, den2, bmat)
!
  use mod_bas, only : mval
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact
  use mod_cc, only : fock, int2x, norb1

  implicit none
  complex(kind(0d0)), intent(in) :: int1e(1:nact, 1:nact)
  complex(kind(0d0)), intent(in) :: int2e(1:nact, 1:nact, 1:nact, 1:nact)
  complex(kind(0d0)), intent(in) :: den1(1:nact, 1:nact)
  complex(kind(0d0)), intent(in) :: den2(1:nact, 1:nact, 1:nact, 1:nact)
  complex(kind(0d0)), intent(out) :: bmat(1:nact, 1:nact)
  integer :: iact, aact, pact, qact, ract
  complex(kind(0d0)) :: tmp1, tmp2

  do iact = 1, norb1
     do aact = norb1 + 1, nact
        if (mval(ncore+iact).ne.mval(ncore+aact)) cycle
        tmp1 = czero
        tmp2 = czero
        do pact = 1, nact
           if (mval(ncore+iact).ne.mval(ncore+pact)) cycle
           tmp1 = tmp1 + int1e(aact, pact) * den1(pact, iact)
           tmp2 = tmp2 + int1e(pact, iact) * den1(aact, pact)
        end do
        do pact = 1, nact
           do qact = 1, nact
              do ract = 1, nact
                 if (mval(ncore+iact)+mval(ncore+qact).ne.&
                     mval(ncore+pact)+mval(ncore+ract)) cycle
                 tmp1 = tmp1 + int2e(aact, pact, qact, ract) * den2(pact, iact, ract, qact)
                 tmp2 = tmp2 + int2e(pact, iact, ract, qact) * den2(aact, pact, qact, ract)
              end do
           end do
        end do
        bmat(iact, aact) = + tmp1 - tmp2
     end do
  end do

end subroutine tdcc_mkbmat_1
!######################################################################
subroutine tdcc_mkbmat_2(cic, hcic, bmat)

  use mod_const, only : czero, runit
  use mod_ormas, only : nact, iprint
  use mod_cc, only : cc_rank, norb1, optcc, bcc, optbcc
  use mod_cc, only : g1inp,t2inp,g2inp,dt2inp,t3inp,g3inp,dt3inp

  implicit none
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(in) :: hcic(1:*)
  complex(kind(0d0)), intent(out) :: bmat(1:nact,1:nact)

  bmat = czero
  if (cc_rank<=2 .and. optcc) return

  if (optbcc) then
     call tdcc_getgcc1(cic, g1inp)
  end if
  if (cc_rank >= 2) then
     call tdcc_gettcc2( cic, t2inp)
     call tdcc_getgcc2( cic, g2inp)
     call tdcc_gettcc2(hcic,dt2inp)
  end if
  if (cc_rank >= 3) then
     call tdcc_gettcc3( cic, t3inp)
     call tdcc_getgcc3( cic, g3inp)
     call tdcc_gettcc3(hcic,dt3inp)
  end if

  if (cc_rank == 2 .and. optbcc) then
     call bccd_bmat2p_main(bmat)
  else if (cc_rank == 3 .and. optcc) then
     call ccdt_bmat2p_main(bmat)
  else if (cc_rank == 3 .and. optbcc) then
     call bccdt_bmat2p_main(bmat)
  else
     stop 'tdcc_mkbmat_2 supports OBCCD,OCCDT,OBCCDT only.'
  end if

end subroutine tdcc_mkbmat_2
!######################################################################
subroutine tdcc_mkbmat_3(cic, hcic, bmat)

  use mod_const, only : czero, runit
  use mod_ormas, only : iprint, nact
  use mod_cc, only : cc_rank, norb1, optcc, bcc, optbcc
  use mod_cc, only : t1inp,t2inp,dg2inp,t3inp,dg3inp

  implicit none
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(in) :: hcic(1:*)
  complex(kind(0d0)), intent(out) :: bmat(1:nact,1:nact)

  bmat = czero
  if (cc_rank<=2) return

  if (optbcc) then
     t1inp = 0d0
  end if
  if (cc_rank >= 2) then
     call tdcc_gettcc2(cic, t2inp)
     call tdcc_getgcc2(hcic, dg2inp)
  end if
  if (cc_rank >= 3) then
     call tdcc_gettcc3(cic, t3inp)
     call tdcc_getgcc3(hcic, dg3inp)
  end if

  if (cc_rank == 3 .and. optcc) then
     call ccdt_bmat3p_main(bmat)
  else if (cc_rank == 3 .and. optbcc) then
     call ccdt_bmat3p_main(bmat)
  else
     stop 'tdcc_mkbmat_3 supports only OCCDT,OBCCDT.'     
  end if

end subroutine tdcc_mkbmat_3
!######################################################################
