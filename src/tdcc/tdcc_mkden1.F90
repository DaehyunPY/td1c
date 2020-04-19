!######################################################################
subroutine tdcc_mkden1(cic, den1)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nact
  use mod_cc, only : norb1,cc_rank,optcc,bcc,optbcc,den1s

  implicit none
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(out) :: den1(1:*)

  call tdcc_getcc(cic)
  den1s(1:nact, 1:nact, 1:2) = 0d0

  if (cc_rank == 2 .and. optcc) then
     call ccd_den1p_main(den1s)
  else if (cc_rank == 2 .and. bcc) then
     call bccd_den1p_main(den1s)
  else if (cc_rank == 2 .and. optbcc) then
     call bccd_den1p_main(den1s)
  else if (cc_rank == 2) then
     call ccsd_den1_main(den1s)
  else if (cc_rank == 3 .and. optcc) then
     call ccdt_den1p_main(den1s)
  else if (cc_rank == 3 .and. bcc) then
     call bccdt_den1p_main(den1s)
  else if (cc_rank == 3 .and. optbcc) then
     call bccdt_den1p_main(den1s)
  else
     stop 'tdcc_den1 nyi.'
  end if

  call tdcc_mkden1_ref1(den1s)
  call tdcc_mkden1_spac1(den1s, den1)

end subroutine tdcc_mkden1
!######################################################################
subroutine tdcc_mkden1_ref1(den1s)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_const, only : runit
  use mod_ormas, only : nact,nelact
  use mod_cc, only : norb1

  implicit none
  complex(kind(0d0)), intent(inout) :: den1s(1:nact, 1:nact, 1:2)

  integer(c_int) :: iact

  do iact = 1, norb1
     den1s(iact, iact, 1) = runit + den1s(iact, iact, 1)
  end do
  if (smul==1.and.nelact(1)==nelact(2)) then
     den1s(:,:,2) = den1s(:,:,1)
  else if (nelact(2)==0) then
     den1s(:,:,2) = 0d0
  else
     stop 'tdcc_mkden1_ref1: A=B.or.B=0 only'
  end if

end subroutine tdcc_mkden1_ref1
!######################################################################
subroutine tdcc_mkden1_spac1(den1s, den1)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, ctwo, chalf
  use mod_bas, only : mval
  use mod_ormas, only : nact, ncore
  use mod_cc, only : norb1, nbiort

  implicit none
  complex(kind(0d0)), intent(in) :: den1s(1:nact, 1:nact, 1:2)
  complex(kind(0d0)), intent(out) :: den1(1:nact, 1:nact)

  !DEBUG
  complex(kind(0d0)) :: ntot
  complex(kind(0d0)), allocatable :: dtmp(:,:)
  complex(kind(0d0)), allocatable :: utmp(:,:)
  !DEBUG
  integer(c_int) :: iact, jact, aact, bact, pact, qact

  ! den1 = 1/2 (den1s( + den1s^+)
  ! spin summation (2) and a factor (1/2) cancel

!  do pact = 1, nact
!  do qact = 1, nact
!     den1s(pact, qact) = den1s(pact, qact) * ctwo
!  end do
!  end do

  if (nbiort == 1) then
     do pact = 1, nact
     do qact = 1, nact
        if (mval(ncore+pact).ne.mval(ncore+qact)) cycle
        den1(pact,qact) = (den1s(pact,qact,1) + conjg(den1s(qact,pact,1))) * chalf &
                        + (den1s(pact,qact,2) + conjg(den1s(qact,pact,2))) * chalf
     end do
     end do
  else
     do pact = 1, nact
     do qact = 1, nact
        if (mval(ncore+pact).ne.mval(ncore+qact)) cycle
        den1(pact, qact) = den1s(pact, qact, 1) + den1s(pact, qact, 2)
     end do
     end do
  end if
!  !DEBUG
!  ntot = czero
!  write(6, "('den1cc-R:')")
!  do iact = 1, nact
!     ntot = ntot + den1(iact, iact)
!     do jact = 1, nact
!        write(6, "(f10.5)", advance = 'no') dble(den1(iact, jact))
!     end do
!     write(6, *)
!  end do
!  write(6, "('den1cc-I:')")
!  do iact = 1, nact
!     ntot = ntot + den1(iact, iact)
!     do jact = 1, nact
!        write(6, "(f10.5)", advance = 'no') aimag(den1(iact, jact))
!     end do
!     write(6, *)
!  end do
!  write(6, "('ntot = ', 2f20.10)") ntot
!  allocate(dtmp(1:nact, 1:nact))
!  allocate(utmp(1:nact, 1:nact))
!  dtmp(1:nact, 1:nact) = den1(1:nact, 1:nact)
!  utmp(1:nact, 1:nact) = czero
!  call util_diag_comp(.true., nact, dtmp, utmp)
!  do iact = 1, nact
!     write(6, "(2f20.10)") dtmp(iact, iact)
!  end do
!  deallocate(utmp)
!  deallocate(dtmp)
!  !DEBUG

end subroutine tdcc_mkden1_spac1
!######################################################################
