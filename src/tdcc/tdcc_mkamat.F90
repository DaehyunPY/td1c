!######################################################################
subroutine tdcc_mkamat(cic, amat)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_ormas, only : nact
  use mod_cc, only : cc_rank,optcc,norb1,t2inp,g2inp,t3inp,g3inp

  implicit none
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(out) :: amat(1:nact**4)
  complex(kind(0d0)),allocatable :: amat_spin(:,:)

  if (cc_rank < 3 .or. .not.optcc) stop 'tdcc_mkamat supports OCCDT only.'

  call tdcc_gettcc2(cic, t2inp)
  call tdcc_getgcc2(cic, g2inp)
  call tdcc_gettcc3(cic, t3inp)
  call tdcc_getgcc3(cic, g3inp)

  allocate(amat_spin(nact**4,1:2))
  amat_spin = 0d0
  call ccdt_amat_main(amat_spin)

  ! RHF only
  amat(:) = amat_spin(:,1) + amat_spin(:,2)
  deallocate(amat_spin)

end subroutine tdcc_mkamat
!######################################################################
