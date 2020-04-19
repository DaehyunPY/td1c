!######################################################################
subroutine tdcc_mkden2x(cic, den2x)

  use mod_const, only : czero
  use mod_ormas, only : nact
  use mod_cc, only : cc_rank,optcc,norb1,t2inp,g2inp,t3inp,g3inp

  implicit none
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(out) :: den2x(1:*)
  complex(kind(0d0)), allocatable :: d2x(:,:,:,:,:,:)

  if (cc_rank < 3 .or. .not.optcc) stop 'tdcc_mkden2x supports OCCDT only.'
  allocate(d2x(1:nact,1:nact,1:nact,1:nact,1:2,1:2))
  d2x = czero

  call tdcc_gettcc2(cic, t2inp)
  call tdcc_getgcc2(cic, g2inp)
  call tdcc_gettcc3(cic, t3inp)
  call tdcc_getgcc3(cic, g3inp)

  call ccdt_den2xp_main(d2x)
!  call ccdt_den2x_main(d2x)
!  call ccdt_den2x_pathak_main(d2x)
!  call ccdt_den2x_old_main(d2x)

  call tdcc_mkden2x_spac2(d2x, den2x)
  deallocate(d2x)

end subroutine tdcc_mkden2x
!######################################################################
subroutine tdcc_mkden2x_spac2(d2x, den2x)

  use mod_const, only : ctwo, chalf
  use mod_ormas, only : nact
  use mod_cc, only : norb1

  implicit none
  complex(kind(0d0)), intent(inout) :: d2x(1:nact, 1:nact, 1:nact, 1:nact, 1:2, 1:2)
  complex(kind(0d0)), intent(out) :: den2x(1:nact, 1:nact, 1:nact, 1:nact, 1:2)

  integer :: pact, qact, ract, sact
  integer :: iact, jact, kact, lact
  integer :: aact, bact, cact, dact
  complex(kind(0d0)) :: tmp

  do iact = 1, norb1
  do jact = 1, norb1
  do aact = norb1 + 1, nact
  do bact = norb1 + 1, nact
!    den2x(iact, aact, jact, bact) = (d2x(aact, bact, iact, jact, 1) &
!                                  +  d2x(aact, bact, iact, jact, 2)) * ctwo
    den2x(iact, aact, jact, bact, 1) =  d2x(aact, bact, iact, jact, 1, 1) &
                                     +  d2x(aact, bact, iact, jact, 2, 1)
!    den2x(iact, aact, jact, bact) = (d2x(aact, bact, iact, jact, 1) &
!                                  +  d2x(aact, bact, iact, jact, 2)) * 100.d+0
  end do
  end do
  end do
  end do

end subroutine tdcc_mkden2x_spac2
!######################################################################
