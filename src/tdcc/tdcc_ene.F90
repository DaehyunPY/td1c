!######################################################################
real(kind(0d0)) function tdcc_enec(int1e, int1x)

  use mod_ormas, only : nact
  use mod_cc, only : norb1

  implicit none
  complex(kind(0d0)), intent(in) :: int1e(1:nact, 1:nact)
  complex(kind(0d0)), intent(in) :: int1x(1:nact, 1:nact)

  integer :: iact
  complex(kind(0d0)) :: ene

  ene = 0d0
  do iact = 1, norb1
     ene = ene + int1e(iact, iact) + int1x(iact, iact)
  end do

  tdcc_enec = dble(ene)
  return

end function tdcc_enec
!######################################################################
complex(kind(0d0)) function tdcc_ene2()
!
! fock,int2x,t1inp,t2inp should be loaded before calling me
!
  use mod_cc, only : cc_rank,optcc,bcc,optbcc

  implicit none
  complex(kind(0d0)) :: ene2
!  complex(kind(0d0)), external :: ccs_ene_main
  complex(kind(0d0)), external :: ccd_ene_main
  complex(kind(0d0)), external :: ccsd_ene_main

  if (optcc) then
     ene2 = ccd_ene_main()
  else if (bcc) then
     ene2 = ccd_ene_main()
  else if (optbcc) then
     ene2 = ccd_ene_main()
  else
     ene2 = ccsd_ene_main()
  end if

  tdcc_ene2 = ene2

end function tdcc_ene2
!######################################################################
