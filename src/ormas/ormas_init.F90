!################################################################################
subroutine ormas_init(idebug, mtot)

  use, intrinsic :: iso_c_binding
  use mod_control, only : name
  use mod_ormas, only : nact, nstr_alph, nstr_beta
  use mod_ormas, only : fab_den2, iprint, cic_old, lcic, ndetx, tdcc

  implicit none
  integer(c_long), intent(in) :: idebug, mtot
!debug  namelist /debug/ cic_old

  iprint = idebug
  cic_old = .false.
!debug  open(unit=99,file=trim(trim(name)//".inp"),status='old',form='formatted')
!debug  read(99,nml=debug)
!debug  write(6,nml=debug)
!debug  close(99)

  ! get info
  call ormas_info()

  if (nact == 0) then
     lcic = 1
     ndetx = 1
     nstr_alph = 1
     nstr_beta = 1
     return
  else
     ! spin-occupation boundaries
     call ormas_occbc()

     ! spin-distributions
     call ormas_dist()
     call ormas_nstr()

     ! spin-strings
     call ormas_arcwgt()
     call ormas_str()

     ! single excitation map
     call ormas_int1x()
!OLD     call ormas_int1xr()

     ! non-redundant orbital rotations
     call ormas_rotoo()
     !call ormas_rotov()

     ! set up arrays for mkden2 version 4 by fabian
     if (fab_den2) call ormas_fab()

     ! mval-adaptation
     call ormas_madapt(mtot);

     ! tdcc setting
     if (.not. tdcc) then
        lcic = ndetx
     else
        lcic = 2*ndetx
        call tdcc_init()
     end if
  end if

!  stop 'for debug @ ormas_init'
end subroutine ormas_init
!################################################################################
