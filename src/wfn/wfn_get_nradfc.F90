!######################################################################
subroutine wfn_get_nradfc(orb)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc, xrad, mapf, mapb, ndvr, ecs_flag, irad_ecs
  use mod_ormas, only : nfcore, thradfc
  use mod_const, only : zero
  use mod_sph, only : lmax1
  use mod_bas, only : mval

  implicit none
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:*)

  real(c_double) :: vtail
  real(c_double), parameter :: tiny = 1.0D-18
! real(c_double), parameter :: thrtail = 1.0D-20
! real(c_double), parameter :: thrtail = 1.0D-15
! real(c_double), parameter :: thrtail = 1.0D-12
! real(c_double), parameter :: thrtail = 1.0D-10

  integer(c_long) :: iorb, irad, l, nradmax, lradmax, nradtmp, max_irad

! Sato_ECS
  if (ecs_flag == 0) then
     max_irad = nrad - 1
  else
     max_irad = irad_ecs - 1
  end if
  nradfc = max_irad
! Sato_ECS

  if (nfcore == 0 .or. thradfc < tiny) then
     continue
  else
     do iorb = 1, nfcore
        nradmax = 1
        lradmax = abs(mval(iorb))
        do l = abs(mval(iorb)), lmax1
           nradtmp = max_irad
           do irad = 1, max_irad
              vtail = sqrt(conjg(orb(irad, l, iorb)) * orb(irad, l, iorb))
              !debug
              write(6, "(3i5, f20.10, 2e20.10)") iorb, l, irad, xrad(irad), vtail, thradfc
              !debug
              if (vtail < thradfc) then
                 nradtmp = irad
                 exit
              end if
           end do
           if (nradtmp > nradmax) then
              nradmax = nradtmp
              lradmax = l
           end if
        end do
        nradfc = min(max_irad, nradmax)
     end do
     nradfc = min(max_irad, mapf(mapb(nradfc)) + ndvr)
  end if
  write(6, "('wfn_get_nradfc: radfc = ', i5, F20.10, E20.10)") nradfc, xrad(nradfc), thradfc
  
end subroutine wfn_get_nradfc
!######################################################################
