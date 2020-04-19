!######################################################################
subroutine bas_pp_gen_trailneeds_num()

  use, intrinsic :: iso_c_binding
  use mod_control, only : psp_type
  use mod_rad, only : nrad, xrad, wrad
  use mod_bas, only : znuc, psp_label, pp_znuc, pp_rloc, pp_cloc, pp_rproj, pp_hproj, &
       pp_vloc, pp_vlocHF, pp_pproj, pp_fproj, pp_gproj, pp_maxl, pp_nump, pp_irmax
  use mod_ppatom

  implicit none
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  integer(c_int), parameter :: ic = 1
  integer(c_int) :: irad, label
  logical :: found
  real(c_double), parameter :: thresh = 1D-20

  label = znuc

  if (psp_type == 4) then
     do irad = 1, nrad - 1
        pp_vlocHF(irad,0) = ppatom_getvs(ic,xrad(irad))
        pp_vlocHF(irad,1) = ppatom_getvp(ic,xrad(irad))
        pp_vlocHF(irad,2) = ppatom_getvloc(ic,xrad(irad))
     end do
     !debug
     do irad = 1, nrad - 1
        write(6,"('vloc: ',4f20.10)") xrad(irad),pp_vlocHF(irad,0:2)
     end do
     !debug
  else if (psp_type == 5) then
     ! local part
     do irad = 1, nrad - 1
        pp_vloc(irad) = ppatom_getvloc(ic,xrad(irad))
     end do
     !debug
     do irad = 1, nrad - 1
        write(6,"('vloc: ',2f20.10)") xrad(irad),pp_vloc(irad)
     end do
     !debug

     ! nonlocal part
     pp_maxl = 1
     pp_nump = 1

     pp_hproj = 0d0
     pp_hproj(1,1,0) = 1d0/ppatom_seig(ic)
     pp_hproj(1,1,1) = 1d0/ppatom_peig(ic)

     pp_irmax = 1
     pp_pproj = 0d0
     pp_fproj = 0d0
     pp_gproj = 0d0
     ! s-channel
     found = .false.
     do irad = 1, nrad - 1
        pp_pproj(irad,1,0) = ppatom_getsfun(ic,xrad(irad))
        pp_fproj(irad,1,0) = pp_pproj(irad,1,0) * xrad(irad) * sqrt(wrad(irad))
        pp_gproj(irad,1,0) = pp_hproj(1,1,0)*pp_fproj(irad,1,0)
        if (.not.found .and. &
             abs(pp_pproj(irad,1,0)) < thresh .and. &
             abs(pp_fproj(irad,1,0)) < thresh .and. &
             abs(pp_gproj(irad,1,0)) < thresh) then
           found = .true.
           pp_irmax(1,0) = irad
        end if
     end do

     ! p-channel
     found = .false.
     do irad = 1, nrad - 1
        pp_pproj(irad,1,1) = ppatom_getpfun(ic,xrad(irad))
        pp_fproj(irad,1,1) = pp_pproj(irad,1,1) * xrad(irad) * sqrt(wrad(irad))
        pp_gproj(irad,1,1) = pp_hproj(1,1,1)*pp_fproj(irad,1,1)
        if (.not.found .and. &
             abs(pp_pproj(irad,1,1)) < thresh .and. &
             abs(pp_fproj(irad,1,1)) < thresh .and. &
             abs(pp_gproj(irad,1,1)) < thresh) then
           found = .true.
           pp_irmax(1,1) = irad
        end if
     end do
  else
     stop 'bas_pp_gen_trailneeds_num: bad psp_type.'
  end if

!  stop 'for debug @ bas_pp_trailneeds_num'

!  contains
  !######################################################################
  !######################################################################
end subroutine bas_pp_gen_trailneeds_num
!######################################################################
