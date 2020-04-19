!######################################################################
subroutine bas_pp_gen_tmpsp()

  use, intrinsic :: iso_c_binding
  use mod_const, only : pi
  use mod_control, only : psp_type
  use mod_rad, only : nrad, xrad, wrad
  use mod_bas, only : pp_vloc, pp_fproj, pp_gproj, pp_maxl, pp_nump, pp_irmax
  use mod_tmpsp

  implicit none
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  integer(c_int) :: irad,l,iproj
  logical :: found
  real(c_double), parameter :: thresh = 1D-20

  ! local part
  do irad = 1, nrad - 1
     pp_vloc(irad) = tmpsp_getvloc(xrad(irad))
  end do

  ! nonlocal part
  pp_maxl = tmpsp_maxl()
  do l = 0, pp_maxl
     pp_nump(l) = 1
     do iproj = 1, pp_nump(l)
        pp_irmax(iproj,l) = 0
        pp_fproj(:,iproj,l) = 0d0
        pp_gproj(:,iproj,l) = 0d0
        found = .false.
        do irad = 1, nrad - 1
           !pp_fproj(irad,iproj,l) = tmpsp_getproj(l,xrad(irad)) * xrad(irad) * sqrt(wrad(irad))
           pp_fproj(irad,iproj,l) = tmpsp_getproj(l,xrad(irad)) * sqrt(wrad(irad))
           pp_gproj(irad,iproj,l) = pp_fproj(irad,iproj,l) * tmpsp_ekb(l)
           if (.not.found .and. &
                abs(pp_fproj(irad,iproj,l)) < thresh .and. &
                abs(pp_gproj(irad,iproj,l)) < thresh) then
              found = .true.
              pp_irmax(iproj,l) = irad
           end if
        end do
     end do
  end do

  !debug
  write(6,"('pp_maxl: ',i5)") pp_maxl
  do l = 0, pp_maxl
     write(6,"('pp_nump: ',2i5)") l,pp_nump(l)
  end do
  do l = 0, pp_maxl
     do iproj = 1, pp_nump(l)
        write(6,"('pp_rmax: ',3i5,f15.10)") l,iproj,pp_irmax(iproj,l),xrad(pp_irmax(iproj,l))
     end do
  end do
!  do irad = 1, nrad - 1
!     write(6,"('pp: ',f10.5)",advance='no') xrad(irad)
!     write(6,"(f10.5)",advance='no') wrad(irad)
!     do l = 0, pp_maxl
!        do iproj = 1, pp_nump(l)
!           write(6,"(2f15.10)",advance='no') pp_fproj(irad,iproj,l), pp_gproj(irad,iproj,l)
!        end do
!     end do
!     write(6,"(2f15.10)") pp_vloc(irad),tmpsp_getd1vloc(xrad(irad))
!  end do
!  stop 'for debug @ bas_pp_tmpsp'
  !debug

!  contains
  !######################################################################
  !######################################################################
end subroutine bas_pp_gen_tmpsp
!######################################################################
