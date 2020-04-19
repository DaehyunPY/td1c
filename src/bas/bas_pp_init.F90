!######################################################################
subroutine bas_pp_init()
  use, intrinsic :: iso_c_binding
  use mod_control, only : PSP, psp_type
  use mod_rad, only : nrad, xrad, wrad
  use mod_bas, only : znuc, pp_vloc, pp_vlocHF, pp_pproj, pp_fproj, pp_gproj, pp_maxl, pp_nump, pp_irmax
       
  implicit none
  integer(c_long) :: irad, l, i, j

   if (.not. PSP) return

  allocate(pp_vloc(1:(nrad-1)))
  allocate(pp_vlocHF(1:(nrad-1), 0:2))
  allocate(pp_pproj(1:(nrad-1), 1:3, 0:2))
  allocate(pp_fproj(1:(nrad-1), 1:3, 0:2))
  allocate(pp_gproj(1:(nrad-1), 1:3, 0:2))

  pp_vloc = 0d0
  pp_vlocHF = 0d0
  pp_pproj = 0d0
  pp_fproj = 0d0
  pp_gproj = 0d0

  if (psp_type == 1) then
     ! C. Hartwigsen, S. Goedecker, and J. Hutter, Phys. Rev. B, 58, 3641 (1998).
     ! ``Relativistic separable dual-space Gaussian pseudopotentials from H to Rn''
     call bas_pp_gen_hgh
  else if (psp_type == 2) then
     ! J. R. Trail and R. J. Needs, J. Chem. Phys, 122, 014112 (2004).
     ! ``Norm-conserving Hartree-Fock pseudopotentials and their asymptotic behavior''
     call bas_pp_gen_trailneeds
  else if (psp_type == 3) then
     ! J. R. Trail and R. J. Needs, J. Chem. Phys, 122, 174109 (2005).
     ! ``Smooth relativistic Hartree-Fock pseudopotentials for H to Ba and Lu to Hg''
     call bas_pp_gen_trailneeds
  else if (psp_type == 4) then
     ! Read and interpolate tabulated semilocal pseudopotentials
     call bas_pp_gen_trailneeds_num
  else if (psp_type == 5) then
     ! Read and interpolate tabulated semilocal pseudopotentials and pseudowavefunctions, and
     ! construct nonlocal separable pseudopotential of the Kleinman-Bylander form. See
     ! L. Kleinman and D. M. Bylander, Phys. Rev. Lett, 48, 1425 (1982). 
     ! ``Efficacious form for model pseudopotentials''
     call bas_pp_gen_trailneeds_num
  else
     stop 'bas_pp_init: bad psp_type.'
  end if

end subroutine bas_pp_init
!######################################################################
