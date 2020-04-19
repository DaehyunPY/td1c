!######################################################################
subroutine bas_pp_init()
  use, intrinsic :: iso_c_binding
  use mod_control, only : name, PSP, psp_type
  use mod_rad, only : nrad, xrad, wrad
  use mod_bas, only : znuc, pp_vloc, pp_vlocHF, pp_pproj, pp_fproj, pp_gproj, pp_maxl, pp_nump, pp_irmax
  use mod_ppatom
  use mod_oncvpsp
  use mod_tmpsp
       
  implicit none
  integer(c_int) :: irad, l, i, j
  integer(c_int), parameter :: ic = 1
  character(len=256) :: ppfile

  if (.not. PSP) return

  allocate(pp_vloc(1:(nrad-1)))
  allocate(pp_vlocHF(1:(nrad-1), 0:3))
  allocate(pp_pproj(1:(nrad-1), 1:3, 0:3))
  allocate(pp_fproj(1:(nrad-1), 1:3, 0:3))
  allocate(pp_gproj(1:(nrad-1), 1:3, 0:3))
!nyi  allocate(pp_fprojx(1:(nrad-1), 1:3, 0:3))
!nyi  allocate(pp_gprojx(1:(nrad-1), 1:3, 0:3))

  pp_vloc = 0d0
  pp_vlocHF = 0d0
  pp_pproj = 0d0
  pp_fproj = 0d0
  pp_gproj = 0d0

  ppfile = trim(trim(name)//".pp")

  if (psp_type == 1) then
     ! C. Hartwigsen, S. Goedecker, and J. Hutter, Phys. Rev. B, 58, 3641 (1998).
     ! ``Relativistic separable dual-space Gaussian pseudopotentials from H to Rn''
     call bas_pp_gen_hgh
!nyi  else if (psp_type == 2) then
!nyi     ! J. R. Trail and R. J. Needs, J. Chem. Phys, 122, 014112 (2004).
!nyi     ! ``Norm-conserving Hartree-Fock pseudopotentials and their asymptotic behavior''
!nyi     call bas_pp_gen_trailneeds
!nyi  else if (psp_type == 3) then
!nyi     ! J. R. Trail and R. J. Needs, J. Chem. Phys, 122, 174109 (2005).
!nyi     ! ``Smooth relativistic Hartree-Fock pseudopotentials for H to Ba and Lu to Hg''
!nyi     call bas_pp_gen_trailneeds
!nyi  else if (psp_type == 4) then
!nyi     ! Read and interpolate tabulated semilocal pseudopotentials
!nyi     call ppatom_init(ic,ppfile)
!nyi     call bas_pp_gen_trailneeds_num
!nyi  else if (psp_type == 5) then
!nyi     ! Read and interpolate tabulated semilocal pseudopotentials and pseudowavefunctions, and
!nyi     ! construct nonlocal separable pseudopotential of the Kleinman-Bylander form. See
!nyi     ! L. Kleinman and D. M. Bylander, Phys. Rev. Lett, 48, 1425 (1982). 
!nyi     ! ``Efficacious form for model pseudopotentials''
!nyi     call ppatom_init(ic,ppfile)
!nyi     call bas_pp_gen_trailneeds_num
  else if (psp_type == 6) then
     ! D. R. Hamann, Phys. Rev. B, 88, 085117 (2013).
     ! ``Optimized norm-conserving Vanderbilt pseudopotentials''
     call oncvpsp_init(ppfile)
     call bas_pp_gen_oncvpsp
  else if (psp_type == 7) then
     ! LDA Troullier-Martins pseudopotentials
     call tmpsp_init(ppfile)
     call bas_pp_gen_tmpsp
  else
     stop 'bas_pp_init: bad psp_type.'
  end if

end subroutine bas_pp_init
!######################################################################
