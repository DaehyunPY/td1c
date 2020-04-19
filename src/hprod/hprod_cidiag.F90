!#######################################################################
subroutine hprod_cidiag(wfn, cic)

  use, intrinsic :: iso_c_binding
  use mod_control, only : exact3j, PSP
  use mod_const, only : zero, czero, runit, ctrue
  use mod_ormas, only : neltot, ndcore, nact, nelact, lcic, tdcc
  use mod_hprod, only : rho2, v2sph
  use mod_hprod, only : orb, orbg, h0orb, h1orb, gorb, gorbg
  use mod_hprod, only : ene_fcore, ene_dcore, ene_core
  use mod_hprod, only : int1e, int2e

  implicit none
  complex(c_double_complex), intent(in) :: wfn(*)
  complex(c_double_complex), intent(inout) :: cic(1:lcic)

  real(c_double) :: lfield(1:9) = zero
  real(c_double), external :: hprod_ene_fcore
  real(c_double), external :: hprod_ene_dcore
  complex(c_double_complex), external :: util_zdotc

  integer(c_int) :: idav, jdav, isub, jsub
  integer(c_int), parameter :: max_dav = 30
  real(c_double), parameter :: thr_dav = 1.D-10
!debug  integer(c_int), parameter :: max_dav = 20
!debug  real(c_double), parameter :: thr_dav = 1.D-8
  real(c_double) :: norm
  complex(c_double_complex) :: ovlp
  complex(c_double_complex), allocatable :: hcic(:), xcic(:,:), ycic(:,:)
  complex(c_double_complex), allocatable :: hmat(:,:), hsub(:,:), usub(:,:)

  if (nact==0 .or. nelact(3)<=1) return
!  call hprod_semicanonical(wfn)

  call hprod_htot_init
  call hprod_orbin(lfield, wfn, orb, orbg)
  call hprod_tprod_dyn(orb, h0orb)
  if (PSP) call hprod_projpp(runit, lfield, orb, h0orb)
  if (neltot(3) >= 2) then
     if (.not.exact3j) then
        call hprod_mkrho2_dyn
        call hprod_mkv2mf_dyn
        call hprod_mfprod2_ene
        call bas_ang2sph1_dyn(gorbg, gorb);
     else
        call hprod_mkrho2_x3j_dyn
        call hprod_mkv2mf_x3j_dyn
        call hprod_mfprodx3j_ene
     end if
  end if
  call hprod_mkint1_sph(ctrue, orb, h0orb, h1orb, gorb)
  if (.not.exact3j) then
     call hprod_mkint2_sph(rho2, v2sph)
  else
     call hprod_mkint2_x3j(rho2, v2sph)
  end if

  ene_dcore = zero
  ene_core = zero
  if (ndcore > 0) ene_dcore = hprod_ene_dcore(ctrue, orb, h0orb, h1orb, gorb)
  ene_core = ene_fcore + ene_dcore

  if (.not.tdcc) then
     call hprod_cidiagx(cic)
  else
     call tdcc_ccsolve(int1e,int2e,cic)
  end if

!  stop "STOP for debug @ hprod_cidiag."

end subroutine hprod_cidiag
!#######################################################################
