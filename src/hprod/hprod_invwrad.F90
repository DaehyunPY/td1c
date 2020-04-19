!######################################################################
subroutine hprod_invwrad_all(wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_ormas, only : nfun

  implicit none
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:*)

  call hprod_invwrad(wfn, 1, nfun, 1, nrad - 1)

end subroutine hprod_invwrad_all
!######################################################################
subroutine hprod_invwrad_dyn(wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_ormas, only : nfcore, nfun

  implicit none
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:*)

  call hprod_invwrad(wfn, nfcore + 1, nfun, 1, nrad - 1)

end subroutine hprod_invwrad_dyn
!######################################################################
subroutine hprod_invwrad_fc(wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore

  implicit none
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:*)

  if (nfcore == 0) return
  call hprod_invwrad(wfn, 1, nfcore, 1, nradfc)

end subroutine hprod_invwrad_fc
!######################################################################
subroutine hprod_invwrad_fc1(wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore

  implicit none
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:*)

  if (nfcore1 == 0) return
  call hprod_invwrad(wfn, nfcore2+1, nfcore, 1, nradfc)

end subroutine hprod_invwrad_fc1
!######################################################################
subroutine hprod_invwrad_fc2(wfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore

  implicit none
  complex(c_double_complex), intent(inout) :: wfn(1:nbas, 1:*)

  if (nfcore2 == 0) return
  call hprod_invwrad(wfn, 1, nfcore2, 1, nradfc)

end subroutine hprod_invwrad_fc2
!######################################################################
subroutine hprod_invwrad(wfn, llfun, ulfun, llr, ulr)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun, llr, ulr
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  integer(c_long) :: llrp, ulrp

  !$omp parallel default(shared) private(llrp, ulrp)
  !###########################
  call util_omp_disp(llr, ulr, llrp, ulrp)
  call hprod_invwradp(wfn, llfun, ulfun, llrp, ulrp)
  !###########################
  !$omp end parallel

end subroutine hprod_invwrad
!######################################################################
subroutine hprod_invwradp(wfn, llfun, ulfun, llrp, ulrp)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, wrad

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun, llrp, ulrp
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_long) :: ifun, l, m, irad

  do ifun = llfun, ulfun
     m = mval(ifun)
     do l = abs(m), lmax1 - 1
        do irad = llrp, ulrp
           wfn(irad, l, ifun) = wfn(irad, l, ifun) / wrad(irad)
        end do
     end do
  end do

end subroutine hprod_invwradp
!######################################################################
