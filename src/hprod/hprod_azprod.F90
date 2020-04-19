!######################################################################
subroutine hprod_azprod_all(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_ormas, only : nfun, nfcore

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

! call hprod_azprod(zfac, wfn, hwfn, 1, nfun, 1, nrad - 1)
 if (nfcore > 0) call hprod_azprod_fc(zfac, wfn, hwfn)
 call hprod_azprod_dyn(zfac, wfn, hwfn)

end subroutine hprod_azprod_all
!######################################################################
subroutine hprod_azprod_dyn(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nrad
  use mod_ormas, only : nfcore, nfun

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  call hprod_azprod(zfac, wfn, hwfn, nfcore + 1, nfun, 1, nrad - 1)

end subroutine hprod_azprod_dyn
!######################################################################
subroutine hprod_azprod_fc(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore == 0) return
  call hprod_azprod(zfac, wfn, hwfn, 1, nfcore, 1, nradfc)

end subroutine hprod_azprod_fc
!######################################################################
subroutine hprod_azprod_fc1(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore1 == 0) return
  call hprod_azprod(zfac, wfn, hwfn, nfcore2+1, nfcore, 1, nradfc)

end subroutine hprod_azprod_fc1
!######################################################################
subroutine hprod_azprod_fc2(zfac, wfn, hwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_rad, only : nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore

  implicit none
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:nbas, 1:*)

  if (nfcore2 == 0) return
  call hprod_azprod(zfac, wfn, hwfn, 1, nfcore2, 1, nradfc)

end subroutine hprod_azprod_fc2
!######################################################################
subroutine hprod_azprod(zfac, wfn, hwfn, llfun, ulfun, llr, ulr)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun, llr, ulr
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)
  integer(c_long) :: llrp, ulrp

  !$omp parallel default(shared) private(llrp, ulrp)
  !###########################
  call util_omp_disp(llr, ulr, llrp, ulrp)
  call hprod_azprodp(zfac, wfn, hwfn, llfun, ulfun, llrp, ulrp)
  !###########################
  !$omp end parallel

end subroutine hprod_azprod
!######################################################################
subroutine hprod_azprodp(zfac, wfn, hwfn, llfun, ulfun, llrp, ulrp)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_bas, only : mval, bas_azfac

  implicit none
  integer(c_long), intent(in) :: llfun, ulfun, llrp, ulrp
  complex(c_double_complex), intent(in) :: zfac
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 0:lmax1, 1:*)

  complex(c_double_complex) :: tmp
  integer(c_long) :: ifun, l, m, irad

  do ifun = llfun, ulfun
     m = mval(ifun)
     do l = abs(m), lmax1 - 1
        do irad = llrp, ulrp
           tmp = zfac * bas_azfac(irad, l, m)
           hwfn(irad, l,     ifun) = hwfn(irad, l,     ifun) + wfn(irad, l + 1, ifun) * tmp
           hwfn(irad, l + 1, ifun) = hwfn(irad, l + 1, ifun) + wfn(irad, l,     ifun) * tmp
        end do
     end do
  end do

end subroutine hprod_azprodp
!######################################################################
