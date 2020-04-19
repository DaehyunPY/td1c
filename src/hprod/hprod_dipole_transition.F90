!######################################################################
subroutine hprod_orbin(lfield, wfn, orb, orbg)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas, ngrid
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : orbe, orbo

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(in) :: wfn(1:nbas, 1:*)
  complex(c_double_complex), intent(out) :: orb(1:nbas, 1:*)
  complex(c_double_complex), intent(out) :: orbg(1:ngrid, 1:*)

  integer(c_int) :: ifun
  integer(c_int) :: llr, ulr, irad, ilat
  complex(c_double_complex) :: faclv

!  call zcopy_omp(nbas*(nfun-nfcore), wfn(1,nfcore+1), orb(1,nfcore+1))
  call zcopy(nbas*(nfun-nfcore), wfn(1,nfcore+1), 1, orb(1,nfcore+1), 1)
  call bas_sph2ang1_dyn(orb, orbg)
  call hprod_fcorb(lfield, orb, orbg)

! ##### 3j selection rule #####
  !if (exact3j) then
  !   call bas_sph2ang1_fc3j(orb, orbe, orbo)
  !   call bas_sph2ang1_dyn3j(orb, orbe, orbo)  
  !end if
! ##### 3j selection rule #####

end subroutine hprod_orbin
!######################################################################
subroutine hprod_fcorb(lfield, orb, orbg)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : nbas
  use mod_control, only : igauge
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : orb0

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(out) :: orb(1:*)
  complex(c_double_complex), intent(out) :: orbg(1:*)

  if (nfcore == 0) return

  if (igauge == 0) then
!     call zcopy_omp(nbas*nfcore, orb0, orb)
     call zcopy(nbas*nfcore, orb0, 1, orb, 1)
     call bas_sph2ang1_fc(orb, orbg)
  else
     call hprod_fcorbg(lfield, orbg)
     call bas_ang2sph1_fc(orbg, orb)
  end if

end subroutine hprod_fcorb
!######################################################################
subroutine hprod_op1trans(lfield, orbg)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : xrad, nrad, nradfc
  use mod_sph, only : lmax1, nlat, cost
  use mod_control, only : igauge
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : orb0
  use mod_const, only : iunit

  implicit none
  real(c_double), intent(in) :: lfield(1:3, 1:3)
  complex(c_double_complex), intent(out) :: horb(1:(nrad-1), 0:lmax1, 1:*)

  integer(c_int) :: ifun
  integer(c_int) :: llr, ulr, irad, ilat
  complex(c_double_complex) :: faclv

  if (nfcore == 0) return

  call bas_sph2ang1_fc(horb, op1trans_orbg)
  if (igauge == 1) then
     !$omp parallel default(shared) private(faclv, llr, ulr)
     !###########################
     call util_omp_disp(1, nradfc, llr, ulr)
     do ilat = 1, nlat
        do irad = llr, ulr
           faclv = exp(-iunit * lfield(3, 3) * xrad(irad) * cost(ilat))
           do ifun = 1, nfcore
              orbg(irad, ilat, ifun) = orbg(irad, ilat, ifun) * faclv
           end do
        end do
     end do
     !###########################
     !$omp end parallel
  end if

end subroutine hprod_fcorbg
!######################################################################
