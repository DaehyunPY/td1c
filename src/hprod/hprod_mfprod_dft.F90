!######################################################################
subroutine hprod_mfprod_dftj(v2j, wfn, hwfn, edftj)

  ! direct Coulomb

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfcore, nfun
  use mod_bas, only : ngrid, wgt
  use mod_const, only : zero, czero

  implicit none
  complex(c_double_complex), intent(in) :: v2j(1:ngrid)
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: hwfn(1:ngrid, 1:nfun)
  real(c_double), intent(out) :: edftj

  real(c_double) :: ene
  complex(c_double_complex) :: tmp1, tmp2
  integer(c_int) :: ifun, igrid, llg, ulg

  ene = zero

  !$omp parallel default(shared) private(ifun, igrid, llg, ulg, tmp1, tmp2) reduction(+:ene)
  !###########################
  call util_omp_disp(1, ngrid, llg, ulg)
  do ifun = 1, nfcore
     do igrid = llg, ulg
        tmp1 = v2j(igrid) * wfn(igrid, ifun)
        hwfn(igrid, ifun) = hwfn(igrid, ifun) + tmp1
     end do
  end do

  tmp2 = czero
  do ifun = nfcore + 1, nfun
     do igrid = llg, ulg
        tmp1 = v2j(igrid) * wfn(igrid, ifun)
        hwfn(igrid, ifun) = hwfn(igrid, ifun) + tmp1
        tmp2 = tmp2 + conjg(wfn(igrid, ifun)) * tmp1 * wgt(igrid)
     end do
  end do
  ene = ene + dble(tmp2)
  !###########################
  !$omp end parallel

  edftj = ene

end subroutine hprod_mfprod_dftj
!######################################################################
subroutine hprod_mfprod_dftx(rhoxc, wfn, hwfn, edftxc)

  ! eXchange and correlation

  use, intrinsic :: iso_c_binding
  use xc_f03_lib_m
  use mod_ormas, only : nfun
  use mod_bas, only : ngrid
  use mod_sph, only : nlat, wlat
  use mod_control, only : dft_type
  use mod_rad, only : nrad, wrad, xrad
  use mod_const, only : zero, two, half, pi, czero, third, fac_lda
  use mod_control, only : xc_funcx, xc_funcc, xc_infox, xc_infoc, func_idx, func_idc

  implicit none
  real(c_double), intent(in) :: rhoxc(1:(nrad-1), 1:nlat)
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: hwfn(1:(nrad-1), 1:nlat, 1:nfun)
  real(c_double), intent(out) :: edftxc

  real(c_double) :: ene, edftx, edftc
  complex(c_double_complex) :: tmp1
  integer(c_int) :: ifun, ilat, irad, llr, ulr
  real(c_double), allocatable :: rho(:,:)
  real(c_double), allocatable :: exc(:,:)
  real(c_double), allocatable :: vxc(:,:)
!libxc
  integer(c_int) :: npt
  npt = ngrid
!libxc

  allocate(rho(1:(nrad-1), 1:nlat))
  allocate(exc(1:(nrad-1), 1:nlat))
  allocate(vxc(1:(nrad-1), 1:nlat))
! rho(1:(nrad-1), 1:nlat) = rhoxc(1:(nrad-1), 1:nlat) * half
  rho(1:(nrad-1), 1:nlat) = rhoxc(1:(nrad-1), 1:nlat)
  exc(1:(nrad-1), 1:nlat) = zero
  vxc(1:(nrad-1), 1:nlat) = zero

  ene = zero
  if (mod(dft_type, 10) .ne. 0) then

     !MANUAL !$omp parallel default(shared) private(llr, ulr, tmp1) reduction(+:ene)
     !MANUAL !###########################
     !MANUAL call util_omp_disp(1, nrad-1, llr, ulr)
     !MANUAL tmp1 = czero
     !MANUAL do ilat = 1, nlat
     !MANUAL    do irad = llr, ulr
     !MANUAL       vxc(irad, ilat) = -fac_lda * rhoxc(irad, ilat) ** third
     !MANUAL       tmp1 = tmp1 + rhoxc(irad, ilat) * vxc(irad, ilat) * wrad(irad) * wlat(ilat) * xrad(irad)**two * (two * pi)
     !MANUAL    end do
     !MANUAL end do
     !MANUAL
     !MANUAL do ifun = 1, nfun
     !MANUAL    do ilat = 1, nlat
     !MANUAL       do irad = llr, ulr
     !MANUAL          hwfn(irad, ilat, ifun) = hwfn(irad, ilat, ifun) + vxc(irad, ilat) * wfn(irad, ilat, ifun)
     !MANUAL       end do
     !MANUAL    end do
     !MANUAL end do
     !MANUAL ene = ene + dble(tmp1) * 3.d+0/4.d+0
     !MANUAL !###########################
     !MANUAL !$omp end parallel

     !dfrep exc(1:(nrad-1), 1:nlat) = zero
     !dfrep vxc(1:(nrad-1), 1:nlat) = zero
     !dfrep call hprod_xc_ldax(1, ngrid, rho, exc, vxc)
     !dfrep write(6, "('exc_ldax from dfrep:')")
     !dfrep do ilat = 1, nlat
     !dfrep    do irad = 1, nrad - 1
     !dfrep       write(6, "(2i10, 2f20.10)") ilat, irad, exc(irad, ilat), vxc(irad, ilat)
     !dfrep    end do
     !dfrep end do

     exc(1:(nrad-1), 1:nlat) = zero
     vxc(1:(nrad-1), 1:nlat) = zero
     call xc_f03_lda_exc_vxc(xc_funcx, npt, rho, exc, vxc)
     exc(1:(nrad-1), 1:nlat) = exc(1:(nrad-1), 1:nlat) * rho(1:(nrad-1), 1:nlat)
     !debug write(6, "('exc_ldax from libxc:')")
     !debug do ilat = 1, nlat
     !debug    do irad = 1, nrad - 1
     !debug       write(6, "(2i10, 2f20.10)") ilat, irad, exc(irad, ilat), vxc(irad, ilat)
     !debug    end do
     !debug end do

     !$omp parallel default(shared) private(llr, ulr, tmp1) reduction(+:ene)
     !###########################
     call util_omp_disp(1, nrad-1, llr, ulr)
     tmp1 = czero
     do ilat = 1, nlat
        do irad = llr, ulr
           tmp1 = tmp1 + exc(irad, ilat) * wrad(irad) * wlat(ilat) * xrad(irad)**two * (two * pi)
        end do
     end do
     do ifun = 1, nfun
        do ilat = 1, nlat
           do irad = llr, ulr
              hwfn(irad, ilat, ifun) = hwfn(irad, ilat, ifun) + vxc(irad, ilat) * wfn(irad, ilat, ifun)
           end do
        end do
     end do
     ene = ene + dble(tmp1)
     !###########################
     !$omp end parallel
  end if
  edftx = ene

  ene = zero
  if (dft_type / 10 .ne. 0) then

     !dfrep exc(1:(nrad-1), 1:nlat) = zero
     !dfrep vxc(1:(nrad-1), 1:nlat) = zero
     !dfrep call hprod_xc_ldac(1, ngrid, rho, exc, vxc)
     !dfrep write(6, "('exc_ldac from dfrep:')")
     !dfrep do ilat = 1, nlat
     !dfrep    do irad = 1, nrad - 1
     !dfrep       write(6, "(2i10, 2f20.10)") ilat, irad, exc(irad, ilat), vxc(irad, ilat)
     !dfrep    end do
     !dfrep end do

     exc(1:(nrad-1), 1:nlat) = zero
     vxc(1:(nrad-1), 1:nlat) = zero
     call xc_f03_lda_exc_vxc(xc_funcc, npt, rho, exc, vxc)
     exc(1:(nrad-1), 1:nlat) = exc(1:(nrad-1), 1:nlat) * rho(1:(nrad-1), 1:nlat)
     !debug write(6, "('exc_ldac from libxc:')")
     !debug do ilat = 1, nlat
     !debug    do irad = 1, nrad - 1
     !debug       write(6, "(2i10, 2f20.10)") ilat, irad, exc(irad, ilat), vxc(irad, ilat)
     !debug    end do
     !debug end do

     !$omp parallel default(shared) private(llr, ulr, tmp1) reduction(+:ene)
     !###########################
     call util_omp_disp(1, nrad-1, llr, ulr)
     tmp1 = czero
     do ilat = 1, nlat
        do irad = llr, ulr
           tmp1 = tmp1 + exc(irad, ilat) * wrad(irad) * wlat(ilat) * xrad(irad)**two * (two * pi)
        end do
     end do
     do ifun = 1, nfun
        do ilat = 1, nlat
           do irad = llr, ulr
              hwfn(irad, ilat, ifun) = hwfn(irad, ilat, ifun) + vxc(irad, ilat) * wfn(irad, ilat, ifun)
           end do
        end do
     end do
     ene = ene + dble(tmp1)
     !###########################
     !$omp end parallel
  end if
  edftc = ene

  edftxc = edftx + edftc
  deallocate(vxc)
  deallocate(exc)
  deallocate(rho)

end subroutine hprod_mfprod_dftx
!######################################################################
