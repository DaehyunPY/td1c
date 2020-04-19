!######################################################################
subroutine hprod_mkrrad(rho_type, wfn, cic, rrad, rradpw)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun
  use mod_hprod, only : orb, orbg, den1

  implicit none
  integer(c_long), intent(in) :: rho_type
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: rrad(0:nrad)
  complex(c_double_complex), intent(out) :: rradpw(0:nrad, 0:lmax1)

  integer(c_long) :: sizel

  sizel = nbas * nfun;
  call ormas_mkden1(cic, den1)
  call hprod_mkrradp(rho_type, wfn, rrad, rradpw)

!old  call zcopy_omp(sizel, wfn, orb)
!old  if (rho_type == 0) then
!old     call bas_sph2ang1(orb, orbg)
!old  else if (rho_type == -1) then
!old     call bas_sph2ang1_fc(orb, orbg)
!old  else if (rho_type == +1) then
!old     call bas_sph2ang1_dyn(orb, orbg)
!old  end if
!old  call hprod_mkrradp_old(rho_type, orbg, rrad)

end subroutine hprod_mkrrad
!######################################################################
subroutine hprod_mkrrad0(rho_type, rad_ipx, wfn, cic, rrad, rradpw)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun
  use mod_hprod, only : orb, orbg, ovlp, den1

  implicit none
  integer(c_long), intent(in) :: rho_type
  real(c_double), intent(in) :: rad_ipx
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: rrad(0:nrad)
  complex(c_double_complex), intent(out) :: rradpw(0:nrad, 0:lmax1)

  integer(c_long) :: sizel, ifun

  sizel = nbas * nfun;
  call hprod_mkovlp(rad_ipx, wfn, wfn, ovlp);
  call ormas_mkden1_smat(ovlp, cic, den1)
  call hprod_mkrradp(rho_type, wfn, rrad, rradpw)

end subroutine hprod_mkrrad0
!######################################################################
subroutine hprod_mkrrad1(rho_type, rad_ipx, wfn, cic, rrad, rradpw)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun
  use mod_hprod, only : orb, orbg, ovlp, den1

  implicit none
  integer(c_long), intent(in) :: rho_type
  real(c_double), intent(in) :: rad_ipx
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: rrad(0:nrad)
  complex(c_double_complex), intent(out) :: rradpw(0:nrad, 0:lmax1)

  integer(c_long) :: sizel, ifun

  sizel = nbas * nfun;
  call hprod_mkovlp(rad_ipx, wfn, wfn, ovlp);
  ovlp = -ovlp
  do ifun = 1, nfun
     ovlp(ifun,ifun) = 1d0+ovlp(ifun,ifun)
  end do

  call ormas_mkden1_smat(ovlp, cic, den1)
  call hprod_mkrradp(rho_type, wfn, rrad, rradpw)

end subroutine hprod_mkrrad1
!######################################################################
subroutine hprod_mkrrad2(rho_type, rion, r2in, r2out, wfn, cic, rrad, rradpw)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, nelact, nact
  use mod_hprod, only : orb, orbg, ovlp, den1

  implicit none
  integer(c_long), intent(in) :: rho_type
  real(c_double), intent(in) :: rion, r2in, r2out
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: rrad(0:nrad)
  complex(c_double_complex), intent(out) :: rradpw(0:nrad, 0:lmax1)

  integer(c_long) :: ifun, iact
  real(c_double) :: neeff, ene1, lfield(1:3,1:3)
  real(c_double), external :: hprod_energy1

  call hprod_mkovlp_shell(r2in, r2out, wfn, wfn, ovlp);
  call ormas_mkden1_smat(ovlp, cic, den1)
  call hprod_mkrradp(rho_type, wfn, rrad, rradpw)

  if (nelact(3) == 2) then
     lfield = 0d0
     neeff = 0d0
     do iact = 1, nact
        neeff = neeff + den1(iact,iact)
     end do
     den1 = den1*nelact(3)/neeff
     neeff = 0d0
     do iact = 1, nact
        neeff = neeff + den1(iact,iact)
     end do
     ovlp = den1
     ene1 = hprod_energy1(lfield, rion, wfn, cic)
     write(6,"('hprod_mkrrad3: ene1 = ', 2f20.10)") neeff, ene1*0.5d0
  end if

end subroutine hprod_mkrrad2
!######################################################################
subroutine hprod_mkrrad3(rho_type, rion, k2in, k2out, wfn, cic, rrad, rradpw)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_bas, only : nbas
  use mod_sph, only : lmax1
  use mod_ormas, only : nfun, nelact, nact
  use mod_hprod, only : orb, orbg, ovlp, den1

  implicit none
  integer(c_long), intent(in) :: rho_type
  real(c_double), intent(in) :: rion, k2in, k2out
  complex(c_double_complex), intent(in) :: wfn(1:*)
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: rrad(0:nrad)
  complex(c_double_complex), intent(out) :: rradpw(0:nrad, 0:lmax1)

  integer(c_long) :: ifun, iact
  real(c_double) :: neeff, ene1, lfield(1:3,1:3)
  real(c_double), external :: hprod_energy1

  call hprod_mkovlp_kshell(k2in, k2out, wfn, wfn, ovlp);
  call ormas_mkden1_smat(ovlp, cic, den1)
  call hprod_mkrradp(rho_type, wfn, rrad, rradpw)

  if (nelact(3) == 2) then
     lfield = 0d0
     neeff = 0d0
     do iact = 1, nact
        neeff = neeff + den1(iact,iact)
     end do
     den1 = den1*nelact(3)/neeff
     neeff = 0d0
     do iact = 1, nact
        neeff = neeff + den1(iact,iact)
     end do
     ovlp = den1
     ene1 = hprod_energy1(lfield, rion, wfn, cic)
     write(6,"('hprod_mkrrad3: ene1 = ', 2f20.10)") neeff, ene1*0.5d0
  end if

end subroutine hprod_mkrrad3
!######################################################################
subroutine hprod_mkrradp(rho_type, orb, rrad, rradpw)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun, nfcore, ncore, nact
  use mod_bas, only : mval
  use mod_rad, only : nrad, nradfc, wrad
  use mod_sph, only : lmax1
  use mod_const, only : one, czero, ctwo
  use mod_hprod, only : den1

  implicit none
  integer(c_long), intent(in) :: rho_type
  complex(c_double_complex), intent(in) :: orb(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: rrad(0:nrad)
  complex(c_double_complex), intent(out) :: rradpw(0:nrad, 0:lmax1)

  integer(c_long) :: llr, ulr, ifun, jfun, iact, jact, l, irad
  complex(c_double_complex) :: num_elec, tmp, numep, fac, facp, facm

  num_elec = czero

  !$omp parallel default(shared) private(tmp, numep, fac, facp, facm, iact, jact, llr, ulr) reduction(+:num_elec)
  !###########################  
  call util_omp_disp(1, nrad - 1, llr, ulr)

  numep = czero
  rrad(llr:ulr) = czero
  rradpw(llr:ulr, 0:lmax1) = czero

  if (rho_type <= 0) then
     do ifun = 1, nfcore
        do l = 0, lmax1
           do irad = llr, ulr
              tmp = conjg(orb(irad, l, ifun)) &
                        * orb(irad, l, ifun) * ctwo
              rradpw(irad, l) = rradpw(irad, l) + tmp / wrad(irad)
              numep = numep + tmp
           end do
        end do
     end do
  end if

  if (rho_type >= 0) then
     do ifun = nfcore + 1, ncore
        do l = 0, lmax1
           do irad = llr, ulr
              tmp = conjg(orb(irad, l, ifun)) &
                        * orb(irad, l, ifun) * ctwo
              rradpw(irad, l) = rradpw(irad, l) + tmp / wrad(irad)
              numep = numep + tmp
           end do
        end do
     end do

     do ifun = ncore + 1, nfun
        do jfun = ncore + 1, nfun
           if (mval(ifun) == mval(jfun)) then
              iact = ifun - ncore
              jact = jfun - ncore
              fac = den1(iact, jact)
              do l = 0, lmax1
                 do irad = llr, ulr
                    tmp = conjg(orb(irad, l, jfun)) &
                              * orb(irad, l, ifun) * fac
                    rradpw(irad, l) = rradpw(irad, l) + tmp / wrad(irad)
                    numep = numep + tmp
                 end do
              end do
           end if
        end do
     end do
  end if

  do l = 0, lmax1
     do irad = llr, ulr
        rrad(irad) = rrad(irad) + rradpw(irad, l)
     end do
  end do

  num_elec = num_elec + numep
  !###########################
  !$omp end parallel

  !DEBUG
  write(6, "('# hprod_mkrrad: num_elec = ', 2f20.10)") num_elec
  !DEBUG

end subroutine hprod_mkrradp
!######################################################################
subroutine hprod_mkrradp_old(rho_type, orbg, rrad)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun, nfcore, ncore, nact
  use mod_bas, only : mval
  use mod_rad, only : nrad, nradfc, wrad
  use mod_sph, only : nlat, wang
  use mod_const, only : czero, ctwo
  use mod_hprod, only : den1

  implicit none
  integer(c_long), intent(in) :: rho_type
  complex(c_double_complex), intent(in) :: orbg(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(out) :: rrad(1:(nrad-1))

  integer(c_long) :: ifun_ll, ifun_ul, llr, ulr, ifun, jfun, iact, jact, ilat, irad
  complex(c_double_complex) :: num_elec, tmp, numep, fac

  if (rho_type == 0) then
     ifun_ll = 1
     ifun_ul = nfun
  else if (rho_type == -1) then
     ifun_ll = 1
     ifun_ul = nfcore
  else if (rho_type == +1) then
     ifun_ll = nfcore + 1
     ifun_ul = nfun
  end if

  !$omp parallel default(shared) private(tmp, numep, fac, iact, jact, llr, ulr) reduction(+:num_elec)
  !###########################  
  numep = czero
  call util_omp_disp(1, nrad - 1, llr, ulr)

  num_elec = czero
  rrad(llr:ulr) = czero

  do ifun = ifun_ll, min(ncore, ifun_ul)
     do ilat = 1, nlat
        fac = ctwo * wang(ilat)
        do irad = llr, ulr
           tmp = conjg(orbg(irad, ilat, ifun)) &
                     * orbg(irad, ilat, ifun) * fac
           rrad(irad) = rrad(irad) + tmp / wrad(irad)
           numep = numep + tmp
        end do
     end do
  end do

  do ifun = max(ncore + 1, ifun_ll), ifun_ul
  do jfun = max(ncore + 1, ifun_ll), ifun_ul
     if (mval(ifun) == mval(jfun)) then
        iact = ifun - ncore
        jact = jfun - ncore
        do ilat = 1, nlat
           fac = den1(iact, jact) * wang(ilat)
           do irad = llr, ulr
              tmp = conjg(orbg(irad, ilat, jfun)) &
                        * orbg(irad, ilat, ifun) * fac
              rrad(irad) = rrad(irad) + tmp / wrad(irad)
              numep = numep + tmp
           end do
        end do
     end if
  end do
  end do

  num_elec = num_elec + numep
  !###########################
  !$omp end parallel

  !DEBUG
  write(6, "('# hprod_mkrrad: num_elec = ', 2f20.10)") num_elec
  !DEBUG

end subroutine hprod_mkrradp_old
!######################################################################
