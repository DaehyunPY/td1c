!######################################################################
subroutine hprod_mkv2mf2()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : rho2, v2sph, v2ang

  implicit none
  character(len=6) :: v2_type
  v2_type = 'tot'
  call hprod_mkv2mf2_poisson(v2_type, rho2, v2sph)
  call bas_sph2ang2(v2sph, v2ang);
  call hprod_mkv2mf2_herm(v2_type, v2sph, v2ang)

end subroutine hprod_mkv2mf2
!######################################################################
subroutine hprod_mkv2mf2_dyn()

  use, intrinsic :: iso_c_binding
  use mod_hprod, only : rho2, v2sph, v2ang

  implicit none
  character(len=6) :: v2_type
  v2_type = 'dyn'
  call hprod_mkv2mf2_poisson(v2_type, rho2, v2sph)
  call bas_sph2ang2_dyn(v2sph, v2ang);
  call hprod_mkv2mf2_herm(v2_type, v2sph, v2ang)
  v2_type = 'fc1dyn'
  call hprod_mkv2mf2_poisson(v2_type, rho2, v2sph)
  call bas_sph2ang2_fc1dyn(v2sph, v2ang);
  call hprod_mkv2mf2_herm(v2_type, v2sph, v2ang)

!  write(6, "('make v2mf for fc2dyn for debug...')")
!  v2_type = 'fc2dyn'
!  call hprod_mkv2mf2_poisson(v2_type, rho2, v2sph)
!  call bas_sph2ang2_fc2dyn(v2sph, v2ang);
!  call hprod_mkv2mf2_herm(v2_type, v2sph, v2ang)

end subroutine hprod_mkv2mf2_dyn
!######################################################################
subroutine hprod_mkv2mf2_poisson(v2_type, rho2, v2sph)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : nfun, nfcore1, nfcore2, nfcore
  use mod_rad, only : nrad, nradfc
  use mod_const, only : czero
  use mod_sph, only : lmax2
  use mod_bas, only : mval

  implicit none
  character(len=*), intent(in) :: v2_type
  complex(c_double_complex), intent(in) ::     rho2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int) :: dim, ifun_ll, ifun_ul, jfun_ll, jfun_ul
  integer(c_int) :: ifun, jfun, lji, mji, ijl, nijl, mapijl(1:3, 1:nfun*nfun*(lmax2+1))

  if (trim(v2_type) == 'tot') then
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfun
     jfun_ll = 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fcfc') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore
     jfun_ll = 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc1fc1') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore2 + 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc2fc2') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = 1
     jfun_ul = nfcore2
  else if (trim(v2_type) == 'fc1dyn') then
     dim = nradfc
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fc2dyn') then
     dim = nradfc
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'dyn') then
     dim = nrad - 1
     ifun_ll = nfcore + 1
     ifun_ul = nfun
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else
     stop "bad v2_type in hprod_mkv2mf2_poisson"
  end if

  nijl = 0
  do ifun = ifun_ll, ifun_ul
     do jfun = max(jfun_ll, ifun), max(jfun_ul, ifun)
        do lji = 0, lmax2
           mji = -mval(ifun) + mval(jfun)
           if (lji < abs(mji)) then
              v2sph(1:(nrad-1), lji, jfun, ifun) = czero
           else
              nijl = nijl + 1
              mapijl(1, nijl) = ifun
              mapijl(2, nijl) = jfun
              mapijl(3, nijl) = lji
           end if
        end do
     end do
  end do

!debug  !$omp parallel default(shared) private(ifun,jfun,lji,mji)
!debug  !$omp do
!debug  do ijl = 1, nijl
!debug     ifun = mapijl(1, ijl)
!debug     jfun = mapijl(2, ijl)
!debug     lji  = mapijl(3, ijl)
!debug     mji = -mval(ifun) + mval(jfun)
!debug     call hprod_mkv2mf2_poissonp(rho2(1,lji,jfun,ifun), v2sph(1,lji,jfun,ifun), dim, lji, mji)
!debug  end do
!debug  !$omp end do
!debug  !$omp end parallel
     call hprod_mkv2mf2_poissonp2(rho2, v2sph, dim, nijl, mapijl)

end subroutine hprod_mkv2mf2_poisson
!######################################################################
subroutine hprod_mkv2mf2_poissonp(rho2, v2sph, dim, l, m)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : ndvr, nrad, xrad
  use mod_const, only : zero, one, czero, iunit
  use mod_bas, only : d2ll, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1

  implicit none
  integer(c_int), intent(in) :: dim, l, m
  complex(c_double_complex), intent(in) ::     rho2(1:(nrad-1))
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1))

  real(c_double) :: tmpr, tmpi, d2fac1
  real(c_double), allocatable :: rrho2(:,:)
  integer(c_int) :: irad, ld, info

  ld = ndvr + 1
  d2fac1 = one / xrad(dim + 1) ** (2 * l + 1)
  allocate(rrho2(1:dim, 1:2))

  ! r.h.s and BC
  tmpr = zero
  tmpi = zero
  do irad = 1, dim
     rrho2(irad, 1) = dble( rho2(irad))
     rrho2(irad, 2) = aimag(rho2(irad))
     tmpr = tmpr + rrho2(irad, 1) * bas_d2rpl0(irad, l)
     tmpi = tmpi + rrho2(irad, 2) * bas_d2rpl0(irad, l)
     rrho2(irad, 1) = rrho2(irad, 1) * bas_d2invr(irad, l)
     rrho2(irad, 2) = rrho2(irad, 2) * bas_d2invr(irad, l)
  end do
  tmpr = tmpr * d2fac1
  tmpi = tmpi * d2fac1

  ! Poisson's equation with vanishing BCs
  call DPBTRS('L', dim, ndvr, 2, d2ll(1,1,l), ld, rrho2, dim, info)
  if (info /= 0) then
     write(6, "('hprod_poisson: dpbtrs bad info. ', i5)") info
     stop
  end if

  ! solution with the correct BC
  do irad = 1, dim
     v2sph(irad) = (rrho2(irad, 1) + tmpr * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) &
                 + (rrho2(irad, 2) + tmpi * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) * iunit
  end do

  deallocate(rrho2)

end subroutine hprod_mkv2mf2_poissonp
!######################################################################
subroutine hprod_mkv2mf2_poissonp2(rho2, v2sph, dim, nijl, mapijl)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : lmax2
  use mod_ormas, only : nfun
  use mod_rad, only : ndvr, nrad, xrad
  use mod_const, only : zero, one, czero, iunit
  use mod_bas, only : d2ll, bas_d2fac2, bas_d2invr, bas_d2rpl0, bas_d2rpl1

  implicit none
  integer(c_int), intent(in) :: dim
  complex(c_double_complex), intent(in) ::     rho2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  integer(c_int), intent(in) :: nijl, mapijl(1:3, 1:nfun*nfun*(lmax2+1))

  real(c_double) :: tmpr, tmpi, d2fac1(0:lmax2)
  real(c_double), allocatable :: rrho2(:,:,:)
  integer(c_int) :: irad, ld, info, iproc, nproc
  integer(c_int) :: ifun, jfun, l, m, ijl
  integer(c_int), external :: util_omp_nproc
  integer(c_int), external :: util_omp_iproc

  ld = ndvr + 1
  nproc = util_omp_nproc()
  allocate(rrho2(1:dim, 1:2, 0:(nproc-1)))
  do l = 0, lmax2
     d2fac1(l) = one / xrad(dim + 1) ** (2 * l + 1)
  end do

  !$omp parallel default(shared) private(ifun,jfun,l,m,tmpr,tmpi,info,iproc)
  !$omp do
  do ijl = 1, nijl
     ifun = mapijl(1, ijl)
     jfun = mapijl(2, ijl)
     l =    mapijl(3, ijl)
     m = -mval(ifun) + mval(jfun)
     iproc = util_omp_iproc()  

     ! r.h.s and BC
     tmpr = zero
     tmpi = zero
     do irad = 1, dim
        rrho2(irad, 1, iproc) = dble( rho2(irad,l,jfun,ifun))
        rrho2(irad, 2, iproc) = aimag(rho2(irad,l,jfun,ifun))
        tmpr = tmpr + rrho2(irad, 1, iproc) * bas_d2rpl0(irad, l)
        tmpi = tmpi + rrho2(irad, 2, iproc) * bas_d2rpl0(irad, l)
        rrho2(irad, 1, iproc) = rrho2(irad, 1, iproc) * bas_d2invr(irad, l)
        rrho2(irad, 2, iproc) = rrho2(irad, 2, iproc) * bas_d2invr(irad, l)
     end do
     tmpr = tmpr * d2fac1(l)
     tmpi = tmpi * d2fac1(l)

     ! Poisson's equation with vanishing BCs
     call DPBTRS('L', dim, ndvr, 2, d2ll(1,1,l), ld, rrho2(1,1,iproc), dim, info)
     if (info /= 0) then
        write(6, "('hprod_poisson: dpbtrs bad info. ', i5)") info
        stop
     end if

     ! solution with the correct BC
     do irad = 1, dim
        v2sph(irad,l,jfun,ifun) = (rrho2(irad, 1, iproc) + tmpr * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) &
                                + (rrho2(irad, 2, iproc) + tmpi * bas_d2rpl1(irad, l)) * bas_d2fac2(irad, l) * iunit
     end do
  end do
  !$omp end do
  !$omp end parallel

  deallocate(rrho2)

end subroutine hprod_mkv2mf2_poissonp2
!######################################################################
subroutine hprod_mkv2mf2_herm(v2_type, v2sph, v2ang)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : lmax2, nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore, nfcore1, nfcore2, nfun

  implicit none
  character(len=6), intent(in) :: v2_type
  complex(c_double_complex), intent(inout) :: v2sph(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: v2ang(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)

  integer(c_int) :: l, ilat, ifun, jfun, irad
  integer(c_int) :: dim, ifun_ll, ifun_ul, jfun_ll, jfun_ul, llr, ulr

  if (trim(v2_type) == 'tot') then
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfun
     jfun_ll = 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fcfc') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore
     jfun_ll = 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc1fc1') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore2 + 1
     jfun_ul = nfcore
  else if (trim(v2_type) == 'fc2fc2') then
!BUG dim = nradfc
     dim = nrad - 1
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = 1
     jfun_ul = nfcore2
  else if (trim(v2_type) == 'fc1dyn') then
     dim = nradfc
     ifun_ll = nfcore2 + 1
     ifun_ul = nfcore
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'fc2dyn') then
     dim = nradfc
     ifun_ll = 1
     ifun_ul = nfcore2
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else if (trim(v2_type) == 'dyn') then
     dim = nrad - 1
     ifun_ll = nfcore + 1
     ifun_ul = nfun
     jfun_ll = nfcore + 1
     jfun_ul = nfun
  else
     stop "bad v2_type in hprod_mkv2mf2_herm"
  end if

  !$omp parallel default(shared) private(llr, ulr)
  !###########################
  call util_omp_disp(1, dim, llr, ulr)
  do ifun = ifun_ll, ifun_ul
     do jfun = max(jfun_ll, ifun + 1), jfun_ul
        do l = abs(-mval(ifun)+mval(jfun)), lmax2
           do irad = llr, ulr
              v2sph(irad, l, ifun, jfun) = conjg(v2sph(irad, l, jfun, ifun))
           end do
        end do
        do ilat = 1, nlat
           do irad = llr, ulr
              v2ang(irad, ilat, ifun, jfun) = conjg(v2ang(irad, ilat, jfun, ifun))
           end do
        end do
     end do
  end do
  !###########################
  !$omp end parallel

end subroutine hprod_mkv2mf2_herm
!######################################################################
