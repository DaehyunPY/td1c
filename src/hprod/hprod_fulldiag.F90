!######################################################################
!
! Construct and diagonalize RHF hamiltonian
!
subroutine hprod_fulldiag(wfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_ormas, only : neltot
  use mod_const, only : czero

  implicit none
  ! ### args ###
  complex(c_double_complex), intent(inout) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  ! ### local ###
  integer(c_int) :: dim
  complex(c_double_complex), allocatable :: fock(:,:,:)

  dim = nrad - 1
  allocate(fock(1:dim, 1:dim, 0:lmax1))
  fock(1:dim, 1:dim, 0:lmax1) = czero

  call hprod_fulldiag_fock1(fock)
  if (neltot(3) > 1) then
     call hprod_fulldiag_fock2(wfn, fock)
  end if
  call hprod_fulldiag_diag(fock)

  deallocate(fock)

end subroutine hprod_fulldiag
!######################################################################
subroutine hprod_fulldiag_fock1(fock)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : znuc
  use mod_rad, only : nrad, ndvr, xrad, radk
  use mod_sph, only : lmax1
  use mod_const, only : one, half, czero

  implicit none

  ! ### args ###
  complex(c_double_complex), intent(inout) :: fock(1:(nrad-1), 1:(nrad-1), 0:lmax1)

  ! ### local ###
  real(c_double) :: oor, zor
  integer(c_int) :: dim, nsub, irad, jrad, l, l2, jll, jul, ib1

  dim = nrad - 1
  nsub = ndvr

  do l = 0, lmax1
     l2 = l * (l + 1)
     do irad = 1, nrad - 1
        oor = one / xrad(irad)
        zor = oor * (-dble(znuc) + dble(l2) * oor * half);

        fock(irad, irad, l) = fock(irad, irad, l) + zor

        jll = max(1,   irad - nsub)
        jul = min(dim, irad)
        do jrad = jll, jul
           ib1 = nsub + 1 + irad - jrad
           fock(irad, jrad, l) = fock(irad, jrad, l) + radk(ib1, jrad)
           if (jrad < irad) fock(jrad, irad, l) = fock(irad, jrad, l)
        end do
     end do
  end do

end subroutine hprod_fulldiag_fock1
!######################################################################
subroutine hprod_fulldiag_fock2(wfn, fock)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: fock(1:(nrad-1), 1:(nrad-1), 0:lmax1)

  call hprod_fulldiag_fock2j(wfn, fock)
  call hprod_fulldiag_fock2k(wfn, fock)

end subroutine hprod_fulldiag_fock2
!######################################################################
subroutine hprod_fulldiag_fock2j(wfn, fock)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : d2ll
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nocc
  use mod_const, only : PI, czero, two, four

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: fock(1:(nrad-1), 1:(nrad-1), 0:lmax1)

  integer(c_int) :: ifun, irad, l
  complex(c_double_complex), allocatable :: rho(:)

  allocate(rho(1:(nrad-1)))

  ! Assume J is spherical
  rho(1:(nrad-1)) = czero
  do ifun = 1, nocc
     do l = 0, lmax1
        do irad = 1, nrad - 1
           rho(irad) = rho(irad) + wfn(irad, l, ifun) * conjg(wfn(irad, l, ifun))
        end do
     end do
  end do
  call hprod_fulldiag_poisson1(0, d2ll, rho)
  ! Assume J is spherical

  !debug
  ! do irad = 1, nrad - 1
  !    write(6, "('hprod_fulldiag_fock2j: ', i5, f10.5, 3f20.10)") &
  !         irad, xrad(irad), wrad(irad), dble(rho(irad)), dble(v2sph(irad))
  ! end do
  !debug

  do l = 0, lmax1
     do irad = 1, nrad - 1
        ! ##### 1 / 4 PI ##############################################
        fock(irad, irad, l) = &
        fock(irad, irad, l) + rho(irad) / (four * PI) * two
        ! #############################################################
     end do
  end do

  deallocate(rho)

end subroutine hprod_fulldiag_fock2j
!######################################################################
subroutine hprod_fulldiag_fock2k(wfn, fock)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : d2ll
  use mod_rad, only : nrad
  use mod_sph, only : lmax1
  use mod_ormas, only : nocc
  use mod_const, only : PI, czero, two, four

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:*)
  complex(c_double_complex), intent(inout) :: fock(1:(nrad-1), 1:(nrad-1), 0:lmax1)

  integer(c_int) :: jfun, irad, jrad, l, ll, jtype
  real(c_double) :: fac
  real(c_double), external :: util_3j000
  real(c_double), parameter :: small = 1.D-15
  complex(c_double_complex), allocatable :: rho(:)

  allocate(rho(1:(nrad-1)))

  do l = 0, lmax1
     do jfun = 1, nocc
        jtype = 0 ! s-only
        do jrad = 1, nrad - 1
           do ll = abs(l - jtype), l + jtype

              fac = util_3j000(jtype, ll, l)
!              fac = (-1.d+0) ** l / sqrt(2.d+0 * l + 1.d+0)

              if (abs(fac) > small) then
                 fac = fac * fac * (2 * jtype + 1)
                 !do irad = 1, nrad - 1
                 !if (irad == jrad) then
                 irad = jrad
                 rho(1:(nrad-1)) = czero
                 rho(irad) = fac * conjg(wfn(irad, 0, jfun))
                 !end if
                 !end do
!                call hprod_fulldiag_poisson1(ll, d2ll(1,1,ll), rho)
                 call hprod_fulldiag_poisson1(ll, d2ll(1:,1:,ll), rho)

                 do irad = 1, nrad - 1
                    ! ##### 1 / 4 PI #############################################
                    fock(irad, jrad, l) = &
                    fock(irad, jrad, l) - wfn(irad, 0, jfun) * rho(irad) &
                    / (PI * PI) !* (2 * ll + 1) ???????????????????????????????
                    ! ############################################################
                 end do
              end if
           end do
        end do
     end do
  end do

  deallocate(rho)

end subroutine hprod_fulldiag_fock2k
!######################################################################
subroutine hprod_fulldiag_poisson1(lval, d2ll, rho)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : ndvr, nrad, xrad, wrad
  use mod_const, only : PI, zero, one, half, four, czero, iunit

  implicit none

  ! ### args ###
  integer(c_int), intent(in) :: lval
  real(c_double), intent(in) :: d2ll(1:(ndvr+1), 1:(nrad-1))
  complex(c_double_complex), intent(inout) :: rho(1:(nrad-1))

  ! ### local ###
  integer(c_int) :: irad, dim, nsub, ld, info
  real(c_double) :: r2pi, tmpr, tmpi
  real(c_double), allocatable :: rho2(:,:)
  real(c_double), allocatable :: fac1(:)   ! 4*pi / (2*l+1)
  real(c_double)              :: fac2      ! 1 / R**(2*l+1)
  real(c_double), allocatable :: invr(:)   ! (2*l+1) / r / sqrt(w)
  real(c_double), allocatable :: rpowl(:)  ! r ** l
  real(c_double), allocatable :: rpowl1(:) ! r ** (l+1) * sqrt(w)

  dim = nrad - 1
  nsub = ndvr
  ld = nsub + 1
  r2pi = half / PI

  allocate(rho2(1:dim, 1:2))
  allocate(fac1(1:dim))
  allocate(invr(1:dim))
  allocate(rpowl(1:dim))
  allocate(rpowl1(1:dim))

  fac2 = one / xrad(nrad) ** (2 * lval + 1)

  do irad = 1, dim
     rpowl(irad) = xrad(irad) ** lval
     invr(irad) = (2 * lval + 1) / (xrad(irad) * sqrt(wrad(irad)))
     rpowl1(irad) = xrad(irad) ** (lval + 1) * sqrt(wrad(irad))
     fac1(irad) = four * pi / (2 * lval + 1) / (xrad(irad) * sqrt(wrad(irad)))
  end do

  tmpr = zero
  tmpi = zero
  do irad = 1, dim
!old rho2(irad, 1) = dble(rho(irad))
!old rho2(irad, 2) = aimag(rho(irad))
!new factor r2pi comes from separately normalized lat and phi functions
     rho2(irad, 1) = dble(rho(irad))  * r2pi
     rho2(irad, 2) = aimag(rho(irad)) * r2pi

     tmpr = tmpr + rho2(irad, 1) * rpowl(irad)
     tmpi = tmpi + rho2(irad, 2) * rpowl(irad)
     rho2(irad, 1) = rho2(irad, 1) * invr(irad)
     rho2(irad, 2) = rho2(irad, 2) * invr(irad)
  end do
  tmpr = tmpr * fac2
  tmpi = tmpi * fac2

  call DPBTRS('L', dim, nsub, 2, d2ll, ld, rho2, dim, info)
  if (info /= 0) then
     write(6, "('hprod_fulldiag_poisson1: dpbtrs bad info. ', i5)") info
     stop
  end if

  ! solution with boundary condition
  do irad = 1, dim
     rho(irad) = (rho2(irad, 1) + tmpr * rpowl1(irad)) * fac1(irad) &
               + (rho2(irad, 2) + tmpi * rpowl1(irad)) * fac1(irad) * iunit
  end do

  deallocate(rpowl1)
  deallocate(rpowl)
  deallocate(invr)
  deallocate(fac1)
  deallocate(rho2)

end subroutine hprod_fulldiag_poisson1
!######################################################################
subroutine hprod_fulldiag_diag(fock)
!
  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit, iunit
  use mod_rad, only : nrad
  use mod_sph, only : lmax1

  implicit none
  ! ### args ###
  complex(c_double_complex), intent(in) :: fock(1:(nrad-1), 1:(nrad-1), 0:lmax1)

  ! ### local ###
  character(len = 9), parameter :: fname = "HFOrb.dat"
  integer(c_int), parameter :: iow = 10
  integer(c_int) :: dim, l, ifun, irad, maxfun
  complex(c_double_complex), allocatable :: ftmp(:,:), utmp(:,:)

  dim = nrad - 1
  maxfun = lmax1
  allocate(ftmp(1:dim, 1:dim))
  allocate(utmp(1:dim, 1:dim))

  open(unit = iow, file = trim(fname), status = 'unknown', form = 'formatted')
  write(iow, "(I25)") nrad
  write(iow, "(I25)") maxfun
  write(iow, "(I25)") lmax1

  do l = 0, lmax1
     utmp(1:dim, 1:dim) = czero
     ftmp(1:dim, 1:dim) = fock(1:dim, 1:dim, l)

     call lapack_zheev(dim, ftmp, utmp)
     write(6, "('Eigenvalues: l = ', i5)") l
     do ifun = 1, maxfun
        write(6, "(i5, f20.10)") ifun, dble(ftmp(ifun, ifun))
     end do

     do ifun = 1, maxfun
        do irad = 1, nrad - 1
           write(iow, "(2E25.15)") utmp(irad, ifun)
        end do
     end do
  end do

  close(unit = iow)

  deallocate(utmp)
  deallocate(ftmp)

end subroutine hprod_fulldiag_diag
!######################################################################
