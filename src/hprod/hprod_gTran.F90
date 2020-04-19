!######################################################################
subroutine hprod_gtran(L2V, lfield, max_irmax, wfn, twfn)

  use, intrinsic :: iso_c_binding
  use mod_rad, only : xrad
  use mod_sph, only : lmax1, nlat, cost
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : iunit
  use mod_bas, only : mval

  implicit none
  logical(c_bool), intent(in) :: L2V
  real(c_double), intent(in) :: lfield(1:3,1:3)
  integer(c_int), intent(in) :: max_irmax
  complex(c_double_complex), intent(in) :: wfn(1:max_irmax,0:lmax1,1:nfun)
  complex(c_double_complex), intent(out) :: twfn(1:max_irmax,0:lmax1,1:nfun)

  integer(c_int) :: ifun
  integer(c_int) :: l,i,irad,ilat
  complex(c_double_complex) :: faclv
  complex(c_double_complex), allocatable :: wang(:,:,:)

  allocate(wang(1:max_irmax,1:nlat,1:nfun))  

  !call bas_sph2ang1_dyn(wfn, wang)
  call hprod_gtran_sph2ang(max_irmax, wfn, wang)

  if (L2V) then
     !###########################
     !$omp parallel default(shared) private(faclv)
     !$omp do
     do ilat = 1, nlat
        do irad = 1, max_irmax
           faclv = exp(-iunit * lfield(3,3) * xrad(irad) * cost(ilat))
           do ifun = nfcore+1, nfun
              wang(irad, ilat, ifun) = wang(irad, ilat, ifun) * faclv
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     !###########################
  else
     !###########################
     !$omp parallel default(shared) private(faclv)
     !$omp do
     do ilat = 1, nlat
        do irad = 1, max_irmax
           faclv = exp(+iunit * lfield(3,3) * xrad(irad) * cost(ilat))
           do ifun = nfcore+1, nfun
              wang(irad, ilat, ifun) = wang(irad, ilat, ifun) * faclv
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     !###########################
  end if

  twfn = 0D0
  !call bas_ang2sph1_dyn(wang, twfn)
  call hprod_gtran_ang2sph(max_irmax, wang, twfn)

  deallocate(wang)

end subroutine hprod_gtran
!######################################################################
subroutine hprod_gtran_sph2ang(max_irmax, orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1, nlat, legf1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  integer(c_int), intent(in) :: max_irmax
  !complex(c_double_complex), intent(in) :: orb1(1:max_irmax, 0:lmax1, 1:*)
  complex(c_double_complex), intent(in) :: orb1(1:max_irmax, 0:lmax1, 1:*)
  complex(c_double_complex), intent(out) :: orb2(1:max_irmax, 1:nlat, 1:*)

  integer(c_int) :: ifun, irad, ilat, l, mi
  integer(c_int) :: il, nil, mapil(1:2, 1:nfun*nlat)

  nil = 0
  do ifun = nfcore + 1, nfun
     do ilat = 1, nlat
        nil = nil + 1
        mapil(1, nil) = ifun
        mapil(2, nil) = ilat
     end do
  end do

  !$omp parallel default(shared) private(mi,ifun,ilat)
  !$omp do
  do il = 1, nil
     ifun = mapil(1, il)
     ilat = mapil(2, il)
     mi = mval(ifun)
     orb2(1:max_irmax, ilat, ifun) = czero
     do l = abs(mi), lmax1
        do irad = 1, max_irmax
           orb2(irad, ilat, ifun) = orb2(irad, ilat, ifun) &
                                  + orb1(irad, l,    ifun) * legf1(l, ilat, mi)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine hprod_gtran_sph2ang
!######################################################################
subroutine hprod_gtran_ang2sph(max_irmax, orb1, orb2)

  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1, nlat, legb1
  use mod_bas, only : mval
  use mod_ormas, only : nfcore, nfun
  use mod_const, only : czero

  implicit none
  integer(c_int), intent(in) :: max_irmax
  complex(c_double_complex), intent(in) :: orb1(1:max_irmax, 1:nlat, 1:*)
  complex(c_double_complex), intent(out) :: orb2(1:max_irmax, 0:lmax1, 1:*)

  integer(c_int) :: ifun, irad, ilat, l, mi
  integer(c_int) :: il, nil, mapil(1:2, 1:nfun*(lmax1+1))

  nil = 0
  do ifun = nfcore + 1, nfun
     mi = mval(ifun)
     !do l = abs(mi), min(lmax1,pp_maxl)
     do l = abs(mi), lmax1
        nil = nil + 1
        mapil(1, nil) = ifun
        mapil(2, nil) = l
     end do
  end do

  !$omp parallel default(shared) private(mi,ifun,l)
  !$omp do
  do il = 1, nil
     ifun = mapil(1, il)
     l    = mapil(2, il)
     mi = mval(ifun)
     orb2(1:max_irmax, l, ifun) = czero
     do ilat = 1, nlat
        do irad = 1, max_irmax
           orb2(irad, l, ifun) = orb2(irad, l,    ifun) &
                               + orb1(irad, ilat, ifun) * legb1(ilat, l, mi)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

end subroutine hprod_gtran_ang2sph
!######################################################################
