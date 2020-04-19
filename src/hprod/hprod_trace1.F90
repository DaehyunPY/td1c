!######################################################################
complex(c_double_complex) function hprod_trace1(dofc, den1, wfn, pwfn)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, nradfc, ecs_flag, irad_ecs
  use mod_const, only : zero, czero, ctwo
  use mod_ormas, only : nact, nfcore, ncore, nfun

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: den1(1:nact, 1:nact)
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: pwfn(1:(nrad-1), 0:lmax1, 1:nfun)

  complex(c_double_complex) :: ttmp, tval
  complex(c_double_complex), allocatable :: twfn(:,:)
  integer(c_int) :: llrf, ulrf, llrd, ulrd, ifun, jfun, iact, jact, l, irad
  integer(c_int) :: rdim

  allocate(twfn(1:(nrad-1), 0:lmax1))
! Orimo_ECS
  if (ecs_flag == 1) then
     rdim = irad_ecs - 1
  else 
     rdim = nrad - 1
  end if
! Orimo_ECS

  tval = czero
  !$omp parallel default(shared) private(ifun, jfun, llrf, ulrf, llrd, ulrd, ttmp) reduction(+:tval)
  !###########################
  call util_omp_disp(1, nradfc,   llrf, ulrf)
! Orimo_ECS
  !call util_omp_disp(1, nrad - 1, llrd, ulrd)
  call util_omp_disp(1, rdim, llrd, ulrd)
! Orimo_ECS

  ttmp = czero
  if (dofc) then
     do ifun = 1, nfcore
        do l = abs(mval(ifun)), lmax1
           do irad = llrf, ulrf
              ttmp = ttmp + conjg(wfn(irad, l, ifun)) * pwfn(irad, l, ifun)
           end do
        end do
     end do
  end if

  do ifun = nfcore + 1, ncore
     do l = abs(mval(ifun)), lmax1
        do irad = llrd, ulrd
           ttmp = ttmp + conjg(wfn(irad, l, ifun)) * pwfn(irad, l, ifun)
        end do
     end do
  end do

  ttmp = ttmp * ctwo

  do iact = 1, nact
     ifun = ncore + iact
     twfn(llrd:ulrd, 0:lmax1) = czero
     do jact = 1, nact
        jfun = ncore + jact
        if (mval(ifun) == mval(jfun)) then
           do l = abs(mval(ifun)), lmax1
              do irad = llrd, ulrd
                 twfn(irad, l) = twfn(irad, l) + wfn(irad, l, jfun) * den1(jact, iact)
              end do
           end do
        end if
     end do
     do l = 0, lmax1
        do irad = llrd, ulrd
           ttmp = ttmp + conjg(twfn(irad, l)) * pwfn(irad, l, ifun)
        end do
     end do
  end do
  tval = tval + ttmp
  !###########################
  !$omp end parallel

  hprod_trace1 = tval
  deallocate(twfn)
  return

end function hprod_trace1
!######################################################################
