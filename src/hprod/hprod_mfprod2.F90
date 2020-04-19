!######################################################################
subroutine hprod_mfprod2()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_hprod, only : orbg, v2ang, torbg, gorbg, v2orbg
  use mod_rad, only : ecs_flag

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(1:ngrid*max(1,nact*nact)))

  if (.not. jfc_implicit) call hprod_mfprod2_fcj(cfalse, orbg, gorbg)
  if (.not. xfc_implicit) call hprod_mfprod2_fcx2(cfalse, orbg, gorbg)
  call hprod_mfprod2_fcx1(cfalse, orbg, v2ang, gorbg)
  call hprod_mfprod2_dcj(cfalse, orbg, v2ang, gorbg, work)
  call hprod_mfprod2_dcx(cfalse, orbg, v2ang, gorbg)
  call hprod_mfprod2_actj(cfalse, orbg, v2ang, v2orbg, work)
  call hprod_mfprod2_actx(cfalse, orbg, v2ang, v2orbg, work)
! Orimo_ECS
  if (ecs_flag == 0) then
     call hprod_mfprod2_act2(orbg, v2ang, torbg, v2orbg, work)
  else
     call hprod_mfprod2_act2_ecs(orbg, v2ang, torbg, v2orbg, work)
  end if
! Orimo_ECS
!  call hprod_mfprod2_act2v2(orbg, v2ang, torbg, v2orbg, work)

  deallocate(work)

end subroutine hprod_mfprod2
!######################################################################
subroutine hprod_mfprod2_ene()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_hprod, only : orbg, v2ang, torbg, gorbg, v2orbg

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(1:ngrid))

  if (.not. jfc_implicit) call hprod_mfprod2_fcj(cfalse, orbg, gorbg)
  if (.not. xfc_implicit) call hprod_mfprod2_fcx2(cfalse, orbg, gorbg)
  call hprod_mfprod2_fcx1(cfalse, orbg, v2ang, gorbg)
  call hprod_mfprod2_dcj(cfalse, orbg, v2ang, gorbg, work)
  call hprod_mfprod2_dcx(cfalse, orbg, v2ang, gorbg)

  deallocate(work)

end subroutine hprod_mfprod2_ene
!######################################################################
subroutine hprod_mfprod2_gfock()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_hprod, only : orbg, v2ang, torbg, gorbg, v2orbg

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(ngrid*max(1,nact*nact)))

  if (.not. jfc_implicit) call hprod_mfprod2_fcj(ctrue, orbg, gorbg)
  if (.not. xfc_implicit) call hprod_mfprod2_fcx2(ctrue, orbg, gorbg)
  call hprod_mfprod2_fcx1(ctrue, orbg, v2ang, gorbg)
  call hprod_mfprod2_dcj(ctrue, orbg, v2ang, gorbg, work)
  call hprod_mfprod2_dcx(ctrue, orbg, v2ang, gorbg)
  call hprod_mfprod2_actj(ctrue, orbg, v2ang, v2orbg, work)
  call hprod_mfprod2_actx(ctrue, orbg, v2ang, v2orbg, work)
  call hprod_mfprod2_act2(orbg, v2ang, torbg, v2orbg, work)
!  call hprod_mfprod2_act2v2(orbg, v2ang, torbg, v2orbg, work)

  deallocate(work)

end subroutine hprod_mfprod2_gfock
!######################################################################
subroutine hprod_mfprod2_mp()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_hprod, only : gorbg
  use mod_hprod, only : orbg, v2ang, torbg, gorbg
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(ngrid*max(1,nact)))

  if (.not. jfc_implicit) call hprod_mfprod2_fcj(ctrue, orbg, gorbg)
  if (.not. xfc_implicit) call hprod_mfprod2_fcx2(ctrue, orbg, gorbg)
  call hprod_mfprod2_fcx1(ctrue, orbg, v2ang, gorbg)
  call hprod_mfprod2_dcj(ctrue, orbg, v2ang, gorbg, work)
  call hprod_mfprod2_dcx(ctrue, orbg, v2ang, gorbg)
  call hprod_mfprod2_actjall(ctrue, orbg, v2ang, gorbg, work)
  call hprod_mfprod2_actxall(ctrue, orbg, v2ang, gorbg, work)

  deallocate(work)

end subroutine hprod_mfprod2_mp
!######################################################################
subroutine hprod_mfprod2_fcj(dofc, wfn, g1wfn)

  ! frozen-core-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : v2jfc

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_long) :: ifun, ilat, irad, llr, ulr

  if (nfcore == 0) return

  ! acting on FC
  if (dofc) then
     !$omp parallel default(shared) private(llr,ulr)
     call util_omp_disp(1, nradfc, llr, ulr)
     do ifun = 1, nfcore
        do ilat = 1, nlat
           do irad = llr, ulr
              g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + v2jfc(irad) * wfn(irad, ilat, ifun)
           end do
        end do
     end do
     !$omp end parallel
  end if

  ! acting on DC and active
  !$omp parallel default(shared) private(llr,ulr)
  call util_omp_disp(1, nrad-1, llr, ulr)
  do ifun = nfcore + 1, nfun
     do ilat = 1, nlat
        do irad = llr, ulr
           g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + v2jfc(irad) * wfn(irad, ilat, ifun)
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod2_fcj
!######################################################################
subroutine hprod_mfprod2_fcx1(dofc, wfn, v2, g1wfn)

  ! frozen-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore, nfun

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_long) :: ifun, jfun, ilat, irad, llr, ulr

  if (nfcore1 == 0) return

  !$omp parallel default(shared) private(llr,ulr)
  call util_omp_disp(1, nradfc, llr, ulr)
  if (dofc) then
     do ifun = 1, nfcore
        do jfun = nfcore2 + 1, nfcore
           do ilat = 1, nlat
              do irad = llr, ulr
                 g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) - v2(irad, ilat, jfun, ifun) * wfn(irad, ilat, jfun)
              end do
           end do
        end do
     end do
  end if
  do ifun = nfcore + 1, nfun
     do jfun = nfcore2 + 1, nfcore
        do ilat = 1, nlat
           do irad = llr, ulr
              g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) - v2(irad, ilat, jfun, ifun) * wfn(irad, ilat, jfun)
           end do
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod2_fcx1
!######################################################################
subroutine hprod_mfprod2_fcx2(dofc, wfn, g1wfn)

  ! frozen-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore, nfun
  use mod_hprod, only : v2xfc

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_long) :: ifun, jfun, ilat, irad, llr, ulr

  if (nfcore2 == 0) return

  !$omp parallel default(shared) private(llr,ulr)
  call util_omp_disp(1, nradfc, llr, ulr)
  if (dofc) then
     do ifun = 1, nfcore
        do ilat = 1, nlat
           do irad = llr, ulr
              g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + v2xfc(irad) * wfn(irad, ilat, ifun)
           end do
        end do
     end do
  end if
  do ifun = nfcore + 1, nfun
     do ilat = 1, nlat
        do irad = llr, ulr
           g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + v2xfc(irad) * wfn(irad, ilat, ifun)
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod2_fcx2
!######################################################################
subroutine hprod_mfprod2_dcj(dofc, wfn, v2, g1wfn, d2v)

  ! dynamical-core-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_const, only : czero, ctwo
  use mod_ormas, only : nfcore, ndcore, ncore, nfun

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid)
  integer(c_long) :: ifun, jfun, igrid, llg, ulg

  if (ndcore == 0) return

  !$omp parallel default(shared) private(llg,ulg,ifun,jfun)
  call util_omp_disp(1, ngrid, llg, ulg)
  ! generate potential
  d2v(llg:ulg) = czero
  do jfun = nfcore + 1, ncore
     do igrid = llg, ulg
        d2v(igrid) = d2v(igrid) + v2(igrid, jfun, jfun) * ctwo
     end do
  end do

  ! acting on DC and active
  call util_omp_disp(1, ngrid, llg, ulg)
  do ifun = nfcore + 1, nfun
     do igrid = llg, ulg
        g1wfn(igrid, ifun) = g1wfn(igrid, ifun) + d2v(igrid) * wfn(igrid, ifun)
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) call hprod_mfprod2_dcj_onFC(wfn, d2v, g1wfn)

  contains
  !=======
  subroutine hprod_mfprod2_dcj_onFC(wfn, d2v, g1wfn)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: d2v(1:(nrad-1), 1:nlat)
    complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_long) :: irad, ilat, llr, ulr
    !$omp parallel default(shared) private(llr,ulr)
    call util_omp_disp(1, nradfc, llr, ulr)
    do ifun = 1, nfcore
       do ilat = 1, nlat
          do irad = llr, ulr
             g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + d2v(irad, ilat) * wfn(irad, ilat, ifun)
          end do
       end do
    end do
    !$omp end parallel
  end subroutine hprod_mfprod2_dcj_onFC
  !=======
end subroutine hprod_mfprod2_dcj
!######################################################################
subroutine hprod_mfprod2_dcx(dofc, wfn, v2, g1wfn)

  ! dynamical-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nfcore, ndcore, ncore, nfun
  use mod_const, only : zero, two, ctwo, pi

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfn(1:ngrid, 1:nfun)
  integer(c_long) :: ifun, jfun, igrid, llg, ulg

  if (ndcore == 0) return

  !$omp parallel default(shared) private(llg,ulg,ifun,jfun)
  call util_omp_disp(1, ngrid, llg, ulg)
  ! acting on DC and active
  do ifun = nfcore + 1, nfun
     do jfun = nfcore + 1, ncore
        do igrid = llg, ulg
           g1wfn(igrid, ifun) = g1wfn(igrid, ifun) - v2(igrid, jfun, ifun) * wfn(igrid, jfun)
        end do
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) call hprod_mfprod2_dcx_onFC(wfn, v2, g1wfn)

  contains
  !=======
  subroutine hprod_mfprod2_dcx_onFC(wfn, v2, g1wfn)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
    complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_long) :: irad, ilat, llr, ulr
    !$omp parallel default(shared) private(llr,ulr)
    call util_omp_disp(1, nradfc, llr, ulr)
    do ifun = 1, nfcore
       do jfun = nfcore + 1, ncore
          do ilat = 1, nlat
             do irad = llr, ulr
                g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) - v2(irad, ilat, jfun, ifun) * wfn(irad, ilat, jfun)
             end do
          end do
       end do
    end do
    !$omp end parallel
  end subroutine hprod_mfprod2_dcx_onFC
  !=======  
end subroutine hprod_mfprod2_dcx
!######################################################################
subroutine hprod_mfprod2_actj(dofc, wfn, v2, g2wfn, d2v)
  
  ! active-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_const, only : czero, chalf
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid)
  integer(c_long) :: ifun, jfun, kfun, jact, kact, igrid, llg, ulg
  complex(c_double_complex) :: dens
  
  if (ncore == 0 .or. nact == 0) return
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,dens)
  call util_omp_disp(1, ngrid, llg, ulg)
  ! generate potential
  d2v(llg:ulg) = czero
  do jact = 1, nact
     jfun = ncore + jact
     do kact = 1, nact
        kfun = ncore + kact
        dens = den1(kact, jact)
        do igrid = llg, ulg
           d2v(igrid) = d2v(igrid) + v2(igrid, jfun, kfun) * dens
        end do
     end do
  end do
  ! acting on DC
  do ifun = nfcore + 1, ncore
     do igrid = llg, ulg
        g2wfn(igrid, ifun) = g2wfn(igrid, ifun) + d2v(igrid) * wfn(igrid, ifun)
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) call hprod_mfprod2_actj_onFC(wfn, d2v, g2wfn)

  contains
  !=======
  subroutine hprod_mfprod2_actj_onFC(wfn, d2v, g2wfn)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: d2v(1:(nrad-1), 1:nlat)
    complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_long) :: irad, ilat, llr, ulr
    !$omp parallel default(shared) private(llr,ulr)
    call util_omp_disp(1, nradfc, llr, ulr)
    do ifun = 1, nfcore
       do ilat = 1, nlat
          do irad = llr, ulr
             g2wfn(irad, ilat, ifun) = g2wfn(irad, ilat, ifun) + d2v(irad, ilat) * wfn(irad, ilat, ifun)
          end do
       end do
    end do
    !$omp end parallel
  end subroutine hprod_mfprod2_actj_onFC
  !=======
end subroutine hprod_mfprod2_actj
!######################################################################
subroutine hprod_mfprod2_actx(dofc, wfn, v2, g2wfn, wd1)
  
  ! active-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, ngrid
  use mod_const, only : czero, chalf
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: wd1(1:ngrid, 1:nact)
  integer(c_long) :: ifun, jfun, kfun, jact, kact, igrid, llg, ulg
  complex(c_double_complex) :: dens
  
  if (ncore == 0 .or. nact == 0) return
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,dens)
  call util_omp_disp(1, ngrid, llg, ulg)
  do jact = 1, nact
     jfun = ncore + jact
     wd1(llg:ulg, jact) = czero
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(jfun) == mval(kfun)) then
           dens = - den1(kact, jact) * chalf
           do igrid = llg, ulg
              wd1(igrid, jact) = wd1(igrid, jact) + wfn(igrid, kfun) * dens
           end do
        end if
     end do
     ! acting on DC
     do ifun = nfcore + 1, ncore
        do igrid = llg, ulg
           g2wfn(igrid, ifun) = g2wfn(igrid, ifun) + wd1(igrid, jact) * v2(igrid, jfun, ifun)
        end do
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) call hprod_mfprod2_actx_onFC(wfn, v2, wd1, g2wfn)

  contains
  !=======
  subroutine hprod_mfprod2_actx_onFC(wfn, v2, wd1, g2wfn)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
    complex(c_double_complex), intent(in) :: wd1(1:(nrad-1), 1:nlat, 1:nact)
    complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_long) :: irad, ilat, llr, ulr
    !$omp parallel default(shared) private(llr,ulr,jfun)
    call util_omp_disp(1, nradfc, llr, ulr)
    do jact = 1, nact
       jfun = ncore + jact
       do ifun = 1, nfcore
          do ilat = 1, nlat
             do irad = llr, ulr
                g2wfn(irad, ilat, ifun) = &
                g2wfn(irad, ilat, ifun) + wd1(irad, ilat, jact) * v2(irad, ilat, jfun, ifun)
             end do
          end do
       end do
    end do
    !$omp end parallel
  end subroutine hprod_mfprod2_actx_onFC
  !=======
end subroutine hprod_mfprod2_actx
!######################################################################
subroutine hprod_mfprod2_actjall(dofc, wfn, v2, g2wfn, d2v)
  
  ! active-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_const, only : czero
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid)
  integer(c_long) :: ifun, jfun, kfun, jact, kact, igrid, llg, ulg
  complex(c_double_complex) :: dens
  
  if (nact == 0) return
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,dens)
  call util_omp_disp(1, ngrid, llg, ulg)
  ! generate potential
  d2v(llg:ulg) = czero
  do jact = 1, nact
     jfun = ncore + jact
     do kact = 1, nact
        kfun = ncore + kact
        dens = den1(kact, jact)
        do igrid = llg, ulg
           d2v(igrid) = d2v(igrid) + v2(igrid, jfun, kfun) * dens
        end do
     end do
  end do
  do ifun = 1, ncore
     do igrid = llg, ulg
        g2wfn(igrid, ifun) = g2wfn(igrid, ifun) + d2v(igrid) * wfn(igrid, ifun)
     end do
  end do
  if (nact >= 2) then
     do ifun = ncore + 1, nfun
        do igrid = llg, ulg
           g2wfn(igrid, ifun) = g2wfn(igrid, ifun) + d2v(igrid) * wfn(igrid, ifun)
        end do
     end do
  end if
  !$omp end parallel

end subroutine hprod_mfprod2_actjall
!######################################################################
subroutine hprod_mfprod2_actxall(dofc, wfn, v2, g2wfn, d2v)
  
  ! active-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, ngrid
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid)
  integer(c_long) :: ifun, jfun, kfun, jact, kact, igrid, llg, ulg
  complex(c_double_complex) :: dens
  
  if (nact == 0) return
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,dens)
  call util_omp_disp(1, ngrid, llg, ulg)
  do jact = 1, nact
     jfun = ncore + jact
     d2v(llg:ulg) = czero
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(jfun) == mval(kfun)) then
           dens = - den1(kact, jact) * chalf
           do igrid = llg, ulg
              d2v(igrid) = d2v(igrid) + wfn(igrid, kfun) * dens
           end do
        end if
     end do
     do ifun = 1, ncore
        do igrid = llg, ulg
           g2wfn(igrid, ifun) = g2wfn(igrid, ifun) + d2v(igrid) * v2(igrid, jfun, ifun)
        end do
     end do
     if (nact >= 2) then
        do ifun = ncore + 1, nfun
           do igrid = llg, ulg
              g2wfn(igrid, ifun) = g2wfn(igrid, ifun) + d2v(igrid) * v2(igrid, jfun, ifun)
           end do
        end do
     end if
  end do
  !$omp end parallel

end subroutine hprod_mfprod2_actxall
!######################################################################
subroutine hprod_mfprod2_act2(wfn, v2, twfn, g2wfn, d2v)
  
  ! active 2RDM contributions

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, ngrid
  use mod_const, only : czero, thrwfn
  use mod_ormas, only : ncore, nact, nfun, nelact
  use mod_hprod, only : den2, rden

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: twfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid)
  integer(c_long) :: ifun, jfun, kfun, lfun, iact, jact, kact, lact, igrid, llg, ulg
  complex(c_double_complex) :: dens, tmp
  
  if (nact == 0) return
  if (nelact(3) < 2) return
  twfn(1:ngrid, 1:nfun) = czero
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,lfun,dens,tmp)
  !###########################
  call util_omp_disp(1, ngrid, llg, ulg)
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, iact
!     do jact = 1, nact
        jfun = ncore + jact
        d2v(llg:ulg) = czero
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              if (mval(ifun) + mval(kfun) == mval(jfun) + mval(lfun)) then
                 dens = den2(jact, iact, lact, kact)
                 if (abs(dens) > thrwfn) then
                    call zaxpy(ulg-llg+1, dens, v2(llg,kfun,lfun), 1, d2v(llg), 1)
!                    do igrid = llg, ulg
!                       d2v(igrid) = d2v(igrid) + v2(igrid, kfun, lfun) * dens
!                    end do
                 end if
              end if
           end do
        end do
  
        do igrid = llg, ulg
           twfn(igrid, ifun) = twfn(igrid, ifun) + wfn(igrid, jfun) * d2v(igrid)
        end do
        if (iact > jact) then
           do igrid = llg, ulg
              twfn(igrid, jfun) = twfn(igrid, jfun) + wfn(igrid, ifun) * conjg(d2v(igrid))
!              twfn(igrid, jfun) = twfn(igrid, jfun) + wfn(igrid, ifun) * conjg(d2v(igrid)) * (-1)**(mval(ifun)-mval(jfun))
           end do
        end if
     end do
  end do

  do iact = 1, nact
     ifun = ncore + iact
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(ifun) == mval(kfun)) then
           tmp = rden(iact, kact)
           call zaxpy(ulg-llg+1, tmp, twfn(llg,ifun), 1, g2wfn(llg,kfun), 1)
!           do igrid = llg, ulg
!              g2wfn(igrid, kfun) = g2wfn(igrid, kfun) + twfn(igrid, ifun) * tmp
!           end do
        end if
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod2_act2
!######################################################################
! Orimo_ECS
subroutine hprod_mfprod2_act2_ecs(wfn, v2, twfn, g2wfn, d2v)
  
  ! active 2RDM contributions

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, ngrid
  use mod_const, only : czero, thrwfn
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : den2, rden

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: twfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid)
  integer(c_long) :: ifun, jfun, kfun, lfun, iact, jact, kact, lact, igrid, llg, ulg
  complex(c_double_complex) :: dens, tmp
  
  if (nact == 0) return
  twfn(1:ngrid, 1:nfun) = czero
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,lfun,dens,tmp)
  !###########################
  call util_omp_disp(1, ngrid, llg, ulg)
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact  !! <<<<<< changed from iact
        jfun = ncore + jact
        d2v(llg:ulg) = czero
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              if (mval(ifun) + mval(kfun) == mval(jfun) + mval(lfun)) then
                 dens = den2(jact, iact, lact, kact)
                 if (abs(dens) > thrwfn) then
                    call zaxpy(ulg-llg+1, dens, v2(llg,kfun,lfun), 1, d2v(llg), 1)
!                    do igrid = llg, ulg
!                       d2v(igrid) = d2v(igrid) + v2(igrid, kfun, lfun) * dens
!                    end do
                 end if
              end if
           end do
        end do
  
        do igrid = llg, ulg
           twfn(igrid, ifun) = twfn(igrid, ifun) + wfn(igrid, jfun) * d2v(igrid)
        end do
!not hermitian        if (iact > jact) then
!not hermitian           do igrid = llg, ulg
!not hermitian              twfn(igrid, jfun) = twfn(igrid, jfun) + wfn(igrid, ifun) * conjg(d2v(igrid))
!not hermitian           end do
!not hermitian        end if
     end do
  end do

  do iact = 1, nact
     ifun = ncore + iact
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(ifun) == mval(kfun)) then
           tmp = rden(iact, kact)
           call zaxpy(ulg-llg+1, tmp, twfn(llg,ifun), 1, g2wfn(llg,kfun), 1)
!           do igrid = llg, ulg
!              g2wfn(igrid, kfun) = g2wfn(igrid, kfun) + twfn(igrid, ifun) * tmp
!           end do
        end if
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod2_act2_ecs
! Orimo_ECS
!######################################################################
subroutine hprod_mfprod2_act2v2(wfn, v2, twfn, g2wfn, d2v)
  
  ! active 2RDM contributions

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, ngrid
  use mod_const, only : czero, thrwfn
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : den2, rden

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: twfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid, 1:nact, 1:nact)
  integer(c_long) :: ifun, jfun, kfun, lfun, iact, jact, kact, lact, igrid, llg, ulg
  complex(c_double_complex) :: dens, tmp
  integer(c_long) :: ij, nij, mapij(1:2, 1:nact*nact)
  complex(c_double_complex) :: denijkl(1:nact**4)
  
  if (nact == 0) return
  twfn(1:ngrid, 1:nfun) = czero

  nij = 0
  do iact = 1, nact
     do jact = 1, iact ! note here
        nij = nij + 1
        mapij(1, nij) = iact
        mapij(2, nij) = jact
     end do
  end do
  
  !$omp parallel default(shared) private(iact,jact,ifun,jfun,kfun,lfun,dens)
  !$omp do
  do ij = 1, nij
     iact = mapij(1, ij); ifun = ncore + iact
     jact = mapij(2, ij); jfun = ncore + jact
     d2v(1:ngrid, jact, iact) = czero
     do kact = 1, nact
        kfun = ncore + kact
        do lact = 1, nact
           lfun = ncore + lact
           if (mval(ifun) + mval(kfun) == mval(jfun) + mval(lfun)) then
              dens = den2(jact, iact, lact, kact)
              if (abs(dens) > thrwfn) then
                 call zaxpy(ngrid, dens, v2(1,kfun,lfun), 1, d2v(1,jact,iact), 1)
              end if
           end if
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,tmp)
  !###########################
  call util_omp_disp(1, ngrid, llg, ulg)
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, iact
        jfun = ncore + jact
        do igrid = llg, ulg
           twfn(igrid, ifun) = twfn(igrid, ifun) + wfn(igrid, jfun) * d2v(igrid,jact,iact)
        end do
        if (iact > jact) then
           do igrid = llg, ulg
              twfn(igrid, jfun) = twfn(igrid, jfun) + wfn(igrid, ifun) * conjg(d2v(igrid,jact,iact))
           end do
        end if
     end do
  end do

  do iact = 1, nact
     ifun = ncore + iact
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(ifun) == mval(kfun)) then
           tmp = rden(iact, kact)
           call zaxpy(ulg-llg+1, tmp, twfn(llg,ifun), 1, g2wfn(llg,kfun), 1)
!           do igrid = llg, ulg
!              g2wfn(igrid, kfun) = g2wfn(igrid, kfun) + twfn(igrid, ifun) * tmp
!           end do
        end if
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod2_act2v2
!######################################################################
