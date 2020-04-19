!######################################################################
subroutine hprod_mfprod()

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_hprod, only : gorbg, v2orbg
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_rad, only : ecs_flag

  implicit none
  integer(c_int) :: iproc, lll, ull
  integer(c_int), external :: util_omp_iproc
  
  !$omp parallel default(shared) private(iproc, lll, ull)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nlat, lll, ull)
  if (.not. jfc_implicit) call hprod_mfprod_fcj(ctrue, gorbg, lll, ull)
  if (.not. xfc_implicit) call hprod_mfprod_fcx2(ctrue, gorbg, lll, ull)
  call hprod_mfprod_fcx1(ctrue, gorbg, lll, ull)
  call hprod_mfprod_dcj(ctrue, gorbg, lll, ull)
  call hprod_mfprod_dcx(ctrue, gorbg, lll, ull)
  call hprod_mfprod_actj(ctrue, v2orbg, lll, ull)
  call hprod_mfprod_actx(ctrue, v2orbg, lll, ull)
! Orimo_ECS
  if (ecs_flag == 0) then
     call hprod_mfprod_act2(v2orbg, lll, ull)
  else
     call hprod_mfprod_act2_ecs(v2orbg, lll, ull)
  end if
! Orimo_ECS
  !###########################
  !$omp end parallel

end subroutine hprod_mfprod
!######################################################################
subroutine hprod_mfprod_ene()

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_hprod, only : gorbg
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit

  implicit none
  integer(c_int) :: iproc, lll, ull
  integer(c_int), external :: util_omp_iproc
  
  !$omp parallel default(shared) private(iproc, lll, ull)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nlat, lll, ull)
!New_ene_FC
!  call hprod_mfprod_fcj(ctrue, gorbg, lll, ull)
!  call hprod_mfprod_fcx2(ctrue, gorbg, lll, ull)
!  call hprod_mfprod_fcx1(ctrue, gorbg, lll, ull)
!  call hprod_mfprod_dcj(ctrue, gorbg, lll, ull)
!  call hprod_mfprod_dcx(ctrue, gorbg, lll, ull)
  if (.not. jfc_implicit) call hprod_mfprod_fcj(cfalse, gorbg, lll, ull)
  if (.not. xfc_implicit) call hprod_mfprod_fcx2(cfalse, gorbg, lll, ull)
  call hprod_mfprod_fcx1(cfalse, gorbg, lll, ull)
  call hprod_mfprod_dcj(cfalse, gorbg, lll, ull)
  call hprod_mfprod_dcx(cfalse, gorbg, lll, ull)
!New_ene_FC
  !###########################
  !$omp end parallel

end subroutine hprod_mfprod_ene
!######################################################################
subroutine hprod_mfprod_gfock()

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_hprod, only : gorbg, v2orbg
  use mod_const, only : ctrue, cfalse

  implicit none
  integer(c_int) :: iproc, lll, ull
  integer(c_int), external :: util_omp_iproc

  stop 'hprod_mfprod_gfock is obsolete.'
  
  !$omp parallel default(shared) private(iproc, lll, ull)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nlat, lll, ull)
  call hprod_mfprod_fcj(ctrue, gorbg, lll, ull)
  call hprod_mfprod_fcx2(ctrue, gorbg, lll, ull)
  call hprod_mfprod_fcx1(ctrue, gorbg, lll, ull)
  call hprod_mfprod_dcj(ctrue, gorbg, lll, ull)
  call hprod_mfprod_dcx(ctrue, gorbg, lll, ull)
  call hprod_mfprod_actj(ctrue, v2orbg, lll, ull)
  call hprod_mfprod_actx(ctrue, v2orbg, lll, ull)
  call hprod_mfprod_act2(v2orbg, lll, ull)
  !###########################
  !$omp end parallel

end subroutine hprod_mfprod_gfock
!######################################################################
subroutine hprod_mfprod_gfock_nofcx()

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_hprod, only : gorbg, v2orbg
  use mod_const, only : ctrue, cfalse

  implicit none
  integer(c_int) :: iproc, lll, ull
  integer(c_int), external :: util_omp_iproc
  
  !$omp parallel default(shared) private(iproc, lll, ull)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nlat, lll, ull)
  call hprod_mfprod_fcj(cfalse, gorbg, lll, ull)
  call hprod_mfprod_fcx2(cfalse, gorbg, lll, ull)
  call hprod_mfprod_fcx1(cfalse, gorbg, lll, ull)
  call hprod_mfprod_dcj(ctrue, gorbg, lll, ull)
  call hprod_mfprod_dcx(cfalse, gorbg, lll, ull)
  call hprod_mfprod_actj(ctrue, v2orbg, lll, ull)
  call hprod_mfprod_actx(cfalse, v2orbg, lll, ull)
  call hprod_mfprod_act2(v2orbg, lll, ull)
!debug  call hprod_mfprod_fcj(cfalse, gorbg, lll, ull)
!debug  call hprod_mfprod_fcx2(cfalse, gorbg, lll, ull)
!debug  call hprod_mfprod_fcx1(cfalse, gorbg, lll, ull)
!debug  call hprod_mfprod_dcj(ctrue, gorbg, lll, ull)
!debug  call hprod_mfprod_dcx(ctrue, gorbg, lll, ull)
!debug  call hprod_mfprod_actj(ctrue, v2orbg, lll, ull)
!debug  call hprod_mfprod_actx(ctrue, v2orbg, lll, ull)
!debug  call hprod_mfprod_act2(v2orbg, lll, ull)
  !###########################
  !$omp end parallel

end subroutine hprod_mfprod_gfock_nofcx
!######################################################################
subroutine hprod_mfprod_mp()

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_hprod, only : gorbg
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit

  implicit none
  integer(c_int) :: iproc, lll, ull
  integer(c_int), external :: util_omp_iproc
  
  !$omp parallel default(shared) private(iproc, lll, ull)
  !###########################
  iproc = util_omp_iproc()
  call util_omp_disp(1, nlat, lll, ull)
  if (.not. jfc_implicit) call hprod_mfprod_fcj(ctrue, gorbg, lll, ull)
  if (.not. xfc_implicit) call hprod_mfprod_fcx2(ctrue, gorbg, lll, ull)
  call hprod_mfprod_fcx1(ctrue, gorbg, lll, ull)
  call hprod_mfprod_dcj(ctrue, gorbg, lll, ull)
  call hprod_mfprod_dcx(ctrue, gorbg, lll, ull)
  call hprod_mfprod_actjall(ctrue, gorbg, lll, ull)
  call hprod_mfprod_actxall(ctrue, gorbg, lll, ull)
  !###########################
  !$omp end parallel

end subroutine hprod_mfprod_mp
!######################################################################
subroutine hprod_mfprod_fcj(dofc, g1wfn, lll, ull)

  ! frozen-core-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : orbg, v2jfc

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, ilat, irad

  if (nfcore == 0) return

  ! acting on FC
  if (dofc) then
     do ifun = 1, nfcore
        do ilat = lll, ull
           do irad = 1, nradfc
              g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + v2jfc(irad) * orbg(irad, ilat, ifun)
           end do
        end do
     end do
  end if
  ! acting on DC and active
  do ifun = nfcore + 1, nfun
     do ilat = lll, ull
        do irad = 1, nrad - 1
           g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + v2jfc(irad) * orbg(irad, ilat, ifun)
        end do
     end do
  end do

end subroutine hprod_mfprod_fcj
!######################################################################
subroutine hprod_mfprod_fcx1(dofc, g1wfn, lll, ull)

  ! frozen-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore, nfun
  use mod_hprod, only : orbg, v2ang

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, ilat, irad

  if (nfcore1 == 0) return

  if (dofc) then
     do ifun = 1, nfcore
        do jfun = nfcore2 + 1, nfcore
           do ilat = lll, ull
              do irad = 1, nradfc
                 g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) - v2ang(irad, ilat, jfun, ifun) * orbg(irad, ilat, jfun)
              end do
           end do
        end do
     end do
  end if
  do ifun = nfcore + 1, nfun
     do jfun = nfcore2 + 1, nfcore
        do ilat = lll, ull
           do irad = 1, nradfc
              g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) - v2ang(irad, ilat, jfun, ifun) * orbg(irad, ilat, jfun)
           end do
        end do
     end do
  end do

end subroutine hprod_mfprod_fcx1
!######################################################################
subroutine hprod_mfprod_fcx2(dofc, g1wfn, lll, ull)

  ! frozen-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore, nfun
  use mod_hprod, only : orbg, v2ang, v2xfc

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, ilat, irad

  if (nfcore2 == 0) return

  if (dofc) then
     do ifun = 1, nfcore
        do ilat = lll, ull
           do irad = 1, nradfc
              g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + v2xfc(irad) * orbg(irad, ilat, ifun)
           end do
        end do
     end do
!exx     do ifun = 1, nfcore
!exx        do jfun = 1, nfcore2
!exx           do ilat = lll, ull
!exx              do irad = 1, nradfc
!exx                 g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) - v2ang(irad, ilat, jfun, ifun) * orbg(irad, ilat, jfun)
!exx              end do
!exx           end do
!exx        end do
!exx     end do
  end if
  do ifun = nfcore + 1, nfun
     do ilat = lll, ull
        do irad = 1, nradfc
           g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + v2xfc(irad) * orbg(irad, ilat, ifun)
        end do
     end do
  end do
!exx  do ifun = nfcore + 1, nfun
!exx     do jfun = 1, nfcore2
!exx        do ilat = lll, ull
!exx           do irad = 1, nradfc
!exx              g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) - v2ang(irad, ilat, jfun, ifun) * orbg(irad, ilat, jfun)
!exx           end do
!exx        end do
!exx     end do
!exx  end do

end subroutine hprod_mfprod_fcx2
!######################################################################
subroutine hprod_mfprod_dcj(dofc, g1wfn, lll, ull)

  ! dynamical-core-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_const, only : czero, ctwo
  use mod_ormas, only : nfcore, ndcore, ncore, nfun
  use mod_hprod, only : orbg, v2ang

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, ilat, irad
  complex(c_double_complex), allocatable :: d2v(:,:)

  if (ndcore == 0) return
  allocate(d2v(1:(nrad-1), lll:ull))

  ! generate potential
  d2v(1:(nrad-1), lll:ull) = czero
  do jfun = nfcore + 1, ncore
     do ilat = lll, ull
        do irad = 1, nrad - 1
           d2v(irad, ilat) = d2v(irad, ilat) + v2ang(irad, ilat, jfun, jfun) * ctwo
        end do
     end do
  end do
  ! acting on FC
  if (dofc) then
     do ifun = 1, nfcore
        do ilat = lll, ull
           do irad = 1, nradfc
              g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + d2v(irad, ilat) * orbg(irad, ilat, ifun)
           end do
        end do
     end do
  end if
  ! acting on DC and active
  do ifun = nfcore + 1, nfun
     do ilat = lll, ull
        do irad = 1, nrad - 1
           g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) + d2v(irad, ilat) * orbg(irad, ilat, ifun)
        end do
     end do
  end do

  deallocate(d2v)
  
end subroutine hprod_mfprod_dcj
!######################################################################
subroutine hprod_mfprod_dcx(dofc, g1wfn, lll, ull)

  ! dynamical-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_ormas, only : nfcore, ndcore, ncore, nfun
  use mod_rad, only : nrad, nradfc, xrad, wrad
  use mod_hprod, only : orbg, v2ang
  use mod_const, only : zero, two, ctwo, pi

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, ilat, irad

  if (ndcore == 0) return

  ! acting on FC
  if (dofc) then
     do ifun = 1, nfcore
        do jfun = nfcore + 1, ncore
           do ilat = lll, ull
              do irad = 1, nradfc
                 g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) - v2ang(irad, ilat, jfun, ifun) * orbg(irad, ilat, jfun)
              end do
           end do
        end do
     end do
  end if
  ! acting on DC and active
  do ifun = nfcore + 1, nfun
     do jfun = nfcore + 1, ncore
        do ilat = lll, ull
           do irad = 1, nrad - 1
              g1wfn(irad, ilat, ifun) = g1wfn(irad, ilat, ifun) - v2ang(irad, ilat, jfun, ifun) * orbg(irad, ilat, jfun)
           end do
        end do
     end do
  end do
  
end subroutine hprod_mfprod_dcx
!######################################################################
subroutine hprod_mfprod_actj(dofc, g2wfn, lll, ull)
  
  ! active-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_const, only : czero, chalf
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : orbg, v2ang, den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, kfun, jact, kact, ilat, irad
  complex(c_double_complex) :: dens
  complex(c_double_complex), allocatable :: d2v(:,:)
  
  if (ncore == 0 .or. nact == 0) return
  allocate(d2v(1:(nrad-1), lll:ull))
  
  ! generate potential
  d2v(1:(nrad-1), lll:ull) = czero
  do jact = 1, nact
     jfun = ncore + jact
     do kact = 1, nact
        kfun = ncore + kact
        dens = den1(kact, jact)
        do ilat = lll, ull
           do irad = 1, nrad - 1
              d2v(irad, ilat) = d2v(irad, ilat) + v2ang(irad, ilat, jfun, kfun) * dens
           end do
        end do
     end do
  end do
  ! acting on FC
  if (dofc) then
     do ifun = 1, nfcore
        do ilat = lll, ull
           do irad = 1, nradfc
              g2wfn(irad, ilat, ifun) = g2wfn(irad, ilat, ifun) + d2v(irad, ilat) * orbg(irad, ilat, ifun)
           end do
        end do
     end do
  end if
  ! acting on DC
  do ifun = nfcore + 1, ncore
     do ilat = lll, ull
        do irad = 1, nrad - 1
           g2wfn(irad, ilat, ifun) = g2wfn(irad, ilat, ifun) + d2v(irad, ilat) * orbg(irad, ilat, ifun)
        end do
     end do
  end do

  deallocate(d2v)
  
end subroutine hprod_mfprod_actj
!######################################################################
subroutine hprod_mfprod_actx(dofc, g2wfn, lll, ull)
  
  ! active-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_const, only : czero, chalf
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : orbg, v2ang, den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, kfun, jact, kact, ilat, irad
  complex(c_double_complex) :: dens
  complex(c_double_complex), allocatable :: d2v(:,:)
  
  if (ncore == 0 .or. nact == 0) return
  allocate(d2v(1:(nrad-1), lll:ull))
  
  do jact = 1, nact
     jfun = ncore + jact
     d2v(1:(nrad-1), lll:ull) = czero
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(jfun) == mval(kfun)) then
           dens = - den1(kact, jact) * chalf
           do ilat = lll, ull
              do irad = 1, nrad - 1
                 d2v(irad, ilat) = d2v(irad, ilat) + orbg(irad, ilat, kfun) * dens
              end do
           end do
        end if
     end do
     ! acting on FC
     if (dofc) then
        do ifun = 1, nfcore
           do ilat = lll, ull
              do irad = 1, nradfc
                 g2wfn(irad, ilat, ifun) = g2wfn(irad, ilat, ifun) + d2v(irad, ilat) * v2ang(irad, ilat, jfun, ifun)
              end do
           end do
        end do
     end if
     ! acting on DC
     do ifun = nfcore + 1, ncore
        do ilat = lll, ull
           do irad = 1, nrad - 1
              g2wfn(irad, ilat, ifun) = g2wfn(irad, ilat, ifun) + d2v(irad, ilat) * v2ang(irad, ilat, jfun, ifun)
           end do
        end do
     end do
  end do

  deallocate(d2v)
  
end subroutine hprod_mfprod_actx
!######################################################################
subroutine hprod_mfprod_actjall(dofc, g2wfn, lll, ull)
  
  ! active-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : orbg, v2ang, den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, kfun, jact, kact, ilat, irad
  complex(c_double_complex) :: dens
  complex(c_double_complex), allocatable :: d2v(:,:)
  
  if (nact == 0) return
  allocate(d2v(1:(nrad-1), lll:ull))
  
  ! generate potential
  d2v(1:(nrad-1), lll:ull) = czero
  do jact = 1, nact
     jfun = ncore + jact
     do kact = 1, nact
        kfun = ncore + kact
        dens = den1(kact, jact)
        do ilat = lll, ull
           do irad = 1, nrad - 1
              d2v(irad, ilat) = d2v(irad, ilat) + v2ang(irad, ilat, jfun, kfun) * dens
           end do
        end do
     end do
  end do
  do ifun = 1, nfun
     do ilat = lll, ull
        do irad = 1, nrad - 1
           g2wfn(irad, ilat, ifun) = g2wfn(irad, ilat, ifun) + d2v(irad, ilat) * orbg(irad, ilat, ifun)
        end do
     end do
  end do

  deallocate(d2v)
  
end subroutine hprod_mfprod_actjall
!######################################################################
subroutine hprod_mfprod_actxall(dofc, g2wfn, lll, ull)
  
  ! active-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : orbg, v2ang, den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, kfun, jact, kact, ilat, irad
  complex(c_double_complex) :: dens
  complex(c_double_complex), allocatable :: d2v(:,:)
  
  if (nact == 0) return
  allocate(d2v(1:(nrad-1), lll:ull))
  
  do jact = 1, nact
     jfun = ncore + jact
     d2v(1:(nrad-1), lll:ull) = czero
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(jfun) == mval(kfun)) then
           dens = - den1(kact, jact) * chalf
           do ilat = lll, ull
              do irad = 1, nrad - 1
                 d2v(irad, ilat) = d2v(irad, ilat) + orbg(irad, ilat, kfun) * dens
              end do
           end do
        end if
     end do
     do ifun = 1, nfun
        do ilat = lll, ull
           do irad = 1, nrad - 1
              g2wfn(irad, ilat, ifun) = g2wfn(irad, ilat, ifun) + d2v(irad, ilat) * v2ang(irad, ilat, jfun, ifun)
           end do
        end do
     end do
  end do

  deallocate(d2v)
  
end subroutine hprod_mfprod_actxall
!######################################################################
subroutine hprod_mfprod_act2(g2wfn, lll, ull)
  
  ! active 2RDM contributions

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_const, only : czero, thrwfn
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : orbg, v2ang, torbg, den2, rden

  implicit none
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, kfun, lfun, iact, jact, kact, lact, irad, ilat
  complex(c_double_complex) :: dens, tmp
  complex(c_double_complex), allocatable :: d2v(:,:)
  
  if (nact == 0) return
  allocate(d2v(1:(nrad-1), lll:ull))
  torbg(1:(nrad-1), lll:ull, 1:nfun) = czero
  
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, iact
        jfun = ncore + jact

        d2v(1:(nrad-1), lll:ull) = czero
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              if (mval(ifun) + mval(kfun) == mval(jfun) + mval(lfun)) then
                 dens = den2(jact, iact, lact, kact)
                 if (abs(dens) > thrwfn) then
                    do ilat = lll, ull
                       do irad = 1, nrad - 1
                          d2v(irad, ilat) = d2v(irad, ilat) + v2ang(irad, ilat, kfun, lfun) * dens
                       end do
                    end do
                 end if
              end if
           end do
        end do
  
        do ilat = lll, ull
           do irad = 1, nrad - 1
              torbg(irad, ilat, ifun) = torbg(irad, ilat, ifun) + orbg(irad, ilat, jfun) * d2v(irad, ilat)
           end do
        end do
        if (iact > jact) then
           do ilat = lll, ull
              do irad = 1, nrad - 1
                 torbg(irad, ilat, jfun) = torbg(irad, ilat, jfun) + orbg(irad, ilat, ifun) * conjg(d2v(irad, ilat))
              end do
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
           do ilat = lll, ull
              do irad = 1, nrad - 1
                 g2wfn(irad, ilat, kfun) = g2wfn(irad, ilat, kfun) + torbg(irad, ilat, ifun) * tmp
              end do
           end do
        end if
     end do
  end do

  deallocate(d2v)
  
end subroutine hprod_mfprod_act2
!######################################################################
subroutine hprod_mfprod_act2_ecs(g2wfn, lll, ull)
  
  ! active 2RDM contributions

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_const, only : czero, thrwfn
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : orbg, v2ang, torbg, den2, rden

  implicit none
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int), intent(in) :: lll, ull
  integer(c_int) :: ifun, jfun, kfun, lfun, iact, jact, kact, lact, irad, ilat
  complex(c_double_complex) :: dens, tmp
  complex(c_double_complex), allocatable :: d2v(:,:)
  
  if (nact == 0) return
  allocate(d2v(1:(nrad-1), lll:ull))
  torbg(1:(nrad-1), lll:ull, 1:nfun) = czero
  
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact !! <<<<<< changed from iact
        jfun = ncore + jact

        d2v(1:(nrad-1), lll:ull) = czero
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              if (mval(ifun) + mval(kfun) == mval(jfun) + mval(lfun)) then
                 dens = den2(jact, iact, lact, kact)
                 if (abs(dens) > thrwfn) then
                    do ilat = lll, ull
                       do irad = 1, nrad - 1
                          d2v(irad, ilat) = d2v(irad, ilat) + v2ang(irad, ilat, kfun, lfun) * dens
                       end do
                    end do
                 end if
              end if
           end do
        end do
  
        do ilat = lll, ull
           do irad = 1, nrad - 1
              torbg(irad, ilat, ifun) = torbg(irad, ilat, ifun) + orbg(irad, ilat, jfun) * d2v(irad, ilat)
           end do
        end do
!not hermitian        if (iact > jact) then
!not hermitian           do ilat = lll, ull
!not hermitian              do irad = 1, nrad - 1
!not hermitian                 torbg(irad, ilat, jfun) = torbg(irad, ilat, jfun) + orbg(irad, ilat, ifun) * conjg(d2v(irad, ilat))
!not hermitian              end do
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
           do ilat = lll, ull
              do irad = 1, nrad - 1
                 g2wfn(irad, ilat, kfun) = g2wfn(irad, ilat, kfun) + torbg(irad, ilat, ifun) * tmp
              end do
           end do
        end if
     end do
  end do

  deallocate(d2v)
  
end subroutine hprod_mfprod_act2_ecs
!######################################################################
