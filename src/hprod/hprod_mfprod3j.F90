!######################################################################
subroutine hprod_mfprod3j()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_hprod, only : orbg, v2ang, torbg, gorbg, v2orbg
  use mod_hprod, only : orbe,orbo,v2ange,v2ango,torb3j,gorbe,gorbo,v2orbe,v2orbo

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(1:2*ngrid*max(1,nact*nact)))

  if (.not. jfc_implicit) call hprod_mfprod3j_fcj(cfalse, orbe, orbo, gorbe, gorbo)
  if (.not. xfc_implicit) call hprod_mfprod3j_fcx2(cfalse, orbe, orbo, gorbe, gorbo)
  call hprod_mfprod3j_fcx1(cfalse, orbe, orbo, v2ange, v2ango, gorbe, gorbo)
  call hprod_mfprod3j_dcj(cfalse, orbe, orbo, v2ange, v2ango, gorbe, gorbo, work)
  call hprod_mfprod3j_dcx(cfalse, orbe, orbo, v2ange, v2ango, gorbe, gorbo)
  call hprod_mfprod3j_actj(cfalse, orbe,orbo,v2ange,v2ango,v2orbe,v2orbo, work)
  call hprod_mfprod3j_actx(cfalse, orbe,orbo,v2ange,v2ango,v2orbe,v2orbo, work)
  call hprod_mfprod3j_act2(orbe,orbo, v2ange,v2ango, torb3j, v2orbe,v2orbo, work)

  deallocate(work)

end subroutine hprod_mfprod3j
!######################################################################
subroutine hprod_mfprod3j_ene()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_hprod, only : orbe,orbo,v2ange,v2ango,torb3j,gorbe,gorbo,v2orbe,v2orbo

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(1:2*ngrid))

  if (.not. jfc_implicit) call hprod_mfprod3j_fcj(cfalse, orbe, orbo, gorbe, gorbo)
  if (.not. xfc_implicit) call hprod_mfprod3j_fcx2(cfalse, orbe, orbo, gorbe, gorbo)
  call hprod_mfprod3j_fcx1(cfalse, orbe, orbo, v2ange, v2ango, gorbe, gorbo)
  call hprod_mfprod3j_dcj(cfalse, orbe, orbo, v2ange, v2ango, gorbe, gorbo, work)
  call hprod_mfprod3j_dcx(cfalse, orbe, orbo, v2ange, v2ango, gorbe, gorbo)

  deallocate(work)

end subroutine hprod_mfprod3j_ene
!######################################################################
subroutine hprod_mfprod3j_gfock()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_hprod, only : orbe,orbo,v2ange,v2ango,torb3j,gorbe,gorbo,v2orbe,v2orbo

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(2*ngrid*max(1,nact*nact)))

  if (.not. jfc_implicit) call hprod_mfprod3j_fcj(ctrue, orbe, orbo, gorbe, gorbo)
  if (.not. xfc_implicit) call hprod_mfprod3j_fcx2(ctrue, orbe, orbo, gorbe, gorbo)
  call hprod_mfprod3j_fcx1(ctrue, orbe, orbo, v2ange, v2ango, gorbe, gorbo)
  call hprod_mfprod3j_dcj(ctrue, orbe, orbo, v2ange, v2ango, gorbe, gorbo, work)
  call hprod_mfprod3j_dcx(ctrue, orbe, orbo, v2ange, v2ango, gorbe, gorbo)
  call hprod_mfprod3j_actj(ctrue, orbe,orbo,v2ange,v2ango,v2orbe,v2orbo, work)
  call hprod_mfprod3j_actx(ctrue, orbe,orbo,v2ange,v2ango,v2orbe,v2orbo, work)
  call hprod_mfprod3j_act2(orbe,orbo, v2ange,v2ango, torb3j, v2orbe,v2orbo, work)

  deallocate(work)

end subroutine hprod_mfprod3j_gfock
!######################################################################
subroutine hprod_mfprod3j_mp()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_hprod, only : orbe,orbo,v2ange,v2ango,torb3j,gorbe,gorbo,v2orbe,v2orbo
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(2*ngrid*max(1,nact)))

  if (.not. jfc_implicit) call hprod_mfprod3j_fcj(ctrue, orbe, orbo, gorbe, gorbo)
  if (.not. xfc_implicit) call hprod_mfprod3j_fcx2(ctrue, orbe, orbo, gorbe, gorbo)
  call hprod_mfprod3j_fcx1(ctrue, orbe, orbo, v2ange, v2ango, gorbe, gorbo)
  call hprod_mfprod3j_dcj(ctrue, orbe, orbo, v2ange, v2ango, gorbe, gorbo, work)
  call hprod_mfprod3j_dcx(ctrue, orbe, orbo, v2ange, v2ango, gorbe, gorbo)
  call hprod_mfprod3j_actjall(ctrue, orbe, orbo, v2ange, v2ango, gorbe, gorbo, work)
  call hprod_mfprod3j_actxall(ctrue, orbe, orbo, v2ange, v2ango, gorbe, gorbo, work)

  deallocate(work)

end subroutine hprod_mfprod3j_mp
!######################################################################
subroutine hprod_mfprod3j_fcj(dofc, wfne, wfno, g1wfne, g1wfno)

  ! frozen-core-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : v2jfc

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfne(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfne(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfno(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_long) :: ifun, ilat, irad, llr, ulr

  if (nfcore == 0) return

  ! acting on FC
  if (dofc) then
     !$omp parallel default(shared) private(llr,ulr)
     call util_omp_disp(1, nradfc, llr, ulr)
     do ifun = 1, nfcore
        do ilat = 1, nlat
           do irad = llr, ulr
              g1wfne(irad, ilat, ifun) = g1wfne(irad, ilat, ifun) + v2jfc(irad) * wfne(irad, ilat, ifun)
              g1wfno(irad, ilat, ifun) = g1wfno(irad, ilat, ifun) + v2jfc(irad) * wfno(irad, ilat, ifun)
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
           g1wfne(irad, ilat, ifun) = g1wfne(irad, ilat, ifun) + v2jfc(irad) * wfne(irad, ilat, ifun)
           g1wfno(irad, ilat, ifun) = g1wfno(irad, ilat, ifun) + v2jfc(irad) * wfno(irad, ilat, ifun)
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod3j_fcj
!######################################################################
subroutine hprod_mfprod3j_fcx1(dofc, wfne,wfno, v2e,v2o, g1wfne,g1wfno)

  ! frozen-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore, nfun

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfne(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(in) :: v2e(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2o(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfne(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfno(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_long) :: ifun, jfun, ilat, irad, llr, ulr

  if (nfcore1 == 0) return

  !$omp parallel default(shared) private(llr,ulr)
  call util_omp_disp(1, nradfc, llr, ulr)
  if (dofc) then
     do ifun = 1, nfcore
        do jfun = nfcore2 + 1, nfcore
           do ilat = 1, nlat
              do irad = llr, ulr
                 g1wfne(irad, ilat, ifun) = g1wfne(irad, ilat, ifun) &
                      - v2e(irad, ilat, jfun, ifun) * wfne(irad, ilat, jfun) &
                      - v2o(irad, ilat, jfun, ifun) * wfno(irad, ilat, jfun)
                 g1wfno(irad, ilat, ifun) = g1wfno(irad, ilat, ifun) &
                      - v2e(irad, ilat, jfun, ifun) * wfno(irad, ilat, jfun) &
                      - v2o(irad, ilat, jfun, ifun) * wfne(irad, ilat, jfun)
              end do
           end do
        end do
     end do
  end if
  do ifun = nfcore + 1, nfun
     do jfun = nfcore2 + 1, nfcore
        do ilat = 1, nlat
           do irad = llr, ulr
              g1wfne(irad, ilat, ifun) = g1wfne(irad, ilat, ifun) &
                   - v2e(irad, ilat, jfun, ifun) * wfne(irad, ilat, jfun) &
                   - v2o(irad, ilat, jfun, ifun) * wfno(irad, ilat, jfun)
              g1wfno(irad, ilat, ifun) = g1wfno(irad, ilat, ifun) &
                   - v2e(irad, ilat, jfun, ifun) * wfno(irad, ilat, jfun) &
                   - v2o(irad, ilat, jfun, ifun) * wfne(irad, ilat, jfun)
           end do
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod3j_fcx1
!######################################################################
subroutine hprod_mfprod3j_fcx2(dofc, wfne,wfno, g1wfne,g1wfno)

  ! frozen-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore1, nfcore2, nfcore, nfun
  use mod_hprod, only : v2xfc

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfne(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfne(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfno(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_long) :: ifun, jfun, ilat, irad, llr, ulr

  if (nfcore2 == 0) return

  !$omp parallel default(shared) private(llr,ulr)
  call util_omp_disp(1, nradfc, llr, ulr)
  if (dofc) then
     do ifun = 1, nfcore
        do ilat = 1, nlat
           do irad = llr, ulr
              g1wfne(irad, ilat, ifun) = g1wfne(irad, ilat, ifun) + v2xfc(irad) * wfne(irad, ilat, ifun)
              g1wfno(irad, ilat, ifun) = g1wfno(irad, ilat, ifun) + v2xfc(irad) * wfno(irad, ilat, ifun)
           end do
        end do
     end do
  end if
  do ifun = nfcore + 1, nfun
     do ilat = 1, nlat
        do irad = llr, ulr
           g1wfne(irad, ilat, ifun) = g1wfne(irad, ilat, ifun) + v2xfc(irad) * wfne(irad, ilat, ifun)
           g1wfno(irad, ilat, ifun) = g1wfno(irad, ilat, ifun) + v2xfc(irad) * wfno(irad, ilat, ifun)
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod3j_fcx2
!######################################################################
subroutine hprod_mfprod3j_dcj(dofc, wfne,wfno, v2e,v2o, g1wfne,g1wfno, d2v)

  ! dynamical-core-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_const, only : czero, ctwo
  use mod_ormas, only : nfcore, ndcore, ncore, nfun

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2e(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2o(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid, 1:2)
  integer(c_long) :: ifun, jfun, igrid, llg, ulg

  if (ndcore == 0) return

  !$omp parallel default(shared) private(llg,ulg,ifun,jfun)
  call util_omp_disp(1, ngrid, llg, ulg)
  ! generate potential
  d2v(llg:ulg, 1:2) = czero
  do jfun = nfcore + 1, ncore
     do igrid = llg, ulg
        d2v(igrid,1) = d2v(igrid,1) + v2e(igrid, jfun, jfun) * ctwo
        d2v(igrid,2) = d2v(igrid,2) + v2o(igrid, jfun, jfun) * ctwo
     end do
  end do
  ! acting on DC and active
  do ifun = nfcore + 1, nfun
     do igrid = llg, ulg
        g1wfne(igrid, ifun) = g1wfne(igrid, ifun) &
             + d2v(igrid,1) * wfne(igrid, ifun) &
             + d2v(igrid,2) * wfno(igrid, ifun)
        g1wfno(igrid, ifun) = g1wfno(igrid, ifun) &
             + d2v(igrid,1) * wfno(igrid, ifun) &
             + d2v(igrid,2) * wfne(igrid, ifun)
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) call hprod_mfprod3j_dcj_onFC(wfne,wfno, d2v, g1wfne,g1wfno)

  contains
  !=======
  subroutine hprod_mfprod3j_dcj_onFC(wfne,wfno, d2v, g1wfne,g1wfno)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfne(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: wfno(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: d2v(1:(nrad-1), 1:nlat, 1:2)
    complex(c_double_complex), intent(inout) :: g1wfne(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(inout) :: g1wfno(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_long) :: irad, ilat, llr, ulr
    !$omp parallel default(shared) private(llr,ulr)
    call util_omp_disp(1, nradfc, llr, ulr)
    do ifun = 1, nfcore
       do ilat = 1, nlat
          do irad = llr, ulr
             g1wfne(irad,ilat,ifun) = g1wfne(irad,ilat,ifun) + d2v(irad,ilat,1) * wfne(irad,ilat,ifun)
             g1wfno(irad,ilat,ifun) = g1wfno(irad,ilat,ifun) + d2v(irad,ilat,2) * wfno(irad,ilat,ifun)
          end do
       end do
    end do
    !$omp end parallel
  end subroutine hprod_mfprod3j_dcj_onFC
  !=======
end subroutine hprod_mfprod3j_dcj
!######################################################################
subroutine hprod_mfprod3j_dcx(dofc, wfne,wfno, v2e,v2o, g1wfne,g1wfno)

  ! dynamical-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nfcore, ndcore, ncore, nfun
  use mod_const, only : zero, two, ctwo, pi

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2e(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2o(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfno(1:ngrid, 1:nfun)
  integer(c_long) :: ifun, jfun, igrid, llg, ulg

  if (ndcore == 0) return

  !$omp parallel default(shared) private(llg,ulg,ifun,jfun)
  call util_omp_disp(1, ngrid, llg, ulg)
  ! acting on DC and active
  do ifun = nfcore + 1, nfun
     do jfun = nfcore + 1, ncore
        do igrid = llg, ulg
           g1wfne(igrid, ifun) = g1wfne(igrid, ifun) &
                - v2e(igrid, jfun, ifun) * wfne(igrid, jfun) &
                - v2o(igrid, jfun, ifun) * wfno(igrid, jfun)
           g1wfno(igrid, ifun) = g1wfno(igrid, ifun) &
                - v2e(igrid, jfun, ifun) * wfno(igrid, jfun) &
                - v2o(igrid, jfun, ifun) * wfne(igrid, jfun)
        end do
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) call hprod_mfprod3j_dcx_onFC(wfne,wfno, v2e,v2o, g1wfne,g1wfno)

  contains
  !=======
  subroutine hprod_mfprod3j_dcx_onFC(wfne,wfno, v2e,v2o, g1wfne,g1wfno)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfne(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: wfno(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: v2e(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
    complex(c_double_complex), intent(in) :: v2o(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
    complex(c_double_complex), intent(inout) :: g1wfne(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(inout) :: g1wfno(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_long) :: irad, ilat, llr, ulr
    !$omp parallel default(shared) private(llr,ulr)
    call util_omp_disp(1, nradfc, llr, ulr)
    do ifun = 1, nfcore
       do jfun = nfcore + 1, ncore
          do ilat = 1, nlat
             do irad = llr, ulr
                g1wfne(irad, ilat, ifun) = g1wfne(irad, ilat, ifun) &
                     - v2e(irad, ilat, jfun, ifun) * wfne(irad, ilat, jfun) &
                     - v2o(irad, ilat, jfun, ifun) * wfno(irad, ilat, jfun)
                g1wfno(irad, ilat, ifun) = g1wfno(irad, ilat, ifun) &
                     - v2e(irad, ilat, jfun, ifun) * wfno(irad, ilat, jfun) &
                     - v2o(irad, ilat, jfun, ifun) * wfne(irad, ilat, jfun)
             end do
          end do
       end do
    end do
    !$omp end parallel
  end subroutine hprod_mfprod3j_dcx_onFC
  !=======  
end subroutine hprod_mfprod3j_dcx
!######################################################################
subroutine hprod_mfprod3j_actj(dofc, wfne,wfno, v2e,v2o, g2wfne,g2wfno, d2v)
  
  ! active-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_const, only : czero, chalf
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2e(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2o(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid, 1:2)
  integer(c_long) :: ifun, jfun, kfun, jact, kact, igrid, llg, ulg
  complex(c_double_complex) :: dens
  
  if (ncore == 0 .or. nact == 0) return
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,dens)
  call util_omp_disp(1, ngrid, llg, ulg)
  ! generate potential
  d2v(llg:ulg, 1:2) = czero
  do jact = 1, nact
     jfun = ncore + jact
     do kact = 1, nact
        kfun = ncore + kact
        dens = den1(kact, jact)
        do igrid = llg, ulg
           d2v(igrid,1) = d2v(igrid,1) + v2e(igrid,jfun,kfun) * dens
           d2v(igrid,2) = d2v(igrid,2) + v2o(igrid,jfun,kfun) * dens
        end do
     end do
  end do
  ! acting on DC
  do ifun = nfcore + 1, ncore
     do igrid = llg, ulg
        g2wfne(igrid,ifun) = g2wfne(igrid,ifun) &
             + d2v(igrid,1) * wfne(igrid,ifun) &
             + d2v(igrid,2) * wfno(igrid,ifun)
        g2wfno(igrid,ifun) = g2wfno(igrid,ifun) &
             + d2v(igrid,1) * wfno(igrid,ifun) &
             + d2v(igrid,2) * wfne(igrid,ifun)
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) call hprod_mfprod3j_actj_onFC(wfne,wfno, d2v, g2wfne,g2wfno)

  contains
  !=======
  subroutine hprod_mfprod3j_actj_onFC(wfne,wfno, d2v, g2wfne,g2wfno)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfne(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: wfno(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: d2v(1:(nrad-1), 1:nlat, 1:2)
    complex(c_double_complex), intent(inout) :: g2wfne(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(inout) :: g2wfno(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_long) :: irad, ilat, llr, ulr
    !$omp parallel default(shared) private(llr,ulr)
    call util_omp_disp(1, nradfc, llr, ulr)
    do ifun = 1, nfcore
       do ilat = 1, nlat
          do irad = llr, ulr
             g2wfne(irad,ilat,ifun) = g2wfne(irad,ilat,ifun) &
                  + d2v(irad,ilat,1) * wfne(irad,ilat,ifun) &
                  + d2v(irad,ilat,2) * wfno(irad,ilat,ifun)
             g2wfno(irad,ilat,ifun) = g2wfno(irad,ilat,ifun) &
                  + d2v(irad,ilat,1) * wfno(irad,ilat,ifun) &
                  + d2v(irad,ilat,2) * wfne(irad,ilat,ifun)
          end do
       end do
    end do
    !$omp end parallel
  end subroutine hprod_mfprod3j_actj_onFC
  !=======
end subroutine hprod_mfprod3j_actj
!######################################################################
subroutine hprod_mfprod3j_actx(dofc, wfne,wfno, v2e,v2o, g2wfne,g2wfno, wd1)
  
  ! active-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, ngrid
  use mod_const, only : czero, chalf
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2e(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2o(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: wd1(1:ngrid, 1:nact, 1:2)
  integer(c_long) :: ifun, jfun, kfun, jact, kact, igrid, llg, ulg
  complex(c_double_complex) :: dens
  
  if (ncore == 0 .or. nact == 0) return
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,dens)
  call util_omp_disp(1, ngrid, llg, ulg)
  do jact = 1, nact
     jfun = ncore + jact
     wd1(llg:ulg, jact, 1:2) = czero
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(jfun) == mval(kfun)) then
           dens = - den1(kact, jact) * chalf
           do igrid = llg, ulg
              wd1(igrid, jact,1) = wd1(igrid, jact,1) + wfne(igrid, kfun) * dens
              wd1(igrid, jact,2) = wd1(igrid, jact,2) + wfno(igrid, kfun) * dens
           end do
        end if
     end do
     ! acting on DC
     do ifun = nfcore + 1, ncore
        do igrid = llg, ulg
           g2wfne(igrid,ifun) = g2wfne(igrid,ifun) &
                + wd1(igrid,jact,1) * v2e(igrid,jfun,ifun) &
                + wd1(igrid,jact,2) * v2o(igrid,jfun,ifun)
           g2wfno(igrid,ifun) = g2wfno(igrid,ifun) &
                + wd1(igrid,jact,1) * v2o(igrid,jfun,ifun) &
                + wd1(igrid,jact,2) * v2e(igrid,jfun,ifun)
        end do
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) call hprod_mfprod3j_actx_onFC(wfne,wfno, v2e,v2o, wd1, g2wfne,g2wfno)

  contains
  !=======
  subroutine hprod_mfprod3j_actx_onFC(wfne,wfno, v2e,v2o, wd1, g2wfne,g2wfno)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfne(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: wfno(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: v2e(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
    complex(c_double_complex), intent(in) :: v2o(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
    complex(c_double_complex), intent(in) :: wd1(1:(nrad-1), 1:nlat, 1:nact, 1:2)
    complex(c_double_complex), intent(inout) :: g2wfne(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(inout) :: g2wfno(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_long) :: irad, ilat, llr, ulr
    !$omp parallel default(shared) private(llr,ulr,jfun)
    call util_omp_disp(1, nradfc, llr, ulr)
    do jact = 1, nact
       jfun = ncore + jact
       do ifun = 1, nfcore
          do ilat = 1, nlat
             do irad = llr, ulr
                g2wfne(irad,ilat,ifun) = g2wfne(irad,ilat,ifun) &
                     + wd1(irad,ilat,jact,1) * v2e(irad,ilat,jfun,ifun) &
                     + wd1(irad,ilat,jact,2) * v2o(irad,ilat,jfun,ifun)
                g2wfno(irad,ilat,ifun) = g2wfno(irad,ilat,ifun) &
                     + wd1(irad,ilat,jact,1) * v2o(irad,ilat,jfun,ifun) &
                     + wd1(irad,ilat,jact,2) * v2e(irad,ilat,jfun,ifun)
             end do
          end do
       end do
    end do
    !$omp end parallel
  end subroutine hprod_mfprod3j_actx_onFC
  !=======
end subroutine hprod_mfprod3j_actx
!######################################################################
subroutine hprod_mfprod3j_actjall(dofc, wfne,wfno, v2e,v2o, g2wfne,g2wfno, d2v)
  
  ! active-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2e(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2o(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid, 1:2)
  integer(c_long) :: ifun, jfun, kfun, jact, kact, igrid, llg, ulg
  complex(c_double_complex) :: dens
  
  if (nact == 0) return
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,dens)
  call util_omp_disp(1, ngrid, llg, ulg)
  ! generate potential
  d2v(llg:ulg,1:2) = czero
  do jact = 1, nact
     jfun = ncore + jact
     do kact = 1, nact
        kfun = ncore + kact
        dens = den1(kact, jact)
        do igrid = llg, ulg
           d2v(igrid,1) = d2v(igrid,1) + v2e(igrid, jfun, kfun) * dens
           d2v(igrid,2) = d2v(igrid,2) + v2o(igrid, jfun, kfun) * dens
        end do
     end do
  end do
  do ifun = 1, nfun
     do igrid = llg, ulg
        g2wfne(igrid, ifun) = g2wfne(igrid, ifun) &
             + d2v(igrid,1) * wfne(igrid, ifun) &
             + d2v(igrid,2) * wfno(igrid, ifun)
        g2wfno(igrid, ifun) = g2wfno(igrid, ifun) &
             + d2v(igrid,1) * wfno(igrid, ifun) &
             + d2v(igrid,2) * wfne(igrid, ifun)
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod3j_actjall
!######################################################################
subroutine hprod_mfprod3j_actxall(dofc, wfne,wfno, v2e,v2o, g2wfne,g2wfno, wd1)
  
  ! active-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, ngrid
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2e(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2o(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: wd1(1:ngrid, 1:2)
  integer(c_long) :: ifun, jfun, kfun, jact, kact, igrid, llg, ulg
  complex(c_double_complex) :: dens
  
  if (nact == 0) return
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,dens)
  call util_omp_disp(1, ngrid, llg, ulg)
  do jact = 1, nact
     jfun = ncore + jact
     wd1(llg:ulg,1:2) = czero
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(jfun) == mval(kfun)) then
           dens = - den1(kact, jact) * chalf
           do igrid = llg, ulg
              wd1(igrid,1) = wd1(igrid,1) + wfne(igrid, kfun) * dens
              wd1(igrid,2) = wd1(igrid,2) + wfno(igrid, kfun) * dens
           end do
        end if
     end do
     do ifun = 1, nfun
        do igrid = llg, ulg
           g2wfne(igrid, ifun) = g2wfne(igrid, ifun) &
                + wd1(igrid,1) * v2e(igrid, jfun, ifun) &
                + wd1(igrid,2) * v2o(igrid, jfun, ifun)
           g2wfno(igrid, ifun) = g2wfno(igrid, ifun) &
                + wd1(igrid,1) * v2o(igrid, jfun, ifun) &
                + wd1(igrid,2) * v2e(igrid, jfun, ifun)
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod3j_actxall
!######################################################################
subroutine hprod_mfprod3j_act2(wfne,wfno, v2e,v2o, twfn, g2wfne,g2wfno, d2v)
  
  ! active 2RDM contributions

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, ngrid
  use mod_const, only : czero, thrwfn
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : den2, rden

  implicit none
  complex(c_double_complex), intent(in) :: wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(in) :: v2e(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: v2o(1:ngrid, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: twfn(1:ngrid, 1:nfun, 1:2)
  complex(c_double_complex), intent(inout) :: g2wfne(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfno(1:ngrid, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:ngrid, 1:2)
  integer(c_long) :: ifun, jfun, kfun, lfun, iact, jact, kact, lact, igrid, llg, ulg
  complex(c_double_complex) :: dens, tmp
  
  if (nact == 0) return
  twfn(1:ngrid, 1:nfun, 1:2) = czero
  
  !$omp parallel default(shared) private(llg,ulg,ifun,jfun,kfun,lfun,dens,tmp)
  !###########################
  call util_omp_disp(1, ngrid, llg, ulg)
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, iact
        jfun = ncore + jact
        d2v(llg:ulg,1:2) = czero
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              if (mval(ifun) + mval(kfun) == mval(jfun) + mval(lfun)) then
                 dens = den2(jact, iact, lact, kact)
                 if (abs(dens) > thrwfn) then
                    call zaxpy(ulg-llg+1, dens, v2e(llg,kfun,lfun), 1, d2v(llg,1), 1)
                    call zaxpy(ulg-llg+1, dens, v2o(llg,kfun,lfun), 1, d2v(llg,2), 1)
                 end if
              end if
           end do
        end do
  
        do igrid = llg, ulg
           twfn(igrid,ifun,1) = twfn(igrid,ifun,1) + wfne(igrid,jfun) * d2v(igrid,1) &
                                                   + wfno(igrid,jfun) * d2v(igrid,2)
           twfn(igrid,ifun,2) = twfn(igrid,ifun,2) + wfne(igrid,jfun) * d2v(igrid,2) &
                                                   + wfno(igrid,jfun) * d2v(igrid,1)
        end do
        if (iact > jact) then
           do igrid = llg, ulg
              twfn(igrid,jfun,1) = twfn(igrid,jfun,1) + wfne(igrid,ifun) * conjg(d2v(igrid,1)) &
                                                      + wfno(igrid,ifun) * conjg(d2v(igrid,2))
              twfn(igrid,jfun,2) = twfn(igrid,jfun,2) + wfne(igrid,ifun) * conjg(d2v(igrid,2)) &
                                                      + wfno(igrid,ifun) * conjg(d2v(igrid,1))
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
           call zaxpy(ulg-llg+1, tmp, twfn(llg,ifun,1), 1, g2wfne(llg,kfun), 1)
           call zaxpy(ulg-llg+1, tmp, twfn(llg,ifun,2), 1, g2wfno(llg,kfun), 1)
        end if
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprod3j_act2
!######################################################################
