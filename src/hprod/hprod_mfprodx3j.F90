!######################################################################
subroutine hprod_mfprodx3j()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_hprod, only : orb, v2sph, torb, gorb, v2orb

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(1:ngrid*max(1,nact*nact)))

  if (.not. jfc_implicit) call hprod_mfprodx3j_fcj(cfalse, orb, gorb)
  if (.not. xfc_implicit) call hprod_mfprodx3j_fcx2(cfalse, orb, gorb)
  call hprod_mfprodx3j_fcx1(cfalse, orb, v2sph, gorb)
  call hprod_mfprodx3j_dcj(cfalse, orb, v2sph, gorb, work)
  call hprod_mfprodx3j_dcx(cfalse, orb, v2sph, gorb)
  call hprod_mfprodx3j_actj(cfalse, orb, v2sph, v2orb, work)
  call hprod_mfprodx3j_actx(cfalse, orb, v2sph, v2orb, work)
  call hprod_mfprodx3j_act2(orb, v2sph, torb, v2orb, work)
!  call hprod_mfprodx3j_act2v2(orb, v2sph, torb, v2orb, work)

  deallocate(work)

end subroutine hprod_mfprodx3j
!######################################################################
subroutine hprod_mfprodx3j_ene()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_hprod, only : orb, v2sph, gorb

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(1:ngrid))

  if (.not. jfc_implicit) call hprod_mfprodx3j_fcj(cfalse, orb, gorb)
  if (.not. xfc_implicit) call hprod_mfprodx3j_fcx2(cfalse, orb, gorb)
  call hprod_mfprodx3j_fcx1(cfalse, orb, v2sph, gorb)
  call hprod_mfprodx3j_dcj(cfalse, orb, v2sph, gorb, work)
  call hprod_mfprodx3j_dcx(cfalse, orb, v2sph, gorb)

  deallocate(work)

end subroutine hprod_mfprodx3j_ene
!######################################################################
subroutine hprod_mfprodx3j_gfock()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit
  use mod_hprod, only : orb, v2sph, torb, gorb, v2orb

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(ngrid*max(1,nact*nact)))

  if (.not. jfc_implicit) call hprod_mfprodx3j_fcj(ctrue, orb, gorb)
  if (.not. xfc_implicit) call hprod_mfprodx3j_fcx2(ctrue, orb, gorb)
  call hprod_mfprodx3j_fcx1(ctrue, orb, v2sph, gorb)
  call hprod_mfprodx3j_dcj(ctrue, orb, v2sph, gorb, work)
  call hprod_mfprodx3j_dcx(ctrue, orb, v2sph, gorb)
  call hprod_mfprodx3j_actj(ctrue, orb, v2sph, v2orb, work)
  call hprod_mfprodx3j_actx(ctrue, orb, v2sph, v2orb, work)
  call hprod_mfprodx3j_act2(orb, v2sph, torb, v2orb, work)
!  call hprod_mfprodx3j_act2v2(orb, v2sph, torb, v2orb, work)

  deallocate(work)

end subroutine hprod_mfprodx3j_gfock
!######################################################################
subroutine hprod_mfprodx3j_mp()

  use, intrinsic :: iso_c_binding
  use mod_bas, only : ngrid
  use mod_ormas, only : nact
  use mod_hprod, only : gorb
  use mod_hprod, only : orb, v2sph, gorb
  use mod_const, only : ctrue, cfalse
  use mod_control, only : jfc_implicit, xfc_implicit

  implicit none
  complex(c_double_complex), allocatable :: work(:)
  allocate(work(ngrid*max(1,nact)))

  if (.not. jfc_implicit) call hprod_mfprodx3j_fcj(ctrue, orb, gorb)
  if (.not. xfc_implicit) call hprod_mfprodx3j_fcx2(ctrue, orb, gorb)
  call hprod_mfprodx3j_fcx1(ctrue, orb, v2sph, gorb)
  call hprod_mfprodx3j_dcj(ctrue, orb, v2sph, gorb, work)
  call hprod_mfprodx3j_dcx(ctrue, orb, v2sph, gorb)
  call hprod_mfprodx3j_actjall(orb, v2sph, gorb, work)
  call hprod_mfprodx3j_actxall(orb, v2sph, gorb, work)

  deallocate(work)

end subroutine hprod_mfprodx3j_mp
!######################################################################
subroutine hprod_mfprodx3j_fcj(dofc, wfn, g1wfn)

  ! frozen-core-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, lval
  use mod_sph, only : lmax1
  use mod_rad, only : nrad, nradfc
  use mod_const, only : pi
  use mod_ormas, only : nfcore, nfun
  use mod_hprod, only : v2jfc

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  integer(c_int) :: ifun, l1, irad, llr, ulr

! real(c_double), parameter :: y00 = 1d0/sqrt(4d0*pi)
  real(c_double), parameter :: y00 = 1d0/(2d0*pi)

  if (nfcore == 0) return

  ! acting on FC
  if (dofc) then
     !$omp parallel default(shared) private(llr,ulr)
     call util_omp_disp(1, nradfc, llr, ulr)
     do ifun = 1, nfcore
        do irad = llr, ulr
!           g1wfn(irad,lval(ifun),ifun) = g1wfn(irad,lval(ifun),ifun) + v2jfc(irad) * wfn(irad,lval(ifun),ifun)
           g1wfn(irad,lval(ifun),ifun) = g1wfn(irad,lval(ifun),ifun) + v2jfc(irad) * wfn(irad,lval(ifun),ifun) * y00
        end do
     end do
     !$omp end parallel
  end if

  ! acting on DC and active
  !$omp parallel default(shared) private(llr,ulr)
  call util_omp_disp(1, nrad-1, llr, ulr)
  do ifun = nfcore + 1, nfun
     do l1 = abs(mval(ifun)), lmax1
        do irad = llr, ulr
!           g1wfn(irad,l1,ifun) = g1wfn(irad,l1,ifun) + v2jfc(irad) * wfn(irad,l1,ifun)
           g1wfn(irad,l1,ifun) = g1wfn(irad,l1,ifun) + v2jfc(irad) * wfn(irad,l1,ifun) * y00
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprodx3j_fcj
!######################################################################
subroutine hprod_mfprodx3j_fcx1(dofc, wfn, v2, g1wfn)

  ! frozen-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, nradfc
  use mod_bas, only : mval, lval
  use mod_sph, only : lmax1,lmax2, sph_gaunt
  use mod_ormas, only : nfcore1, nfcore2, nfcore, nfun

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  integer(c_int) :: ifun,jfun,l1,l12,irad,llr,ulr,mij,lll12,ull12

  if (nfcore1 == 0) return

  !$omp parallel default(shared) private(llr,ulr,mij,lll12,ull12)
  call util_omp_disp(1, nradfc, llr, ulr)
  if (dofc) then
     do ifun = 1, nfcore
        do jfun = nfcore2 + 1, nfcore
           mij = -mval(jfun)+mval(ifun)
           lll12 = max(abs(mij), abs(lval(ifun)-lval(jfun)))
           ull12 = min(lmax1, lval(ifun)+lval(jfun))
           do l12 = lll12, ull12
              do irad = llr, ulr
                 g1wfn(irad,lval(ifun),ifun) = g1wfn(irad,lval(ifun),ifun) &
                      - sph_gaunt(l12,lval(ifun),lval(jfun),mval(ifun),mval(jfun)) &
                      * wfn(irad,lval(jfun),jfun) * v2(irad,l12,jfun,ifun)
              end do
           end do
        end do
     end do
  end if
  do ifun = nfcore + 1, nfun
     do jfun = nfcore2 + 1, nfcore
        mij = -mval(jfun)+mval(ifun)
        do l1 = abs(mval(ifun)), lmax1
           lll12 = max(abs(mij), abs(l1-lval(jfun)))
           ull12 = min(lmax1, l1+lval(jfun))
           do l12 = lll12, ull12
              do irad = llr, ulr
                 g1wfn(irad,l1,ifun) = g1wfn(irad,l1,ifun) &
                      - sph_gaunt(l12,l1,lval(jfun),mval(ifun),mval(jfun)) &
                      * wfn(irad,lval(jfun),jfun) * v2(irad,l12,jfun,ifun)
              end do
           end do
        end do
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprodx3j_fcx1
!######################################################################
subroutine hprod_mfprodx3j_fcx2(dofc, wfn, g1wfn)

  ! frozen-core-eXchange

  use, intrinsic :: iso_c_binding
  use mod_sph, only : nlat
  use mod_rad, only : nrad, nradfc
  use mod_ormas, only : nfcore2, nfcore, nfun
  use mod_hprod, only : v2xfc

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
  complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
  integer(c_int) :: ifun, ilat, irad, llr, ulr

  if (nfcore2 == 0) return
  if (nfcore2 > 0) stop 'hprod_mfprodx3j_fcx2: nyi.'

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

end subroutine hprod_mfprodx3j_fcx2
!######################################################################
subroutine hprod_mfprodx3j_dcj(dofc, wfn, v2, g1wfn, d2v)

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
  integer(c_int) :: ifun, jfun, igrid, llg, ulg

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
  if (dofc) call hprod_mfprodx3j_dcj_onFC(wfn, d2v, g1wfn)

  contains
  !=======
  subroutine hprod_mfprodx3j_dcj_onFC(wfn, d2v, g1wfn)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: d2v(1:(nrad-1), 1:nlat)
    complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_int) :: irad, ilat, llr, ulr
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
  end subroutine hprod_mfprodx3j_dcj_onFC
  !=======
end subroutine hprod_mfprodx3j_dcj
!######################################################################
subroutine hprod_mfprodx3j_dcx(dofc, wfn, v2, g1wfn)

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
  integer(c_int) :: ifun, jfun, igrid, llg, ulg

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
  if (dofc) call hprod_mfprodx3j_dcx_onFC(wfn, v2, g1wfn)

  contains
  !=======
  subroutine hprod_mfprodx3j_dcx_onFC(wfn, v2, g1wfn)
    use, intrinsic :: iso_c_binding
    use mod_sph, only : nlat
    use mod_rad, only : nrad, nradfc
    implicit none
    complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 1:nlat, 1:nfun)
    complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 1:nlat, 1:nfun, 1:nfun)
    complex(c_double_complex), intent(inout) :: g1wfn(1:(nrad-1), 1:nlat, 1:nfun)
    integer(c_int) :: irad, ilat, llr, ulr
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
  end subroutine hprod_mfprodx3j_dcx_onFC
  !=======  
end subroutine hprod_mfprodx3j_dcx
!######################################################################
subroutine hprod_mfprodx3j_actj(dofc, wfn, v2, g2wfn, d2v)
  
  ! active-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, lval
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax1,lmax2,sph_gaunt
  use mod_const, only : czero, chalf
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:(nrad-1), 0:lmax2)
  integer(c_int) :: ifun,jfun,iact,jact,irad,llr,ulr,l1,l2,l12,lll1,ull1
  complex(c_double_complex) :: d2vo(1:(nrad-1))
  
  if (ncore == 0 .or. nact == 0) return
  d2v = 0d0
  
  !$omp parallel default(shared) private(llr,ulr,ifun,jfun,lll1,ull1)
  call util_omp_disp(1, nrad-1, llr, ulr)
  ! generate potential
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        if (mval(ifun).ne.mval(jfun)) cycle
        do l12 = 0, lmax2
           do irad = llr, ulr
              d2v(irad,l12) = d2v(irad,l12) + v2(irad,l12,jfun,ifun) * den1(iact,jact)
           end do
        end do
     end do
  end do

  ! acting on DC
  do ifun = nfcore+1, ncore
     do l12 = 0, lmax2
        do l2 = abs(mval(ifun)), lmax1
           do irad = llr, ulr
              d2vo(irad) = wfn(irad,l2,ifun) * d2v(irad,l12)
           end do
           lll1 = max(abs(mval(ifun)),abs(l12-l2))
           ull1 = min(lmax1,l12+l2)
           do l1 = lll1, ull1
              do irad = llr, ulr
                 g2wfn(irad,l1,ifun) = g2wfn(irad,l1,ifun) + sph_gaunt(l12,l1,l2,mval(ifun),mval(ifun)) * d2vo(irad)
              end do
           end do
        end do
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) then
     !$omp parallel default(shared) private(llr,ulr,ifun,jfun,lll1,ull1)
     call util_omp_disp(1, nradfc, llr, ulr)
     do ifun = 1, nfcore
        do l12 = 0, lmax2
           do irad = llr, ulr
              d2vo(irad) = wfn(irad,lval(ifun),ifun) * d2v(irad,l12)
           end do
           lll1 = max(abs(mval(ifun)),abs(l12-lval(ifun)))
           ull1 = min(lmax1,l12+lval(ifun))
           do l1 = lll1, ull1
              do irad = llr, ulr
                 g2wfn(irad,l1,ifun) = g2wfn(irad,l1,ifun) + sph_gaunt(l12,l1,lval(ifun),mval(ifun),mval(ifun)) * d2vo(irad)
              end do
           end do
        end do
     end do
     !$omp end parallel
  end if

end subroutine hprod_mfprodx3j_actj
!######################################################################
subroutine hprod_mfprodx3j_actx(dofc, wfn, v2, g2wfn, wd1)
  
  ! active-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval, lval
  use mod_rad, only : nrad, nradfc
  use mod_sph, only : lmax1,lmax2,sph_gaunt
  use mod_const, only : czero, chalf
  use mod_ormas, only : nfcore, ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  logical(c_bool), intent(in) :: dofc
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: wd1(1:(nrad-1), 0:lmax1, 1:nact)
  integer(c_int) :: ifun,jfun,iact,jact,irad,llr,ulr,l1,l2,l12,lll1,ull1,mij
  complex(c_double_complex) :: wd1v(1:(nrad-1))
  complex(c_double_complex) :: dens
  
  if (ncore == 0 .or. nact == 0) return
  wd1 = czero
  
  !$omp parallel default(shared) private(llr,ulr,ifun,jfun,dens,mij,lll1,ull1)
  call util_omp_disp(1, nrad-1, llr, ulr)
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        if (mval(ifun)==mval(jfun)) then
           dens = - den1(jact, iact) * chalf
           do l1 = abs(mval(ifun)), lmax1
              do irad = llr, ulr
                 wd1(irad,l1,iact) = wd1(irad,l1,iact) + wfn(irad,l1,jfun) * dens
              end do
           end do
        end if
     end do
  end do

  ! acting on DC
  do ifun = nfcore + 1, ncore
     do jact = 1, nact
        jfun = ncore + jact
        mij = mval(ifun)-mval(jfun)
        do l12 = abs(mij), lmax2
           do l2 = abs(mval(jfun)), lmax1
              do irad = llr, ulr
                 wd1v(irad) = wd1(irad,l2,jact)*v2(irad,l12,jfun,ifun)
              end do
              lll1 = max(abs(mval(ifun)),abs(l12-l2))
              ull1 = min(lmax1,l12+l2)
              do l1 = lll1, ull1
                 do irad = llr, ulr
                    g2wfn(irad,l1,ifun) = g2wfn(irad,l1,ifun) &
                         + sph_gaunt(l12,l1,l2,mval(ifun),mval(jfun)) * wd1v(irad)
                 end do
              end do
           end do
        end do
     end do
  end do
  !$omp end parallel

  ! acting on FC
  if (dofc) then
     !$omp parallel default(shared) private(llr,ulr,ifun,jfun,dens,mij,lll1,ull1)
     call util_omp_disp(1, nradfc, llr, ulr)
     do ifun = 1, nfcore
        do jact = 1, nact
           jfun = ncore + jact
           mij = mval(ifun)-mval(jfun)
           do l12 = abs(mij), min(lval(ifun)+lmax1, lmax2)
              do l2 = abs(mval(jfun)), lmax1
                 do irad = llr, ulr
                    wd1v(irad) = wd1(irad,l2,jact)*v2(irad,l12,jfun,ifun)
                 end do
                 lll1 = max(abs(mval(ifun)),abs(l12-l2))
                 ull1 = min(lmax1,l12+l2)
                 do l1 = lll1, ull1
                    do irad = llr, ulr
                       g2wfn(irad,l1,ifun) = g2wfn(irad,l1,ifun) &
                            + sph_gaunt(l12,l1,l2,mval(ifun),mval(jfun)) * wd1v(irad)
                    end do
                 end do
              end do
           end do
        end do
     end do
     !$omp end parallel
  end if

end subroutine hprod_mfprodx3j_actx
!######################################################################
subroutine hprod_mfprodx3j_actjall(wfn, v2, g2wfn, d2v)
  
  ! active-Coulomb

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1,lmax2,sph_gaunt
  use mod_bas, only : mval
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:(nrad-1), 0:lmax2)
  integer(c_int) :: ifun,jfun,iact,jact,irad,llr,ulr,l1,l2,l12,lll1,ull1
  complex(c_double_complex) :: d2vo(1:(nrad-1))
  
  if (nact == 0) return
  d2v = 0d0
  
  !$omp parallel default(shared) private(llr,ulr,ifun,jfun,lll1,ull1)
  call util_omp_disp(1, nrad-1, llr, ulr)
  ! generate potential
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        if (mval(ifun).ne.mval(jfun)) cycle
        do l12 = 0, lmax2
           do irad = llr, ulr
              d2v(irad,l12) = d2v(irad,l12) + v2(irad,l12,jfun,ifun) * den1(iact,jact)
           end do
        end do
     end do
  end do

  do ifun = 1, ncore
     do l12 = 0, lmax2
        do l2 = abs(mval(ifun)), lmax1
           do irad = llr, ulr
              d2vo(irad) = wfn(irad,l2,ifun) * d2v(irad,l12)
           end do
           lll1 = max(abs(mval(ifun)),abs(l12-l2))
           ull1 = min(lmax1,l12+l2)
           do l1 = lll1, ull1
              do irad = llr, ulr
                 g2wfn(irad,l1,ifun) = g2wfn(irad,l1,ifun) + sph_gaunt(l12,l1,l2,mval(ifun),mval(ifun)) * d2vo(irad)
              end do
           end do
        end do
     end do
  end do
  if (nact >= 2) then
     do ifun = ncore + 1, nfun
        do l12 = 0, lmax2
           do l2 = abs(mval(ifun)), lmax1
              do irad = llr, ulr
                 d2vo(irad) = wfn(irad,l2,ifun) * d2v(irad,l12)
              end do
              lll1 = max(abs(mval(ifun)),abs(l12-l2))
              ull1 = min(lmax1,l12+l2)
              do l1 = lll1, ull1
                 do irad = llr, ulr
                    g2wfn(irad,l1,ifun) = g2wfn(irad,l1,ifun) + sph_gaunt(l12,l1,l2,mval(ifun),mval(ifun)) * d2vo(irad)
                 end do
              end do
           end do
        end do
     end do
  end if
  !$omp end parallel

end subroutine hprod_mfprodx3j_actjall
!######################################################################
subroutine hprod_mfprodx3j_actxall(wfn, v2, g2wfn, wd1)
  
  ! active-eXchange

  use, intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_rad, only : nrad
  use mod_sph, only : lmax1,lmax2,sph_gaunt
  use mod_const, only : czero, chalf
  use mod_ormas, only : ncore, nact, nfun
  use mod_hprod, only : den1

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: wd1(1:(nrad-1), 0:lmax1, 1:nact)
  integer(c_int) :: ifun,jfun,iact,jact,irad,llr,ulr,l1,l2,l12,lll1,ull1,mij
  complex(c_double_complex) :: wd1v(1:(nrad-1))
  complex(c_double_complex) :: dens
  
  if (nact == 0) return
  wd1 = czero

  !$omp parallel default(shared) private(llr,ulr,ifun,jfun,dens,mij,lll1,ull1)
  call util_omp_disp(1, nrad-1, llr, ulr)
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        if (mval(ifun)==mval(jfun)) then
           dens = - den1(jact, iact) * chalf
           do l1 = abs(mval(ifun)), lmax1
              do irad = llr, ulr
                 wd1(irad,l1,iact) = wd1(irad,l1,iact) + wfn(irad,l1,jfun) * dens
              end do
           end do
        end if
     end do
  end do

  ! acting on all orbitals
  do ifun = 1, nfun
     if (nact <= 1 .and. ifun > ncore) cycle
     do jact = 1, nact
        jfun = ncore + jact
        mij = mval(ifun)-mval(jfun)
        do l12 = abs(mij), lmax2
           do l2 = abs(mval(jfun)), lmax1
              do irad = llr, ulr
                 wd1v(irad) = wd1(irad,l2,jact)*v2(irad,l12,jfun,ifun)
              end do
              lll1 = max(abs(mval(ifun)),abs(l12-l2))
              ull1 = min(lmax1,l12+l2)
              do l1 = lll1, ull1
                 do irad = llr, ulr
                    g2wfn(irad,l1,ifun) = g2wfn(irad,l1,ifun) &
                         + sph_gaunt(l12,l1,l2,mval(ifun),mval(jfun)) * wd1v(irad)
                 end do
              end do
           end do
        end do
     end do
  end do
  !$omp end parallel
  
end subroutine hprod_mfprodx3j_actxall
!######################################################################
subroutine hprod_mfprodx3j_act2(wfn, v2, twfn, g2wfn, d2v)
  
  ! active 2RDM contributions

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad
  use mod_sph, only : lmax1,lmax2,sph_gaunt
  use mod_bas, only : mval
  use mod_const, only : czero, thrwfn, pi
  use mod_ormas, only : ncore, nact, nfun, nelact
  use mod_hprod, only : den2, rden

  implicit none
  complex(c_double_complex), intent(in) :: wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(in) :: v2(1:(nrad-1), 0:lmax2, 1:nfun, 1:nfun)
  complex(c_double_complex), intent(out) :: twfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(inout) :: g2wfn(1:(nrad-1), 0:lmax1, 1:nfun)
  complex(c_double_complex), intent(out) :: d2v(1:(nrad-1), 0:lmax2, 1:nact, 1:nact)
  integer(c_int) :: ifun,jfun,kfun,lfun,iact,jact,kact,lact,irad,l1,l2,l12,mij,mfac,llr,ulr,lll1,ull1
  complex(c_double_complex) :: d2vo(1:(nrad-1))
  
  if (nact == 0) return
  if (nelact(3) < 2) return
  d2v = 0d0
  twfn = 0d0
  
  !$omp parallel default(shared) private(llr,ulr,ifun,jfun,kfun,lfun,mij,mfac,lll1,ull1)
  !###########################
  call util_omp_disp(1, nrad-1, llr, ulr)
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, iact
!     do jact = 1, nact
        jfun = ncore + jact
        mij = mval(ifun)-mval(jfun)
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
              if (mval(ifun)+mval(kfun) == mval(jfun)+mval(lfun) .and. &
                   abs(den2(jact,iact,lact,kact)) > thrwfn) then
                 do l12 = abs(mij), lmax2
                    do irad = llr, ulr
                       d2v(irad,l12,jact,iact) = d2v(irad,l12,jact,iact) + v2(irad,l12,kfun,lfun)*den2(jact,iact,lact,kact)
                    end do
                 end do
              end if
           end do
        end do
     end do
  end do

  do iact = 1, nact
     ifun = ncore + iact
     do jact = iact + 1, nact
        jfun = ncore + jact
        mij = mval(ifun)-mval(jfun)
        mfac = (-1)**mij
        do l12 = abs(mij), lmax2
           do irad = llr, ulr
!              d2v(irad,l12,jact,iact) = conjg(d2v(irad,l12,iact,jact))
              d2v(irad,l12,jact,iact) = conjg(d2v(irad,l12,iact,jact))*mfac
           end do
        end do
     end do
  end do

!NOTE HERE. SEE ALSO hprod_mkrho2_x3j and bas_gen_d2fac
!  do iact = 1, nact
!     do jact = 1, nact
!        mij = mval(ifun)-mval(jfun)
!        do l12 = abs(mij), lmax2
!           do irad = llr, ulr
!              d2v(irad,l12,jact,iact) = d2v(irad,l12,jact,iact)*sqrt(2d0*pi)
!           end do
!        end do
!     end do
!  end do
!NOTE HERE. SEE ALSO hprod_mkrho2_x3j and bas_gen_d2fac

  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        mij = mval(ifun)-mval(jfun)        
        do l12 = abs(mij), lmax2
           do l2 = abs(mval(jfun)), lmax1
              do irad = llr, ulr
                 d2vo(irad) = wfn(irad,l2,jfun)*d2v(irad,l12,jact,iact)
              end do
              lll1 = max(abs(mval(ifun)),abs(l12-l2))
              ull1 = min(lmax1,l12+l2)
              do l1 = lll1, ull1
                 do irad = llr, ulr
                    twfn(irad,l1,ifun) = twfn(irad,l1,ifun) &
                         + sph_gaunt(l12,l1,l2,mval(ifun),mval(jfun)) * d2vo(irad)
                 end do
              end do
           end do
        end do
     end do
  end do

  do iact = 1, nact
     ifun = ncore + iact
     do kact = 1, nact
        kfun = ncore + kact
        if (mval(ifun) == mval(kfun)) then
           do l1 = abs(mval(ifun)), lmax1
              do irad = llr, ulr
                 g2wfn(irad,l1,ifun) = g2wfn(irad,l1,ifun) + twfn(irad,l1,kfun) * rden(kact,iact)
              end do
           end do
        end if
     end do
  end do
  !$omp end parallel

end subroutine hprod_mfprodx3j_act2
!######################################################################
