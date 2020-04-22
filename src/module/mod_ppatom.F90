!###############
module mod_ppatom
  use iso_c_binding
  use fgsl
  implicit none

  public ppatom_init
  public ppatom_final
  public ppatom_getvs
  public ppatom_getvp
  public ppatom_getvloc
  public ppatom_getsfun
  public ppatom_getpfun
  public ppatom_getd1vloc
  public ppatom_seig
  public ppatom_peig
  public ppatom_rcut
  public ppatom_rmin
  public ppatom_rmax

  private
  integer(c_int) :: ncenter
  integer(c_int), allocatable :: zeff(:)
  real(c_double), allocatable :: rmin(:),rmax(:),rcut(:),seig(:),peig(:)
  type(fgsl_interp_accel), allocatable :: vs_acc(:),vp_acc(:),vloc_acc(:),sfun_acc(:),pfun_acc(:)
  type(fgsl_spline), allocatable :: vs_spline(:),vp_spline(:),vloc_spline(:),sfun_spline(:),pfun_spline(:)

  contains
  !##########
  subroutine ppatom_init(natom, ppfile)
    !-----
    integer(c_int), intent(in) :: natom
    character(len=256), intent(in) :: ppfile(1:natom)
    !-----
    integer(c_int) :: A

#ifdef FGSL_OLD
    stop 'ppatom_init nyi @ beagle, change interface to fgsl_spline_init.'
#endif

    ncenter = natom
    allocate(zeff(1:ncenter))
    allocate(rmin(1:ncenter))
    allocate(rmax(1:ncenter))
    allocate(rcut(1:ncenter))
    allocate(seig(1:ncenter))
    allocate(peig(1:ncenter))
    allocate(vs_acc(1:ncenter))
    allocate(vp_acc(1:ncenter))
    allocate(vloc_acc(1:ncenter))
    allocate(sfun_acc(1:ncenter))
    allocate(pfun_acc(1:ncenter))
    allocate(vs_spline(1:ncenter))
    allocate(vp_spline(1:ncenter))
    allocate(vloc_spline(1:ncenter))
    allocate(sfun_spline(1:ncenter))
    allocate(pfun_spline(1:ncenter))
    do A = 1, ncenter
       call ppatom_init1(A,ppfile(A),zeff(A),rmin(A),rmax(A),rcut(A),seig(A),peig(A), &
            vs_acc(A),vp_acc(A),vloc_acc(A),sfun_acc(A),pfun_acc(A),vs_spline(A), &
            vp_spline(A),vloc_spline(A),sfun_spline(A), pfun_spline(A))
    end do
  end subroutine ppatom_init
  !##########
  subroutine ppatom_init1(A,ppfile1,zeff1,rmin1,rmax1,rcut1,seig1,peig1, &
       vs_acc1,vp_acc1,vloc_acc1,sfun_acc1,pfun_acc1,vs_spline1, &
       vp_spline1,vloc_spline1,sfun_spline1,pfun_spline1)
    !-----
    integer(c_int), intent(in) :: A
    character(len=256), intent(in) :: ppfile1
    integer(c_int), intent(out) :: zeff1
    real(c_double), intent(out) :: rmin1,rmax1,rcut1,seig1,peig1
    type(fgsl_interp_accel), intent(out) :: vs_acc1,vp_acc1,vloc_acc1,sfun_acc1,pfun_acc1
    type(fgsl_spline), intent(out) :: vs_spline1,vp_spline1,vloc_spline1,sfun_spline1,pfun_spline1
    !-----
    character(len=256) :: sdum
    integer(c_int) :: irad, nlast
    integer(fgsl_size_t) :: nrad
    real(c_double), allocatable :: rad(:)
    real(c_double), allocatable :: ppot(:,:)
    real(c_double), allocatable :: pwfn(:,:)
    real(c_double), allocatable :: vwfn(:,:)
    real(c_double), allocatable :: pvwfn(:,:)

    type(fgsl_interp_accel) :: tmp_acc
    type(fgsl_spline) :: tmp_spline
    integer(fgsl_int) :: tmp_status

    ! read data
    open(unit=1, file=trim(ppfile1), status='old', form='formatted')
    read(1,*) sdum,zeff1
    read(1,*) sdum,nrad
    read(1,*)
    nlast = nrad - 1 ! count as 0:nlast
    allocate(rad(0:nlast))
    allocate(ppot(0:nlast,0:2))
    allocate(pwfn(0:nlast,0:1))
    allocate(vwfn(0:nlast,0:1))
    allocate(pvwfn(0:nlast,0:1))

    do irad = 0, nlast
       read(1,*) rad(irad),ppot(irad,0),ppot(irad,1),ppot(irad,2), &
            pwfn(irad,0),pwfn(irad,1)!,pwfn(irad,2)
    end do
    close(unit=1)

    vwfn = 0d0
    pvwfn = 0d0

    !>> 1st run: with radial factor 
    do irad = 1, nlast
       ppot(irad,0) = ppot(irad,0)*0.5d0 ! Rydberg to Hartree units
       ppot(irad,1) = ppot(irad,1)*0.5d0 ! Rydberg to Hartree units
       ppot(irad,2) = ppot(irad,2)*0.5d0 ! Rydberg to Hartree units
       vwfn(irad,0) = (ppot(irad,0)-ppot(irad,2))*pwfn(irad,0)/rad(irad)
       vwfn(irad,1) = (ppot(irad,1)-ppot(irad,2))*pwfn(irad,1)/rad(irad)
       pvwfn(irad,0) = pwfn(irad,0)*vwfn(irad,0)
       pvwfn(irad,1) = pwfn(irad,1)*vwfn(irad,1)
    end do

    rcut1 = -1d0
    do irad = 1, nlast
       if (max(abs(vwfn(irad,0)/rad(irad)),abs(vwfn(irad,1)/rad(irad))) < 1d-20) then
          rcut1 = rad(irad)
          exit
       end if
    end do
    !%%%%% safety factor %%%%%
    rcut1 = rcut1 + 1d0
    !%&%%%%%%%%%%%%%%%%%%%%%%%
    rmin1 = rcut1*1d-10
    rmax1 = rad(nlast)
    write(6,"('# A = ',i5,' zeff = ',i25)") A,zeff1
    write(6,"('# A = ',i5,' rcut = ',e25.10e2)") A,rcut1
    write(6,"('# A = ',i5,' rmin = ',e25.10e2)") A,rmin1
    write(6,"('# A = ',i5,' rmax = ',e25.10e2)") A,rmax1

    !>> cubic spline coefficients
    vs_acc1 = fgsl_interp_accel_alloc()
    vp_acc1 = fgsl_interp_accel_alloc()
    vloc_acc1 = fgsl_interp_accel_alloc()
    sfun_acc1 = fgsl_interp_accel_alloc()
    pfun_acc1 = fgsl_interp_accel_alloc()

    vs_spline1 = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
    vp_spline1 = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
    vloc_spline1 = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
    sfun_spline1 = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
    pfun_spline1 = fgsl_spline_alloc(fgsl_interp_cspline, nrad)

#ifndef FGSL_OLD
    tmp_status = fgsl_spline_init(vs_spline1, rad, ppot(:,0))
    tmp_status = fgsl_spline_init(vp_spline1, rad, ppot(:,1))
    tmp_status = fgsl_spline_init(vloc_spline1, rad, ppot(:,2))
    tmp_status = fgsl_spline_init(sfun_spline1, rad, vwfn(:,0))
    tmp_status = fgsl_spline_init(pfun_spline1, rad, vwfn(:,1))
#endif
    !<< cubic spline coefficients

    !>> denominator of Kleinman-Bylander transform
    tmp_acc = fgsl_interp_accel_alloc()
    tmp_spline = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
#ifndef FGSL_OLD
    tmp_status = fgsl_spline_init(tmp_spline, rad, pvwfn(:,0))
#endif
    seig1 = fgsl_spline_eval_integ(tmp_spline, rad(0), rad(nlast), tmp_acc)
#ifndef FGSL_OLD
    tmp_status = fgsl_spline_init(tmp_spline, rad, pvwfn(:,1))
#endif
    peig1 = fgsl_spline_eval_integ(tmp_spline, rad(0), rad(nlast), tmp_acc)
    call fgsl_spline_free(tmp_spline)
    call fgsl_interp_accel_free(tmp_acc)
    write(6,"('# A = ',i5,' seig = ',e25.10e2)") A,seig1
    write(6,"('# A = ',i5,' peig = ',e25.10e2)") A,peig1
    !<< denominator of Kleinman-Bylander transform

    !>> 2nd run: without radial factor 
    ppot(0,0) = ppatom_getvs(A,rmin(A))/rmin(A)
    ppot(0,1) = ppatom_getvp(A,rmin(A))/rmin(A)
    ppot(0,2) = ppatom_getvloc(A,rmin(A))/rmin(A)
    vwfn(0,0) = ppatom_getsfun(A,rmin(A))/rmin(A)
    !vwfn(0,1) = ppatom_getpfun(rmin)/rmin
    vwfn(0,1) = 0d0
    do irad = 1, nlast
       ppot(irad,0) = ppot(irad,0)/rad(irad)
       ppot(irad,1) = ppot(irad,1)/rad(irad)
       ppot(irad,2) = ppot(irad,2)/rad(irad)
       vwfn(irad,0) = vwfn(irad,0)/rad(irad)
       vwfn(irad,1) = vwfn(irad,1)/rad(irad)
    end do

    !>> cubic spline coefficients
#ifndef FGSL_OLD
    tmp_status = fgsl_spline_init(vs_spline1, rad, ppot(:,0))
    tmp_status = fgsl_spline_init(vp_spline1, rad, ppot(:,1))
    tmp_status = fgsl_spline_init(vloc_spline1, rad, ppot(:,2))
    tmp_status = fgsl_spline_init(sfun_spline1, rad, vwfn(:,0))
    tmp_status = fgsl_spline_init(pfun_spline1, rad, vwfn(:,1))
#endif
    !<< cubic spline coefficients

    deallocate(rad)
    deallocate(ppot)
    deallocate(pwfn)
    deallocate(vwfn)
    deallocate(pvwfn)
  end subroutine ppatom_init1
  !##########
  subroutine ppatom_final
    integer(c_int) :: A
    do A = 1, ncenter
       call fgsl_spline_free(vs_spline(A))
       call fgsl_spline_free(vp_spline(A))
       call fgsl_spline_free(vloc_spline(A))
       call fgsl_spline_free(sfun_spline(A))
       call fgsl_spline_free(pfun_spline(A))
       call fgsl_interp_accel_free(vs_acc(A))
       call fgsl_interp_accel_free(vp_acc(A))
       call fgsl_interp_accel_free(vloc_acc(A))
       call fgsl_interp_accel_free(sfun_acc(A))
       call fgsl_interp_accel_free(pfun_acc(A))
    end do
    deallocate(rmin)
    deallocate(rmax)
    deallocate(rcut)
    deallocate(zeff)
    deallocate(seig)
    deallocate(peig)
    deallocate(vs_acc)
    deallocate(vp_acc)
    deallocate(vloc_acc)
    deallocate(sfun_acc)
    deallocate(pfun_acc)
    deallocate(vs_spline)
    deallocate(vp_spline)
    deallocate(vloc_spline)
    deallocate(sfun_spline)
    deallocate(pfun_spline)    
  end subroutine ppatom_final
  !##########
  real(c_double) function ppatom_getvs(A,r)
    !----------
    integer(c_int), intent(in) :: A
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,yi
    if (abs(r) > rcut(A)) then
       ppatom_getvs = -zeff(A)/r
    else
       xi = r
       yi = fgsl_spline_eval(vs_spline(A), xi, vs_acc(A))
       ppatom_getvs = yi
    end if
  end function ppatom_getvs
  !##########
  real(c_double) function ppatom_getvp(A,r)
    !----------
    integer(c_int), intent(in) :: A
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,yi
    if (abs(r) > rcut(A)) then
       ppatom_getvp = -zeff(A)/r
    else
       xi = r
       yi = fgsl_spline_eval(vp_spline(A), xi, vp_acc(A))
       ppatom_getvp = yi
    end if
  end function ppatom_getvp
  !##########
  real(c_double) function ppatom_getvloc(A,r)
    !----------
    integer(c_int), intent(in) :: A
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,yi
    if (abs(r) > rcut(A)) then
       ppatom_getvloc = -zeff(A)/r
    else
       xi = r
       yi = fgsl_spline_eval(vloc_spline(A), xi, vloc_acc(A))
       ppatom_getvloc = yi
    end if
  end function ppatom_getvloc
  !##########
  real(c_double) function ppatom_getd1vloc(A,r)
    !----------
    integer(c_int), intent(in) :: A
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,dyi
    xi = r
    dyi = fgsl_spline_eval_deriv(vloc_spline(A), xi, vloc_acc(A))
    ppatom_getd1vloc = dyi
  end function ppatom_getd1vloc
  !##########
  real(c_double) function ppatom_getsfun(A,r)
    !----------
    integer(c_int), intent(in) :: A
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,yi
    if (abs(r) > rcut(A)) then
       ppatom_getsfun = 0d0
    else
       xi = r
       yi = fgsl_spline_eval(sfun_spline(A), xi, sfun_acc(A))
       ppatom_getsfun = yi
    end if
  end function ppatom_getsfun
  !##########
  real(c_double) function ppatom_getpfun(A,r)
    !----------
    integer(c_int), intent(in) :: A
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,yi
    if (abs(r) > rcut(A)) then
       ppatom_getpfun = 0d0
    else
       xi = r
       yi = fgsl_spline_eval(pfun_spline(A), xi, pfun_acc(A))
       ppatom_getpfun = yi
    end if
  end function ppatom_getpfun
  !##########
  real(c_double) function ppatom_seig(A)
    !----------
    integer(c_int), intent(in) :: A
    !----------
    ppatom_seig = seig(A)
  end function ppatom_seig
  !##########
  real(c_double) function ppatom_peig(A)
    !----------
    integer(c_int), intent(in) :: A
    !----------
    ppatom_peig = peig(A)
  end function ppatom_peig
  !##########
  real(c_double) function ppatom_rcut(A)
    !----------
    integer(c_int), intent(in) :: A
    !----------
    ppatom_rcut = rcut(A)
  end function ppatom_rcut
  !##########
  real(c_double) function ppatom_rmax(A)
    !----------
    integer(c_int), intent(in) :: A
    !----------
    ppatom_rmax = rmax(A)
  end function ppatom_rmax
  !##########
  real(c_double) function ppatom_rmin(A)
    !----------
    integer(c_int), intent(in) :: A
    !----------
    ppatom_rmin = rmin(A)
  end function ppatom_rmin
  !##########
!  subroutine ppatom_test(ppfile)
!    !-----
!    character(len=256), intent(in) :: ppfile
!    !-----
!    character(len=256) :: sdum
!    integer(c_int) :: irad, nrad, nlast
!    real(c_double), allocatable :: rad(:)
!    real(c_double), allocatable :: ppot(:,:)
!    real(c_double), allocatable :: pwfn(:,:)
!    real(c_double), allocatable :: vwfn(:,:)
!    real(c_double) :: norm(0:2), ovlp(0:2), val0,val1,val2,fac0,fac1,fac2,radp
!
!    ! read data
!    open(unit=1, file=trim(ppfile), status='old', form='formatted')
!    read(1,*) sdum,nrad
!    read(1,*)
!    nlast = nrad - 1 ! count as 0:nrad
!    allocate(rad(0:nlast))
!    allocate(ppot(0:nlast,0:2))
!    allocate(pwfn(0:nlast,0:2))
!    allocate(vwfn(0:nlast,0:2))
!    do irad = 0, nlast
!       read(1,*) rad(irad), ppot(irad,0:2), pwfn(irad,0:2)
!    end do
!    close(unit=1)
!
!    vwfn = 0d0
!    do irad = 1, nlast
!       ppot(irad,0) = ppot(irad,0)*0.5d0/rad(irad) ! Rydberg to Hartree units
!       ppot(irad,1) = ppot(irad,1)*0.5d0/rad(irad) ! Rydberg to Hartree units
!       ppot(irad,2) = ppot(irad,2)*0.5d0/rad(irad) ! Rydberg to Hartree units
!       vwfn(irad,0) = (ppot(irad,0)-ppot(irad,2))*pwfn(irad,0)
!       vwfn(irad,1) = (ppot(irad,1)-ppot(irad,2))*pwfn(irad,1)
!       vwfn(irad,2) =  ppot(irad,2)              *pwfn(irad,2)
!    end do
!
!    norm = 0d0
!    ovlp = 0d0
!
!    ! Simpson 1/3. See http://www5.synapse.ne.jp/windandwing/sub9.html
!    do irad = 2, nlast, 2
!       fac0 = (rad(irad)-rad(irad-2))*(-rad(irad)-2d0*rad(irad-2)+3d0*rad(irad-1))/(rad(irad-1)-rad(irad-2))
!       fac1 = (rad(irad)-rad(irad-2))**3/(rad(irad-1)-rad(irad-2))/(rad(irad)-rad(irad-1))
!       fac2 = (rad(irad)-rad(irad-2))*(2d0*rad(irad)+rad(irad-2)-3d0*rad(irad-1))/(rad(irad)-rad(irad-1))
!    
!       val0 = pwfn(irad-2,0)**2
!       val1 = pwfn(irad-1,0)**2
!       val2 = pwfn(irad,0)**2
!       norm(0) = norm(0) + val0*fac0+val1*fac1+val2*fac2
!    
!       val0 = pwfn(irad-2,1)**2
!       val1 = pwfn(irad-1,1)**2
!       val2 = pwfn(irad,1)**2
!       norm(1) = norm(1) + val0*fac0+val1*fac1+val2*fac2
!    
!       val0 = pwfn(irad-2,2)**2
!       val1 = pwfn(irad-1,2)**2
!       val2 = pwfn(irad,2)**2
!       norm(2) = norm(2) + val0*fac0+val1*fac1+val2*fac2
!    
!       val0 = pwfn(irad-2,0)*vwfn(irad-2,0)
!       val1 = pwfn(irad-1,0)*vwfn(irad-1,0)
!       val2 = pwfn(irad,0)  *vwfn(irad,0)
!       ovlp(0) = ovlp(0) + val0*fac0+val1*fac1+val2*fac2
!    
!       val0 = pwfn(irad-2,1)*vwfn(irad-2,1)
!       val1 = pwfn(irad-1,1)*vwfn(irad-1,1)
!       val2 = pwfn(irad,1)  *vwfn(irad,1)
!       ovlp(1) = ovlp(1) + val0*fac0+val1*fac1+val2*fac2
!    
!       val0 = pwfn(irad-2,2)*vwfn(irad-2,2)
!       val1 = pwfn(irad-1,2)*vwfn(irad-1,2)
!       val2 = pwfn(irad,2)  *vwfn(irad,2)  
!       ovlp(2) = ovlp(2) + val0*fac0+val1*fac1+val2*fac2
!    end do
!    norm = norm/6d0
!    ovlp = ovlp/6d0
!    
!    write(6,"('# norm:',3f20.10)") norm(0:2)
!    write(6,"('# ovlp:',3f20.10)") ovlp(0:2)
!    write(6,"('# ovlp:',3f20.10)") sign(1d0,ovlp(0)),sign(1d0,ovlp(1)),sign(1d0,ovlp(2))
!
!    write(6,"(6E25.15)") rad(0), ppot(1,2), vwfn(1,0:1)/rad(1), vwfn(1,0:1)/rad(1)/sqrt(abs(ovlp(0:1)))
!    radp = rad(0)
!    do irad = 1, nrad
!       if (rad(irad)-radp > 1D-4) then
!          vwfn(irad,0:1) = vwfn(irad,0:1)/rad(irad)
!          write(6,"(6E25.15)") rad(irad), ppot(irad,2), vwfn(irad,0:1), vwfn(irad,0:1)/sqrt(abs(ovlp(0:1)))
!          radp = rad(irad)
!       end if
!!       if (max(abs(vwfn(irad,0)),abs(vwfn(irad,0))) < 1D-20) exit
!    end do
!    
!    deallocate(rad)
!    deallocate(ppot)
!    deallocate(pwfn)
!    deallocate(vwfn)
!  end subroutine ppatom_test
  !##########
end module mod_ppatom
!###############
