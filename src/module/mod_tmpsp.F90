!###############
module mod_tmpsp
  use iso_c_binding
  use fgsl
  implicit none

  public tmpsp_init
  public tmpsp_final
  public tmpsp_getproj
  public tmpsp_getvloc
  public tmpsp_getd1vloc
  public tmpsp_ekb
  public tmpsp_maxl
  public tmpsp_rcut
  public tmpsp_rmax

  private
  integer(c_int) :: maxl,lloc
  real(c_double) :: zeff,rmin,rmax,rcut,ekb(0:4)
  type(fgsl_interp_accel) :: proj_acc(0:4),vloc_acc
  type(fgsl_spline) :: proj_spline(0:4),vloc_spline

  contains
  !##########
  subroutine tmpsp_init(ppfile)
    use mod_const, only : pi
    implicit none
    !-----
    character(len=256), intent(in) :: ppfile
    !-----
    integer(c_int) :: A

    character(len=256) :: sdum
    integer(c_int) :: irad, l, nrad, nlast, lmaxt
    real(c_double) :: norm1,norm2,delv,delr
    real(c_double), allocatable :: rad(:),pwfn(:,:),vloc(:,:),proj(:,:),proj2(:)
    integer(fgsl_int) :: tmp_status
    type(fgsl_interp_accel) :: proj2_acc
    type(fgsl_spline) :: proj2_spline

#ifdef FGSL_OLD
    stop 'tmpsp_init nyi @ beagle, change interface to fgsl_spline_init.'
#endif

    ! read data
    open(unit=1, file=trim(ppfile), status='old', form='formatted')
    read(1,*)
    read(1,*) sdum,zeff,sdum,sdum
    read(1,*) sdum,sdum,maxl,lloc,nrad,sdum,sdum
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)

!debug
    write(6,"('zeff  = ',f20.10)") zeff
    write(6,"('maxl  = ',i20)") maxl
    write(6,"('lloc  = ',i20)") lloc
    write(6,"('nrad  = ',i20)") nrad
!debug

    nlast = nrad    ! count as 0:nlast
    nrad = nrad + 1 ! count as 0:nlast
    allocate(rad(0:nlast))
    allocate(vloc(0:nlast,0:4))
    allocate(pwfn(0:nlast,0:4))
    allocate(proj(0:nlast,0:4))

    do l = 0, maxl
       read(1,*) sdum,sdum
       rad(0) = 0.0
       pwfn(0,l) = 0.0
       do irad = 1, nlast
          read(1,*) sdum,rad(irad),pwfn(irad,l),vloc(irad,l)
       end do
       vloc(0,l) = vloc(1,l) - (vloc(2,l)-vloc(1,l))/(rad(2)-rad(1))*rad(1)
    end do
    close(unit=1)

    proj = 0d0
    do l = 0, maxl
       do irad = 1, nlast
          proj(irad,l) = (vloc(irad,l)-vloc(irad,lloc))*pwfn(irad,l)
       end do
    end do

    rcut = rad(nlast)

    allocate(proj2(0:nlast))
    proj2_acc = fgsl_interp_accel_alloc()
    proj2_spline = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
    do l = 0, maxl
       proj_acc(l) = fgsl_interp_accel_alloc()
       proj_spline(l) = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
#ifndef FGSL_OLD
       tmp_status = fgsl_spline_init(proj_spline(l), rad, proj(:,l))
#endif
       ekb(l) = 0d0;
       if (l .ne. lloc) then
          proj2(0) = 0d0
          do irad = 1, nlast
             proj2(irad) = pwfn(irad,l) * proj(irad,l)
          end do
          tmp_status = fgsl_spline_init(proj2_spline, rad, proj2(:))
          ekb(l) = 1d0/fgsl_spline_eval_integ(proj2_spline, rad(1), rad(nlast), proj2_acc)
          write(6,"('l = ',i5,': ekb = ',f20.10)") l, ekb(l)
       end if
    end do
    call fgsl_spline_free(proj2_spline)
    call fgsl_interp_accel_free(proj2_acc)
    deallocate(proj2)

    vloc_acc = fgsl_interp_accel_alloc()
    vloc_spline = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
#ifndef FGSL_OLD
    tmp_status = fgsl_spline_init(vloc_spline, rad, vloc(:,lloc))
#endif

    deallocate(proj)
    deallocate(pwfn)
    deallocate(vloc)
    deallocate(rad)
  end subroutine tmpsp_init
  !##########
  subroutine tmpsp_final
    integer(c_int) :: l
    call fgsl_spline_free(vloc_spline)
    call fgsl_interp_accel_free(vloc_acc)
    do l = 0, maxl
       call fgsl_spline_free(proj_spline(l))
       call fgsl_interp_accel_free(proj_acc(l))
    end do
  end subroutine tmpsp_final
  !##########
  real(c_double) function tmpsp_getproj(l,r)
    !----------
    integer(c_int), intent(in) :: l
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,yi
    if (l == lloc) then
       tmpsp_getproj = 0.0
    else if (abs(r) > rcut) then
       tmpsp_getproj = 0.0
    else
       xi = r
       yi = fgsl_spline_eval(proj_spline(l), xi, proj_acc(l))
       tmpsp_getproj = yi
    end if
  end function tmpsp_getproj
  !##########
  real(c_double) function tmpsp_getvloc(r)
    !----------
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,yi
    if (abs(r) > rcut) then
       tmpsp_getvloc = -zeff/r
    else
       xi = r
       yi = fgsl_spline_eval(vloc_spline, xi, vloc_acc)
       tmpsp_getvloc = yi
    end if
  end function tmpsp_getvloc
  !##########
  real(c_double) function tmpsp_getd1vloc(r)
    !----------
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,dyi
    if (abs(r) > rcut) then
       tmpsp_getd1vloc = zeff/(r*r)
    else
       xi = r
       dyi = fgsl_spline_eval_deriv(vloc_spline, xi, vloc_acc)
       tmpsp_getd1vloc = dyi
    end if
  end function tmpsp_getd1vloc
  !##########
  real(c_double) function tmpsp_ekb(l)
    !----------
    integer(c_int), intent(in) :: l
    !----------
    tmpsp_ekb = ekb(l)
  end function tmpsp_ekb
  !##########
  integer(c_int) function tmpsp_maxl()
    tmpsp_maxl = maxl
  end function tmpsp_maxl
  !##########
  real(c_double) function tmpsp_rcut()
    tmpsp_rcut = rcut
  end function tmpsp_rcut
  !##########
  real(c_double) function tmpsp_rmax()
    tmpsp_rmax = rmax
  end function tmpsp_rmax
  !##########
end module mod_tmpsp
!###############
