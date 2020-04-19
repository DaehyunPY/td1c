!###############
module mod_oncvpsp
  use iso_c_binding
  use fgsl
  implicit none

  public oncvpsp_init
  public oncvpsp_final
  public oncvpsp_getproj
  public oncvpsp_getvloc
  public oncvpsp_getd1vloc
  public oncvpsp_ekb
  public oncvpsp_maxl
  public oncvpsp_nump
  public oncvpsp_rcut
  public oncvpsp_rmax

  private
  integer(c_int) :: maxl,nproj(0:4)
  real(c_double) :: zeff,rmin,rmax,rcut,ekb(1:2,0:4)
  type(fgsl_interp_accel) :: proj_acc(1:2,0:4),vloc_acc
  type(fgsl_spline) :: proj_spline(1:2,0:4),vloc_spline

  contains
  !##########
  subroutine oncvpsp_init(ppfile)
    implicit none
    !-----
    character(len=256), intent(in) :: ppfile
    !-----
    integer(c_int) :: A

    character(len=256) :: line
    character(len=256) :: sdum
    integer(c_int) :: irad, l, iproj, lloc, nrad, nlast
    real(c_double), allocatable :: rad(:),proj(:,:,:),vloc(:)
    integer(fgsl_int) :: tmp_status

#ifdef FGSL_OLD
    stop 'oncvpsp_init nyi @ beagle, change interface to fgsl_spline_init.'
#endif

    ! read data
    open(unit=1, file=trim(ppfile), status='old', form='formatted')
    read(1,*)
    read(1,*) sdum,zeff,sdum,sdum
    read(1,*) sdum,sdum,maxl,lloc,nrad,sdum,sdum
    read(1,*)
    read(1,*) nproj(0:4),sdum
    read(1,*)

!debug
    write(6,"('zeff  = ',f20.10)") zeff
    write(6,"('maxl  = ',i20)") maxl
    write(6,"('lloc  = ',i20)") lloc
    write(6,"('nrad  = ',i20)") nrad
    do l = 0, maxl
       write(6,"('nproj = ',2i20)") l, nproj(l)
    end do
!    stop
!debug

    nlast = nrad - 1 ! count as 0:nlast
    allocate(rad(0:nlast))
    allocate(proj(0:nlast,1:2,0:4))
    allocate(vloc(0:nlast))

    do l = 0, maxl
       read(1,*) sdum,ekb(1:nproj(l),l)
       do irad = 0, nlast
          read(1,*) sdum,rad(irad),proj(irad,1:nproj(l),l)
       end do
    end do
    read(1,*) sdum
    do irad = 0, nlast
       read(1,*) sdum,rad(irad),vloc(irad)
    end do

    close(unit=1)

!debug
!    do l = 0, maxl
!       do iproj = 1, nproj(l)
!          write(6,"('ekb: ',2i5,f15.10)") l,iproj,ekb(iproj,l)
!       end do
!    end do
!    do irad = 0, nlast
!       write(6,"('prj: ',f10.5)",advance='no') rad(irad)
!       do l = 0, maxl
!          do iproj = 1, nproj(l)
!             write(6,"(f15.10)",advance='no') proj(irad,iproj,l)
!          end do
!       end do
!       write(6,"(f15.10)") vloc(irad)
!    end do
!    stop
!debug

    rcut = rad(nlast)

    do l = 0, maxl
       do iproj = 1, nproj(l)
          proj_acc(iproj,l) = fgsl_interp_accel_alloc()
          proj_spline(iproj,l) = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
#ifndef FGSL_OLD
          tmp_status = fgsl_spline_init(proj_spline(iproj,l), rad, proj(:,iproj,l))
#endif
       end do
    end do
    vloc_acc = fgsl_interp_accel_alloc()
    vloc_spline = fgsl_spline_alloc(fgsl_interp_cspline, nrad)
#ifndef FGSL_OLD
    tmp_status = fgsl_spline_init(vloc_spline, rad, vloc(:))
#endif

    deallocate(vloc)
    deallocate(proj)
    deallocate(rad)
  end subroutine oncvpsp_init
  !##########
  subroutine oncvpsp_final
    integer(c_int) :: iproj, l
    call fgsl_spline_free(vloc_spline)
    call fgsl_interp_accel_free(vloc_acc)
    do l = 0, maxl
       do iproj = 1, nproj(l)
          call fgsl_spline_free(proj_spline(iproj,l))
          call fgsl_interp_accel_free(proj_acc(iproj,l))
       end do
    end do
  end subroutine oncvpsp_final
  !##########
  real(c_double) function oncvpsp_getproj(i,l,r)
    !----------
    integer(c_int), intent(in) :: i,l
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,yi
    if (abs(r) > rcut) then
       oncvpsp_getproj = 0.0
    else
       xi = r
       yi = fgsl_spline_eval(proj_spline(i,l), xi, proj_acc(i,l))
       oncvpsp_getproj = yi
    end if
  end function oncvpsp_getproj
  !##########
  real(c_double) function oncvpsp_getvloc(r)
    !----------
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,yi
    if (abs(r) > rcut) then
       oncvpsp_getvloc = -zeff/r
    else
       xi = r
       yi = fgsl_spline_eval(vloc_spline, xi, vloc_acc)
       oncvpsp_getvloc = yi
    end if
  end function oncvpsp_getvloc
  !##########
  real(c_double) function oncvpsp_getd1vloc(r)
    !----------
    real(c_double), intent(in) :: r
    !----------
    real(fgsl_double) :: xi,dyi
    if (abs(r) > rcut) then
       oncvpsp_getd1vloc = zeff/(r*r)
    else
       xi = r
       dyi = fgsl_spline_eval_deriv(vloc_spline, xi, vloc_acc)
       oncvpsp_getd1vloc = dyi
    end if
  end function oncvpsp_getd1vloc
  !##########
  real(c_double) function oncvpsp_ekb(i,l)
    !----------
    integer(c_int), intent(in) :: i,l
    !----------
    oncvpsp_ekb = ekb(i,l)
  end function oncvpsp_ekb
  !##########
  integer(c_int) function oncvpsp_nump(l)
    !----------
    integer(c_int), intent(in) :: l
    !----------
    oncvpsp_nump = nproj(l)
  end function oncvpsp_nump
  !##########
  integer(c_int) function oncvpsp_maxl()
    oncvpsp_maxl = maxl
  end function oncvpsp_maxl
  !##########
  real(c_double) function oncvpsp_rcut()
    oncvpsp_rcut = rcut
  end function oncvpsp_rcut
  !##########
  real(c_double) function oncvpsp_rmax()
    oncvpsp_rmax = rmax
  end function oncvpsp_rmax
  !##########
end module mod_oncvpsp
!###############
