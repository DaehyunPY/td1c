!################################################################################
subroutine bas_pp_gen_hgh()

  use, intrinsic :: iso_c_binding
  use mod_rad, only : nrad, xrad, wrad
  use mod_bas, only : znuc, psp_label, pp_znuc, pp_rloc, pp_cloc, pp_rproj, pp_hproj, &
       pp_vloc, pp_pproj, pp_fproj, pp_gproj, pp_maxl, pp_nump, pp_irmax
  implicit none
  !--------------------------------------------------------------------
  !--------------------------------------------------------------------
  logical(c_bool) :: found
  integer(c_long) :: irad,l,i,j
  real(c_double) :: thresh = 1.D-6
  real(c_double) :: thresh2 = 1.D-20

!  real(c_double) :: rloc, cloc(1:4), rproj(0:2), hproj(1:3,1:3,0:2)

  call hgh_get_znuc(psp_label,pp_znuc)
  call hgh_get_rloc(psp_label,pp_rloc)
  call hgh_get_cloc(psp_label,pp_cloc)
  call hgh_get_rproj(psp_label,pp_rproj)
  call hgh_get_hproj(psp_label,pp_hproj)

  pp_maxl = -1
  pp_nump = 0
  if (pp_rproj(2) > thresh) then
     pp_maxl = 2
  else if (pp_rproj(1) > thresh) then
     pp_maxl = 1
  else if (pp_rproj(0) > thresh) then
     pp_maxl = 0
  end if
  do l = 0, pp_maxl
     if (abs(pp_hproj(3,3,l)) > thresh) then
        pp_nump(l) = 3
     else if (abs(pp_hproj(2,2,l)) > thresh) then
        pp_nump(l) = 2
     else if (abs(pp_hproj(1,1,l)) > thresh) then
        pp_nump(l) = 1
     end if
  end do

  ! local part
  do irad = 1, nrad - 1
     pp_vloc(irad) = &
          - pp_znuc/xrad(irad)*erf(xrad(irad)/sqrt(2.d+0)/pp_rloc) &
          + exp(-0.5d+0*(xrad(irad)/pp_rloc)**2.d+0) &
          * (pp_cloc(1) &
           + pp_cloc(2)*(xrad(irad)/pp_rloc)**2.d+0 &
           + pp_cloc(3)*(xrad(irad)/pp_rloc)**4.d+0 &
           + pp_cloc(4)*(xrad(irad)/pp_rloc)**6.d+0)
  end do

  ! nonlocal part
  pp_irmax = 1
  pp_pproj = 0.d+0
  pp_fproj = 0.d+0
  pp_gproj = 0.d+0
  do l = 0, pp_maxl
     do i = 1, pp_nump(l)
        found = .false.
        do irad = 1, nrad - 1
           pp_pproj(irad,i,l) = &
                sqrt(2.d+0) * xrad(irad)**(l+2*(i-1)) * exp(-(xrad(irad)/pp_rproj(l))**2/2.d+0) &
                / (pp_rproj(l)**(l+(4*i-1)/2.d+0) * sqrt(gamma(l+(4*i-1)/2.d+0)))
           pp_fproj(irad,i,l) = pp_pproj(irad,i,l) * xrad(irad) * sqrt(wrad(irad))
           do j = 1, pp_nump(l)
              pp_gproj(irad,i,l) = pp_gproj(irad,i,l) + pp_hproj(i,j,l)*pp_fproj(irad,j,l)
           end do
           if (.not.found .and. &
                abs(pp_pproj(irad,i,l)) < thresh2 .and. &
                abs(pp_fproj(irad,i,l)) < thresh2 .and. &
                abs(pp_gproj(irad,i,l)) < thresh2) then
              found = .true.
              pp_irmax(i,l) = irad
           end if
        end do
     end do
  end do

!debug
  do irad = 1, nrad - 1
     write(6,"('vloc: ',2f20.10)") xrad(irad),pp_vloc(irad)
  end do
!  write(6,"('pp_maxl  = ', i10)") pp_maxl
!  write(6,"('pp_nump  = ', 3i10)") pp_nump(0:2)
!  do l = 0, pp_maxl
!     do i = 1, 3!pp_nump(l)
!        write(6,"('pp_rmax = ',2i5,i10,4e12.5)") i,l,pp_irmax(i,l),xrad(pp_irmax(i,l)), &
!             pp_pproj(pp_irmax(i,l),i,l), &
!             pp_fproj(pp_irmax(i,l),i,l), &
!             pp_gproj(pp_irmax(i,l),i,l)
!     end do
!     write(6,"('pp_hproj: l = ', i5)") l
!     do i = 1, 3!pp_nump(l)
!     do j = 1, 3!pp_nump(l)
!        write(6, "(e12.5)", advance='no') pp_hproj(j,i,l)
!     end do
!     write(6,*)
!     end do
!     write(6,*)
!  end do
!
!  open(unit=99,file=trim(trim(name)//".pp"),status='unknown',form='formatted')
!  do irad = 1, nrad - 1
!     write(99,"(f10.5, 2e20.8)",advance='no') xrad(irad), -znuc/xrad(irad), pp_vloc(irad)
!     do l = 0, pp_maxl
!        do i = 1, pp_nump(l)
!           write(99,"(e20.8)",advance='no') pp_pproj(irad,i,l)
!        end do
!     end do
!     write(99,*)
!  end do
!  close(99)
!debug

  contains
  !######################################################################
  subroutine hgh_get_znuc(label, znuc)
    implicit none
    !--------------------------------------------------------------------
    integer(c_long), intent(in) :: label
    real(c_double), intent(out) :: znuc
    !--------------------------------------------------------------------
    if (label == 1) then
       znuc = 1d0
    else if (label == 2) then
       znuc = 2d0
    else if (label == 3) then
       znuc = 1d0
    else if (label == -3) then
       znuc = 3d0
    else if (label == 4) then
       znuc = 2d0
    else if (label == -4) then
       znuc = 4d0
    else if (label == 5) then
       znuc = 3d0
    else if (label == 6) then
       znuc = 4d0
    else if (label == 7) then
       znuc = 5d0
    else if (label == 8) then
       znuc = 6d0
    else if (label == 9) then
       znuc = 7d0
    else if (label == 10) then
       znuc = 8d0
    else if (label == 18) then
       znuc = 8d0
    else if (label == 36) then
       znuc = 8d0
    else
       stop 'hgh_get_znuc: bad label'
    end if
  end subroutine hgh_get_znuc
  !######################################################################
  subroutine hgh_get_rloc(label, rloc)
    implicit none
    !--------------------------------------------------------------------
    integer(c_long), intent(in) :: label
    real(c_double), intent(out) :: rloc
    !--------------------------------------------------------------------
    if (label == 1) then
       rloc = 0.2d0
    else if (label == 2) then
       rloc = 0.2d0
    else if (label == 3) then
       rloc = 0.787553d0
    else if (label == -3) then
       rloc = 0.4d0
    else if (label == 4) then
       rloc = 0.739009d0
    else if (label == -4) then
       rloc = 0.325d0
    else if (label == 5) then
       rloc = 0.433930d0
    else if (label == 6) then
       rloc = 0.348830d0
    else if (label == 7) then
       rloc = 0.289179d0
    else if (label == 8) then
       rloc = 0.247621d0
    else if (label == 9) then
       rloc = 0.218525d0
    else if (label == 10) then
       rloc = 0.19d0
    else if (label == 18) then
       rloc = 0.4d0
    else if (label == 36) then
       rloc = 0.5d0
    else
       stop 'hgh_get_rloc: bad label'
    end if
  end subroutine hgh_get_rloc
  !######################################################################
  subroutine hgh_get_cloc(label, cloc)
    implicit none
    !--------------------------------------------------------------------
    integer(c_long), intent(in) :: label
    real(c_double), intent(out) :: cloc(1:4)
    !--------------------------------------------------------------------
    cloc = 0d0
    if (label == 1) then
       cloc(1) =  -4.180237d0; cloc(2) = 0.725075d0
    else if (label == 2) then
       cloc(1) =  -9.112023d0; cloc(2) = 1.698368d0
    else if (label == 3) then
       cloc(1) = -1.892612d0; cloc(2) = 0.286060d0
    else if (label == -3) then
       cloc(1) = -14.034868d0; cloc(2) = 9.553476d0
       cloc(3) =  -1.766488d0; cloc(4) = 0.084370d0
    else if (label == 4) then
       cloc(1) = -2.592951d0; cloc(2) = 0.354839d0
    else if (label == -4) then
       cloc(1) = -24.015041d0; cloc(2) = 17.204014d0
       cloc(3) =  -3.326390d0; cloc(4) = 0.165419d0
    else if (label == 5) then
       cloc(1) = -5.578642d0; cloc(2) = 0.804251d0
    else if (label == 6) then
       cloc(1) = -8.513771d0; cloc(2) = 1.228432d0
    else if (label == 7) then
       cloc(1) = -12.234820d0; cloc(2) = 1.766407d0
    else if (label == 8) then
       cloc(1) = -16.580318d0; cloc(2) = 2.395701d0
    else if (label == 9) then
       cloc(1) = -21.307361d0; cloc(2) = 3.072869d0
    else if (label == 10) then
       cloc(1) = -27.692852d0; cloc(2) = 4.005906d0
    else if (label == 18) then
       cloc(1) = -7.1d0
    else if (label == 36) then
    else
       stop 'hgh_get_cloc: bad label'
    end if
  end subroutine hgh_get_cloc
  !######################################################################
  subroutine hgh_get_rproj(label, rproj)
    implicit none
    !--------------------------------------------------------------------
    integer(c_long), intent(in) :: label
    real(c_double), intent(out) :: rproj(0:2)
    !--------------------------------------------------------------------
    rproj = 0d0
    if (label == 1) then
    else if (label == 2) then
    else if (label == 3) then
       rproj(0) = 0.666375d0; rproj(1) = 1.079306d0
    else if (label == -3) then
    else if (label == 4) then
       rproj(0) = 0.528797d0; rproj(1) = 0.658153d0
    else if (label == -4) then
    else if (label == 5) then
       rproj(0) = 0.373843d0; rproj(1) = 0.360393d0
    else if (label == 6) then
       rproj(0) = 0.304553d0; rproj(1) = 0.232677d0
    else if (label == 7) then
       rproj(0) = 0.256605d0; rproj(1) = 0.270134d0
    else if (label == 8) then
       rproj(0) = 0.221786d0; rproj(1) = 0.256829d0
    else if (label == 9) then
       rproj(0) = 0.195567d0; rproj(1) = 0.174268d0
    else if (label == 10) then
       rproj(0) = 0.179488d0; rproj(1) = 0.214913d0
    else if (label == 18) then
       rproj(0) = 0.317381d0; rproj(1) = 0.351619d0
    else if (label == 36) then
       rproj(0) = 0.410759d0; rproj(1) = 0.430256d0; rproj(2) = 0.517120d0
    else
       stop 'hgh_get_rproj: bad label'
    end if
  end subroutine hgh_get_rproj
  !######################################################################
  subroutine hgh_get_hproj(label, hproj)
    implicit none
    !--------------------------------------------------------------------
    integer(c_long), intent(in) :: label
    real(c_double), intent(out) :: hproj(1:3,1:3,0:2)
    !--------------------------------------------------------------------
    hproj = 0d0
    if (label == 1) then
    else if (label == 2) then
    else if (label == 3) then
       hproj(1,1,0) = 1.858811d0
       hproj(1,1,1) = -0.005895d0
    else if (label == -3) then
    else if (label == 4) then
       hproj(1,1,0) = 3.061666d0
       hproj(1,1,1) = 0.092462d0
    else if (label == -4) then
    else if (label == 5) then
       hproj(1,1,0) = 6.233928d0
    else if (label == 6) then
       hproj(1,1,0) = 9.522842d0
    else if (label == 7) then
       hproj(1,1,0) = 13.552243d0
    else if (label == 8) then
       hproj(1,1,0) = 18.266917d0
    else if (label == 9) then
       hproj(1,1,0) = 23.584942d0
    else if (label == 10) then
       hproj(1,1,0) = 28.506098d0
       hproj(2,2,0) = -1.076245d0
       hproj(1,1,1) = -0.000090d0
    else if (label == 18) then
       hproj(1,1,0) = 10.249487d0
       hproj(2,2,0) =  5.602516d0
       hproj(1,1,1) =  4.978801d0
    else if (label == 36) then
       hproj(1,1,0) =  5.911194d0
       hproj(2,2,0) =  1.967372d0
       hproj(3,3,0) = -1.458069d0
       hproj(1,1,1) =  3.524357d0
       hproj(2,2,1) = -0.691198d0
       hproj(1,1,2) =  0.629228d0
    else
       stop 'hgh_get_hproj: bad label'
    end if

    hproj(1,2,0) = -0.5d0*sqrt(3d0/5d0)*hproj(2,2,0)
    hproj(1,3,0) = +0.5d0*sqrt(5d0/21d0)*hproj(3,3,0)
    hproj(2,3,0) = -0.5d0*sqrt(100d0/63d0)*hproj(3,3,0)
    hproj(1,2,1) = -0.5d0*sqrt(5d0/7d0)*hproj(2,2,1)
    hproj(1,3,1) = +1d0/6d0*sqrt(35d0/11d0)*hproj(3,3,1)
    hproj(2,3,1) = -14d0/6d0/sqrt(11d0)*hproj(3,3,1)
    hproj(1,2,2) = -0.5d0*sqrt(7d0/9d0)*hproj(2,2,2)
    hproj(1,3,2) = +0.5d0*sqrt(63d0/143d0)*hproj(3,3,2)
    hproj(2,3,2) = -9d0/sqrt(143d0)*hproj(3,3,2)

    hproj(2,1,0) = hproj(1,2,0)
    hproj(3,1,0) = hproj(1,3,0)
    hproj(3,2,0) = hproj(2,3,0)
    hproj(2,1,1) = hproj(1,2,1)
    hproj(3,1,1) = hproj(1,3,1)
    hproj(3,2,1) = hproj(2,3,1)
    hproj(2,1,2) = hproj(1,2,2)
    hproj(3,1,2) = hproj(1,3,2)
    hproj(3,2,2) = hproj(2,3,2)
  end subroutine hgh_get_hproj
  !######################################################################
end subroutine bas_pp_gen_hgh
!################################################################################
