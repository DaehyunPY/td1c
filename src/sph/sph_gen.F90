!#######################################################################
subroutine sph_gen(lmax1, lmax2, mmax1, mmax2, nlat, nphi, &
     wlat, wphi, wang, cost, sint, legf1, legb1, legf2, legb2)

  use, intrinsic :: iso_c_binding
  use mod_control, only : exact3j
  use mod_const, only : PI
  use shtools

  implicit none
  integer(c_long), intent(in) :: lmax1, lmax2, mmax1, mmax2, nlat, nphi
  real(c_double),  intent(out) :: wlat(1:nlat), wphi, wang(1:nlat, 1:nphi), cost(1:nlat), sint(1:nlat)
  real(c_double),  intent(out) :: legf1(0:lmax1, 1:nlat, -mmax1:mmax1)
  real(c_double),  intent(out) :: legb1(1:nlat, 0:lmax1, -mmax1:mmax1)
  real(c_double),  intent(out) :: legf2(0:lmax2, 1:nlat, -mmax2:mmax2)
  real(c_double),  intent(out) :: legb2(1:nlat, 0:lmax2, -mmax2:mmax2)
  real(c_double), allocatable :: SHT_p1x(:,:)
  integer(c_long) :: ilat, iphi, l, m, lm
  !hack for SHTOOLS
  integer(4) :: lmax2_int4, i0, i1, i4
  !hack for SHTOOLS

!  if (exact3j) then
!     call sph_gen_x3j(lmax1, lmax2, mmax1, mmax2, nlat, nphi, &
!     wlat, wphi, wang, cost, sint, legf1, legb1, legf2, legb2)
!     return
!  end if

  allocate(SHT_p1x(1:nlat, 1:((lmax2+1)*(lmax2+2))/2))

  i0 = 0
  i1 = 1
  i4 = 4
  lmax2_int4 = lmax2
  wphi = 2.D+0 * PI / nphi;
!  call SHGLQ(lmax2_int4, cost, wlat, SHT_p1x, i4, -i1, i0)
  call SHGLQ(lmax2_int4, cost, wlat, SHT_p1x, i4, -i1, i1)

!debug
!  do l = 0, lmax1
!     do m = -l, l
!        lm = (l*(l+1))/2 + m + 1
!        do ilat = 1, nlat
!           write(6,"(4i5,2f20.10)") l,m,lm,ilat,cost(ilat),SHT_p1x(ilat,lm)
!        end do
!     end do
!  end do
!  stop
!debug

  do iphi = 1, nphi
     do ilat = 1, nlat
        wang(ilat, iphi) = wlat(ilat)
!       wang(ilat, iphi) = wlat(ilat) * wphi
     end do
  end do

  do ilat = 1, nlat
     sint(ilat) = sqrt(1-cost(ilat)**2.D+0)
  end do

  ! for positive magnetic quantum numbers
  do m = 0, mmax1
     do l = abs(m), lmax1
        lm = (l*(l+1))/2 + m + 1
        do ilat = 1, nlat
           legf1(l, ilat, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI)
           legb1(ilat, l, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * wlat(ilat)
        end do
     end do
  end do
  do m = 0, mmax2
     do l = abs(m), lmax2
        lm = (l*(l+1))/2 + m + 1
        do ilat = 1, nlat
! 2016/9/5
           legf2(l, ilat, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI)
           legb2(ilat, l, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * wlat(ilat)
!this doesn't work...
!           legf2(l, ilat, m) = SHT_p1x(ilat, lm) * (2.D+0 * PI)
!           legb2(ilat, l, m) = SHT_p1x(ilat, lm) * wlat(ilat)
!this doesn't work...
! 2016/9/5
        end do
     end do
  end do

  ! for negative magnetic quantum numbers
  do m = 1, mmax1
     do l = abs(m), lmax1
        do ilat = 1, nlat
!           legf1(l, ilat, -m) = -legf1(l, ilat, m)
!           legb1(ilat, l, -m) = -legb1(ilat, l, m)
           legf1(l, ilat, -m) = legf1(l, ilat, m)
           legb1(ilat, l, -m) = legb1(ilat, l, m)
        end do
     end do
  end do
  do m = 1, mmax2
     do l = abs(m), lmax2
        do ilat = 1, nlat
!           legf2(l, ilat, -m) = -legf2(l, ilat, m)
!           legb2(ilat, l, -m) = -legb2(ilat, l, m)
           legf2(l, ilat, -m) = legf2(l, ilat, m)
           legb2(ilat, l, -m) = legb2(ilat, l, m)
        end do
     end do
  end do

!DEBUG
!  write(6, "('sph_gen: nodes and weights')")
!  do ilat = 1, nlat
!     write(6, "(i10, 4f20.10)") ilat, acos(cost(ilat))/PI*180.D+0, cost(ilat), sint(ilat), wlat(ilat)
!  end do
!
!  write(6, "('sph_gen: Associated Legendre polynomials')")
!  do ilat = 1, nlat
!     write(6, "('ilat:', i10)", advance='no') ilat
!     lm = 0
!     do l = 0, lmax2
!        do m = 0, +l
!           lm = lm + 1
!           write(6, "(f20.10)", advance='no') SHT_p1x(ilat, lm)
!        end do
!     end do
!     write(6, *)
!  end do
!
!  write(6, "('sph_gen: forward Legendre transform')")
!! do m = 0, lmax2
!  do m = 0, mmax2
!     write(6, "('m =  ', i10)", advance = 'no') m
!     do ilat = 1, nlat
!        write(6, "(i20)", advance='no') ilat
!     end do
!     write(6, *)
!     do l = abs(m), lmax2
!        lm = (l*(l+1))/2 + m + 1
!        write(6, "('l:   ', i10)", advance='no') l
!        do ilat = 1, nlat
!           write(6, "(f20.10)", advance='no') SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI)
!        end do
!        write(6, *)
!     end do
!  end do
!
!  write(6, "('sph_gen: backward Legendre transform')")
!! do m = 0, lmax2
!  do m = 0, mmax2
!     write(6, "('m =  ', i10)", advance = 'no') m
!     do l = abs(m), lmax2
!        write(6, "(i20)", advance='no') l
!     end do
!     write(6, *)
!     do ilat = 1, nlat
!        write(6, "('ilat:', i10)", advance='no') ilat
!        do l = abs(m), lmax2
!           lm = (l*(l+1))/2 + m + 1
!           write(6, "(f20.10)", advance='no') SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * wlat(ilat)
!        end do
!        write(6, *)
!     end do
!  end do
!DEBUG

  deallocate(SHT_p1x)

end subroutine sph_gen
!#######################################################################
subroutine sph_gen_x3j(lmax1, lmax2, mmax1, mmax2, nlat, nphi, &
     wlat, wphi, wang, cost, sint, legf1, legb1, legf2, legb2)

  use, intrinsic :: iso_c_binding
  use mod_const, only : PI
  use shtools

  implicit none
  integer(c_long), intent(in) :: lmax1, lmax2, mmax1, mmax2, nlat, nphi
  real(c_double),  intent(out) :: wlat(1:nlat), wphi, wang(1:nlat, 1:nphi), cost(1:nlat), sint(1:nlat)
  real(c_double),  intent(out) :: legf1(0:lmax1, 1:nlat, -mmax1:mmax1)
  real(c_double),  intent(out) :: legb1(1:nlat, 0:lmax1, -mmax1:mmax1)
  real(c_double),  intent(out) :: legf2(0:lmax2, 1:nlat, -mmax2:mmax2)
  real(c_double),  intent(out) :: legb2(1:nlat, 0:lmax2, -mmax2:mmax2)
  real(c_double), allocatable :: SHT_p1x(:,:)
  integer(c_long) :: ilat, iphi, l, m, lm
  !hack for SHTOOLS
  integer(4) :: lmax2_int4, i0, i1, i4
  !hack for SHTOOLS

!DEBUG
!  write(6, "('size(cost) = ', i5)") size(cost)
!DEBUG

  allocate(SHT_p1x(1:nlat, 1:((lmax2+1)*(lmax2+2))/2))

  i0 = 0
  i1 = 1
  i4 = 4
  lmax2_int4 = lmax2
  wphi = 2.D+0 * PI / nphi;
!  call SHGLQ(lmax2_int4, cost, wlat, SHT_p1x, i4, -i1, i0)
!  call SHGLQ(lmax2_int4, cost, wlat, SHT_p1x, i4, -i1, i1)
  call SHGLQ(lmax2_int4, cost, wlat, SHT_p1x, i4, i1, i1)

!debug
!  do l = 0, lmax1
!     do m = -l, l
!        lm = (l*(l+1))/2 + m + 1
!        do ilat = 1, nlat
!           write(6,"(4i5,2f20.10)") l,m,lm,ilat,cost(ilat),SHT_p1x(ilat,lm)
!        end do
!     end do
!  end do
!  stop
!debug

  do iphi = 1, nphi
     do ilat = 1, nlat
        wang(ilat, iphi) = wlat(ilat)
!       wang(ilat, iphi) = wlat(ilat) * wphi
     end do
  end do

  do ilat = 1, nlat
     sint(ilat) = sqrt(1-cost(ilat)**2.D+0)
  end do

  ! for positive magnetic quantum numbers
  do m = 0, mmax1
     do l = abs(m), lmax1
        lm = (l*(l+1))/2 + m + 1
        do ilat = 1, nlat
!           legf1(l, ilat, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI)
!           legb1(ilat, l, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * wlat(ilat)
           legf1(l, ilat, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * (-1)**m
           legb1(ilat, l, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * wlat(ilat) * (-1)**m
        end do
     end do
  end do
  do m = 0, mmax2
     do l = abs(m), lmax2
        lm = (l*(l+1))/2 + m + 1
        do ilat = 1, nlat
! 2016/9/5
!           legf2(l, ilat, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI)
!           legb2(ilat, l, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * wlat(ilat)
           legf2(l, ilat, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * (-1)**m
           legb2(ilat, l, m) = SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * wlat(ilat) * (-1)**m
!this doesn't work...
!           legf2(l, ilat, m) = SHT_p1x(ilat, lm) * (2.D+0 * PI)
!           legb2(ilat, l, m) = SHT_p1x(ilat, lm) * wlat(ilat)
!this doesn't work...
! 2016/9/5
        end do
     end do
  end do

  ! for negative magnetic quantum numbers
  do m = 1, mmax1
     do l = abs(m), lmax1
        do ilat = 1, nlat
           legf1(l, ilat, -m) = -legf1(l, ilat, m)
           legb1(ilat, l, -m) = -legb1(ilat, l, m)
!           legf1(l, ilat, -m) = legf1(l, ilat, m)
!           legb1(ilat, l, -m) = legb1(ilat, l, m)
        end do
     end do
  end do
  do m = 1, mmax2
     do l = abs(m), lmax2
        do ilat = 1, nlat
           legf2(l, ilat, -m) = -legf2(l, ilat, m)
           legb2(ilat, l, -m) = -legb2(ilat, l, m)
!           legf2(l, ilat, -m) = legf2(l, ilat, m)
!           legb2(ilat, l, -m) = legb2(ilat, l, m)
        end do
     end do
  end do

!DEBUG
!  write(6, "('sph_gen: nodes and weights')")
!  do ilat = 1, nlat
!     write(6, "(i10, 4f20.10)") ilat, acos(cost(ilat))/PI*180.D+0, cost(ilat), sint(ilat), wlat(ilat)
!  end do
!
!  write(6, "('sph_gen: Associated Legendre polynomials')")
!  do ilat = 1, nlat
!     write(6, "('ilat:', i10)", advance='no') ilat
!     lm = 0
!     do l = 0, lmax2
!        do m = 0, +l
!           lm = lm + 1
!           write(6, "(f20.10)", advance='no') SHT_p1x(ilat, lm)
!        end do
!     end do
!     write(6, *)
!  end do
!
!  write(6, "('sph_gen: forward Legendre transform')")
!! do m = 0, lmax2
!  do m = 0, mmax2
!     write(6, "('m =  ', i10)", advance = 'no') m
!     do ilat = 1, nlat
!        write(6, "(i20)", advance='no') ilat
!     end do
!     write(6, *)
!     do l = abs(m), lmax2
!        lm = (l*(l+1))/2 + m + 1
!        write(6, "('l:   ', i10)", advance='no') l
!        do ilat = 1, nlat
!           write(6, "(f20.10)", advance='no') SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI)
!        end do
!        write(6, *)
!     end do
!  end do
!
!  write(6, "('sph_gen: backward Legendre transform')")
!! do m = 0, lmax2
!  do m = 0, mmax2
!     write(6, "('m =  ', i10)", advance = 'no') m
!     do l = abs(m), lmax2
!        write(6, "(i20)", advance='no') l
!     end do
!     write(6, *)
!     do ilat = 1, nlat
!        write(6, "('ilat:', i10)", advance='no') ilat
!        do l = abs(m), lmax2
!           lm = (l*(l+1))/2 + m + 1
!           write(6, "(f20.10)", advance='no') SHT_p1x(ilat, lm) * sqrt(2.D+0 * PI) * wlat(ilat)
!        end do
!        write(6, *)
!     end do
!  end do
!DEBUG

  deallocate(SHT_p1x)

end subroutine sph_gen_x3j
!#######################################################################
