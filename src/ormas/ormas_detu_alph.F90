!################################################################################
subroutine ormas_detu(iipx, istr, jstr, ncore, nel, nela, orb, s0, work, dmat)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet, nfun
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: iipx, istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-iipx), 1:(nel-iipx))
  complex(c_double_complex), intent(out) :: dmat(1:nfun, 1:nfun)
  !--------------------------------------------------------------------

  if (iipx == nel-1) then
     call ormas_detu_null(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)
  else if (iipx == 0) then
     call ormas_detu0(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)
  else if (iipx == 1) then
     call ormas_detu1(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)
  else if (iipx == 2) then
     call ormas_detu2(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)
  else if (iipx == 3) then
     call ormas_detu3(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)
  else
     stop 'ormas_detu: bad iipx.'
  end if

  return

end subroutine ormas_detu
!################################################################################
subroutine ormas_detu0(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet, nfun
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-1), 1:(nel-1))
  complex(c_double_complex), intent(out) :: dmat(1:nfun, 1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,j0,sgn0
  complex(c_double_complex), external :: util_det

  dmat = czero
  do i0 = 1, nel
  do j0 = 1, nel
     if (i0 <= ncore) then
        ifun = i0
     else
        ifun = ncore + orb(i0-ncore, istr)
     end if
     if (j0 <= ncore) then
        jfun = j0
     else
        jfun = ncore + orb(j0-ncore, jstr)
     end if
     sgn0 = (-1)**(i0+j0)
     call util_window1(nel, i0, j0, s0, work)
     dmat(ifun, jfun) = dmat(ifun, jfun) + util_det(nel-1, thrdet, work) * sgn0
  end do
  end do

end subroutine ormas_detu0
!################################################################################
subroutine ormas_detu1(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet, nfun
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-2), 1:(nel-2))
  complex(c_double_complex), intent(out) :: dmat(1:nfun, 1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,i1,j0,j1,i0x,i1x,j0x,j1x,sgn0,sgn01
  complex(c_double_complex), external :: util_det

  dmat = czero
  do i0 = 1, nel
  do j0 = 1, nel
     if (i0 <= ncore) then
        ifun = i0
     else
        ifun = ncore + orb(i0-ncore, istr)
     end if
     if (j0 <= ncore) then
        jfun = j0
     else
        jfun = ncore + orb(j0-ncore, jstr)
     end if
     sgn0 = (-1)**(i0+j0)

     do i1 = 1, nel
     do j1 = 1, nel
        if (i1 == i0 .or. j1 == j0) cycle
        if (i1 <= ncore .and. i1.ne.j1) cycle
        if (i1 > ncore .and. j1 <= ncore) cycle
        if (i1 > ncore .and. orb(i1-ncore,istr) /= orb(j1-ncore,jstr)) cycle
        sgn01 = sgn0 * (-1)**(i1+j1)
        i0x = min(i0,i1)
        i1x = max(i0,i1)
        j0x = min(j0,j1)
        j1x = max(j0,j1)
        call util_window2(nel, i0x, i1x, j0x, j1x, s0, work)
        dmat(ifun, jfun) = dmat(ifun, jfun) + util_det(nel-2, thrdet, work) * sgn01
     end do
     end do
  end do
  end do

end subroutine ormas_detu1
!################################################################################
subroutine ormas_detu2(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet, nfun
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-3), 1:(nel-3))
  complex(c_double_complex), intent(out) :: dmat(1:nfun, 1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,i1,i2,j0,j1,j2,i0x,i1x,i2x,j0x,j1x,j2x, &
       sgn0,sgn01,sgn012
  complex(c_double_complex), external :: util_det

  dmat = czero
  do i0 = 1, nel
  do j0 = 1, nel
     if (i0 <= ncore) then
        ifun = i0
     else
        ifun = ncore + orb(i0-ncore, istr)
     end if
     if (j0 <= ncore) then
        jfun = j0
     else
        jfun = ncore + orb(j0-ncore, jstr)
     end if
     sgn0 = (-1)**(i0+j0)

     do i1 = 1, nel
     do j1 = 1, nel
        if (i1 == i0 .or. j1 == j0) cycle
        if (i1 <= ncore .and. i1.ne.j1) cycle
        if (i1 > ncore .and. j1 <= ncore) cycle
        if (i1 > ncore .and. orb(i1-ncore,istr) /= orb(j1-ncore,jstr)) cycle
        sgn01 = sgn0 * (-1)**(i1+j1)
        do i2 = i1 + 1, nel
        do j2 = j1 + 1, nel
           if (i2 == i0 .or. j2 == j0) cycle
           if (i2 <= ncore .and. i2.ne.j2) cycle
           if (i2 > ncore .and. j2 <= ncore) cycle
           if (i2 > ncore .and. orb(i2-ncore,istr) /= orb(j2-ncore,jstr)) cycle
           sgn012 = sgn01 * (-1)**(i2+j2)
           i0x = min(i0,i1)
           i1x = min(max(i0,i1),min(i0,i2))
           i2x = max(i0,i2)
           j0x = min(j0,j1)
           j1x = min(max(j0,j1),min(j0,j2))
           j2x = max(j0,j2)
           call util_window3(nel, i0x, i1x, i2x, j0x, j1x, j2x, s0, work)
           dmat(ifun, jfun) = dmat(ifun, jfun) + util_det(nel-3, thrdet, work) * sgn012
        end do
        end do
     end do
     end do
  end do
  end do

end subroutine ormas_detu2
!################################################################################
subroutine ormas_detu3(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet, nfun
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-4), 1:(nel-4))
  complex(c_double_complex), intent(out) :: dmat(1:nfun, 1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,i1,i2,i3,j0,j1,j2,j3,i0x,i1x,i2x,i3x, &
       j0x,j1x,j2x,j3x,sgn0, sgn01, sgn012, sgn0123
  complex(c_double_complex), external :: util_det

  dmat = czero
  do i0 = 1, nel
  do j0 = 1, nel
     if (i0 <= ncore) then
        ifun = i0
     else
        ifun = ncore + orb(i0-ncore, istr)
     end if
     if (j0 <= ncore) then
        jfun = j0
     else
        jfun = ncore + orb(j0-ncore, jstr)
     end if
     sgn0 = (-1)**(i0+j0)

     do i1 = 1, nel
     do j1 = 1, nel
        if (i1 == i0 .or. j1 == j0) cycle
        if (i1 <= ncore .and. i1.ne.j1) cycle
        if (i1 > ncore .and. j1 <= ncore) cycle
        if (i1 > ncore .and. orb(i1-ncore,istr) /= orb(j1-ncore,jstr)) cycle
        sgn01 = sgn0 * (-1)**(i1+j1)
        do i2 = i1 + 1, nel
        do j2 = j1 + 1, nel
           if (i2 == i0 .or. j2 == j0) cycle
           if (i2 <= ncore .and. i2.ne.j2) cycle
           if (i2 > ncore .and. j2 <= ncore) cycle
           if (i2 > ncore .and. orb(i2-ncore,istr) /= orb(j2-ncore,jstr)) cycle
           sgn012 = sgn01 * (-1)**(i2+j2)
           do i3 = i2 + 1, nel
           do j3 = j2 + 1, nel
              if (i3 == i0 .or. j3 == j0) cycle
              if (i3 <= ncore .and. i3.ne.j3) cycle
              if (i3 > ncore .and. j3 <= ncore) cycle
              if (i3 > ncore .and. orb(i3-ncore,istr) /= orb(j3-ncore,jstr)) cycle
              sgn0123 = sgn012 * (-1)**(i3+j3)

              i0x = min(i0,i1)
              i1x = min(max(i0,i1),min(i0,i2))
              i2x = max(i0,i2)
              j0x = min(j0,j1)
              j1x = min(max(j0,j1),min(j0,j2))
              j2x = max(j0,j2)
              call util_window4(nel,i0x,i1x,i2x,i3x,j0x,j1x,j2x,j3x,s0,work)
              dmat(ifun,jfun) = dmat(ifun,jfun) + util_det(nel-4,thrdet,work) * sgn0123
           end do
           end do
        end do
        end do
     end do
     end do
  end do
  end do

end subroutine ormas_detu3
!################################################################################
subroutine ormas_detu_null(istr, jstr, ncore, nel, nela, orb, s0, work, dmat)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet, nfun
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-1), 1:(nel-1))
  complex(c_double_complex), intent(out) :: dmat(1:nfun, 1:nfun)
  !--------------------------------------------------------------------
  integer(c_int) :: ifun,jfun,i0,j0,sgn0
  complex(c_double_complex), external :: util_det

  dmat = czero
  do i0 = 1, ncore
     dmat(i0,i0) = runit
  end do

  do j0 = 1, nel
     if (i0 <= ncore) then
        ifun = i0
     else
        ifun = ncore + orb(i0-ncore, istr)
     end if
     if (j0 <= ncore) then
        jfun = j0
     else
        jfun = ncore + orb(j0-ncore, jstr)
     end if
     sgn0 = (-1)**(i0+j0)
     call util_window1(nel, i0, j0, s0, work)
     dmat(ifun, jfun) = dmat(ifun, jfun) + util_det(nel-1, thrdet, work) * sgn0
  end do
  end do

end subroutine ormas_detu0
!################################################################################
