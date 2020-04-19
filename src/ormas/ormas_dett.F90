!################################################################################
complex(c_double_complex) function ormas_dett(iipx, istr, jstr, ncore, nel, nela, orb, s0, work)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: iipx, istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-iipx), 1:(nel-iipx))
  !--------------------------------------------------------------------
  complex(c_double_complex), external :: ormas_dett1, ormas_dett2, ormas_dett3
  complex(c_double_complex), external :: ormas_dett2_old

  if (iipx > nel) then
     ormas_dett = czero
     return
  else if (iipx == nel) then
     if (istr == jstr) then
        ormas_dett = runit
     else
        ormas_dett = czero
     end if
     return
  else if (iipx == nel - 1) then
     ! use special routine ormas_dett_1.
     stop 'bad iipx in ormas_dett.'
  end if

! iipx < nel

  if (iipx <= 0 .or. iipx >= 4) then
     stop 'bad iipx in ormas_dett.'
  else if (iipx == 1) then
     ormas_dett = ormas_dett1(istr, jstr, ncore, nel, nela, orb, s0, work)
  else if (iipx == 2) then
     ormas_dett = ormas_dett2(istr, jstr, ncore, nel, nela, orb, s0, work)
!     ormas_dett = ormas_dett2_old(istr, jstr, ncore, nel, nela, orb, s0, work)
  else if (iipx == 3) then
     ormas_dett = ormas_dett3(istr, jstr, ncore, nel, nela, orb, s0, work)
! else if (iipx == 4) then
!    ormas_dett = ormas_dett4(istr, jstr, ncore, nel, nela, orb, s0, work)
  end if

  return

end function ormas_dett
!################################################################################
!################################################################################
complex(c_double_complex) function ormas_dett1(istr, jstr, ncore, nel, nela, orb, s0, work)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-1), 1:(nel-1))
  !--------------------------------------------------------------------
  integer(c_int) :: iel, jel, iela, jela
  complex(c_double_complex) :: sgnij
  complex(c_double_complex), external :: util_det

  ormas_dett1 = czero
  do iel = 1, ncore
     call util_window1(nel, iel, iel, s0, work)
     ormas_dett1 = ormas_dett1 + util_det(nel - 1, thrdet, work)
  end do

  do iel = ncore + 1, nel
     iela = iel - ncore
     do jel = ncore + 1, nel
        jela = jel - ncore
        if (orb(iela, istr) == orb(jela, jstr)) then
           sgnij = (-runit) ** (iela + jela)

           call util_window1(nel, iel, jel, s0, work)
           ormas_dett1 = ormas_dett1 + util_det(nel - 1, thrdet, work) * sgnij
           exit
        end if
     end do
  end do

end function ormas_dett1
!################################################################################
!################################################################################
complex(c_double_complex) function ormas_dett2_old(istr, jstr, ncore, nel, nela, orb, s0, work)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-2), 1:(nel-2))
  !--------------------------------------------------------------------
  integer(c_int) :: iel, jel, kel, lel, iela, jela, kela, lela
  complex(c_double_complex) :: sgnik, sgnjl, sgnijkl
  complex(c_double_complex), external :: util_det

  ormas_dett2_old = czero

  ! i < j <= ncore
  do jel = 1, ncore
     do iel = 1, jel - 1
        call util_window2(nel, iel, jel, iel, jel, s0, work)
        ormas_dett2_old = ormas_dett2_old + util_det(nel - 2, thrdet, work)
     end do
  end do

  do jel = ncore + 1, nel
     jela = jel - ncore
     do lel = ncore + 1, nel
        lela = lel - ncore
        if (orb(jela, istr) == orb(lela, jstr)) then
           sgnjl = (-runit) ** (jela + lela)

           ! i <= ncore < j
           do iel = 1, ncore
              call util_window2(nel, iel, jel, iel, lel, s0, work)
              ormas_dett2_old = ormas_dett2_old + util_det(nel - 2, thrdet, work) * sgnjl
           end do

           ! ncore < i < j
           do iel = ncore + 1, jel - 1
              iela = iel - ncore
              do kel = ncore + 1, lel - 1
                 kela = kel - ncore
                 if (orb(iela, istr) == orb(kela, jstr)) then
                    sgnik = (-runit) ** (iela + kela)
                    sgnijkl = sgnik * sgnjl
                    call util_window2(nel, iel, jel, kel, lel, s0, work)
                    ormas_dett2_old = ormas_dett2_old + util_det(nel - 2, thrdet, work) * sgnijkl
                    exit
                 end if
              end do
           end do
           exit
        end if
     end do
  end do

end function ormas_dett2_old
!################################################################################
complex(c_double_complex) function ormas_dett2(istr, jstr, ncore, nel, nela, orb, s0, work)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-2), 1:(nel-2))
  !--------------------------------------------------------------------
  integer(c_int) :: i1,i2,j1,j2
  complex(c_double_complex) :: sgn2, sgn21
  complex(c_double_complex), external :: util_det

  ormas_dett2 = czero

  ! i1 < i2 <= ncore
  do i2 = 1, ncore
     do i1 = 1, i2 - 1
        call util_window2(nel,i1,i2,i1,i2,s0,work)
        ormas_dett2 = ormas_dett2 + util_det(nel-2, thrdet, work)
     end do
  end do

  ! i1 <= ncore < i2
  do i2 = ncore+1, nel
  do j2 = ncore+1, nel
     if (orb(i2-ncore,istr) /= orb(j2-ncore,jstr)) cycle
     sgn2 = (-1)**(i2+j2)
     do i1 = 1, ncore
        call util_window2(nel,i1,i2,i1,j2,s0,work)
        ormas_dett2 = ormas_dett2 + util_det(nel-2, thrdet, work)*sgn2
     end do
  end do
  end do

  ! ncore < i1 < i2
  do i2 = ncore+1, nel
  do j2 = ncore+1, nel
     if (orb(i2-ncore,istr) /= orb(j2-ncore,jstr)) cycle
     sgn2 = (-1)**(i2+j2)
     do i1 = ncore+1, i2 - 1
     do j1 = ncore+1, j2 - 1
        if (orb(i1-ncore,istr) /= orb(j1-ncore,jstr)) cycle
        sgn21 = sgn2 * (-1)**(i1+j1)
        call util_window2(nel,i1,i2,j1,j2,s0,work)
        ormas_dett2 = ormas_dett2 + util_det(nel-2, thrdet, work)*sgn21
     end do
     end do
  end do
  end do

end function ormas_dett2
!################################################################################
complex(c_double_complex) function ormas_dett3(istr, jstr, ncore, nel, nela, orb, s0, work)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-3), 1:(nel-3))
  !--------------------------------------------------------------------
  integer(c_int) :: i1,i2,i3,j1,j2,j3
  complex(c_double_complex) :: sgn3, sgn32, sgn321
  complex(c_double_complex), external :: util_det

  ormas_dett3 = czero

  ! i1 < i2 < i3 <= ncore
  do i3 = 1, ncore
     do i2 = 1, i3 - 1
        do i1 = 1, i2 - 1
           call util_window3(nel,i1,i2,i3,i1,i2,i3,s0,work)
           ormas_dett3 = ormas_dett3 + util_det(nel-3, thrdet, work)
        end do
     end do
  end do

  ! i1 < i2 <= ncore < i3
  do i3 = ncore+1, nel
  do j3 = ncore+1, nel
     if (orb(i3-ncore,istr) /= orb(j3-ncore,jstr)) cycle
     sgn3 = (-1)**(i3+j3)
     do i2 = 1, ncore
        do i1 = 1, i2 - 1
           call util_window3(nel,i1,i2,i3,i1,i2,j3,s0,work)
           ormas_dett3 = ormas_dett3 + util_det(nel-3, thrdet, work)*sgn3
        end do
     end do
  end do
  end do

  ! i1 <= ncore < i2 < i3
  do i3 = ncore+1, nel
  do j3 = ncore+1, nel
     if (orb(i3-ncore,istr) /= orb(j3-ncore,jstr)) cycle
     sgn3 = (-1)**(i3+j3)
     do i2 = ncore+1, i3 - 1
     do j2 = ncore+1, j3 - 1
        if (orb(i2-ncore,istr) /= orb(j2-ncore,jstr)) cycle
        sgn32 = sgn3 * (-1)**(i2+j2)
        do i1 = 1, ncore
           call util_window3(nel,i1,i2,i3,i1,j2,j3,s0,work)
           ormas_dett3 = ormas_dett3 + util_det(nel-3, thrdet, work)*sgn32
        end do
     end do
     end do
  end do
  end do

  ! ncore < i1 < i2 < i3
  do i3 = ncore+1, nel
  do j3 = ncore+1, nel
     if (orb(i3-ncore,istr) /= orb(j3-ncore,jstr)) cycle
     sgn3 = (-1)**(i3+j3)
     do i2 = ncore+1, i3 - 1
     do j2 = ncore+1, j3 - 1
        if (orb(i2-ncore,istr) /= orb(j2-ncore,jstr)) cycle
        sgn32 = sgn3 * (-1)**(i2+j2)
        do i1 = ncore+1, i2 - 1
        do j1 = ncore+1, j2 - 1
           if (orb(i1-ncore,istr) /= orb(j1-ncore,jstr)) cycle
           sgn321 = sgn32 * (-1)**(i1+j1)
           call util_window3(nel,i1,i2,i3,j1,j2,j3,s0,work)
           ormas_dett3 = ormas_dett3 + util_det(nel-3, thrdet, work)*sgn321
        end do
        end do
     end do
     end do
  end do
  end do

end function ormas_dett3
!################################################################################
