!################################################################################
complex(c_double_complex) function ormas_dett(iipx, istr, jstr, ncore, nel, nela, orb, s0, work)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_long), intent(in) :: iipx, istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-iipx), 1:(nel-iipx))
  !--------------------------------------------------------------------
  complex(c_double_complex), external :: ormas_dett1, ormas_dett2

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
     stop 'bad iipx in ormas_dett.'
  end if

! iipx < nel

  if (iipx <= 0 .or. iipx >= 3) then
     stop 'bad iipx in ormas_dett.'
  else if (iipx == 1) then
     ormas_dett = ormas_dett1(istr, jstr, ncore, nel, nela, orb, s0, work)
  else if (iipx == 2) then
     ormas_dett = ormas_dett2(istr, jstr, ncore, nel, nela, orb, s0, work)
!nyi  else if (iipx == 3) then
!nyi     ormas_dett = ormas_dett3(istr, jstr, ncore, nel, nela, orb, s0, work)
!nyi  else if (iipx == 4) then
!nyi     ormas_dett = ormas_dett4(istr, jstr, ncore, nel, nela, orb, s0, work)
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
  integer(c_long), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-1), 1:(nel-1))
  !--------------------------------------------------------------------
  integer(c_long) :: iel, jel, iela, jela
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
complex(c_double_complex) function ormas_dett2(istr, jstr, ncore, nel, nela, orb, s0, work)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet
  use mod_const, only : czero, runit

  implicit none
  !--------------------------------------------------------------------
  integer(c_long), intent(in) :: istr, jstr, ncore, nel, nela, orb(0:nela, *)
  complex(c_double_complex), intent(in) :: s0(1:nel, 1:nel)
  complex(c_double_complex), intent(out) :: work(1:(nel-2), 1:(nel-2))
  !--------------------------------------------------------------------
  integer(c_long) :: iel, jel, kel, lel, iela, jela, kela, lela
  complex(c_double_complex) :: sgnik, sgnjl, sgnijkl
  complex(c_double_complex), external :: util_det

  ormas_dett2 = czero

  ! i < j <= ncore
  do jel = 1, ncore
     do iel = 1, jel - 1
        call util_window2(nel, iel, jel, iel, jel, s0, work)
        ormas_dett2 = ormas_dett2 + util_det(nel - 2, thrdet, work)
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
              ormas_dett2 = ormas_dett2 + util_det(nel - 2, thrdet, work) * sgnjl
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
                    ormas_dett2 = ormas_dett2 + util_det(nel - 2, thrdet, work) * sgnijkl
                    exit
                 end if
              end do
           end do
           exit
        end if
     end do
  end do

end function ormas_dett2
!################################################################################
