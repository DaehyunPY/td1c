!################################################################################
subroutine ormas_ipx_core(max_ipx, ovlp, ipx)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, one, czero, ctwo
  use mod_ormas, only : thrdet, nelact, neltot, ncore, nact, nfun

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: max_ipx
  complex(c_double_complex), intent(in) :: ovlp(1:nfun, 1:nfun)
  real(c_double), intent(out) :: ipx(0:*)
  !--------------------------------------------------------------------
  complex(c_double_complex) :: dets
  integer(c_int), external :: util_bicoeff
  complex(c_double_complex), external :: util_zdotu, util_zdotup
  complex(c_double_complex), external :: ormas_dets0, ormas_dets1, ormas_dett, ormas_dett_1
  complex(c_double_complex), allocatable :: sa0(:)   ! windowed ovlp for alpha strings
  complex(c_double_complex), allocatable :: sb0(:)   ! windowed ovlp for beta strings
  complex(c_double_complex), allocatable :: sa1(:)   ! windowed ovlp for alpha ionized strings
  complex(c_double_complex), allocatable :: sb1(:)   ! windowed ovlp for beta ionized strings
  complex(c_double_complex), allocatable :: cdeta(:) ! cdeta(j, k) = sum_i conjg(cic(i, j)) * det(sa(i, k))
  complex(c_double_complex), allocatable :: detbc(:) ! detbc(j, k) = sum_i det(sa(j, i)) * cic(k, i) 

  integer(c_int) :: iipx, jipx, iipx_a, iipx_b, ne2a, ne2b
  integer(c_int) :: dum(1)

  if (max_ipx < 0 .or. max_ipx > 8) then
     stop 'bad max_ipx in ormas_ipx_core.'
  else if (max_ipx > neltot(3)) then
     stop 'max_ipx > neltot(3) in ormas_ipx_core.'
  end if

  ne2a = neltot(1) * neltot(1)
  ne2b = neltot(2) * neltot(2)
  allocate(sa0(1:ne2a))
  allocate(sb0(1:ne2b))
  allocate(sa1(1:ne2a))
  allocate(sb1(1:ne2b))
  allocate(cdeta(0:max_ipx))
  allocate(detbc(0:max_ipx))

  cdeta = czero
  detbc = czero

  ! a(0)
  dets = ormas_dets0(1, 1, ncore, neltot(1), nelact(1), nfun, dum, ovlp, sa0)
  if (abs(dets) > thrdet) then
     cdeta(0) = cdeta(0) + dets
  end if
  ! a(i > 0)
  do iipx = 1, min(neltot(1), max_ipx)
     if (iipx == neltot(1) - 1) then
        dets = ormas_dett_1(1, 1, ncore, nact, nfun, dum, dum, dum, dum, dum, ovlp)
     else
        dets = ormas_dett(iipx, 1, 1, ncore, neltot(1), nelact(1), dum, sa0, sa1)
     end if
     if (abs(dets) > thrdet) then
        cdeta(iipx) = cdeta(iipx) + dets
     end if
  end do

  ! b(0)
  dets = ormas_dets0(1, 1, ncore, neltot(2), nelact(2), nfun, dum, ovlp, sb0)
  if (abs(dets) > thrdet) then
     detbc(0) = detbc(0) + dets
  end if
  ! b(i > 0)
  do iipx = 1, min(neltot(2), max_ipx)
     if (iipx == neltot(2) - 1) then
        dets = ormas_dett_1(1, 1, ncore, nact, nfun, dum, dum, dum, dum, dum, ovlp)
     else
        dets = ormas_dett(iipx, 1, 1, ncore, neltot(2), nelact(2), dum, sb0, sb1)
     end if
     if (abs(dets) > thrdet) then
        detbc(iipx) = detbc(iipx) + dets
     end if
  end do

  ! current limitation:
  ! neltot(1) <= 4
  ! neltot(2) <= 4
  ipx(0:max_ipx) = zero
  do iipx = 0, max_ipx
     do iipx_a = 0, min(4, min(neltot(1), max_ipx))
        do iipx_b = 0, min(4, min(neltot(2), max_ipx))
           if (iipx_a + iipx_b == iipx) then
              ipx(iipx) = ipx(iipx) + dble(util_zdotup(1, cdeta(iipx_a), detbc(iipx_b)))
           end if
        end do
     end do
  end do

  !debug
  ! write(6, "('ormas_ipx_core. T:')", advance = 'no')
  ! do iipx = 0, max_ipx
  !    write(6, "(f15.8)", advance = 'no') ipx(iipx)
  ! end do
  ! write(6, *)
  !debug

  do iipx = max_ipx, 1, -1
     do jipx = 1, iipx
        ipx(iipx) = ipx(iipx) &
                & + ipx(iipx - jipx) * util_bicoeff(neltot(3) - iipx + jipx, jipx) * (-one) ** jipx
     end do
  end do

  !debug
  ! write(6, "('ormas_ipx_core. P:')", advance = 'no')
  ! do iipx = 0, max_ipx
  !    write(6, "(f15.8)", advance = 'no') ipx(iipx)
  ! end do
  ! write(6, *)
  !debug

  deallocate(detbc)
  deallocate(cdeta)
  deallocate(sb1)
  deallocate(sa1)
  deallocate(sb0)
  deallocate(sa0)

end subroutine ormas_ipx_core
!################################################################################
