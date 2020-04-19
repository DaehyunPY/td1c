!################################################################################
subroutine ormas_ipx_ras(max_ipx, ovlp, cic, ipx)

  use, intrinsic :: iso_c_binding
  use mod_const, only : zero, one, czero, ctwo
  use mod_ormas, only : thrdet, nelact, neltot, ncore, nact, nfun, nstr_alph, nstr_beta, &
       & orb_alph, orb_beta, n1x_alph, n1x_beta, p1x_alph, p1x_beta, h1x_alph, h1x_beta, &
       & eq1x_alph, eq1x_beta, sgn1x_alph, sgn1x_beta, ndetx, ntot_alph_beta, &
       & nstr_alph_beta, llstr_alph_beta, nstr_beta_alph, llstr_beta_alph

  implicit none
  !--------------------------------------------------------------------
  integer(c_long), intent(in) :: max_ipx
  complex(c_double_complex), intent(in) :: ovlp(1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: cic(1:ndetx)
  real(c_double), intent(out) :: ipx(0:*)
  !--------------------------------------------------------------------
  complex(c_double_complex) :: dets
  integer(c_long), external :: util_bicoeff
  complex(c_double_complex), external :: util_zdotu, util_zdotup
  complex(c_double_complex), external :: ormas_dets0, ormas_dets1, ormas_dett, ormas_dett_1
  complex(c_double_complex), allocatable :: sa0(:,:)     ! windowed ovlp for alpha strings
  complex(c_double_complex), allocatable :: sb0(:,:)     ! windowed ovlp for beta strings
  complex(c_double_complex), allocatable :: sa1(:,:)     ! windowed ovlp for alpha ionized strings
  complex(c_double_complex), allocatable :: sb1(:,:)     ! windowed ovlp for beta ionized strings
!toolarge  complex(c_double_complex), allocatable :: cdeta(:,:,:) ! cdeta(j, k) = sum_i conjg(cic(i, j)) * det(sa(i, k))
!toolarge  complex(c_double_complex), allocatable :: detbc(:,:,:) ! detbc(j, k) = sum_i det(sa(j, i)) * cic(k, i) 
!toolarge  complex(c_double_complex), allocatable :: ccic(:,:)    ! ccic(i, j) = conjg(cic(j, i))
  complex(c_double_complex), allocatable :: cdeta(:,:,:) ! cdeta(j, k) = sum_i conjg(cic(i, j)) * det(sa(i, k))
  complex(c_double_complex), allocatable :: detbc(:,:,:) ! detbc(j, k) = sum_i det(sa(j, i)) * cic(k, i) 

  integer(c_long) :: ndet, iipx, jipx, iipx_a, iipx_b, istr, jstr, ii, ne2a, ne2b, iproc, nproc
  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc

  if (max_ipx < 0 .or. max_ipx > 8) then
     stop 'bad max_ipx in ormas_ipx_ras.'
  else if (max_ipx > neltot(3)) then
     stop 'max_ipx > neltot(3) in ormas_ipx_ras.'
  end if

  nproc = util_omp_nproc()
  ndet = nstr_alph * nstr_beta
  ne2a = neltot(1) * neltot(1)
  ne2b = neltot(2) * neltot(2)
  allocate(sa0(1:ne2a, 0:(nproc-1)))
  allocate(sb0(1:ne2b, 0:(nproc-1)))
  allocate(sa1(1:ne2a, 0:(nproc-1)))
  allocate(sb1(1:ne2b, 0:(nproc-1)))
!toolarge  allocate(cdeta(1:nstr_beta, 1:nstr_alph, 0:max_ipx))
!toolarge  allocate(detbc(1:nstr_beta, 1:nstr_alph, 0:max_ipx))
!toolarge  allocate(ccic(1:nstr_beta, 1:nstr_alph))
  allocate(cdeta(1:nstr_beta, 1:nstr_alph, 0:max_ipx))
  allocate(detbc(1:nstr_beta, 1:nstr_alph, 0:max_ipx))

  call util_zcopy(ndet*(max_ipx+1), czero, 0, cdeta, 1)
  call util_zcopy(ndet*(max_ipx+1), czero, 0, detbc, 1)

!toolarge  ! transposed cic for better vectorization
!toolarge  ccic = czero
!toolarge  do jstr = 1, nstr_beta
!toolarge     do istr = llstr_alph_beta(jstr), llstr_alph_beta(jstr)+nstr_alph_beta(jstr)-1
!toolarge        ccic(jstr,istr) = conjg(cic(ntot_alph_beta(jstr)+istr))
!toolarge     end do
!toolarge  end do

  !$omp parallel default(shared) private(istr,jstr,dets,ii,iipx,iproc)
  iproc = util_omp_iproc()
  !$omp do
  do jstr = 1, nstr_alph
     do istr = 1, nstr_alph
        ! a(0)
        dets = ormas_dets0(istr, jstr, ncore, neltot(1), nelact(1), nfun, orb_alph, ovlp, sa0(1,iproc))
        if (abs(dets) > thrdet) then
!          do ii = 1, nstr_beta
           do ii = llstr_beta_alph(istr), llstr_beta_alph(istr)+nstr_alph_beta(istr)-1
!toolarge     cdeta(ii, jstr, 0) = cdeta(ii, jstr, 0) + ccic(ii,istr) * dets
              cdeta(ii, jstr, 0) = cdeta(ii, jstr, 0) + conjg(cic(ntot_alph_beta(ii)+istr)) * dets
           end do
        end if

        ! a(i > 0)
        do iipx = 1, min(neltot(1), max_ipx)
           if (iipx == neltot(1) - 1) then
              dets = ormas_dett_1(istr, jstr, ncore, nact, nfun, n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph, ovlp)
           else
              dets = ormas_dett(iipx, istr, jstr, ncore, neltot(1), nelact(1), orb_alph, sa0(1,iproc), sa1(1,iproc))
           end if
           if (abs(dets) > thrdet) then
!             do ii = 1, nstr_beta
              do ii = llstr_beta_alph(istr), llstr_beta_alph(istr)+nstr_alph_beta(istr)-1
!toolarge        cdeta(ii, jstr, iipx) = cdeta(ii, jstr, iipx) + ccic(ii,istr) * dets
                 cdeta(ii, jstr, iipx) = cdeta(ii, jstr, iipx) + conjg(cic(ntot_alph_beta(ii)+istr)) * dets
              end do
           end if
        end do
     end do
  end do
  !$omp end do

  !$omp do
  do istr = 1, nstr_beta
     do jstr = 1, nstr_beta
        ! b(0)
        dets = ormas_dets0(istr, jstr, ncore, neltot(2), nelact(2), nfun, orb_beta, ovlp, sb0(1,iproc))
        if (abs(dets) > thrdet) then
!          do ii = 1, nstr_alph
           do ii = llstr_alph_beta(jstr), llstr_alph_beta(jstr)+nstr_alph_beta(jstr)-1
!              detbc(istr, ii, 0) = detbc(istr, ii, 0) + dets * cic(mapf_detx(ii,jstr))
              detbc(istr, ii, 0) = detbc(istr, ii, 0) + dets * cic(ntot_alph_beta(jstr)+ii)
           end do
        end if

        ! b(i > 0)
        do iipx = 1, min(neltot(2), max_ipx)
           if (iipx == neltot(2) - 1) then
              dets = ormas_dett_1(istr, jstr, ncore, nact, nfun, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta, ovlp)
           else
              dets = ormas_dett(iipx, istr, jstr, ncore, neltot(2), nelact(2), orb_beta, sb0(1,iproc), sb1(1,iproc))
           end if
           if (abs(dets) > thrdet) then
!             do ii = 1, nstr_alph
              do ii = llstr_alph_beta(jstr), llstr_alph_beta(jstr)+nstr_alph_beta(jstr)-1
!                 detbc(istr, ii, iipx) = detbc(istr, ii, iipx) + dets * cic(mapf_detx(ii,jstr))
                 detbc(istr, ii, iipx) = detbc(istr, ii, iipx) + dets * cic(ntot_alph_beta(jstr)+ii)
              end do
           end if
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  ! current limitation:
  ! neltot(1) <= 4
  ! neltot(2) <= 4
  ipx(0:max_ipx) = zero
  do iipx = 0, max_ipx
     do iipx_a = 0, min(4, min(neltot(1), max_ipx))
        do iipx_b = 0, min(4, min(neltot(2), max_ipx))
           if (iipx_a + iipx_b == iipx) then
!              ipx(iipx) = ipx(iipx) + dble(util_zdotu(ndet, cdeta(1,1,iipx_a), 1, detbc(1,1,iipx_b), 1))
              ipx(iipx) = ipx(iipx) + dble(util_zdotup(ndet, cdeta(1,1,iipx_a), detbc(1,1,iipx_b)))
           end if
        end do
     end do
  end do

  !debug
  ! write(6, "('ormas_ipx_ras. T:')", advance = 'no')
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
  ! write(6, "('ormas_ipx_ras. P:')", advance = 'no')
  ! do iipx = 0, max_ipx
  !    write(6, "(f15.8)", advance = 'no') ipx(iipx)
  ! end do
  ! write(6, *)
  !debug

!toolarge  deallocate(ccic)
  deallocate(detbc)
  deallocate(cdeta)
  deallocate(sb1)
  deallocate(sa1)
  deallocate(sb0)
  deallocate(sa0)

end subroutine ormas_ipx_ras
!################################################################################
