!################################################################################
subroutine ormas_denipx_ras(max_ipx,ovlp,cic,ipx,denipx)

  use,intrinsic :: iso_c_binding
  use mod_bas, only : mval
  use mod_const,only : zero,one,czero,ctwo
  use mod_ormas,only : thrdet,nelact,neltot,nfcore,ncore,nact,nfun,nstr_alph,nstr_beta,&
       & orb_alph,orb_beta,ndetx,ntot_alph_beta,nstr_alph_beta,llstr_alph_beta,&
       nstr_beta_alph,llstr_beta_alph

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: max_ipx
  complex(c_double_complex),intent(in) :: ovlp(1:nfun,1:nfun)
  complex(c_double_complex),intent(in) :: cic(1:ndetx)
  real(c_double),intent(out) :: ipx(0:max_ipx)
  complex(c_double),intent(out) :: denipx(1:nfun,1:nfun,0:max_ipx)
!  complex(c_double) :: denipx(1:nfun,1:nfun,0:max_ipx)
  !--------------------------------------------------------------------
  complex(c_double_complex) :: dets
  integer(c_int),external :: util_bicoeff
  complex(c_double_complex),external :: util_zdotup
  complex(c_double_complex),allocatable :: sa0(:,:) ! windowed ovlp for alpha strings
  complex(c_double_complex),allocatable :: sb0(:,:) ! windowed ovlp for beta strings
  complex(c_double_complex),allocatable :: sa1(:,:) ! windowed ovlp for alpha ionized strings
  complex(c_double_complex),allocatable :: sb1(:,:) ! windowed ovlp for beta ionized strings
  complex(c_double_complex),allocatable :: cdeta0(:,:,:) ! cdeta(j,k) = sum_i conjg(cic(i,j)) * det(sa(i,k))
  complex(c_double_complex),allocatable :: detbc0(:,:,:) ! detbc(j,k) = sum_i det(sa(j,i)) * cic(k,i) 
  complex(c_double_complex),allocatable :: cdeta1(:,:,:,:,:) ! cdeta(j,k) = sum_i conjg(cic(i,j)) * det(sa(i,k))
  complex(c_double_complex),allocatable :: detbc1(:,:,:,:,:) ! detbc(j,k) = sum_i det(sa(j,i)) * cic(k,i) 
  complex(c_double_complex),allocatable :: dmat(:,:,:)
  integer(c_int) :: ndet,iipx,jipx,iipx_a,iipx_b,istr,jstr,ii,ne2a,ne2b,iproc,nproc,ifun,jfun
  integer(c_int),external :: util_omp_nproc
  integer(c_int),external :: util_omp_iproc
  if (max_ipx < 0 .or. max_ipx > 4) then
     stop 'bad max_ipx in ormas_denipx_ras.'
  else if (max_ipx > neltot(3)) then
     stop 'max_ipx > neltot(3) in ormas_denipx_ras.'
  end if

  nproc = util_omp_nproc()
  ndet = nstr_alph * nstr_beta
  allocate(sa0(1:neltot(1)**2,0:(nproc-1)))
  allocate(sb0(1:neltot(2)**2,0:(nproc-1)))
  allocate(sa1(1:neltot(1)**2,0:(nproc-1)))
  allocate(sb1(1:neltot(2)**2,0:(nproc-1)))
  allocate(cdeta0(1:nstr_beta,1:nstr_alph,0:max_ipx))
  allocate(detbc0(1:nstr_beta,1:nstr_alph,0:max_ipx))
  allocate(cdeta1(1:nstr_beta,1:nstr_alph,1:nfun,1:nfun,0:max_ipx))
  allocate(detbc1(1:nstr_beta,1:nstr_alph,1:nfun,1:nfun,0:max_ipx))
  allocate(dmat(1:nfun,1:nfun,0:(nproc-1)))
  call util_zcopy(ndet*(max_ipx+1),czero,0,cdeta0,1)
  call util_zcopy(ndet*(max_ipx+1),czero,0,detbc0,1)
  call util_zcopy(nfun**2*ndet*(max_ipx+1),czero,0,cdeta1,1)
  call util_zcopy(nfun**2*ndet*(max_ipx+1),czero,0,detbc1,1)

  !$omp parallel default(shared) private(istr,jstr,dets,ii,iipx,iproc)
  iproc = util_omp_iproc()
  !$omp do
  do jstr = 1,nstr_alph
     do istr = 1,nstr_alph
        call ormas_sij(istr,jstr,ncore,neltot(1),nelact(1),nfun,orb_alph,ovlp,sa0(1,iproc))
        do iipx = 0,min(neltot(1),max_ipx)
           call ormas_str_dot_str_0(iipx,istr,jstr,ncore,neltot(1),nelact(1), &
                orb_alph,sa0(1,iproc),sa1(1,iproc),dets)
           call ormas_str_dot_str_1(iipx,istr,jstr,ncore,neltot(1),nelact(1), &
                orb_alph,sa0(1,iproc),sa1(1,iproc),dmat(1,1,iproc))
           if (abs(dets) > thrdet) then
              do ii = llstr_beta_alph(istr),llstr_beta_alph(istr)+nstr_beta_alph(istr)-1
                 cdeta0(ii,jstr,iipx) = cdeta0(ii,jstr,iipx) &
                      + conjg(cic(ntot_alph_beta(ii)+istr)) * dets
              end do
           end if
           do ifun = 1,nfun
              do jfun = 1,nfun
                 if (mval(ifun) .ne. mval(jfun)) cycle
                 if (abs(dmat(jfun,ifun,iproc)) < thrdet) cycle
                 do ii = llstr_beta_alph(istr),llstr_beta_alph(istr)+nstr_beta_alph(istr)-1
                    cdeta1(ii,jstr,jfun,ifun,iipx) = cdeta1(ii,jstr,jfun,ifun,iipx) &
                         + conjg(cic(ntot_alph_beta(ii)+istr)) * dmat(jfun,ifun,iproc)
                 end do
              end do
           end do
        end do
     end do
  end do
  !$omp end do

  !$omp do
  do istr = 1,nstr_beta
     do jstr = 1,nstr_beta
        call ormas_sij(istr,jstr,ncore,neltot(2),nelact(2),nfun,orb_beta,ovlp,sb0(1,iproc))
        do iipx = 0,min(neltot(2),max_ipx)
           call ormas_str_dot_str_0(iipx,istr,jstr,ncore,neltot(2),nelact(2),orb_beta, &
                sb0(1,iproc),sb1(1,iproc),dets)
           call ormas_str_dot_str_1(iipx,istr,jstr,ncore,neltot(2),nelact(2),orb_beta, &
                sb0(1,iproc),sb1(1,iproc),dmat(1,1,iproc))
           if (abs(dets) > thrdet) then
              do ii = llstr_alph_beta(jstr),llstr_alph_beta(jstr)+nstr_alph_beta(jstr)-1
                 detbc0(istr,ii,iipx) = detbc0(istr,ii,iipx) + dets * cic(ntot_alph_beta(jstr)+ii)
              end do
           end if
           do ifun = 1,nfun
              do jfun = 1,nfun
                 if (mval(ifun) .ne. mval(jfun)) cycle
                 if (abs(dmat(jfun,ifun,iproc)) < thrdet) cycle
                 do ii = llstr_alph_beta(jstr),llstr_alph_beta(jstr)+nstr_alph_beta(jstr)-1
                    detbc1(istr,ii,jfun,ifun,iipx) = detbc1(istr,ii,jfun,ifun,iipx) &
                         + dmat(jfun,ifun,iproc) * cic(ntot_alph_beta(jstr)+ii)
                 end do
              end do
           end do
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  denipx = zero
  do iipx = 0,max_ipx
     do iipx_a = 0,min(neltot(1),max_ipx)
        do iipx_b = 0,min(neltot(2),max_ipx)
           if (iipx_a + iipx_b .ne. iipx) cycle
           do ifun = 1,nfun
              do jfun = 1,nfun
                 if (mval(ifun) .ne. mval(jfun)) cycle
                 denipx(jfun,ifun,iipx) = denipx(jfun,ifun,iipx) &
                      + util_zdotup(ndet,cdeta1(1,1,ifun,jfun,iipx_a),detbc0(1,1,iipx_b)) &
                      + util_zdotup(ndet,cdeta0(1,1,iipx_a),detbc1(1,1,ifun,jfun,iipx_b))
              end do
           end do
        end do
     end do
  end do

  ! xi->T(0) formula: See note of 2018/7/17
  ipx(0) = zero
  do ifun = 1,nfun
     do jfun = 1,nfun
        if (mval(ifun) .ne. mval(jfun)) cycle
        ipx(0) = ipx(0) + ovlp(ifun,jfun)*denipx(jfun,ifun,0)
     end do
  end do  
  ipx(0) = ipx(0)/neltot(3)

  ! xi->T(iipx>0) formula: See note of 2018/7/17
  do iipx = 1,max_ipx
     ipx(iipx) = zero
     do ifun = 1, nfun
        ipx(iipx) = ipx(iipx) + denipx(ifun,ifun,iipx-1)
     end do
     ipx(iipx) = ipx(iipx)/iipx
  end do
!  write(6,"('ormas_denipx_ras. T:')")
!  do iipx = 0,max_ipx
!     write(6,"(2i5,f15.8)") max_ipx,iipx,ipx(iipx)
!  end do

  ! T-->P formula: See note of 2015/1/22
  do iipx = max_ipx,1,-1
     do jipx = 1,iipx
        ipx(iipx) = ipx(iipx) &
                & + ipx(iipx-jipx) * util_bicoeff(neltot(3)-iipx+jipx,jipx) * (-one)**jipx
     end do
  end do
!  write(6,"('ormas_denipx_ras P:')")
!  do iipx = 0,max_ipx
!     write(6,"(2i5,f15.8)") max_ipx,iipx,ipx(iipx)
!  end do

  ! xi(iipx)-->rho(iipx) formula: See note of 2018/7/17
  do iipx = max_ipx,1,-1
     do jipx = 1,iipx
        do ifun = 1,nfun
           do jfun = 1,nfun
              if (mval(ifun) .ne. mval(jfun)) cycle
              denipx(jfun,ifun,iipx) &
                   = denipx(jfun,ifun,iipx) &
                   + denipx(jfun,ifun,iipx-jipx)*util_bicoeff(neltot(3)-1-iipx+jipx,jipx)*(-one)**jipx
           end do
        end do
     end do
  end do

  deallocate(dmat)
  deallocate(detbc1)
  deallocate(cdeta1)
  deallocate(detbc0)
  deallocate(cdeta0)
  deallocate(sb1)
  deallocate(sa1)
  deallocate(sb0)
  deallocate(sa0)

end subroutine ormas_denipx_ras
!################################################################################
