!#######################################################################
subroutine hprod_cidiagx(cic)

  use, intrinsic :: iso_c_binding
  use mod_control, only : exact3j, PSP
  use mod_const, only : zero, czero, runit, ctrue
  use mod_ormas, only : neltot, ndcore, nact, nelact, lcic
  use mod_hprod, only : ene_fcore, ene_dcore, ene_core, ene_act
  use mod_hprod, only : int1e, int2e

  implicit none
  complex(c_double_complex), intent(inout) :: cic(1:lcic)
  complex(c_double_complex), external :: util_zdotc

  integer(c_int) :: idav, jdav, isub, jsub
  integer(c_int), parameter :: max_dav = 40
  real(c_double), parameter :: thr_dav = 1.D-15
!debug  integer(c_int), parameter :: max_dav = 20
!debug  real(c_double), parameter :: thr_dav = 1.D-8
  real(c_double) :: norm
  complex(c_double_complex) :: ovlp
  complex(c_double_complex), allocatable :: hcic(:), xcic(:,:), ycic(:,:)
  complex(c_double_complex), allocatable :: hmat(:,:), hsub(:,:), usub(:,:)

  allocate(hcic(1:lcic))
  allocate(xcic(1:lcic, 1:max_dav))
  allocate(ycic(1:lcic, 1:max_dav))
  allocate(hmat(1:max_dav, 1:max_dav))

  xcic(1:lcic,1) = cic(1:lcic)
  call wfn_madapt(xcic(1,1))
!  call wfn_ladaptc(xcic(1,1))
  norm = dble(util_zdotc(lcic,xcic(1,1),1,xcic(1,1),1))
  xcic(:,1) = xcic(:,1)/sqrt(abs(norm))

  do idav = 1, max_dav
     ycic(:,idav) = czero
     call ormas_hcic(int1e, int2e, xcic(1,idav), ycic(1,idav), zero)
     call wfn_madapt(ycic(1,idav))
!     call wfn_ladaptc(ycic(1,idav))
     do jdav = 1, idav
        hmat(jdav, idav) = util_zdotc(lcic,xcic(1,jdav),1,ycic(1,idav),1)
     end do

     allocate(hsub(1:idav, 1:idav))
     allocate(usub(1:idav, 1:idav))
     do isub = 1, idav
        do jsub = 1, isub - 1
           hsub(jsub, isub) = hmat(jsub, isub)
           hsub(isub, jsub) = conjg(hmat(jsub, isub))
        end do
        hsub(isub, isub) = hmat(isub, isub)
     end do
     do isub = 1, idav
        do jsub = 1, idav
           write(6,"(e12.4)",advance='no') dble(hsub(jsub,isub))
        end do
        write(6, *)
     end do
     call util_diag_comp(.false., idav, hsub, usub)
     ene_act = dble(hsub(1, 1))

     cic = czero
     hcic = czero
     do isub = 1, idav
        cic(:) = cic(:) + xcic(:,isub) * usub(isub, 1)
        hcic(:) = hcic(:) + ycic(:,isub) * usub(isub, 1)
     end do
     hcic = hcic - cic * ene_act
     deallocate(usub)
     deallocate(hsub)

     do isub = 1, idav
        ovlp = util_zdotc(lcic,xcic(1,isub),1,hcic,1)
        hcic(:) = hcic(:) - xcic(:,isub) * ovlp
     end do

     norm = dble(util_zdotc(lcic,hcic,1,hcic,1))
     write(6, "('hprod_cidiag: ', i5, 4f20.10)") idav, norm, ene_act, ene_core+ene_act, &
          dble(util_zdotc(lcic,cic,1,cic,1))
     if (idav == max_dav .or. abs(norm) < thr_dav) exit

     xcic(:,idav+1) = hcic(:) / sqrt(abs(norm))
     call wfn_madapt(xcic(1,idav+1))
!     call wfn_ladaptc(xcic(1,idav+1))
  end do

  write(6, "('hprod_cidiag: ene = ', 4f20.10)") ene_fcore, ene_dcore, ene_act, ene_core+ene_act

  deallocate(hmat)
  deallocate(ycic)
  deallocate(xcic)
  deallocate(hcic)
!  stop "STOP for debug @ hprod_cidiag."

end subroutine hprod_cidiagx
!#######################################################################
