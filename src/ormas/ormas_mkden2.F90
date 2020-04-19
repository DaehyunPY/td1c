!######################################################################
subroutine ormas_mkden2(cic, den1, den2)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_ormas, only : den2_type,nact,nelact,fab_den2,cic_old,tdcc

  implicit none
  complex(c_double_complex) , intent(in) :: cic(1:*)
  complex(c_double_complex) , intent(in) :: den1(1:*)
  complex(c_double_complex) , intent(out) :: den2(1:*)

  if (nact == 0) return

  if (tdcc) then
     call tdcc_mkden2(cic, den2)
  else if (den2_type == 0) then
     if (.not. cic_old) then
        call ormas_mkden2_ras(cic, den1, den2)
     else
        call ormas_mkden2_old(cic, den1, den2)
     end if
!  else if (den2_type == 1) then
!     call ormas_mkden2_v1(cic, den1, den2)
!  else if (den2_type == 2) then ! previous default
!     call ormas_mkden2_v2(cic, den1, den2)
!  else if (den2_type == 3) then
!     call ormas_mkden2_v3(cic, den1, den2)
!  else if (den2_type == 4) then
!     if (fab_den2 .and. smul == 1 .and. nelact(1) == nelact(2)) then
!        call ormas_mkden2_v4(cic, den1, den2) ! new version by fabian
!     else
!        call ormas_mkden2_v2(cic, den1, den2)
!     end if
  else
     write(6, "('bad den2_type.')")
     stop
  end if

end subroutine ormas_mkden2
!######################################################################
subroutine ormas_mkden2_full(cic, den2)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : ncore,nact,nfun

  implicit none
  complex(c_double_complex), intent(in) :: cic(1:*)
  complex(c_double_complex), intent(out) :: den2(1:nfun,1:nfun,1:nfun,1:nfun)
  complex(c_double_complex), allocatable :: den1a(:,:),den2a(:,:,:,:)
  integer(c_int) :: i,j,k,l,t,u,v,w
  complex(c_double_complex) :: d1

  allocate(den1a(1:nact,1:nact))
  allocate(den2a(1:nact,1:nact,1:nact,1:nact))

  call ormas_mkden1(cic,den1a)
  call ormas_mkden2(cic,den1a,den2a)

  den2 = 0d0
  do i = 1, ncore
     do j = 1, ncore
        den2(i,i,j,j) = den2(i,i,j,j) + 4d0
        den2(i,j,j,i) = den2(i,j,j,i) - 2d0
     end do
  end do

  do t = ncore+1, nfun
     do u = ncore+1, nfun
        do i = 1, ncore
           d1 = den1a(t-ncore,u-ncore)
           den2(t,u,i,i) = den2(t,u,i,i) + d1*2d0
           den2(i,i,t,u) = den2(i,i,t,u) + d1*2d0
           den2(t,i,i,u) = den2(t,i,i,u) - d1
           den2(i,u,t,i) = den2(i,u,t,i) - d1
        end do
     end do
  end do

  do t = ncore+1, nfun
     do u = ncore+1, nfun
        do v = ncore+1, nfun
           do w = ncore+1, nfun
              den2(t,u,v,w) = den2(t,u,v,w) &
                   + den2a(t-ncore,u-ncore,v-ncore,w-ncore)
           end do
        end do
     end do
  end do

  deallocate(den2a)
  deallocate(den1a)

end subroutine ormas_mkden2_full
!######################################################################
subroutine ormas_print_den2(den2)

  use, intrinsic :: iso_c_binding
  use mod_bas, only : smul
  use mod_ormas, only : den2_type, nact, nelact, fab_den2, cic_old

  implicit none
  complex(c_double_complex) , intent(in) :: den2(1:nact,1:nact,1:nact,1:nact)
  integer(c_int) :: i,j,k,l

  if (nact == 0) return
  do i = 1, nact
  do j = 1, nact
  do k = 1, nact
  do l = 1, nact
     write(6,"(4i5,2f10.5)") i,j,k,l,den2(i,j,k,l)
  end do
  end do
  end do
  end do

end subroutine ormas_print_den2
!######################################################################
