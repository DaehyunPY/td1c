!################################################################################
subroutine ormas_ipd(max_ipd, ovlp, cic, ipd)

  use, intrinsic :: iso_c_binding
  use mod_ormas, only : thrdet, nelact, neltot, ncore, nact, nfun, nstr_alph, nstr_beta

  implicit none
  !--------------------------------------------------------------------
  integer(c_int), intent(in) :: max_ipd
  complex(c_double_complex), intent(in) :: ovlp(1:nfun, 1:nfun)
  complex(c_double_complex), intent(in) :: cic(1:*)
  real(c_double), intent(out) :: ipd(0:*)
  !--------------------------------------------------------------------
  integer(c_int) :: iact, jact, kact, lact, ifun, jfun, kfun, lfun
  real(kind(0d0)) :: p0, p1, p2
  complex(kind(0d0)) :: s0, s1, s2, s1c, s1a, s2cc, s2ca, s2aa, tmp
  complex(kind(0d0)), allocatable :: ovlx(:,:)      ! outer region overlap
  complex(kind(0d0)), allocatable :: den1(:,:)      ! 1RDM
  complex(kind(0d0)), allocatable :: den2(:,:,:,:)  ! 2RDM

  if (max_ipd < 0 .or. max_ipd > 2) then
     stop 'bad max_ipd in ormas_ipd.'
  else if (max_ipd > neltot(3)) then
     stop 'max_ipd > neltot(3) in ormas_ipd.'
  end if

  allocate(ovlx(1:nfun, 1:nfun))
  allocate(den1(1:nact, 1:nact))
  allocate(den2(1:nact, 1:nact, 1:nact, 1:nact))

  call ormas_mkden1(cic, den1)
  call ormas_mkden2(cic, den1, den2)  

  do ifun = 1, nfun
     do jfun = 1, ifun - 1
        ovlx(jfun, ifun) = -ovlp(jfun, ifun)
        ovlx(ifun, jfun) =  conjg(ovlx(jfun, ifun))
     end do
     ovlx(ifun, ifun) = 1d0 - ovlp(ifun, ifun)
  end do

  ! S0: unity by definition
  s0 = 1d0

  ! S1: core contribution
  s1c = 0d0
  do ifun = 1, ncore
     s1c = s1c + 2d0 * ovlx(ifun, ifun)
  end do

  ! S1: active contribution
  s1a = 0d0
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        s1a = s1a + den1(jact, iact) * ovlx(ifun, jfun)
     end do
  end do
  s1 = s1c + s1a

  ! S2: core-core contribution
  s2cc = 0d0
  do ifun = 1, ncore
     do jfun = 1, ncore
        s2cc = s2cc + 2d0 * ovlx(ifun, ifun) * ovlx(jfun, jfun) &
                         & - ovlx(ifun, jfun) * ovlx(jfun, ifun)
     end do
  end do
!2013.03.05  s2cc = 2d0 * s2cc

  ! S2: core-active contribution
  s2ca = 0d0
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact

        tmp = 0d0
        do kfun = 1, ncore
           tmp = tmp + 2d0 * ovlx(ifun, jfun) * ovlx(kfun, kfun) &
                          & - ovlx(ifun, kfun) * ovlx(kfun, jfun)
        end do
        
        s2ca = s2ca + den1(jact, iact) * tmp
     end do
  end do

  ! S2: active-active contribution
!2013.03.05:  s2ca = 0d0
  s2aa = 0d0
  do iact = 1, nact
     ifun = ncore + iact
     do jact = 1, nact
        jfun = ncore + jact
        do kact = 1, nact
           kfun = ncore + kact
           do lact = 1, nact
              lfun = ncore + lact
!bug              s2aa = s2aa + den2(kact, lact, iact, jact) &
!bug                   &  * ovlx(ifun, kfun) * ovlx(jfun, lfun)
              s2aa = s2aa + den2(kact, iact, lact, jact) &
                   &  * ovlx(ifun, kfun) * ovlx(jfun, lfun)
           end do
        end do
     end do
  end do
  s2aa = s2aa * 0.5d0
  s2 = s2cc + s2ca + s2aa

  p0 = dble(s0 - s1 + s2)
  p1 = dble(s1 - s2 * 2d0)
  p2 = dble(s2)

  ipd(0) = p0
  ipd(1) = p1
  ipd(2) = p2

  deallocate(den2)
  deallocate(den1)
  deallocate(ovlx)
  
end subroutine ormas_ipd
!######################################################################
