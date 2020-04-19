!######################################################################
subroutine ormas_mkdden1_old(fac, cic, dcic, dden1)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero, runit, ctwo
  use mod_ormas, only : nelact, ncore, nact, nsub, lorb_sub, mval
  use mod_ormas, only : nstr_alph, n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph
  use mod_ormas, only : nstr_beta, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta
  use mod_ormas, only : n1x_m_alph, map1x_m_alph
  use mod_ormas, only : n1x_m_beta, map1x_m_beta

  implicit none
  complex(c_double_complex), intent(in) :: fac
  complex(c_double_complex), intent(in) :: cic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(in) :: dcic(1:nstr_alph, 1:nstr_beta)
  complex(c_double_complex), intent(inout) :: dden1(1:nact, 1:nact)

  integer(c_long) :: istr, jstr, kstr, i1x_m, i1x, isub, jsub, iact, jact, iproc, nproc
  integer(c_long), external :: util_omp_nproc
  integer(c_long), external :: util_omp_iproc
  complex(c_double_complex), allocatable :: dden1p(:,:,:)

  nproc = util_omp_nproc()
  allocate(dden1p(1:nact, 1:nact, 0:(nproc-1)))
  dden1p(1:nact, 1:nact, 0:(nproc-1)) = czero

!$omp parallel default(shared) private(iact,jact,istr,jstr,kstr,i1x,i1x_m,iproc)
  iproc = util_omp_iproc()
!$omp do 
  do istr = 1, nstr_beta
!M-adapt
!    do i1x = 1, n1x_beta(0, istr)
     do i1x_m = 1, n1x_m_beta(0, istr)
        i1x = map1x_m_beta(i1x_m, 0, istr)
!M-adapt
        iact = h1x_beta (i1x, istr)
        jact = p1x_beta (i1x, istr)
        kstr = eq1x_beta(i1x, istr)
        !TMP
        !if (kstr <= 0) cycle
        !TMP

        do jstr = 1, nstr_alph
           dden1p(jact, iact, iproc) = dden1p(jact, iact, iproc) &
              & + conjg(cic(jstr, istr)) &
              &      * dcic(jstr, kstr) * sgn1x_beta(i1x, istr)
        end do
     end do
  end do
!$omp end do
  if (nelact(1) == nelact(2)) then
     dden1p(1:nact, 1:nact, iproc) = dden1p(1:nact, 1:nact, iproc) * ctwo
  else
!$omp do
     do istr = 1, nstr_alph
!M-adapt
!       do i1x = 1, n1x_alph(0, istr)
        do i1x_m = 1, n1x_m_alph(0, istr)
           i1x = map1x_m_alph(i1x_m, 0, istr)
!M-adapt
           iact = h1x_alph (i1x, istr)
           jact = p1x_alph (i1x, istr)
           kstr = eq1x_alph(i1x, istr)
           !TMP
           !if (kstr <= 0) cycle
           !TMP

           do jstr = 1, nstr_beta
              dden1p(jact, iact, iproc) = dden1p(jact, iact, iproc) &
                 & + conjg(cic(istr, jstr)) &
                 &      * dcic(kstr, jstr) * sgn1x_alph(i1x, istr)
           end do
        end do
     end do
!$omp end do
  end if
!$omp end parallel

  do iproc = 1, nproc - 1
     dden1p(1:nact, 1:nact, 0) = dden1p(1:nact, 1:nact, 0) + dden1p(1:nact, 1:nact, iproc)
  end do

  ! unify the bra-derivative contributions, which are a.h.c. of ket-derivatives
  do isub = 1, nsub
     do jsub = 1, isub - 1
        do iact = lorb_sub(1, isub), lorb_sub(2, isub)
           do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
              dden1p(jact, iact, 0) = dden1p(jact, iact, 0) - conjg(dden1p(iact, jact, 0))
              dden1p(iact, jact, 0) = -conjg(dden1p(jact, iact, 0))
           end do
        end do
        do iact = lorb_sub(1, isub), lorb_sub(2, isub)
           do jact = lorb_sub(1, jsub), lorb_sub(2, jsub)
              if (mval(ncore+iact) == mval(ncore+jact)) then
                 dden1(iact, jact) = dden1(iact, jact) + dden1p(iact, jact, 0) * fac
                 dden1(jact, iact) = dden1(jact, iact) + dden1p(jact, iact, 0) * fac
              end if
           end do
        end do
     end do
  end do

  deallocate(dden1p)

!debug
!  write(6, "('ormas_mkdden1_old input:')")
!  do istr = 1, nstr_beta
!     do jstr = 1, nstr_alph
!        write(6,"(2i5,2f20.10)") jstr,istr,dble(cic(jstr,istr)),dble(dcic(jstr,istr))
!     end do
!  end do
!  write(6, "('ormas_mkdden1_old output:')")
!  do iact = 1, nact
!     do jact = 1, nact
!        write(6, "(f20.10)", advance = 'no') dble(dden1(jact, iact))
!     end do
!     write(6, *)
!  end do
!debug
end subroutine ormas_mkdden1_old
!######################################################################
