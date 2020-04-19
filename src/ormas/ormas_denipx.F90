!######################################################################
subroutine ormas_denipx(max_ipx,ovlp,cic,ipx,denipx)

  use,intrinsic :: iso_c_binding
  use mod_ormas,only : nfun,ndetx,cic_old,tdcc

  implicit none
  !--------------------------------------------------------------------
  integer(c_int),intent(in) :: max_ipx
  complex(c_double_complex),intent(in) :: ovlp(1:nfun,1:nfun)
  complex(c_double_complex),intent(in) :: cic(1:ndetx)
  real(c_double),intent(out) :: ipx(0:max_ipx)
  complex(c_double),intent(out) :: denipx(1:nfun,1:nfun,0:max_ipx)
  !--------------------------------------------------------------------

  if (tdcc) then
     write(6,"('denipx nyi for tdcc.')")
     stop
  else if (cic_old) then
     write(6,"('denipx nyi for cic_old.')")
     stop
  else
     call ormas_denipx_ras(max_ipx,ovlp,cic,ipx,denipx)
  end if

end subroutine ormas_denipx
!######################################################################
