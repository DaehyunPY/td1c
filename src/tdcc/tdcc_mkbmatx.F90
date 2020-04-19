!################################################################################
subroutine tdcc_mkbmatx(cic,dcic,bmat)

  use, intrinsic :: iso_c_binding
  use mod_const, only : czero
  use mod_ormas, only : nfun,nact,iprint
  use mod_cc, only : cc_rank,optcc,bcc,optbcc,norb1

  implicit none
  complex(kind(0d0)), intent(in) :: cic(1:*)
  complex(kind(0d0)), intent(in) :: dcic(1:*)
  complex(kind(0d0)), intent(out) :: bmat(1:nact,1:nact)
  !--------------------------------------------------------------------
  complex(kind(0d0)), allocatable :: bmat2(:,:)
  complex(kind(0d0)), allocatable :: bmat3(:,:)

  !debug
  !integer :: a,i
  !debug

  allocate(bmat2(1:nact,1:nact))
  allocate(bmat3(1:nact,1:nact))

  bmat(1:nact,1:nact) = czero
  bmat2(1:nact,1:nact) = czero
  bmat3(1:nact,1:nact) = czero

  if (cc_rank==2 .and. optcc) then
  else if (cc_rank==2 .and. optbcc) then
     call tdcc_mkbmat_2(cic,dcic,bmat2)
  else if (cc_rank==2) then
     stop 'tdcc_mkbmatx: CCSD nyi.'
  else if (cc_rank==3 .and. optcc) then
     call tdcc_mkbmat_2(cic,dcic,bmat2)
     call tdcc_mkbmat_3(cic,dcic,bmat3)
  else if (cc_rank==3 .and. optbcc) then
     call tdcc_mkbmat_2(cic,dcic,bmat2)
     call tdcc_mkbmat_3(cic,dcic,bmat3)
  else if (cc_rank==3) then
     stop 'tdcc_mkbmatx: CCSDT nyi.'
  else
     stop 'tdcc_mkbmatx supports CCD/CCDT/OBCCD/OBCCDT only'
  end if

  bmat = bmat2 - bmat3
!debug  do a = norb1+1,nact
!debug     do i = 1, norb1
!debug        bmat(a,i) = bmat2(i,a) - bmat3(i,a)
!debug     end do
!debug  end do

  if (iprint > 4) then
     write(6, "('# tdcc_mkbmatx: bmatx-2')")
     call util_matoutc(6,nact,bmat2)
     write(6, "('# tdcc_mkbmatx: bmatx-3')")
     call util_matoutc(6,nact,bmat3)
     write(6, "('# tdcc_mkbmatx: bmatx-total')")
     call util_matoutc(6,nact,bmat)
  end if
  !stop 'for debug @ tdcc_mkbmatx.'

  deallocate(bmat2)
  deallocate(bmat3)

end subroutine tdcc_mkbmatx
!######################################################################
