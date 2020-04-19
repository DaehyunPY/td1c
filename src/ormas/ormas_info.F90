!################################################################################
subroutine ormas_info()

  use, intrinsic :: iso_c_binding  
  use mod_ormas, only : iprint, MAX_NBLOCK, MAX_NSUB
  use mod_ormas, only : nblock, type_block, nfun_block, nelcore, nelact, neltot, &
       & nsub, min_sub, max_sub, norb_sub, lorb_sub, nfcore2, nfcore1, nfcore, &
       & ndcore, ncore, nact, nocc, nvir, nfun

  implicit none
  integer(c_long) :: iblock, isub

  if (nblock <= 0 .or. nblock > max_nblock) stop "bad nblock @ ormas_info."
  do iblock = 1, nblock
     if (type_block(iblock) < -2) stop "error in wfn.type_block."
     if (nfun_block(iblock) < 0) stop "error in wfn.nfun_block."
  end do

  nfcore2 = 0
  nfcore1 = 0
  nfcore = 0
  ndcore = 0
  ncore = 0
  nact = 0
  nocc = 0
  nvir = 0
  nfun = 0
  do iblock = 1, nblock
     if (type_block(iblock) == -2) then
        nfcore2 = nfcore2 + nfun_block(iblock)
     else if (type_block(iblock) == -1) then
        nfcore1 = nfcore1 + nfun_block(iblock)
     else if (type_block(iblock) == 0) then
        ndcore = ndcore + nfun_block(iblock)
     else if (type_block(iblock) == 1) then
        nact = nact + nfun_block(iblock)
     else if (type_block(iblock) == 2) then
        nvir = nvir + nfun_block(iblock)
     else
        stop 'bad type_block @ ormas_info.'
     end if
  end do

  nfcore = nfcore2 + nfcore1
  ncore = nfcore + ndcore
  nocc = ncore + nact
  nfun = nocc + nvir

  nelcore(1) = ncore
  nelcore(2) = ncore
  nelcore(3) = nelcore(1) + nelcore(2)
  nelact(1) = neltot(1) - nelcore(1)
  nelact(2) = neltot(2) - nelcore(2)
  nelact(3) = nelact(1) + nelact(2)

  if (nsub == 1) then
    norb_sub(1) = nact
    min_sub(1) = nelact(3)
    max_sub(1) = nelact(3)
  end if

  write(6, *) 'nblock     =', nblock
  write(6, *) 'type_block =', type_block(1:nblock)
  write(6, *) 'nfun_block =', nfun_block(1:nblock)
  write(6, *) 'nelcore    =', nelcore(1:3)
  write(6, *) 'nelact     =', nelact(1:3)
  write(6, *) 'neltot     =', neltot(1:3)
  write(6, *) 'nsub       =', nsub
  write(6, *) 'min_sub    =', min_sub(1:nsub)
  write(6, *) 'max_sub    =', max_sub(1:nsub)
  write(6, *) 'norb_sub   =', norb_sub(1:nsub)
  write(6, *) 'lorb_sub   =', lorb_sub(1:2, 1:nsub)
  write(6, *) 'nfcore2    =', nfcore2
  write(6, *) 'nfcore1    =', nfcore1
  write(6, *) 'nfcore     =', nfcore
  write(6, *) 'ndcore     =', ndcore
  write(6, *) 'ncore      =', ncore
  write(6, *) 'nact       =', nact
  write(6, *) 'nocc       =', nocc
  write(6, *) 'nvir       =', nvir
  write(6, *) 'nfun       =', nfun

  if (nsub <= 0 .or. nsub > max_nsub) stop "bad nsub @ ormas_info."
  if (sum(norb_sub(1:nsub)) /= nact) stop 'error in wfn.norb_sub (1).'
  if (sum(norb_sub(1:nsub))*2 < nelact(3)) stop 'error in wfn.norb_sub (2).'
  if (sum(min_sub(1:nsub)) > nelact(3)) stop 'error in wfn.min_sub.'
  if (sum(max_sub(1:nsub)) < nelact(3)) stop 'error in wfn.max_sub.'

  do isub = 1, nsub
     if (min_sub(isub) > max_sub(isub)) stop 'error in wfn.min/max_sub.'
     if (max_sub(isub) > min(nelact(3), 2*norb_sub(isub))) stop 'error in wfn.min/max_sub.'
  end do
  lorb_sub(1, 1) = 1
  lorb_sub(2, 1) = norb_sub(1)
  do isub = 2, nsub
     lorb_sub(1, isub) = lorb_sub(2, isub-1) + 1
     lorb_sub(2, isub) = lorb_sub(2, isub-1) + norb_sub(isub)
  end do

  if (iprint > 0) then
     write(6, *) 'nblock     =', nblock
     write(6, *) 'type_block =', type_block(1:nblock)
     write(6, *) 'nfun_block =', nfun_block(1:nblock)
     write(6, *) 'nelcore    =', nelcore(1:3)
     write(6, *) 'nelact     =', nelact(1:3)
     write(6, *) 'neltot     =', neltot(1:3)
     write(6, *) 'nsub       =', nsub
     write(6, *) 'min_sub    =', min_sub(1:nsub)
     write(6, *) 'max_sub    =', max_sub(1:nsub)
     write(6, *) 'norb_sub   =', norb_sub(1:nsub)
     write(6, *) 'lorb_sub   =', lorb_sub(1:2, 1:nsub)
     write(6, *) 'nfcore2    =', nfcore2
     write(6, *) 'nfcore1    =', nfcore1
     write(6, *) 'nfcore     =', nfcore
     write(6, *) 'ndcore     =', ndcore
     write(6, *) 'ncore      =', ncore
     write(6, *) 'nact       =', nact
     write(6, *) 'nocc       =', nocc
     write(6, *) 'nvir       =', nvir
     write(6, *) 'nfun       =', nfun
  end if

end subroutine ormas_info
!################################################################################
