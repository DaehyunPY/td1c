!######################################################################
subroutine tdcc_mkint1x(int1e,int2e,int1x)

  use mod_bas,only : smul
  use mod_ormas,only : nact,nelact
  use mod_cc,only : norb1

  implicit none
  complex(kind(0d0)),intent(in) :: int1e(1:nact,1:nact)
  complex(kind(0d0)),intent(in) :: int2e(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(out) :: int1x(1:nact,1:nact,1:2)

  integer :: iact,jact,kact
  complex(kind(0d0)) :: refjk

  int1x(1:nact,1:nact,1) = 0d0
  if (nelact(1)<1) return

  do iact = 1,nact
     do jact = 1,iact
        int1x(jact,iact,1) = int1e(jact,iact)
        do kact = 1,norb1
           !##### RHF only #####
           if (smul==1.and.nelact(1)==nelact(2)) then
              refjk = 2d0 * int2e(jact,iact,kact,kact) &
                          - int2e(jact,kact,kact,iact)
           else if (nelact(2)==0) then
              refjk = int2e(jact,iact,kact,kact) &
                    - int2e(jact,kact,kact,iact)
           end if
           !##### RHF only #####
           int1x(jact,iact,1) = int1x(jact,iact,1) + refjk
        end do
     end do
  end do
  do iact = 1,nact
     do jact = 1,iact-1
        int1x(iact,jact,1) = conjg(int1x(jact,iact,1))
     end do
  end do

  if (smul==1.and.nelact(1)==nelact(2)) then
     int1x(:,:,2) = int1x(:,:,1)
  else if (nelact(2)==0) then
     int1x(:,:,2) = 0d0
  else
     stop 'tdcc_mkint1x: A=B.or.B=0 only'
  end if

end subroutine tdcc_mkint1x
!######################################################################
