!######################################################################
subroutine tdcc_mkint2x(int2e,int2x)

  use, intrinsic :: iso_c_binding
  use mod_bas,only : smul
  use mod_ormas,only : nact,nelact

  implicit none
  complex(kind(0d0)),intent(in) :: int2e(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(out) :: int2x(1:nact,1:nact,1:nact,1:nact,1:6)

  integer(c_int) :: iact,jact,kact,lact

  int2x = 0d0
  if (nelact(1)>=2) then
     do iact = 1,nact
        do jact = 1,iact-1
           if (iact == jact) cycle
           do kact = 1,nact
              do lact = 1,kact-1
                 ! (++||++)
                 int2x(iact,jact,kact,lact,1) = +int2e(iact,kact,jact,lact) &
                                                -int2e(iact,lact,jact,kact)
                 int2x(iact,jact,lact,kact,1) = -int2x(iact,jact,kact,lact,1)
                 int2x(jact,iact,kact,lact,1) = -int2x(iact,jact,kact,lact,1)
                 int2x(jact,iact,lact,kact,1) = +int2x(iact,jact,kact,lact,1)
              end do
           end do
        end do
     end do
  end if
  if (nelact(2)>=1) then
     do iact = 1,nact
        do jact = 1,nact
           do kact = 1,nact
              do lact = 1,nact
                 int2x(lact,jact,kact,iact,3) = +int2e(lact,kact,jact,iact)   ! (+-||+-)
                 int2x(lact,jact,kact,iact,5) = -int2e(lact,iact,jact,kact)   ! (+-||-+)
              end do
           end do
        end do
     end do
  end if

  if (smul==1.and.nelact(1)==nelact(2)) then
     int2x(:,:,:,:,2) = int2x(:,:,:,:,1) ! (--||--)
     int2x(:,:,:,:,4) = int2x(:,:,:,:,3) ! (-+||-+)
     int2x(:,:,:,:,6) = int2x(:,:,:,:,5) ! (-+||+-)
  else if (nelact(2)==0) then
     int2x(:,:,:,:,2) = 0d0
     int2x(:,:,:,:,4) = 0d0
     int2x(:,:,:,:,6) = 0d0
  else
     stop 'tdcc_mkint2x: A=B.or.B=0 only'
  end if

end subroutine tdcc_mkint2x
!######################################################################
