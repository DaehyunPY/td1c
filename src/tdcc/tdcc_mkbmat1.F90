!######################################################################
subroutine tdcc_mkbmat1(fac,int1e,int2e,den1,den2,bmat)

  use,intrinsic :: iso_c_binding
  use mod_bas,only : mval
  use mod_const,only : czero,chalf
  use mod_ormas,only : ncore,nact,nrotaa,rotaa_mapb
!debug
!  use mod_cc,only : nXai,aiX
!debug

  implicit none
  complex(kind(0d0)),intent(in) :: fac
  complex(kind(0d0)),intent(in) :: int1e(1:nact,1:nact)
  complex(kind(0d0)),intent(in) :: int2e(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(in) :: den1(1:nact,1:nact)
  complex(kind(0d0)),intent(in) :: den2(1:nact,1:nact,1:nact,1:nact)
  complex(kind(0d0)),intent(inout) :: bmat(1:nact,1:nact)
  integer(c_int) :: p1,p2,p,q,r,s,irot
  complex(kind(0d0)) :: tmp1,tmp2

!debug  integer(c_int) :: ai
!debug  do ai = 1, nXai
!debug     irot = aiX(ai)
  do irot = 1,nrotaa
     p2 = rotaa_mapb(irot,1)
     p1 = rotaa_mapb(irot,2)
     tmp1 = czero
     tmp2 = czero
     do p = 1,nact
        if (mval(ncore+p1).ne.mval(ncore+p)) cycle
        tmp1 = tmp1 + int1e(p2,p)*den1(p,p1)
        tmp2 = tmp2 + int1e(p1,p)*den1(p,p2)
     end do
     do p = 1,nact
        do q = 1,nact
           do r = 1,nact
              if (mval(ncore+p1)+mval(ncore+q).ne.&
                  mval(ncore+p)+mval(ncore+r)) cycle
              tmp1 = tmp1 + int2e(p2,p,q,r)*den2(p,p1,r,q)
              tmp2 = tmp2 + int2e(p1,p,q,r)*den2(p,p2,r,q)
           end do
        end do
     end do
     bmat(p2,p1) = bmat(p2,p1) + fac*tmp1
     bmat(p1,p2) = bmat(p1,p2) + fac*tmp2
  end do

end subroutine tdcc_mkbmat1
!######################################################################
