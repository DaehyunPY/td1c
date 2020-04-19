subroutine ormas_fab()

  use, intrinsic :: iso_c_binding
  use mod_ormas , only : nact, nstr_alph, onv_alph, nelact
  use mod_ormas, only : fab_nr1x, fab_n1x,  fab_h1x, fab_p1x, fab_eq1x, fab_sgn1x
  implicit none

  integer(c_long) :: inverse(2**nact)
  integer(c_long) :: Conf_local(nact)

  integer(c_long) :: k,k1,k2,rev,count,i_up,j_down

  !======== Define the control arrays for the mkden2 version 4 by fabian ======

  allocate(fab_eq1x(nstr_alph,nelact(1)*(nact-nelact(1)+1)))
  allocate(fab_sgn1x(nstr_alph,nelact(1)*(nact-nelact(1)+1)))
  allocate(fab_h1x(nstr_alph,nelact(1)*(nact-nelact(1)+1)))
  allocate(fab_p1x(nstr_alph,nelact(1)*(nact-nelact(1)+1)))
  allocate(fab_nr1x(nstr_alph))

  !======== We assume a singlet spin state ====================================

  inverse=0
  do k1=1,nstr_alph

     rev=onv_alph(1,k1)
     do k=2,nact
        rev=rev*2+onv_alph(k,k1)
     enddo

     inverse(rev)=k1
  enddo

  do k1=1,nstr_alph

     count=0
     do j_down=1,nact
        do i_up=j_down,nact

           Conf_local(:)=onv_alph(:,k1)

           Conf_local(j_down)= Conf_local(j_down)-1
           If(any(Conf_local.lt.0)) cycle

           Conf_local(i_up)= Conf_local(i_up)+1
           If(any(Conf_local.gt.1)) cycle

           count=count+1

           rev=Conf_local(1)
           do k=2,nact
              rev=rev*2+Conf_local(k)
           enddo
           fab_eq1x(k1,count)=inverse(rev)

           fab_h1x(k1,count)=i_up
           fab_p1x(k1,count)=j_down

           fab_sgn1x(k1,count) = 1
           do k = min(i_up,j_down) + 1, max(i_up,j_down) - 1
              if (onv_alph(k, k1).ne.0) fab_sgn1x(k1,count)=(-1)*fab_sgn1x(k1,count)
           end do

        enddo
     enddo

     fab_nr1x(k1)=count

     do i_up=1,nact
        do j_down=i_up+1,nact

           Conf_local(:)=onv_alph(:,k1)

           Conf_local(j_down)= Conf_local(j_down)-1
           If(any(Conf_local.lt.0)) cycle

           Conf_local(i_up)= Conf_local(i_up)+1
           If(any(Conf_local.gt.1)) cycle

           count=count+1

           rev=Conf_local(1)
           do k=2,nact
              rev=rev*2+Conf_local(k)
           enddo
           fab_eq1x(k1,count)=inverse(rev)

           fab_h1x(k1,count)=i_up
           fab_p1x(k1,count)=j_down

           fab_sgn1x(k1,count) = 1
           do k = min(i_up,j_down) + 1, max(i_up,j_down) - 1
              if (onv_alph(k, k1).ne.0) fab_sgn1x(k1,count)=(-1)*fab_sgn1x(k1,count)
           end do

        enddo
     enddo

     fab_n1x=nelact(1)*(nact-nelact(1)+1) !equivalent: Ntot=count

  enddo

end subroutine ormas_fab
