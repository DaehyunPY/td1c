!######################################################################
subroutine bas_gen_tmat(tmat)
!
  use, intrinsic :: iso_c_binding
  use mod_sph, only : lmax1, mmax1
  use mod_const, only : zero, one, half, iunit
  use mod_control, only : SAE, PSP, PSP_Type, name
! Orimo_ECS
  use mod_rad, only : nrad, ndvr, xrad, radk, radk0, ecs_flag, theta, recs, wrad, cwrad, cxrad, radkI_ecs, mapf, mapb, irad_ecs
! Orimo_ECS
  use mod_bas, only : znuc, pp_vloc, pp_vlocHF
!
! atomic hamiltonian in banded storage
!
! example: KL = 2, KU = 1:
! <--------- input ---------->     <--------- output --------->
!  *    *    *    +    +    +       *    *    *   u14  u25  u36
!  *    *    +    +    +    +       *    *   u13  u24  u35  u46
!  *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
! a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
! a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
! a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
  implicit none
  complex(c_double_complex), intent(out) :: tmat(1:(2*ndvr+1), 1:(nrad-1), 0:lmax1)
  integer(c_int) :: ifun, irad, jrad, krad, l, l2, jll, jul, jb1, ctdep1, ctdep2
  integer(c_int) :: dim, nsub, ld1
  complex(c_double_complex) :: oor, zor, zor1
! Orimo_ECS
  complex(c_double_complex) :: pot(1:(nrad-1)), grad(1:(nrad-1))
! Orimo_ECS
!DEBUG
!  integer(c_int) :: tmp_ecs_flag, idvr, jb2
!DEBUG
  
  dim = nrad - 1
  nsub = ndvr
  ld1 = 2 * nsub + 1
  ctdep1 = 2 * ndvr + 1
  ctdep2 = ndvr
  tmat (1:ld1, 1:dim, 0:lmax1) = zero

!DEBUG
!  tmp_ecs_flag = ecs_flag
!  ecs_flag = 0
!DEBUG

! Orimo_ECS
!  if(PSP .and. ecs_flag == 1) then
!     stop 'bas_gen_tmat (PSP) : ecs nyi.'
!  end if
  if(ecs_flag == 0) then
     do irad = 1, nrad-1 
        grad(irad) = xrad(irad)
     end do
  else if(ecs_flag == 1) then
     do irad = 1, nrad-1 
        grad(irad) = cxrad(irad)
     end do
  end if
! Orimo_ECS

  do l = 0, lmax1
     l2 = l * (l + 1)
     do irad = 1, nrad - 1
! Orimo_ECS
        !oor = one / xrad(irad)
        oor = one / grad(irad)
! Orimo_ECS
        if (PSP) then
           if (psp_type == 1 .or. psp_type == 5 .or. psp_type == 6 .or. psp_type == 7) then
              ! fully nonlocal separable PP
              zor = oor**2.d+0*dble(l2)*half + pp_vloc(irad)
           else
              ! semilocal PP
              zor = oor**2.d+0*dble(l2)*half + pp_vlocHF(irad,min(l,2))
           end if
        else if (SAE) then
           ! For He: J.Phys.B: At. Mol. Opt. Phys. 38, 2593 (2005)
! Orimo_ECS
           !zor = oor * (-one &
           !     - 1.231D+0 * exp(-0.662D+0*xrad(irad)) &
           !     + 1.325D+0 * exp(-1.236D+0*xrad(irad)) * xrad(irad) &
           !     + 0.231D+0 * exp(-0.480D+0*xrad(irad)) &
           !     + dble(l2) * oor * half);
           zor = oor * (-one &
                - 1.231D+0 * exp(-0.662D+0*xrad(irad)) &
                + 1.325D+0 * exp(-1.236D+0*xrad(irad)) * grad(irad) &
                + 0.231D+0 * exp(-0.480D+0*xrad(irad)) &
                + dble(l2) * oor * half);
! Orimo_ECS
        else
           ! pure Coulomb
           zor = oor*(-dble(znuc) + dble(l2)*oor*half);
        end if
        tmat(1 + nsub, irad, l) = zor
     end do
  end do

  if(ecs_flag == 0) then
     do irad = 1, nrad - 1
        jll = max(1,        irad - nsub)
        jul = min(nrad - 1, irad + nsub)
        do jrad = jll, jul
           jb1 = nsub + 1 + jrad - irad
           tmat(jb1, irad, 0) = tmat(jb1, irad, 0) + radk(jb1, irad)
        end do
     end do
     do l = 1, lmax1
        do irad = 1, nrad - 1
           jll = max(1,        irad - nsub)
           jul = min(nrad - 1, irad + nsub)
           do jrad = jll, jul
              jb1 = nsub + 1 + jrad - irad
              tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk0(jb1, irad)
           end do
        end do
     end do
! Orimo_ECS
  else if (ecs_flag == 1) then
     ! for l >= 0
     do l = 0, 0
        ctdep1 = ndvr * 2 + 1
        ctdep2 = ndvr
        do irad = 1, nrad - 1

           jll = max(1,        irad - nsub)
           jul = min(nrad - 1, irad + nsub)
           do jrad = jll, jul
              jb1 = nsub + 1 + jrad - irad

              if( irad < irad_ecs  .and. irad + ndvr < irad_ecs) then 
                 ! not ECS area
                 tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk(jb1, irad)

              else if( irad < irad_ecs .and. irad + ndvr >= irad_ecs ) then
                 if(mapb(irad) == 0) then !first finite element
                    if(jb1 == ctdep1 - 1) then
                       ! for ecs, when left boundary ife = I-1 (ECS will start from ife = I.)
                       tmat(jb1, irad, l) = tmat(jb1, irad, l) + &
                            & radk(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) )
                       ctdep1 = ctdep1 - 1
                    else 
                       ! not ECS area
                       tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk(jb1, irad)
                    end if
                 else !not first finite element
                    if(jb1 == ctdep1) then
                       ! for ecs, when left boundary ife = I-1 (ECS will start from ife = I.)
                       ! In my note, blue area
                       tmat(jb1, irad, l) = tmat(jb1, irad, l) + &
                            & radk(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) )
                       ctdep1 = ctdep1 - 1
                    else 
                       ! not ECS area
                       tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk(jb1, irad)
                    end if
                 end if

              else if( irad == irad_ecs ) then
                 if(jb1 < ndvr + 1) then
                    ! for ecs, when upper boundary between ife = I - 1 and I  (ECS will start from ife = I.)
                    ! In my note, blue area
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + &
                         & radk(jb1, irad) * sqrt( wrad(irad) / cwrad(irad) )
                 else if(jb1 == ndvr + 1) then
                    ! for ecs, when center of boundary between ife = I - 1 and I  (ECS will start from ife = I.)
                    ! In my note, black area
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + radkI_ecs
                 else if(jb1 > ndvr + 1) then
                    ! for ecs, when lower boundary between ife = I - 1 and I  (ECS will start from ife = I.)
                    ! In my note, red area
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + &
                         & radk(jb1, irad) * sqrt( wrad(irad) / cwrad(irad) ) * EXP( (-1.5) * iunit * theta)
                 end if

              else if( irad > irad_ecs .and. irad - ndvr <= irad_ecs ) then
                 if(jb1 == ctdep2) then
                    ! In my note, red area
! Sato_ECS
!                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + & ! SATO: check here
!                         & + radk(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) ) * &
!                         & EXP( (-1.5) * iunit * theta)
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) &
                         & + radk(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) ) * &
                         & EXP( (-1.5) * iunit * theta)
! Sato_ECS
                    ctdep2 = ctdep2 - 1
                 else
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk(jb1, irad) * EXP( (-2.0) * iunit * theta)
                 end if

              else if( irad - ndvr > irad_ecs  .and. irad - ndvr >  irad_ecs) then
                 ! completely ECS area
                 tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk(jb1, irad) * EXP( (-2.0) * iunit * theta)   
              end if

           end do
        end do
     end do

     ! for l >= 1
     do l = 1, lmax1
        ctdep1 = ndvr * 2 + 1
        ctdep2 = ndvr
        do irad = 1, nrad - 1

           jll = max(1,        irad - nsub)
           jul = min(nrad - 1, irad + nsub)
           do jrad = jll, jul
              jb1 = nsub + 1 + jrad - irad

              if( irad < irad_ecs  .and. irad + ndvr < irad_ecs) then 
                 ! not ECS area
                 tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk0(jb1, irad)

              else if( irad < irad_ecs .and. irad + ndvr >= irad_ecs ) then
                 if(mapb(irad) == 0) then !first finite element
                    if(jb1 == ctdep1 - 1) then
                       ! for ecs, when left boundary ife = I-1 (ECS will start from ife = I.)
                       tmat(jb1, irad, l) = tmat(jb1, irad, l) + &
                            & radk0(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) )
                       ctdep1 = ctdep1 - 1
                    else 
                       ! not ECS area
                       tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk0(jb1, irad)
                    end if
                 else !not first finite element
                    if(jb1 == ctdep1) then
                       ! for ecs, when left boundary ife = I-1 (ECS will start from ife = I.)
                       ! In my note, blue area
                       tmat(jb1, irad, l) = tmat(jb1, irad, l) + &
                            & radk0(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) )
                       ctdep1 = ctdep1 - 1
                    else 
                       ! not ECS area
                       tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk0(jb1, irad)
                    end if
                 end if

              else if( irad == irad_ecs ) then
                 if(jb1 < ndvr + 1) then
                    ! for ecs, when upper boundary between ife = I - 1 and I  (ECS will start from ife = I.)
                    ! In my note, blue area
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + &
                         & radk0(jb1, irad) * sqrt( wrad(irad) / cwrad(irad) )
                 else if(jb1 == ndvr + 1) then
                    ! for ecs, when center of boundary between ife = I - 1 and I  (ECS will start from ife = I.)
                    ! In my note, black area
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + radkI_ecs
                 else if(jb1 > ndvr + 1) then
                    ! for ecs, when lower boundary between ife = I - 1 and I  (ECS will start from ife = I.)
                    ! In my note, red area
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + &
                         & radk0(jb1, irad) * sqrt( wrad(irad) / cwrad(irad) ) * EXP( (-1.5) * iunit * theta)
                 end if
                 
              else if( irad > irad_ecs .and. irad - ndvr <= irad_ecs ) then
                 if(jb1 == ctdep2) then
                    ! In my note, red area
! Sato_ECS
!                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + & ! SATO: check here
!                         & + radk0(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) ) * &
!                         & EXP( (-1.5) * iunit * theta)
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) &
                         & + radk0(jb1, irad) * sqrt( wrad(irad_ecs) / cwrad(irad_ecs) ) * &
                         & EXP( (-1.5) * iunit * theta)
! Sato_ECS
                    ctdep2 = ctdep2 - 1
                 else
                    tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk0(jb1, irad) * EXP( (-2.0) * iunit * theta)
                 end if
                 
              else if( irad - ndvr > irad_ecs  .and. irad - ndvr >  irad_ecs) then
                 ! completely ECS area
                 tmat(jb1, irad, l) = tmat(jb1, irad, l) + radk0(jb1, irad) * EXP( (-2.0) * iunit * theta)   
              end if
              
           end do
        end do
     end do
  end if
! Orimo_ECS

!DEBUG
!  write(6,"('bas_gen_tmat:')")
!  do irad = 1, nrad - 1
!     jll = max(1,        irad - nsub)
!     jul = min(nrad - 1, irad + nsub)
!     do jrad = jll, jul
!        jb1 = nsub + 1 + jrad - irad
!        jb2 = nsub + 1 + irad - jrad
!        write(6, "(4i5,4f20.10)") irad,jrad,jb1,jb2, &
!             radk (jb1,irad),radk (jb2,jrad), &
!             radk0(jb1,irad),radk0(jb2,jrad)
!     end do
!  end do
!  do l = 0, lmax1
!     do irad = 1, nrad - 1
!        jll = max(1,        irad - nsub)
!        jul = min(nrad - 1, irad + nsub)
!        do jrad = jll, jul
!           jb1 = nsub + 1 + jrad - irad
!           jb2 = nsub + 1 + irad - jrad
!           write(6, "(3i5,5f20.10)") l,irad,jrad,&
!                tmat(jb1,irad,l),tmat(jb2,jrad,l),&
!                abs(tmat(jb1,irad,l)-tmat(jb2,jrad,l))
!        end do
!     end do
!  end do
!  stop
!DEBUG
!stop 'stop for debug @ bas_gen_tmat.'

!DEBUG
!  ecs_flag = tmp_ecs_flag
!DEBUG

end subroutine bas_gen_tmat
!######################################################################
