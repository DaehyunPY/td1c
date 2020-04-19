!#######################################################################
subroutine sph_3j_init()

  use, intrinsic :: iso_c_binding
  use mod_const, only : PI
  use mod_sph, only : sph_gaunt,lmax1,lmax2,mmax1,mmax2
  use fgsl

  implicit none
  integer(c_int) :: l1,l2,l12,m1,m2,m12
  integer(fgsl_int), parameter :: fgsl0 = 0
  integer(fgsl_int) :: tl1,tl2,tl12,tm1,tm2,tm12
  real(c_double) :: lfac, lmfac

  allocate(sph_gaunt(0:lmax2,0:lmax1,0:lmax1,-mmax1:mmax1,-mmax1:mmax1))
  sph_gaunt = 0d0

  do l1 = 0, lmax1
     tl1 = l1*2
     do l2 = 0, lmax1
        tl2 = l2*2
        do l12 = abs(l1-l2), min(lmax2,l1+l2)
           if (mod(l1+l2+l12,2).ne.0) cycle
           tl12 = l12*2
           lfac = fgsl_sf_coupling_3j(tl1,tl2,tl12,fgsl0,fgsl0,fgsl0) &
                * sqrt(((2*l1+1)*(2*l2+1)*(2*l12+1))/(4*PI))
           do m1 = -min(l1,mmax1), min(l1,mmax1)
              tm1 = m1*2
!NOTE HERE
              lmfac = lfac*(-1)**m1
!1              lmfac = lfac
!2              if (m1 < 0) then
!2                 lmfac = lfac*(-1)**m1
!2              else
!2                 lmfac = lfac
!2              end if
!NOTE HERE
              do m2 = -min(l2,mmax1), min(l2,mmax1)
                 tm2 = m2*2
                 m12 = m1 - m2
                 if (l12<abs(m12)) cycle

                 tm12 = tm1 - tm2
                 sph_gaunt(l12,l1,l2,m1,m2) = fgsl_sf_coupling_3j(tl1,tl2,tl12,-tm1,tm2,tm12)*lmfac
              end do
           end do
        end do
     end do
  end do

!Hochstuhl & Bonitz  do l1 = 0, lmax1
!Hochstuhl & Bonitz     tl1 = l1*2
!Hochstuhl & Bonitz     do l2 = 0, lmax1
!Hochstuhl & Bonitz        tl2 = l2*2
!Hochstuhl & Bonitz        do l12 = abs(l1-l2), min(lmax2,l1+l2)
!Hochstuhl & Bonitz           if (mod(l1+l2+l12,2).ne.0) cycle
!Hochstuhl & Bonitz           tl12 = l12*2
!Hochstuhl & Bonitz           lfac = fgsl_sf_coupling_3j(tl1,tl2,tl12,fgsl0,fgsl0,fgsl0) &
!Hochstuhl & Bonitz                * sqrt(((2*l1+1)*(2*l2+1)*(2*l12+1))/(4*PI))
!Hochstuhl & Bonitz           do m1 = -min(l1,mmax1), min(l1,mmax1)
!Hochstuhl & Bonitz              tm1 = m1*2
!Hochstuhl & Bonitz              do m2 = -min(l2,mmax1), min(l2,mmax1)
!Hochstuhl & Bonitz                 tm2 = m2*2
!Hochstuhl & Bonitz                 m12 = -m1-m2
!Hochstuhl & Bonitz                 if (l12<abs(m12)) cycle
!Hochstuhl & Bonitz
!Hochstuhl & Bonitz                 tm12 = -tm1-tm2
!Hochstuhl & Bonitz                 sph_gaunt(l12,l1,l2,m1,m2) = fgsl_sf_coupling_3j(tl1,tl2,tl12,tm1,tm2,tm12)*lfac
!Hochstuhl & Bonitz              end do
!Hochstuhl & Bonitz           end do
!Hochstuhl & Bonitz        end do
!Hochstuhl & Bonitz     end do
!Hochstuhl & Bonitz  end do

  !DEBUG
  do l1 = 0, lmax1
     tl1 = l1*2
     do l2 = 0, lmax1
        tl2 = l2*2
        do l12 = abs(l1-l2), min(lmax2,l1+l2)
           tl12 = l12*2
           if (mod(l1+l2+l12,2).ne.0) cycle
!           lfac = fgsl_sf_coupling_3j(tl1,tl2,tl12,fgsl0,fgsl0,fgsl0)
!           write(6,"('3j000: ',3i5,f20.10)") l1,l2,l12,lfac
           do m1 = -min(l1,mmax1), min(l1,mmax1)
              tm1 = m1*2
              do m2 = -min(l2,mmax1), min(l2,mmax1)
                 tm2 = m2*2
                 m12 = m1 - m2
                 if (l12<abs(m12)) cycle

                 tm12 = tm1 - tm2
                 lmfac = fgsl_sf_coupling_3j(tl1,tl2,tl12,-tm1,tm2,tm12)
!                 write(6,"('3j: ',6i5,f20.10)") l1,l2,l12,m1,m2,m12,lmfac
                 write(6,"('gaunt: ',6i5,f20.10)") l1,l2,l12,m1,m2,m12,sph_gaunt(l12,l1,l2,m1,m2)
              end do
           end do
        end do
     end do
  end do
  !DEBUG

end subroutine sph_3j_init
!#######################################################################
