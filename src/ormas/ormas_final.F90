!################################################################################
subroutine ormas_final()

  use mod_ormas, only : nact, cic_old
  use mod_ormas, only : sub_orb, dist, det_allowed, &
       & mapr_detx, rotoo_mapf, rotoo_mapb, rotca_mapf, rotca_mapb, &
       & rotaa_mapf, rotaa_mapb
  use mod_ormas, only : dist_alph, nstr_alph_dist, lstr_alph_dist, &
       & nstr_alph_dist_sub, arc_alph, dist_str_alph, substr_alph, onv_alph, &
       & orb_alph, n1x_alph, p1x_alph, h1x_alph, eq1x_alph, sgn1x_alph, &
       & n1x_m_alph, map1x_m_alph
!OLD       & n1xr_alph, r1xr_alph, l1xr_alph, sgn1xr_alph, 
  use mod_ormas, only : dist_beta, nstr_beta_dist, lstr_beta_dist, &
       & nstr_beta_dist_sub, arc_beta, dist_str_beta, substr_beta, onv_beta, &
       & orb_beta, n1x_beta, p1x_beta, h1x_beta, eq1x_beta, sgn1x_beta, &
       & n1x_m_beta, map1x_m_beta
!OLD       & n1xr_beta, r1xr_beta, l1xr_beta, sgn1xr_beta, 
  use mod_ormas, only : fab_den2, fab_nr1x, fab_n1x,  fab_h1x, fab_p1x, fab_eq1x, fab_sgn1x
  use mod_ormas, only : mval_alph, nstr_alph_beta, llstr_alph_beta, ntot_alph_beta
  use mod_ormas, only : mval_beta, nstr_beta_alph, llstr_beta_alph, ntot_beta_alph
  use mod_ormas, only : nstr_dist_m_alph, llstr_dist_m_alph
  use mod_ormas, only : nstr_dist_m_beta, llstr_dist_m_beta
  use mod_ormas, only : map2to1_alph, map1to2_alph
  use mod_ormas, only : map2to1_beta, map1to2_beta
  use mod_ormas, only : tdcc

  implicit none
  if (nact == 0) return

  if (tdcc) call tdcc_final()

  deallocate(sub_orb)

  deallocate(dist_alph)
  deallocate(dist_beta)
  deallocate(dist)

  deallocate(det_allowed)
  deallocate(nstr_alph_dist_sub)
  deallocate(nstr_beta_dist_sub)
  deallocate(nstr_alph_dist)
  deallocate(nstr_beta_dist)
  deallocate(lstr_alph_dist)
  deallocate(lstr_beta_dist)

  deallocate(mapr_detx)
!  deallocate(mapf_detx)
  if (.not. cic_old) then
     deallocate(nstr_dist_m_beta)
     deallocate(nstr_dist_m_alph)
     deallocate(llstr_dist_m_beta)
     deallocate(llstr_dist_m_alph)
     deallocate(mval_beta)
     deallocate(mval_alph)
     deallocate(llstr_alph_beta)
     deallocate(llstr_beta_alph)
     deallocate(ntot_alph_beta)
     deallocate(ntot_beta_alph)
     deallocate(nstr_alph_beta)
     deallocate(nstr_beta_alph)
     deallocate(map2to1_beta)
     deallocate(map2to1_alph)
     deallocate(map1to2_beta)
     deallocate(map1to2_alph)
  end if

  deallocate(arc_alph)
  deallocate(arc_beta)

  deallocate(dist_str_alph)
  deallocate(dist_str_beta)
  deallocate(substr_alph)
  deallocate(substr_beta)
  deallocate(onv_alph)
  deallocate(onv_beta)
  deallocate(orb_alph)
  deallocate(orb_beta)

  deallocate(n1x_alph)
  deallocate(p1x_alph)
  deallocate(h1x_alph)
  deallocate(eq1x_alph)
  deallocate(sgn1x_alph)
  deallocate(n1x_beta)
  deallocate(p1x_beta)
  deallocate(h1x_beta)
  deallocate(eq1x_beta)
  deallocate(sgn1x_beta)

!OLD  deallocate(n1xr_alph)
!OLD  deallocate(r1xr_alph)
!OLD  deallocate(l1xr_alph)
!OLD  deallocate(sgn1xr_alph)
!OLD  deallocate(n1xr_beta)
!OLD  deallocate(r1xr_beta)
!OLD  deallocate(l1xr_beta)
!OLD  deallocate(sgn1xr_beta)

  deallocate(rotaa_mapf)
  deallocate(rotaa_mapb)
  deallocate(rotca_mapf)
  deallocate(rotca_mapb)
  deallocate(rotoo_mapf)
  deallocate(rotoo_mapb)

  deallocate(n1x_m_alph)
  deallocate(n1x_m_beta)
  deallocate(map1x_m_alph)
  deallocate(map1x_m_beta)

  if (fab_den2) then
     deallocate(fab_eq1x)
     deallocate(fab_sgn1x)
     deallocate(fab_h1x)
     deallocate(fab_p1x)
     deallocate(fab_nr1x)
  end if

end subroutine ormas_final
!################################################################################
