!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SRO_EOS is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SRO_EOS.  If not, see <http://www.gnu.org/licenses/>.
!
MODULE Skyrme_Bulk_Density_Derivatives_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  !USE Define_Operators_Mod
  USE Skyrme_Coefficients_Mod
  USE Fermi_Integrals_Mod
  USE Skyrme_Bulk_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE SKYRME_BULK_DENSITY_DERIVATIVES( &
    log10_n_n, log10_n_p, T, G_n, G_p, &
    Deta_n_Dn_n, Deta_n_Dn_p, Deta_p_Dn_n, Deta_p_Dn_p, &
    Dtau_n_Dn_n, Dtau_n_Dn_p, Dtau_p_Dn_n, Dtau_p_Dn_p, &
    Dmu_n_Dn_n,  Dmu_n_Dn_p,  Dmu_p_Dn_n,  Dmu_p_Dn_p,  &
    Dv_n_Dn_n,   Dv_n_Dn_p,   Dv_p_Dn_n,   Dv_p_Dn_p,   &
    DP_d_n_n,    DP_d_n_p)

    USE Physical_Constants_Mod, &
        ONLY : ONE, TWO, THREE, FOUR, Hbarc_Square, &
               Mass_n, Mass_p, Neut_Prot_Mass_Diff

    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: T, log10_n_n, log10_n_p
    REAL(DP), INTENT(OUT) :: G_n, G_p
    REAL(DP), INTENT(OUT) :: Deta_n_Dn_n, Deta_n_Dn_p, Deta_p_Dn_n, Deta_p_Dn_p
    REAL(DP), INTENT(OUT) :: Dtau_n_Dn_n, Dtau_n_Dn_p, Dtau_p_Dn_n, Dtau_p_Dn_p
    REAL(DP), INTENT(OUT) :: Dmu_n_Dn_n, Dmu_n_Dn_p, Dmu_p_Dn_n, Dmu_p_Dn_p
    REAL(DP), INTENT(OUT) :: Dv_n_Dn_n, Dv_n_Dn_p, Dv_p_Dn_n, Dv_p_Dn_p

    REAL(DP), INTENT(OUT) :: DP_d_n_n, DP_d_n_p
    REAL(DP) :: n_n,n_p,n,y
    REAL(DP) :: Meff_n, Meff_p
    REAL(DP) :: eta_n, eta_p
    REAL(DP) :: tau_n, tau_p
    REAL(DP) :: v_n, v_p
    REAL(DP) :: mu_n, mu_p
    REAL(DP) :: V_II, V_IJ
    REAL(DP) :: W_nn, W_np, W_pn, W_pp
    REAL(DP) :: P, F

    CALL SKYRME_BULK_PROPERTIES(log10_n_n,log10_n_p,T,&
      n_n,n_p,n,y,Meff_n,Meff_p,eta_n,eta_p, &
      tau_n,tau_p,v_n,v_p,mu_n,mu_p,P,F)

  !   obtain G_t = 2*F_{+1/2}(eta_t)/F_{-1/2}(eta_t)
    IF (eta_n<-200.0_DP) THEN
      G_n = ONE
    ELSE
      G_n = TWO*FHALF(eta_n)
    ENDIF

    IF (eta_p<-200.0_DP) THEN
      G_p = ONE
    ELSE
      G_p = TWO*FHALF(eta_p)
    ENDIF

  !   tau derivatives
    Dtau_n_Dn_n = Meff_n/Hbarc_Square*( THREE*T*G_n &
       +(9.0_DP*Meff_n*T*n_n*G_n/Hbarc_Square &
       - 5.0_DP*tau_n)*coeff_alpha1)

    Dtau_n_Dn_p = Meff_n/Hbarc_Square*( &
        (9.0_DP*Meff_n*T*n_n*G_n/Hbarc_Square &
       - 5.0_DP*tau_n)*coeff_alpha2)

    Dtau_p_Dn_n = Meff_p/Hbarc_Square*( &
        (9.0_DP*Meff_p*T*n_p*G_p/Hbarc_Square &
       - 5.0_DP*tau_p)*coeff_alpha2)

    Dtau_p_Dn_p = Meff_p/Hbarc_Square*( THREE*T*G_p &
       +(9.0_DP*Meff_p*T*n_p*G_p/Hbarc_Square &
       - 5.0_DP*tau_p)*coeff_alpha1)

  !   v derivatives
    V_II = TWO*COEFF_A + DOT_PRODUCT(COEFF_C*COEFF_DELTA*(ONE+COEFF_DELTA), &
          n**(COEFF_DELTA-ONE)) +  &
          FOUR*DOT_PRODUCT(COEFF_D*(COEFF_DELTA-ONE)*(COEFF_DELTA-TWO), &
          EXP((COEFF_DELTA-THREE)*LOG(n)+LOG(n_p)+LOG(n_n)))

    V_IJ = V_II + FOUR*COEFF_B +  &
          FOUR*DOT_PRODUCT(COEFF_D*COEFF_DELTA,n**(COEFF_DELTA-ONE))

    W_nn = V_II + 8.0_DP*DOT_PRODUCT(COEFF_D*(COEFF_DELTA-ONE),&
           EXP((COEFF_DELTA-TWO)*LOG(n)+LOG(n_p)))

    W_np = V_IJ

    W_pn = V_IJ

    W_pp = V_II + 8.0_DP*DOT_PRODUCT(COEFF_D*(COEFF_DELTA-ONE),&
           EXP((COEFF_DELTA-TWO)*LOG(n)+LOG(n_n)))

    Dv_n_Dn_n = COEFF_ALPHA1*Dtau_n_Dn_n + COEFF_ALPHA2*Dtau_p_Dn_n + W_nn

    Dv_n_Dn_p = COEFF_ALPHA1*Dtau_n_Dn_p + COEFF_ALPHA2*Dtau_p_Dn_p + W_np

    Dv_p_Dn_n = COEFF_ALPHA1*Dtau_p_Dn_n + COEFF_ALPHA2*Dtau_n_Dn_n + W_pn

    Dv_p_Dn_p = COEFF_ALPHA1*Dtau_p_Dn_p + COEFF_ALPHA2*Dtau_n_Dn_p + W_pp

  !  eta derivatives
    Deta_n_Dn_n = G_n/n_n*(one + three/hbarc_square*n_n*Meff_n*coeff_alpha1)

    Deta_n_Dn_p = G_n/n_n*(      three/hbarc_square*n_n*Meff_n*coeff_alpha2)

    Deta_p_Dn_n = G_p/n_p*(      three/hbarc_square*n_p*Meff_p*coeff_alpha2)

    Deta_p_Dn_p = G_p/n_p*(one + three/hbarc_square*n_p*Meff_p*coeff_alpha1)

  !   mu derivatives
    Dmu_n_Dn_n = T*Deta_n_Dn_n + Dv_n_Dn_n

    Dmu_n_Dn_p = T*Deta_n_Dn_p + Dv_n_Dn_p

    Dmu_p_Dn_n = T*Deta_p_Dn_n + Dv_p_Dn_n

    Dmu_p_Dn_p = T*Deta_p_Dn_p + Dv_p_Dn_p

  !   pressure derivatives
    DP_d_n_n = n_n*Dmu_n_Dn_n + n_p*Dmu_p_Dn_n
    DP_d_n_p = n_n*Dmu_n_Dn_p + n_p*Dmu_p_Dn_p

  END SUBROUTINE SKYRME_BULK_DENSITY_DERIVATIVES

END MODULE Skyrme_Bulk_Density_Derivatives_Mod

! SUBROUTINE SKYRME_BULK_DENSITY_SECOND_DERIVATIVES(&
!   log10_n_n,log10_n_p,T,    &
!   d2_eta_n_d_n_n_n_n, d2_eta_n_d_n_p_n_n, &
!   d2_eta_n_d_n_n_n_p, d2_eta_n_d_n_p_n_p, &
!   d2_eta_p_d_n_n_n_n, d2_eta_p_d_n_p_n_n, &
!   d2_eta_p_d_n_n_n_p, d2_eta_p_d_n_p_n_p)
!
!   USE Physical_Constants_Mod, &
!       ONLY : ZERO, ONE, TWO, THREE, FOUR, TEN, R_3_2, R_5_2, &
!              Hbarc_Square, PI, TWO_PI_SQUARE, &
!              log10_dens_min, log10_dens_max, dens_min, dens_max, &
!              Mass_n, Mass_p, Neut_Prot_Mass_Diff
!
!   REAL(DP), INTENT(IN)  :: T, log10_n_n, log10_n_p
!   REAL(DP) :: Deta_n_Dn_n, Deta_n_Dn_p
!   REAL(DP) :: Deta_p_Dn_n, Deta_p_Dn_p
!   REAL(DP) :: Dtau_n_Dn_n, Dtau_n_Dn_p
!   REAL(DP) :: Dtau_p_Dn_n, Dtau_p_Dn_p
!   REAL(DP) :: Dmu_n_Dn_n, Dmu_n_Dn_p
!   REAL(DP) :: Dmu_p_Dn_n, Dmu_p_Dn_p
!   REAL(DP) :: Dv_n_Dn_n, Dv_n_Dn_p
!   REAL(DP) :: Dv_p_Dn_n, Dv_p_Dn_p
!   REAL(DP) :: DP_d_n_n, DP_d_n_p
!   REAL(DP) :: n_n,n_p,n,y
!   REAL(DP) :: Meff_n, Meff_p
!   REAL(DP) :: eta_n, eta_p
!   REAL(DP) :: tau_n, tau_p
!   REAL(DP) :: v_n, v_p
!   REAL(DP) :: mu_n, mu_p
!   REAL(DP) :: G_n, G_p, V_II, V_IJ
!   REAL(DP) :: W_nn, W_np, W_pn, W_pp
!
!   CALL SKYRME_BULK_PROPERTIES(log10_n_n,log10_n_p,T,&
!     n_n,n_p,n,y,   &
!     Meff_n,Meff_p,eta_n,eta_p, &
!     tau_n,tau_p,v_n,v_p,mu_n,mu_p)
!
!   CALL SKYRME_BULK_DENSITY_DERIVATIVES(       &
!       log10_n_n,log10_n_p,T,    &
!       Deta_n_Dn_n, Deta_n_Dn_p, &
!       Deta_p_Dn_n, Deta_p_Dn_p, &
!       Dtau_n_Dn_n, Dtau_n_Dn_p, &
!       Dtau_p_Dn_n, Dtau_p_Dn_p, &
!       Dmu_n_Dn_n, Dmu_n_Dn_p,   &
!       Dmu_p_Dn_n, Dmu_p_Dn_p,   &
!       Dv_n_Dn_n, Dv_n_Dn_p,     &
!       Dv_p_Dn_n, Dv_p_Dn_p,     &
!       DP_d_n_n, DP_d_n_p)
!
! END SUBROUTINE SKYRME_BULK_DENSITY_SECOND_DERIVATIVES
