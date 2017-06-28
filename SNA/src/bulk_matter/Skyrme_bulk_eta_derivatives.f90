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
MODULE Skyrme_Bulk_Eta_Derivatives_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  !USE Define_Operators_Mod
  USE Skyrme_Coefficients_Mod
  USE Fermi_Integrals_Mod
  USE Skyrme_Bulk_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE SKYRME_BULK_ETA_DERIVATIVES( n, T, Y, n_n, n_p, &
    Meff_n, Meff_p, eta_n, eta_p, tau_n, tau_p, v_n, v_p, mu_n, mu_p, &
    Deta_n_Dn_n, Deta_n_Dn_p, Deta_p_Dn_n, Deta_p_Dn_p, &
    Dtau_n_Dn_n, Dtau_n_Dn_p, Dtau_p_Dn_n, Dtau_p_Dn_p, &
    Dmu_n_Dn_n,  Dmu_n_Dn_p,  Dmu_p_Dn_n,  Dmu_p_Dn_p,  &
    Dv_n_Dn_n,   Dv_n_Dn_p,   Dv_p_Dn_n,   Dv_p_Dn_p,   &
    DP_Dn_n,     DP_Dn_p,     DP_Deta_n,   DP_Deta_p,   &
    Dn_n_Deta_n,   Dn_n_Deta_p,   Dn_p_Deta_n,   Dn_p_Deta_p,   &
    Dmu_n_Deta_n,  Dmu_n_Deta_p,  Dmu_p_Deta_n,  Dmu_p_Deta_p,  &
    Dtau_n_Deta_n, Dtau_n_Deta_p, Dtau_p_Deta_n, Dtau_p_Deta_p )

    USE Physical_Constants_Mod, &
        ONLY : ONE, TWO, THREE, FOUR, Hbarc_Square, &
               Mass_n, Mass_p, Neut_Prot_Mass_Diff

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: T, n, y, n_n, n_p
    REAL(DP), INTENT(IN) :: Meff_n, Meff_p
    REAL(DP), INTENT(IN) :: eta_n, eta_p
    REAL(DP), INTENT(IN) :: tau_n, tau_p
    REAL(DP), INTENT(IN) :: v_n, v_p
    REAL(DP), INTENT(IN) :: mu_n, mu_p
    REAL(DP), INTENT(IN) :: Deta_n_Dn_n, Deta_n_Dn_p, Deta_p_Dn_n, Deta_p_Dn_p
    REAL(DP), INTENT(IN) :: Dtau_n_Dn_n, Dtau_n_Dn_p, Dtau_p_Dn_n, Dtau_p_Dn_p
    REAL(DP), INTENT(IN) :: Dmu_n_Dn_n, Dmu_n_Dn_p, Dmu_p_Dn_n, Dmu_p_Dn_p
    REAL(DP), INTENT(IN) :: Dv_n_Dn_n, Dv_n_Dn_p, Dv_p_Dn_n, Dv_p_Dn_p
    REAL(DP), INTENT(IN) :: DP_Dn_n, DP_Dn_p

    REAL(DP), INTENT(OUT) :: DP_Deta_n,   DP_Deta_p
    REAL(DP), INTENT(OUT) :: Dn_n_Deta_n,   Dn_n_Deta_p,   Dn_p_Deta_n,   Dn_p_Deta_p
    REAL(DP), INTENT(OUT) :: Dmu_n_Deta_n,  Dmu_n_Deta_p,  Dmu_p_Deta_n,  Dmu_p_Deta_p
    REAL(DP), INTENT(OUT) :: Dtau_n_Deta_n, Dtau_n_Deta_p, Dtau_p_Deta_n, Dtau_p_Deta_p
    REAL(DP) :: Dv_n_Deta_n,   Dv_n_Deta_p,   Dv_p_Deta_n,   Dv_p_Deta_p

    REAL(DP) :: G_n, G_p, V_II, V_IJ
    REAL(DP) :: W_nn, W_np, W_pn, W_pp
    REAL(DP) :: aux_n, aux_p, aux_q

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

    aux_n = THREE*Meff_n*n_n/Hbarc_Square
    aux_p = THREE*Meff_p*n_p/Hbarc_Square
    aux_q = ONE + coeff_alpha1*(aux_n+aux_p) &
          + aux_n*aux_p*(coeff_alpha1**TWO-coeff_alpha2**TWO)
!
! derivative of n_t w.r.t. eta_r (fixed T)
!
    Dn_n_Deta_n =   n_n/G_n*(ONE+coeff_alpha1*aux_p)/aux_q
    Dn_n_Deta_p = - n_p/G_p*coeff_alpha2*aux_n/aux_q
    Dn_p_Deta_n = - n_n/G_n*coeff_alpha2*aux_p/aux_q
    Dn_p_Deta_p =   n_p/G_p*(ONE+coeff_alpha1*aux_n)/aux_q
!
! derivative of mu_t w.r.t. eta_r (fixed T)
!
    Dmu_n_Deta_n = T + Dv_n_Dn_n*Dn_n_Deta_n + Dv_n_Dn_p*Dn_p_Deta_n
    Dmu_n_Deta_p =     Dv_n_Dn_n*Dn_n_Deta_p + Dv_n_Dn_p*Dn_p_Deta_p
    Dmu_p_Deta_n =     Dv_p_Dn_n*Dn_n_Deta_n + Dv_p_Dn_p*Dn_p_Deta_n
    Dmu_p_Deta_p = T + Dv_p_Dn_n*Dn_n_Deta_p + Dv_p_Dn_p*Dn_p_Deta_p
!
!  derivative of tau_t w.r.t. eta_r (fixed T)
!
    Dtau_n_Deta_n = Dtau_n_Dn_n*Dn_n_Deta_n + Dtau_n_Dn_p*Dn_p_Deta_n
    Dtau_n_Deta_p = Dtau_n_Dn_n*Dn_n_Deta_p + Dtau_n_Dn_p*Dn_p_Deta_p
    Dtau_p_Deta_n = Dtau_p_Dn_n*Dn_n_Deta_n + Dtau_p_Dn_p*Dn_p_Deta_n
    Dtau_p_Deta_p = Dtau_p_Dn_n*Dn_n_Deta_p + Dtau_p_Dn_p*Dn_p_Deta_p
!
! derivative of V_t w.r.t. eta_r (fixed T)
!
    Dv_n_Deta_n = Dv_n_Dn_n*Dn_n_Deta_n + Dv_n_Dn_p*Dn_p_Deta_n
    Dv_n_Deta_p = Dv_n_Dn_n*Dn_n_Deta_p + Dv_n_Dn_p*Dn_p_Deta_p
    Dv_p_Deta_n = Dv_p_Dn_n*Dn_n_Deta_n + Dv_p_Dn_p*Dn_p_Deta_n
    Dv_p_Deta_p = Dv_p_Dn_n*Dn_n_Deta_p + Dv_p_Dn_p*Dn_p_Deta_p
!
! derivative of P_out w.r.t. eta_r (fixed T)
!
    DP_Deta_n = DP_Dn_n*Dn_n_Deta_n + DP_Dn_p*Dn_p_Deta_n
    DP_Deta_p = DP_Dn_n*Dn_n_Deta_p + DP_Dn_p*Dn_p_Deta_p

  END SUBROUTINE SKYRME_BULK_ETA_DERIVATIVES

END MODULE Skyrme_Bulk_Eta_Derivatives_Mod
