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
MODULE Skyrme_Bulk_Temperature_Derivatives_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  !USE Define_Operators_Mod
  USE Skyrme_Coefficients_Mod, ONLY : A => Coeff_a, B => Coeff_b, C => Coeff_c,&
                  D => Coeff_d, alpha1 => Coeff_alpha1, alpha2 => Coeff_alpha2,&
                  delta => Coeff_delta
  USE Skyrme_Bulk_Mod
  USE Fermi_Integrals_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE SKYRME_BULK_TEMPEARTURE_DERIVATIVES( n_n, n_p, T, &
    Meff_n, Meff_p, eta_n, eta_p, tau_n, tau_p, v_n, v_p, mu_n, mu_p, &
    Dn_n_DT, Dn_p_DT, Dtau_n_DT, Dtau_p_DT, Dmu_n_DT, Dmu_p_DT, DP_DT, fixed )
    ! if fixed = 0 then determine derivatives w.r.t. T with fixed eta_p, eta_n
    ! if fixed = 1 then determine derivatives w.r.t. T with fixed n_n, n_p

    USE Physical_Constants_Mod, &
        ONLY : ZERO, ONE, TWO, THREE, FOUR, R_3_2, R_5_2, R_5_3, R_9_2, &
        Mass_n, Mass_p, hc2 => Hbarc_Square

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: n_n, n_p, T
    REAL(DP), INTENT(IN) :: Meff_n, Meff_p
    REAL(DP), INTENT(IN) :: eta_n, eta_p
    REAL(DP), INTENT(IN) :: tau_n, tau_p
    REAL(DP), INTENT(IN) :: v_n, v_p
    REAL(DP), INTENT(IN) :: mu_n, mu_p

    INTEGER(I4B), INTENT(IN) :: fixed

    REAL(DP), INTENT(OUT) :: Dn_n_DT, Dn_p_DT
    REAL(DP), INTENT(OUT) :: Dtau_n_DT, Dtau_p_DT
    REAL(DP), INTENT(OUT) :: Dmu_n_DT, Dmu_p_DT
    REAL(DP), INTENT(OUT) :: DP_DT

    REAL(DP) :: n, G_n, G_p
    REAL(DP) :: Dv_n_DT, Dv_p_DT, DMeff_n_DT, DMeff_p_DT
    REAL(DP) :: Deta_n_DT, Deta_p_DT
    REAL(DP) :: aux_a, aux_b, aux_c, aux_d, aux_e, aux_f
    REAL(DP) :: log_a, log_b, log_c, log_d, log_e, log_f
    REAL(DP) :: sgn_a, sgn_b, sgn_c, sgn_d, sgn_e, sgn_f
    REAL(DP), DIMENSION(non_local_terms) :: dp1, d00, dm1, dm2, dm3

    IF (fixed == 0) THEN
      n = n_n + n_p

      aux_a = ONE+THREE*n_n*Meff_n/hc2*alpha1
      aux_b = ONE+THREE*n_p*Meff_p/hc2*alpha1
      aux_c = R_3_2*n_n/T
      aux_d = R_3_2*n_p/T
      aux_e =     THREE*n_n*Meff_n/hc2*alpha2
      aux_f =     THREE*n_p*Meff_p/hc2*alpha2
!
! derivative of n_t w.r.t. T (fixed eta)
!
      Dn_n_DT = (aux_b*aux_c-aux_d*aux_e)/(aux_a*aux_b-aux_e*aux_f)
      Dn_p_DT = (aux_a*aux_d-aux_c*aux_f)/(aux_a*aux_b-aux_e*aux_f)
!
! derivative of M^star w.r.t. T (fixed eta)
!
      DMeff_n_DT = - (alpha1*Dn_n_DT + alpha2*Dn_p_DT)*TWO*Meff_n**TWO / hc2
      DMeff_p_DT = - (alpha2*Dn_n_DT + alpha1*Dn_p_DT)*TWO*Meff_p**TWO / hc2
!
! derivative of tau_t w.r.t. T (fixed eta)
!
      Dtau_n_DT = R_5_2*tau_n/T + R_5_2*tau_n*DMeff_n_DT/Meff_n
      Dtau_p_DT = R_5_2*tau_p/T + R_5_2*tau_p*DMeff_p_DT/Meff_p
!
!   some auxiliary quantities
!
      aux_a = Dn_n_DT
      aux_b = Dn_p_DT
      aux_c = aux_a + aux_b
      aux_d = n_n*Dn_p_DT
      aux_e = n_p*Dn_n_DT
      aux_f = aux_d + aux_e


      sgn_a = SIGN(ONE,aux_a)
      sgn_b = SIGN(ONE,aux_b)
      sgn_c = SIGN(ONE,aux_c)
      sgn_f = SIGN(ONE,aux_f)

      log_a = LOG(ABS(aux_a))
      log_b = LOG(ABS(aux_b))
      log_c = LOG(ABS(aux_c))
      log_d = LOG(n_n)
      log_e = LOG(n_p)
      log_f = LOG(ABS(aux_f))

      dp1 = delta + ONE
      d00 = delta
      dm1 = delta - ONE
      dm2 = delta - TWO
      dm3 = delta - THREE
!
! derivative of V_t w.r.t. T (fixed eta)
!
      Dv_n_DT = alpha1*Dtau_n_DT + alpha2*Dtau_p_DT + TWO*a*aux_c + FOUR*b*aux_b &
              + DOT_PRODUCT(c*aux_c*dp1*d00,n**dm1) + FOUR*(  &
                DOT_PRODUCT(d*dm1*dm2*sgn_c,EXP(dm3*LOG(n)+log_c+log_d+log_e)) &
              + DOT_PRODUCT(d*dm1*sgn_f,EXP(dm2*LOG(n)+log_f)) &
              + DOT_PRODUCT(d*dm1*sgn_c,EXP(dm2*LOG(n)+log_e+log_c)) &
              + DOT_PRODUCT(d*sgn_b,EXP(dm1*LOG(n)+log_b)) )

      Dv_p_DT = alpha2*Dtau_n_DT + alpha1*Dtau_p_DT + TWO*a*aux_c + FOUR*b*aux_a &
              + DOT_PRODUCT(c*aux_c*dp1*d00,n**dm1) + FOUR*(  &
              + DOT_PRODUCT(d*dm1*dm2*sgn_c,EXP(dm3*LOG(n)+log_c+log_d+log_e)) &
              + DOT_PRODUCT(d*dm1*sgn_f,EXP(dm2*LOG(n)+log_f)) &
              + DOT_PRODUCT(d*dm1*sgn_c,EXP(dm2*LOG(n)+log_d+log_c)) &
              + DOT_PRODUCT(d*sgn_a,EXP(dm1*LOG(n)+log_a)) )

!
! derivative of mu_t w.r.t. T (fixed eta)
!
      Dmu_n_DT = eta_n + Dv_n_DT
      Dmu_p_DT = eta_p + Dv_p_DT
!
! derivative of P_out w.r.t. T (fixed eta)
!
      DP_DT = hc2/TWO*(R_5_3/Meff_n - ONE/Mass_n)*Dtau_n_DT &
            + hc2/TWO*(R_5_3/Meff_p - ONE/Mass_p)*Dtau_p_DT &
            + R_5_3*tau_n*(alpha1*Dn_n_DT + alpha2*Dn_p_DT) &
            + R_5_3*tau_p*(alpha1*Dn_p_DT + alpha2*Dn_n_DT) &
            + TWO*a*n*aux_c + FOUR*b*aux_f &
            + DOT_PRODUCT(c*d00*dp1*sgn_c,EXP(d00*LOG(n)+log_c)) + FOUR*( &
            + DOT_PRODUCT(d*d00*dm1*sgn_c,EXP(dm2*LOG(n)+log_c+log_d+log_e)) &
            + DOT_PRODUCT(d*d00*sgn_f,EXP(dm1*LOG(n)+log_f)) )
  ENDIF

  IF (fixed == 1) THEN
    n = n_n + n_p
!
! derivative of n_t w.r.t. T (fixed n_i, y_i)
!
    Dn_n_DT = ZERO
    Dn_p_DT = ZERO
!
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
!
! derivative of eta_t w.r.t. T (fixed n_i, x_i)
!
    Deta_n_DT = -R_3_2/T*G_n
    Deta_p_DT = -R_3_2/T*G_p
!
! derivative of tau_t w.r.t. T (fixed n_i, x_i)
!
    Dtau_n_DT = R_5_2*tau_n/T-R_9_2*Meff_n/HC2*n_n*G_n
    Dtau_p_DT = R_5_2*tau_p/T-R_9_2*Meff_p/HC2*n_p*G_p
!
! derivative of V_t w.r.t. T (fixed n_i, x_i)
!
    Dv_n_DT = alpha1*Dtau_n_DT + alpha2*Dtau_p_DT
    Dv_p_DT = alpha1*Dtau_p_DT + alpha2*Dtau_n_DT
!
! derivative of P_in w.r.t. T (fixed n_i, x_i)
!
    DP_DT = HC2/TWO*(R_5_3/Meff_n-ONE/Mass_n)*Dtau_n_DT &
          + HC2/TWO*(R_5_3/Meff_p-ONE/Mass_p)*Dtau_p_DT
!
! derivative of mu_t w.r.t. T (fixed n_i, x_i)
!
    Dmu_n_DT = eta_n + T*Deta_n_DT + Dv_n_DT
    Dmu_p_DT = eta_p + T*Deta_p_DT + Dv_p_DT

  ENDIF

  END SUBROUTINE SKYRME_BULK_TEMPEARTURE_DERIVATIVES

END MODULE Skyrme_Bulk_Temperature_Derivatives_Mod
