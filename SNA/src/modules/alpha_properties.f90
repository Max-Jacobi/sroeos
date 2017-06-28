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
MODULE Alpha_Properties_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B

  IMPLICIT NONE

CONTAINS

  SUBROUTINE ALPHA_PROPERTIES (objective, n_n, n_p, mu_n, mu_p, T, F_o, P_o, &
    n_alpha, mu_alpha, F_alpha, P_alpha)

    USE Physical_Constants_Mod, &
        ONLY : ZERO, ONE, TWO, THREE, FOUR, R_3_2, Hbarc_Square, &
               v_alpha, b_alpha, Hbarc_Square, PI, Mass_n

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: objective
    REAL(DP), INTENT(IN)  :: n_n, n_p, mu_n, mu_p, T, F_o, P_o
    REAL(DP), INTENT(OUT) :: n_alpha, mu_alpha, F_alpha, P_alpha
    REAL(DP) :: mass_Q, mu_over_T

!   alpha particle's properties
    mass_Q = Mass_n*T/TWO/PI/Hbarc_Square
    mass_Q = mass_Q**R_3_2
    mu_alpha = TWO*(mu_n+mu_p) + b_alpha - P_o*v_alpha
    mu_over_T = mu_alpha/T
!   check if alpha particle chemical potential too large/small
    IF (mu_alpha>ZERO) mu_over_T = MIN(mu_over_T, 2.d2)
    IF (mu_alpha<ZERO) mu_over_T = MAX(mu_over_T,-2.d2)
!   density of alpha particles
    n_alpha = 8.D0*mass_Q*EXP(mu_over_T)
!   alpha particle pressure and free energy if no heavy nuclei present
    P_alpha = n_alpha*T
    F_alpha = n_alpha*(mu_alpha-b_alpha-T)

    IF (objective == 0) RETURN

  END SUBROUTINE ALPHA_PROPERTIES

  SUBROUTINE ALPHA_DERIVATIVES ( n_n, n_p, n_alpha, mu_alpha, T,&
      Dmu_n_Dn_n,   Dmu_p_Dn_n,   Dmu_n_Dn_p,   Dmu_p_Dn_p,     &
      Dmu_n_Deta_n, Dmu_p_Deta_n, Dmu_n_Deta_p, Dmu_p_Deta_p,   &
      Dmu_n_DT, Dmu_p_DT, DP_DT,                                &
      Dn_n_Deta_n,  Dn_p_Deta_n,  Dn_n_Deta_p,  Dn_p_Deta_p,    &
      DP_Dn_n, DP_Dn_p, DP_Deta_n, DP_Deta_p,                   &
      Dn_alpha_Dn_n, Dn_alpha_Dn_p, Dmu_alpha_Dn_n, Dmu_alpha_Dn_p, &
      Dn_alpha_Deta_n, Dn_alpha_Deta_p, Dmu_alpha_Deta_n, Dmu_alpha_Deta_p, &
      DP_alpha_Deta_n, DP_alpha_Deta_p, Dn_alpha_DT, DP_alpha_DT, Dmu_alpha_DT )

    USE Physical_Constants_Mod, ONLY : ZERO, TWO, V_ALPHA, R_3_2

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: n_n, n_p, n_alpha, mu_alpha, T
    REAL(DP), INTENT(IN) :: Dmu_n_Dn_n, Dmu_p_Dn_n
    REAL(DP), INTENT(IN) :: Dmu_n_Dn_p, Dmu_p_Dn_p
    REAL(DP), INTENT(IN) :: Dmu_n_Deta_n, Dmu_p_Deta_n
    REAL(DP), INTENT(IN) :: Dmu_n_Deta_p, Dmu_p_Deta_p
    REAL(DP), INTENT(IN) :: DP_Dn_n, DP_Dn_p
    REAL(DP), INTENT(IN) :: Dn_n_Deta_n, Dn_p_Deta_n
    REAL(DP), INTENT(IN) :: Dn_n_Deta_p, Dn_p_Deta_p
    REAL(DP), INTENT(IN) :: DP_Deta_n, DP_Deta_p

    REAL(DP), INTENT(IN) :: Dmu_n_DT, Dmu_p_DT, DP_DT

    REAL(DP), INTENT(OUT) :: Dn_alpha_Dn_n, Dn_alpha_Dn_p
    REAL(DP), INTENT(OUT) :: Dmu_alpha_Dn_n, Dmu_alpha_Dn_p
    REAL(DP), INTENT(OUT) :: Dn_alpha_Deta_n, Dn_alpha_Deta_p
    REAL(DP), INTENT(OUT) :: Dmu_alpha_Deta_n, Dmu_alpha_Deta_p
    REAL(DP), INTENT(OUT) :: DP_alpha_Deta_n, DP_alpha_Deta_p
    REAL(DP), INTENT(OUT) :: Dn_alpha_DT, DP_alpha_DT, Dmu_alpha_DT

    REAL(DP) :: DP_alpha_Dn_n, DP_alpha_Dn_p

    IF (n_alpha == zero) THEN
      Dn_alpha_Dn_n  = ZERO
      Dn_alpha_Dn_p  = ZERO
      Dmu_alpha_Dn_n = ZERO
      Dmu_alpha_Dn_p = ZERO
      DP_alpha_Dn_n  = ZERO
      DP_alpha_Dn_p  = ZERO
      Dn_alpha_Deta_n  = ZERO
      Dn_alpha_Deta_p  = ZERO
      Dmu_alpha_Deta_n = ZERO
      Dmu_alpha_Deta_p = ZERO
      DP_alpha_Deta_n  = ZERO
      DP_alpha_Deta_p  = ZERO
      Dn_alpha_DT  = ZERO
      Dmu_alpha_DT = ZERO
      DP_alpha_DT  = ZERO
      RETURN
    ENDIF

!   derivatives of alpha particle pressure w.r.t. nucleon densities
    Dmu_alpha_Dn_n = TWO*(Dmu_n_Dn_n+Dmu_p_Dn_n) - DP_Dn_n*v_alpha
    Dmu_alpha_Dn_p = TWO*(Dmu_n_Dn_p+Dmu_p_Dn_p) - DP_Dn_p*v_alpha

!   derivatives of alpha particle pressure w.r.t. nucleon densities
    DP_alpha_Dn_n = n_alpha*( (two-n_n*v_alpha)*Dmu_n_Dn_n + &
                              (two-n_p*v_alpha)*Dmu_p_Dn_n )
    DP_alpha_Dn_p = n_alpha*( (two-n_n*v_alpha)*Dmu_n_Dn_p + &
                              (two-n_p*v_alpha)*Dmu_p_Dn_p )

!   derivatives of alpha particle density w.r.t. nucleon densities
    Dn_alpha_Dn_n = DP_alpha_Dn_n/T
    Dn_alpha_Dn_p = DP_alpha_Dn_p/T

!   derivatives of alpha particle density w.r.t. nucleon degeneracy parameter
    Dn_alpha_Deta_n = Dn_alpha_Dn_n*Dn_n_Deta_n + Dn_alpha_Dn_p*Dn_p_Deta_n
    Dn_alpha_Deta_p = Dn_alpha_Dn_n*Dn_n_Deta_p + Dn_alpha_Dn_p*Dn_p_Deta_p

!   derivatives of alpha particle pressure w.r.t. nucleon degeneracy parameter
    DP_alpha_Deta_n = T*Dn_alpha_Deta_n
    DP_alpha_Deta_p = T*Dn_alpha_Deta_p
!   derivatives of alpha particle chem. pot. w.r.t. nucleon degeneracy parameter
    Dmu_alpha_Deta_n = TWO*(Dmu_n_Deta_n+Dmu_p_Deta_n)-v_alpha*DP_Deta_n
    Dmu_alpha_Deta_p = TWO*(Dmu_n_Deta_p+Dmu_p_Deta_p)-v_alpha*DP_Deta_p

!   derivatives of alpha particle pressure w.r.t. nucleon densities
    Dmu_alpha_DT = TWO*(Dmu_n_DT + Dmu_p_DT) - DP_DT*v_alpha

!   derivative of n_alpha w.r.t. T (fixed eta)
    Dn_alpha_DT = n_alpha/T*(R_3_2 - mu_alpha/T + Dmu_alpha_DT)

!   derivative of P_alpha w.r.t. T (fixed eta)
    DP_alpha_DT = n_alpha + T*Dn_alpha_DT

  END SUBROUTINE ALPHA_DERIVATIVES

END MODULE Alpha_Properties_Mod
