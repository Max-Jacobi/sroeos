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
MODULE Skyrme_Bulk_Observables_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Skyrme_Bulk_Mod
  USE Skyrme_Bulk_Density_Derivatives_Mod

  IMPLICIT NONE

CONTAINS

  FUNCTION SKYRME_BULK_INCOMPRESSIBILITY(log10_n_n,log10_n_p,T) &
    RESULT(INCOMPRESSIBILITY)

    USE Physical_Constants_Mod, &
    ONLY : ZERO, ONE, TWO, FOUR, R_5_3, Hbarc_Square, TWO_PI_SQUARE, &
           Mass_n, Mass_p, Neut_Prot_Mass_Diff

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: log10_n_n,log10_n_p,T
    REAL(DP) :: INCOMPRESSIBILITY
    REAL(DP) :: Meff_n, Meff_p
    REAL(DP) :: eta_n, eta_p
    REAL(DP) :: tau_n, tau_p
    REAL(DP) :: v_n, v_p
    REAL(DP) :: G_n, G_p
    REAL(DP) :: mu_n, mu_p
    REAL(DP) :: n_n, n_p, n, y
    REAL(DP) :: Deta_n_Dn_n, Deta_n_Dn_p
    REAL(DP) :: Deta_p_Dn_n, Deta_p_Dn_p
    REAL(DP) :: Dtau_n_Dn_n, Dtau_n_Dn_p
    REAL(DP) :: Dtau_p_Dn_n, Dtau_p_Dn_p
    REAL(DP) :: Dmu_n_Dn_n, Dmu_n_Dn_p
    REAL(DP) :: Dmu_p_Dn_n, Dmu_p_Dn_p
    REAL(DP) :: Dv_n_Dn_n, Dv_n_Dn_p
    REAL(DP) :: Dv_p_Dn_n, Dv_p_Dn_p
    REAL(DP) :: DP_d_n_n, DP_d_n_p
    REAL(DP) :: Press, dPress_dn, d2f_dn2, Free

    CALL SKYRME_BULK_PROPERTIES(log10_n_n,log10_n_p,T, &
      n_n,n_p,n,y,Meff_n,Meff_p,eta_n,eta_p,tau_n,tau_p,v_n,v_p,mu_n,mu_p, &
      Press,Free)

    CALL SKYRME_BULK_DENSITY_DERIVATIVES(    &
        log10_n_n, log10_n_p, T, G_n, G_p, &
        Deta_n_Dn_n, Deta_n_Dn_p, Deta_p_Dn_n, Deta_p_Dn_p, &
        Dtau_n_Dn_n, Dtau_n_Dn_p, Dtau_p_Dn_n, Dtau_p_Dn_p, &
        Dmu_n_Dn_n, Dmu_n_Dn_p,   Dmu_p_Dn_n, Dmu_p_Dn_p,   &
        Dv_n_Dn_n, Dv_n_Dn_p,     Dv_p_Dn_n, Dv_p_Dn_p,     DP_d_n_n, DP_d_n_p)

    ! Press = SKYRME_BULK_PRESSURE(log10_n_n,log10_n_p,T)

    dPress_dn = (one-y)*DP_d_n_n + y*DP_d_n_p

    d2f_dn2 = (- two*press/n + dPress_dn)/n**two

    INCOMPRESSIBILITY = 9.0_DP*n*n*(d2f_dn2)

  END FUNCTION SKYRME_BULK_INCOMPRESSIBILITY

  ! calculate this numerically cause I'm too lazy to do the derivatives
  FUNCTION SKYRME_BULK_SKEWNESS(log10_n_n,log10_n_p,T) RESULT(SKEWNESS)

    USE Physical_Constants_Mod, &
    ONLY : ZERO, ONE, TWO, FOUR, R_5_3, Hbarc_Square, TWO_PI_SQUARE, &
           Mass_n, Mass_p, Neut_Prot_Mass_Diff

    REAL(DP), INTENT(IN) :: log10_n_n,log10_n_p,T
    REAL(DP) :: log10_n_n_up,log10_n_p_up
    REAL(DP) :: log10_n_n_dn,log10_n_p_dn
    REAL(DP) :: Meff_n, Meff_p
    REAL(DP) :: eta_n, eta_p
    REAL(DP) :: tau_n, tau_p
    REAL(DP) :: v_n, v_p
    REAL(DP) :: mu_n, mu_p
    REAL(DP) :: n_n, n_p, n, y
    REAL(DP) :: d2F_dn2_up, d2F_dn2_dn, d3F_dn3
    REAL(DP) :: SKEWNESS, Press, Free
    REAL(DP), PARAMETER :: eps = 1.0D-4

    CALL SKYRME_BULK_PROPERTIES(log10_n_n,log10_n_p,T, &
      n_n,n_p,n,y,Meff_n,Meff_p,eta_n,eta_p,tau_n,tau_p,v_n,v_p,mu_n,mu_p, &
      Press,Free)

    log10_n_n_up = log10(n_n*(one+eps))
    log10_n_p_up = log10(n_p*(one+eps))

    d2F_dn2_up = SKYRME_BULK_INCOMPRESSIBILITY(log10_n_n_up,log10_n_p_up,T)/&
                 (3.0_DP*n*(one+eps))**two

    log10_n_n_dn = log10(n_n*(one-eps))
    log10_n_p_dn = log10(n_p*(one-eps))

    d2F_dn2_dn = SKYRME_BULK_INCOMPRESSIBILITY(log10_n_n_dn,log10_n_p_dn,T)/&
                 (3.0_DP*n*(one-eps))**two

    d3f_dn3 = (d2F_dn2_up - d2F_dn2_dn)/(two*eps*n)

    SKEWNESS = 27.D0*n*n*n*d3f_dn3

  END FUNCTION SKYRME_BULK_SKEWNESS

END MODULE Skyrme_Bulk_Observables_Mod
