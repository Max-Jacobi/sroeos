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
MODULE Meff_bulk_mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  !USE Define_Operators_Mod
  USE Skyrme_Coefficients_Mod
  USE Fermi_Integrals_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE MEFF_BULK(log10_n_n,log10_n_p,Temperature, &
        Meff_n, Meff_p, &
        DMeff_n_Dn_n, DMeff_n_Dn_p, &
        DMeff_p_Dn_n, DMeff_p_Dn_p &
        DMeff_n_D2n_nn, DMeff_n_D2n_np, DMeff_n_D2n_pp &
        DMeff_p_D2n_nn, DMeff_p_D2n_np, DMeff_p_D2n_pp &
      )
    USE Physical_Constants_Mod, &
        ONLY : ZERO, ONE, TWO, FOUR, TEN, R_3_2, R_5_2, &
               R_1_3, R_2_3, R_4_3, R_5_3, &
               Hbarc_Square, PI, TWO_PI_SQUARE, &
               log10_dens_min, log10_dens_max, dens_min, dens_max, &
               Mass_n, Mass_p, Neut_Prot_Mass_Diff
    USE, INTRINSIC :: IEEE_ARITHMETIC

    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: log10_n_n, log10_n_p, Temperature
    REAL(DP), INTENT(OUT) :: Meff_n, Meff_p
    REAL(DP), INTENT(OUT), Optional :: DMeff_n_Dn_n, DMeff_n_Dn_p
    REAL(DP), INTENT(OUT), Optional :: DMeff_p_Dn_n, DMeff_p_Dn_p
    REAL(DP), INTENT(OUT), Optional :: DMeff_n_D2n_nn, DMeff_n_D2n_pp
    REAL(DP), INTENT(OUT), Optional :: DMeff_n_D2n_np
    REAL(DP), INTENT(OUT), Optional :: DMeff_p_D2n_nn, DMeff_p_D2n_pp
    REAL(DP), INTENT(OUT), Optional :: DMeff_p_D2n_np
    REAL(DP) :: n_n, n_p, n, y
    REAL(DP) :: eff_n_n, eff_n_p

    REAL(DP) :: n_n_m2_3, n_n_m1_3
    REAL(DP) :: n_n_1_3, n_n_2_3, n_n_4_3, n_n_5_3
    REAL(DP) :: n_p_1_3, n_p_2_3, n_p_4_3, n_p_5_3

    REAL(DP) :: coeff_term_n, coeff_term_p
    REAL(DP) :: sig_fac, eps_fac, eps_term, delta_eps

    REAL(DP) :: Dcoeff_term_n_Dn_n, Dcoeff_term_n_Dn_p
    REAL(DP) :: Dcoeff_term_p_Dn_n, Dcoeff_term_p_Dn_p
    REAL(DP) :: Dsig_fac_Dn, Deps_fac_Dn
    REAL(DP) :: Deps_term_Dn_n, Deps_term_Dn_p

    REAL(DP) :: Dcoeff_term_n_D2n_nn, Dcoeff_term_n_D2n_pp
    REAL(DP) :: Dcoeff_term_p_D2n_nn, Dcoeff_term_p_D2n_pp
    REAL(DP) :: Dsig_fac_D2n, Deps_fac_D2n
    REAL(DP) :: Deps_term_D2n_nn, Deps_term_D2n_np, Deps_term_D2n_pp

    LOGICAL :: calc_first_derivative, calc_second_derivative

    calc_first_derivative = &
        PRESENT(DMeff_n_Dn_n) .OR. PRESENT(DMeff_n_Dn_p) .OR. &
        PRESENT(DMeff_p_Dn_n) .OR. PRESENT(DMeff_p_Dn_p)
    calc_second_derivative = &
        PRESENT(DMeff_n_D2n_nn) .OR. PRESENT(DMeff_n_D2n_np) .OR. &
        PRESENT(DMeff_n_D2n_pp) .OR. PRESENT(DMeff_p_D2n_nn) .OR. &
        PRESENT(DMeff_p_D2n_np) .OR. PRESENT(DMeff_p_D2n_pp)

!   neutron and proton effective densities
!    used to avoind underflow/overflow in calculations below
    eff_n_n = max(log10_n_n,log10_dens_min)
    eff_n_n = min(eff_n_n,log10_dens_max)
    eff_n_n = ten**eff_n_n

    eff_n_p = max(log10_n_p,log10_dens_min)
    eff_n_p = min(eff_n_p,log10_dens_max)
    eff_n_p = ten**eff_n_p

    n_n = eff_n_n
    n_p = eff_n_p

    n = n_n + n_p
    y = n_p/n

!   fix nuclear density and proton fractions in extreme cases
    IF (log10_n_n>log10_n_p+16.D0) THEN
      n = n_n
      y = n_p/n
      IF (ieee_is_nan(y)) y = zero
      n_p = eff_n_p
    ELSEIF (log10_n_n<log10_n_p-16.D0) THEN
      n = n_p
      y = one
      n_n = eff_n_n
    ENDIF

    n_n_m2_3 = n_n**(-R_2_3)
    n_n_m1_3 = n_n**(-R_1_3)
    n_n_1_3 = n_n**(R_1_3)
    n_n_2_3 = n_n**(R_2_3)
    n_n_4_3 = n_n**(R_4_3)
    n_n_5_3 = n_n**(R_5_3)
    n_p_m2_3 = n_p**(-R_2_3)
    n_p_m1_3 = n_p**(-R_1_3)
    n_p_1_3 = n_p**(R_1_3)
    n_p_2_3 = n_p**(R_2_3)
    n_p_4_3 = n_p**(R_4_3)
    n_p_5_3 = n_p**(R_5_3)

    delta_eps = coeff_eps_p - coeff_eps_n

    sig_fac = (1 + exp(5*n))**(-1)
    eps_fac = (1 - exp(-10*n)) / (1 + exp(-5*(n - coeff_n_off)))
    eps_term = coeff_eps_n + y * delta_eps - 1

    ceff_term_n = coeff_alpha1*n_n + coeff_beta1*n_p &
           + coeff_alpha2*n_n_4_3 + coeff_beta2*n_p_4_3 &
           + coeff_alpha3*n_n_5_3 + coeff_beta3*n_p_5_3

    ceff_term_p = coeff_alpha1*n_p + coeff_beta1*n_n &
           + coeff_alpha2*n_p_4_3 + coeff_beta2*n_n_4_3 &
           + coeff_alpha3*n_p_5_3 + coeff_beta3*n_n_5_3

!   neutron and proton effective masses as in Huth et al. 2021
    Meff_n = 1 + coeff_term_n*sig_fac + eps_term*eps_fac
    Meff_n = Meff_n * Mass_n

    Meff_p = 1 + coeff_term_p*sig_fac + eps_term*eps_fac
    Meff_p = Meff_p * Mass_p

!   1st derivatives
    if (.not. (calc_first_derivative .or. calc_second_derivative)) return

    Dcoeff_term_n_Dn_n = coeff_alpha1 &
          + coeff_alpha2*R_4_3*n_n_1_3 &
          + coeff_alpha3*R_5_3*n_n_2_3

    Dcoeff_term_p_Dn_p = coeff_alpha1 &
          + coeff_alpha2*R_4_3*n_p_1_3 &
          + coeff_alpha3*R_5_3*n_p_2_3

    Dcoeff_term_n_Dn_p = coeff_beta1 &
          + coeff_beta2*R_4_3*n_p_1_3 &
          + coeff_beta3*R_5_3*n_p_2_3

    Dcoeff_term_p_Dn_n = coeff_beta1 &
          + coeff_beta2*R_4_3*n_n_1_3 &
          + coeff_beta3*R_5_3*n_n_2_3

    Dsig_fac_Dn = 5*sig_fac**2 - 5*sig_fac

    Deps_fac_Dn = (&
          5*exp(-5*n) ( &
            exp(5*coeff_n_off) &
            + 2*exp(5*n) + exp(5*(coeff_n_off &
            + 2*n))
          )
        ) &
        / (exp(5*coeff_n_off) + exp(5*n))**2


    Deps_term_Dn_n = -y/n * delta_eps
    Deps_term_Dn_p = (1-y)/n * delta_eps

    if (PRESENT(DMeff_n_Dn_n)) then
        DMeff_n_Dn_n = Dcoeff_term_n_Dn_n * sig_fac &
            + coeff_term_n * Dsig_fac_Dn &
            + Deps_term_Dn_n*eps_fac &
            + eps_term*Deps_fac_Dn
        DMeff_n_Dn_n = DMeff_n_Dn_n * Mass_n
    endif

    if (PRESENT(DMeff_n_Dn_p)) then
       DMeff_n_Dn_p = Dcoeff_term_n_Dn_p * sig_fac &
           + coeff_term_n * Dsig_fac_Dn &
           + Deps_term_Dn_p*eps_fac &
           + eps_term*Deps_fac_Dn
        DMeff_n_Dn_p = DMeff_n_Dn_p * Mass_n
    endif

    if (PRESENT(DMeff_p_Dn_p)) then
        DMeff_p_Dn_n = Dcoeff_term_p_Dn_n * sig_fac &
            + coeff_term_p * Dsig_fac_Dn &
            + Deps_term_Dn_n*eps_fac &
            + eps_term*Deps_fac_Dn
        DMeff_p_Dn_p = DMeff_p_Dn_p * Mass_p
    endif

    if (PRESENT(DMeff_p_Dn_p)) then
        DMeff_p_Dn_p = Dcoeff_term_p_Dn_p * sig_fac &
            + coeff_term_p * Dsig_fac_Dn &
            + Deps_term_Dn_p*eps_fac &
            + eps_term*Deps_fac_Dn
        DMeff_p_Dn_p = DMeff_p_Dn_p * Mass_p
    endif

    ! 2nd derivatives
    if (.not. calc_second_derivative) return

    Dcoeff_term_n_D2n_nn = coeff_alpha2*R_4_3*R_1_3*n_n_m2_3 &
        + coeff_alpha3*R_5_3*R_2_3*n_n_m1_3

    Dcoeff_term_p_D2n_pp = coeff_alpha2*R_4_3*R_1_3*n_p_m2_3 &
        + coeff_alpha3*R_5_3*R_2_3*n_p_m1_3

    Dcoeff_term_n_D2n_pp = coeff_beta2*R_4_3*R_1_3*n_p_m2_3 &
        + coeff_beta3*R_5_3*R_2_3*n_p_m1_3

    Dcoeff_term_p_D2n_nn = coeff_beta2*R_4_3*R_1_3*n_n_m2_3 &
        + coeff_beta3*R_5_3*R_2_3*n_n_m1_3

    Dsig_fac_D2n = (10*sig_fac - 5)*Dsig_fac_Dn

    Deps_fac_D2n = -(25*exp(-5*n)*( &
        exp(10*coeff_n_off) + 4*exp(10*n) &
        + 3*exp(5*(coeff_n_off + n)) &
        - exp(10*(coeff_n_off + n)) &
        + exp(5*(coeff_n_off + 3*n)) &
        )) &
        /(exp(5*coeff_n_off) + exp(5*n))**3

    Deps_term_D2n_nn= delta_eps * 2*y/n**2
    Deps_term_D2n_pp= delta_eps * 2*(y-1)/n**2
    Deps_term_D2n_np= delta_eps * (2y-1)/n**2


    if (PRESENT(DMeff_n_D2n_nn)) then
        DMeff_n_D2n_nn = Dcoeff_term_n_D2n_nn * sig_fac &
            + 2 * Dcoeff_term_n_Dn_n * Dsig_fac_Dn &
            + coeff_term_n * Dsig_fac_D2n &
            + Deps_term_D2n_nn * eps_fac &
            + 2 * Deps_term_Dn_n * Deps_fac_Dn &
            + eps_term * Deps_fac_D2n
        DMeff_n_D2n_nn = DMeff_n_D2n_nn * Mass_n
    endif

    if (PRESENT(DMeff_n_D2n_np)) then
        DMeff_n_D2n_np =  &
            + Dcoeff_term_n_Dn_n * Dsig_fac_Dn &
            + Dcoeff_term_n_Dn_p * Dsig_fac_Dn &
            + coeff_term_n * Dsig_fac_D2n &
            + Deps_term_D2n_np * eps_fac &
            + (Deps_term_Dn_n + Deps_term_Dn_p) * Deps_fac_Dn &
            + eps_term * Deps_fac_D2n
        DMeff_n_D2n_np = DMeff_n_D2n_np * Mass_n
    endif

    if (PRESENT(DMeff_p_D2n_pp)) then
        DMeff_n_D2n_pp = Dcoeff_term_n_D2n_pp * sig_fac &
            + 2 * Dcoeff_term_n_Dn_p * Dsig_fac_Dn &
            + coeff_term_n * Dsig_fac_D2n &
            + Deps_term_D2n_pp * eps_fac &
            + 2 * Deps_term_Dn_p * Deps_fac_Dn &
            + eps_term * Deps_fac_D2n
        DMeff_n_D2n_pp = DMeff_n_D2n_pp * Mass_n
     endif

    if (PRESENT(DMeff_p_D2n_nn)) then
        DMeff_p_D2n_nn = Dcoeff_term_p_D2n_nn * sig_fac &
            + 2 * Dcoeff_term_p_Dn_n * Dsig_fac_Dn &
            + coeff_term_p * Dsig_fac_D2n &
            + Deps_term_D2n_nn * eps_fac &
            + 2 * Deps_term_Dn_n * Deps_fac_Dn &
            + eps_term * Deps_fac_D2n
        DMeff_p_D2n_nn = DMeff_p_D2n_nn * Mass_p
    endif

    if (PRESENT(DMeff_p_D2n_np)) then
        DMeff_p_D2n_np =  &
            + Dcoeff_term_p_Dn_n * Dsig_fac_Dn &
            + Dcoeff_term_p_Dn_p * Dsig_fac_Dn &
            + coeff_term_p * Dsig_fac_D2n &
            + Deps_term_D2n_np * eps_fac &
            + (Deps_term_Dn_n + Deps_term_Dn_p) * Deps_fac_Dn &
            + eps_term * Deps_fac_D2n
        DMeff_p_D2n_np = DMeff_p_D2n_np * Mass_p
     endif

     if (PRESENT(DMeff_p_D2n_pp)) then
        DMeff_p_D2n_pp = Dcoeff_term_p_D2n_pp * sig_fac &
            + 2 * Dcoeff_term_p_Dn_p * Dsig_fac_Dn &
            + coeff_term_p * Dsig_fac_D2n &
            + Deps_term_D2n_pp * eps_fac &
            + 2 * Deps_term_Dn_p * Deps_fac_Dn &
            + eps_term * Deps_fac_D2n
        DMeff_p_D2n_pp = DMeff_p_D2n_pp * Mass_p
      endif

    END SUBROUTINE MEFF_BULK
END MODULE MEFF_BULK_MOD
