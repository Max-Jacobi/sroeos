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
MODULE Skyrme_Bulk_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  !USE Define_Operators_Mod
  USE Skyrme_Coefficients_Mod
  USE Fermi_Integrals_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE SKYRME_BULK_PROPERTIES(log10_n_n,log10_n_p,Temperature, &
    n_n,n_p,n,y,Meff_n,Meff_p,eta_n,eta_p,tau_n,tau_p,v_n,v_p,mu_n,mu_p,P,F)

    USE Physical_Constants_Mod, &
        ONLY : ZERO, ONE, TWO, FOUR, TEN, R_3_2, R_5_2, R_5_3, &
               Hbarc_Square, PI, TWO_PI_SQUARE, &
               log10_dens_min, log10_dens_max, dens_min, dens_max, &
               Mass_n, Mass_p, Neut_Prot_Mass_Diff
    USE, INTRINSIC :: IEEE_ARITHMETIC

    IMPLICIT NONE

    REAL(DP), INTENT(IN)  :: log10_n_n, log10_n_p, Temperature
    REAL(DP), INTENT(OUT) :: n_n, n_p, n, y
    REAL(DP), INTENT(OUT) :: Meff_n, Meff_p
    REAL(DP), INTENT(OUT) :: eta_n, eta_p
    REAL(DP), INTENT(OUT) :: tau_n, tau_p
    REAL(DP), INTENT(OUT) :: v_n, v_p
    REAL(DP), INTENT(OUT) :: mu_n, mu_p
    REAL(DP), INTENT(OUT) :: P, F
    REAL(DP) :: eff_n_n, eff_n_p
    REAL(DP) :: two_mneut_t, two_mprot_t

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

!   neutron and proton effective masses
    Meff_n = (Hbarc_Square/two) / (Hbarc_Square/(two*Mass_n) &
                          + coeff_alpha1*n_n + coeff_alpha2*n_p)
    Meff_p = (Hbarc_Square/two) / (Hbarc_Square/(two*Mass_p) &
                          + coeff_alpha1*n_p + coeff_alpha2*n_n)

!   neutron and proton degeneracy parameters eta
    two_mneut_t = two*Meff_n*Temperature
    two_mprot_t = two*Meff_p*Temperature

    eta_n = TWO_PI_SQUARE*n_n*(hbarc_square/two_mneut_t)**R_3_2
    eta_n = inverse_fermi_one_half(eta_n)
    eta_p = TWO_PI_SQUARE*n_p*(hbarc_square/two_mprot_t)**R_3_2
    eta_p = inverse_fermi_one_half(eta_p)

!   if dens_t < dens_min then
!    correct value for degeneracy parameters
!    using relevant asymptotic formulae
!     F_-1/2(x) = sqrt(pi)*exp(x)
!     F_+1/2(x) = (1/2)*sqrt(pi)*exp(x)
!     F_+3/2(x) = (3/4)*sqrt(pi)*exp(x)
!     F_+5/2(x) = (15/8)*sqrt(pi)*exp(x)
      if (log10_n_n<log10_dens_min) then
        eta_n = ((two_mneut_t/hbarc_square)**R_3_2)/TWO_PI_SQUARE
        eta_n = log10_n_n*log(ten)-log(sqrt(pi)*eta_n/two)
      endif
      if (log10_n_p<log10_dens_min) then
        eta_p = ((two_mprot_t/hbarc_square)**R_3_2)/TWO_PI_SQUARE
        eta_p = log10_n_p*log(ten)-log(sqrt(pi)*eta_p/two)
      endif

!   neutron and proton kinetic energies
!    cap tau_t to a min value of -200.d0
!    to avoid underflow issues
    tau_n = (two_mneut_t/hbarc_square)**R_5_2
    tau_n = tau_n*fermi_three_halves(max(eta_n,-2.d2))/TWO_PI_SQUARE
    tau_p = (two_mprot_t/hbarc_square)**R_5_2
    tau_p = tau_p*fermi_three_halves(max(eta_p,-2.d2))/TWO_PI_SQUARE

!   neutron and proton force potentials
    v_n = coeff_alpha1*tau_n + coeff_alpha2*tau_p &
          + two*coeff_a*n + four*coeff_b*n_p &
          + dot_product(coeff_c*(one+coeff_delta),n**coeff_delta) &
          + dot_product(four*coeff_d*(coeff_delta-one)*n_n*n_p,&
            n**(coeff_delta-two)) &
          + dot_product(four*coeff_d*n_p,n**(coeff_delta-one))

    v_p = coeff_alpha1*tau_p + coeff_alpha2*tau_n &
          + two*coeff_a*n + four*coeff_b*n_n &
          + dot_product(coeff_c*(one+coeff_delta),n**coeff_delta) &
          + dot_product(four*coeff_d*(coeff_delta-one)*n_n*n_p, &
            n**(coeff_delta-two)) &
          + dot_product(four*coeff_d*n_n,n**(coeff_delta-one)) &
          - Neut_Prot_Mass_Diff

!   neutron and proton chemical potentials
    mu_n = eta_n*Temperature+v_n
    mu_p = eta_p*Temperature+v_p

!   pressure of bulk nucleons
    P = Hbarc_Square/TWO*(R_5_3/Meff_n-ONE/Mass_n)*tau_n &
      + Hbarc_Square/TWO*(R_5_3/Meff_p-ONE/Mass_p)*tau_p &
      + coeff_a*n**two + four*coeff_b*n_p*n_n &
      + dot_product(coeff_c*coeff_delta,n**(one+coeff_delta)) &
      +four*dot_product(coeff_d*coeff_delta,n**(coeff_delta-one))*n_n*n_p

!   free energy of bulk nucleons
    F = n_n*mu_n+n_p*mu_p-P

  END SUBROUTINE SKYRME_BULK_PROPERTIES

  FUNCTION SKYRME_BULK_ENERGY(log10_n_n,log10_n_p,Temperature) &
    RESULT(ENERGY)

    USE Physical_Constants_Mod, &
    ONLY : ONE, TWO, FOUR, Hbarc_Square, &
           Mass_n, Mass_p, Neut_Prot_Mass_Diff

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: log10_n_n,log10_n_p,Temperature
    REAL(DP) :: ENERGY, P, F
    REAL(DP) :: Meff_n, Meff_p
    REAL(DP) :: eta_n, eta_p
    REAL(DP) :: tau_n, tau_p
    REAL(DP) :: v_n, v_p
    REAL(DP) :: mu_n, mu_p
    REAL(DP) :: n_n, n_p, n, y

    CALL SKYRME_BULK_PROPERTIES(log10_n_n,log10_n_p,Temperature,&
      n_n,n_p,n,y,Meff_n,Meff_p,eta_n,eta_p,tau_n,tau_p,v_n,v_p,mu_n,mu_p,P,F)

    ENERGY = hbarc_square/two* (tau_n/Meff_n+tau_p/Meff_p) &
    + coeff_a*n**two + four*coeff_b*n_n*n_p &
    + DOT_PRODUCT(coeff_c,n**(coeff_delta+one)) &
    + four*DOT_PRODUCT(coeff_d,n**(coeff_delta-one))*n_n*n_p

  END FUNCTION SKYRME_BULK_ENERGY

  FUNCTION SKYRME_BULK_PRESSURE(log10_n_n,log10_n_p,Temperature) &
    RESULT(PRESSURE)

    USE Physical_Constants_Mod, &
    ONLY : ZERO, ONE, TWO, FOUR, R_5_3, Hbarc_Square, TWO_PI_SQUARE, &
           Mass_n, Mass_p, Neut_Prot_Mass_Diff

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: log10_n_n,log10_n_p,Temperature
    REAL(DP) :: PRESSURE, P, F
    REAL(DP) :: Meff_n, Meff_p
    REAL(DP) :: eta_n, eta_p
    REAL(DP) :: tau_n, tau_p
    REAL(DP) :: v_n, v_p
    REAL(DP) :: mu_n, mu_p
    REAL(DP) :: n_n, n_p, n, y

    CALL SKYRME_BULK_PROPERTIES(log10_n_n,log10_n_p,Temperature,&
      n_n,n_p,n,y,Meff_n,Meff_p,eta_n,eta_p,tau_n,tau_p,v_n,v_p,mu_n,mu_p,P,F)

    PRESSURE = Hbarc_Square/TWO*(R_5_3/Meff_n-ONE/Mass_n)*tau_n &
             + Hbarc_Square/TWO*(R_5_3/Meff_p-ONE/Mass_p)*tau_p &
             + coeff_a*n**two + four*coeff_b*n_p*n_n &
             + dot_product(coeff_c*coeff_delta,n**(one+coeff_delta)) &
             +four*dot_product(coeff_d*coeff_delta,n**(coeff_delta-one))*n_n*n_p

  END FUNCTION SKYRME_BULK_PRESSURE

END MODULE Skyrme_Bulk_Mod
