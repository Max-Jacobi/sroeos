!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it AND/or modIFy
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
!    along with SRO_EOS.  IF not, see <http://www.gnu.org/licenses/>.
!
MODULE Surface_Tension_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, HALF
  USE Skyrme_Coefficients_Mod, ONLY : Coeff_qnn, Coeff_qnp, Coeff_qpn, Coeff_qpp
  USE Skyrme_Bulk_Mod, ONLY : SKYRME_BULK_PROPERTIES, SKYRME_BULK_PRESSURE

  IMPLICIT NONE

CONTAINS

  FUNCTION Woods_Saxon(z,rt,at,n_ti,n_to) RESULT(WS)
!*****************************************************************************80
!
!! g_surf calculates the surface thermodinamic potential
!   considerina a woods-saxon density distribution for
!   neutron and proton densities with radii R_n, R_p
!   and thickness a_n and a_p, i.e.,
!
!   n_t(r) = n_ti + (n_to-n_ti)/(1+exp((r-R_t)/a_t)), where t = n, p
!

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: z,rt,at,n_ti,n_to
    REAL(DP) :: den, WS

    den = exp((z-rt)/at)
    if ((z-rt)/at> 200.d0) den = exp( 200.d0)
    if ((z-rt)/at<-200.d0) den = exp(-200.d0)
    WS = n_to + (n_ti - n_to)/(one+den)

  END FUNCTION Woods_Saxon

  FUNCTION Surface_Tension (r_neut, r_prot, a_neut, a_prot, &
        log10_dens_neut_in, log10_dens_prot_in, log10_dens_neut_out, &
        log10_dens_prot_out, temperature) RESULT(thermo_potential)
!*****************************************************************************80
!
!! g_surf calculates the surface thermodinamic potential
!   and from it and neutron excess obtains the surface free energy.
!   We consider woods-saxon type density distributions for both
!   neutron and proton densities with radii R_n, R_p
!   and thickness a_n and a_p, i.e.,
!   n_t(r) = n_ti + (n_to-n_ti)/(1+exp((r-R_t)/a_t))
!    where t = n, p
!
    IMPLICIT NONE

    INTEGER(I4B) :: i, imax
    REAL(DP), INTENT(IN) :: r_neut, r_prot, a_neut, a_prot
    REAL(DP), INTENT(IN) :: log10_dens_neut_in, log10_dens_prot_in
    REAL(DP), INTENT(IN) :: log10_dens_neut_out, log10_dens_prot_out
    REAL(DP), INTENT(IN) :: temperature
    REAL(DP) :: thermo_potential
    REAL(DP) :: press, press_in, press_out
    REAL(DP) :: free_energy, free_energy_in, free_energy_out
    REAL(DP) :: log10_dens_neut, log10_dens_prot
    REAL(DP) :: mu_neut,     mu_prot,     dens_neut,     dens_prot
    REAL(DP) :: mu_neut_in,  mu_prot_in,  dens_neut_in,  dens_prot_in
    REAL(DP) :: mu_neut_out, mu_prot_out, dens_neut_out, dens_prot_out
    REAL(DP) :: dens_nucl,prot_frac
    REAL(DP) :: effective_mass_neut,effective_mass_prot
    REAL(DP) :: eta_neut,eta_prot,tau_neut,tau_prot,v_neut,v_prot
    REAL(DP) :: d_dens_neut_out_dz, d_dens_prot_out_dz
    REAL(DP) :: z, z0, dz, eps, tol, q_tt

    EPS = 1.D-4
    TOL = 1.D-6

!   OBTAIN OUTSIDE CHEMICAL POTENTIAL AND PRESSURE
    CALL SKYRME_BULK_PROPERTIES(log10_dens_neut_out,log10_dens_prot_out,   &
              Temperature,dens_neut_out,dens_prot_out,dens_nucl,prot_frac, &
              effective_mass_neut,effective_mass_prot,eta_neut,eta_prot,   &
              tau_neut,tau_prot,v_neut,v_prot,mu_neut_out,mu_prot_out,     &
              press_out,free_energy_out)

    ! press_out = SKYRME_BULK_PRESSURE(log10_dens_neut_out, &
    !                                  log10_dens_prot_out,Temperature)
    !
    ! FREE_ENERGY_OUT = mu_neut_out*dens_neut_out &
    !                 + mu_prot_out*dens_prot_out - press_out

!   OBTAIN INSIDE CHEMICAL POTENTIAL AND PRESSURE
    CALL SKYRME_BULK_PROPERTIES(log10_dens_neut_in,log10_dens_prot_in,   &
              Temperature,dens_neut_in,dens_prot_in,dens_nucl,prot_frac, &
              effective_mass_neut,effective_mass_prot,eta_neut,eta_prot, &
              tau_neut,tau_prot,v_neut,v_prot,mu_neut_in,mu_prot_in,     &
              press_in,free_energy_in)

    ! press_in  = SKYRME_BULK_PRESSURE(log10_dens_neut_in, &
    !                                  log10_dens_prot_in,Temperature)
    !
    ! FREE_ENERGY_IN  = mu_neut_in*dens_neut_in &
    !                 + mu_prot_in*dens_prot_in - press_in

!   START INTEGRATING ONCE
    DO i = -100, -1
      z = DBLE(i)
      dens_neut = Woods_Saxon(Z,r_neut,a_neut,dens_neut_in,dens_neut_out)
      IF ((dens_neut_in-dens_neut)/dens_neut_in<tol) Z0 = Z
    ENDDO

!   NUMBER OF POINTS USED IN THE INTEGRAL
    imax = 1000
!   SET TIME STEPS
    dz = -TWO*z0/DBLE(imax)
!   SET LOCATION DETERMINED ABOVE AS FIRST STEP
    i = 0
!   SET DENSITY TO INNER DENSITY TO
!   FIND OUT LAST INTEGRATION POINT
    dens_neut = dens_neut_in

!   SET INTEGRAL TO ZERO
    thermo_potential = ZERO

!   PERFORM INTEGRAL
    DO WHILE ( (dens_neut-dens_neut_out)/dens_neut_out > tol)

!     OBTAIN POSITION
      z = z0 + dz*DBLE(i)
      i = i + 1

!     OBTAIN DENSITIES
      dens_neut = Woods_Saxon(Z,r_neut,a_neut,dens_neut_in,dens_neut_out)
      dens_prot = Woods_Saxon(Z,r_prot,a_prot,dens_prot_in,dens_prot_out)

      dens_nucl = dens_neut + dens_prot

      log10_dens_neut = log10(dens_neut)
      log10_dens_prot = log10(dens_prot)

!     OBTAIN DENSITY DERIVATIVES NUMERICALLY
      d_dens_neut_out_dz = &
            (Woods_Saxon(z+eps,r_neut,a_neut,dens_neut_in,dens_neut_out)  &
         -   Woods_Saxon(z-eps,r_neut,a_neut,dens_neut_in,dens_neut_out)) &
            /(two*eps)
      d_dens_prot_out_dz = &
            (Woods_Saxon(z+eps,r_prot,a_prot,dens_prot_in,dens_prot_out)  &
         -   Woods_Saxon(z-eps,r_prot,a_prot,dens_prot_in,dens_prot_out)) &
            /(two*eps)
!   OBTAIN CHEMICAL POTENTIAL AND PRESSURE AT Z
      CALL SKYRME_BULK_PROPERTIES(log10_dens_neut,log10_dens_prot,          &
                Temperature,dens_neut,dens_prot,dens_nucl,prot_frac,        &
                effective_mass_neut,effective_mass_prot,eta_neut,eta_prot,  &
                tau_neut,tau_prot,v_neut,v_prot,mu_neut,mu_prot,            &
                press,free_energy)

      ! press = SKYRME_BULK_PRESSURE(log10_dens_neut,log10_dens_prot,Temperature)
      !
      ! free_energy = mu_neut*dens_neut + mu_prot*dens_prot - press

      Q_TT = HALF*( Coeff_qnn*d_dens_neut_out_dz*d_dens_neut_out_dz &
                  + Coeff_qnp*d_dens_neut_out_dz*d_dens_prot_out_dz &
                  + Coeff_qpn*d_dens_prot_out_dz*d_dens_neut_out_dz &
                  + Coeff_qpp*d_dens_prot_out_dz*d_dens_prot_out_dz )

      thermo_potential = thermo_potential &
                        + ( free_energy + Q_TT - FREE_ENERGY_OUT &
                        - mu_neut_out*(dens_neut-dens_neut_out)   &
                        - mu_prot_out*(dens_prot-dens_prot_out) )*DZ

    ENDDO

  END FUNCTION Surface_Tension

END MODULE Surface_Tension_Mod
