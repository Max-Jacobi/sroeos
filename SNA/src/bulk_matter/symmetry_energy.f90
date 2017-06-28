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
MODULE Symmetry_Energy_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  !USE Define_Operators_Mod
  USE Skyrme_Coefficients_Mod

CONTAINS

! for now only works for dens_neut=dens_prot=dens_nucl_sat/two
  SUBROUTINE SKYRME_SYMMETRY_ENERGY(J_sym,J_half,L_sym,K_sym,Q_sym)

    USE Physical_Constants_Mod, &
    ONLY : ZERO, ONE, TWO, THREE, FOUR, TEN, &
           R_1_3, R_2_3, R_5_3, R_5_9, R_5_18, R_10_9, &
           Hbarc_Square, TWO_PI_SQUARE, C_TAU, &
           Mass_n, Mass_p, Neut_Prot_Mass_Diff
    USE Nuclear_Matter_Properties_Mod, ONLY : Nuc_Sat_Dens
    USE Make_tables_Mod, ONLY : output_directory

    IMPLICIT NONE

    INTEGER(I4B) :: i
    REAL(DP), INTENT(OUT) :: J_sym,J_half,L_sym,K_sym,Q_sym
    REAL(DP) :: log10_dens_neut_up,log10_dens_prot_up
    REAL(DP) :: log10_dens_neut_dn,log10_dens_prot_dn
    REAL(DP) :: effective_Mass_n, effective_Mass_p
    REAL(DP) :: eta_neut, eta_prot
    REAL(DP) :: tau_neut, tau_prot
    REAL(DP) :: v_neut, v_prot
    REAL(DP) :: mu_neut, mu_prot
    REAL(DP) :: dens_neut, dens_prot, dens_nucl, prot_frac
    REAL(DP) :: dummy_coeff, x0
    REAL(DP), DIMENSION(:), ALLOCATABLE :: dummy, expt
    REAL(DP), PARAMETER :: eps = 1.0D-4
    CHARACTER(LEN=64) :: sym_file

    sym_file = trim(output_directory) // '/symmetry.dat'

    ALLOCATE(expt(size(coeff_delta)))
    ALLOCATE(dummy(size(coeff_d)))

    dummy_coeff = (two*coeff_alpha1-coeff_alpha2)

    open (11,file=sym_file)
    do i = 1, 10000
      dens_nucl = dble(i)*Nuc_Sat_Dens/1.d3
      dens_prot = dens_nucl/two
      dens_neut = dens_nucl/two

      dummy = coeff_d
      expt  = coeff_delta

      J_sym = + R_5_18*hbarc_square*c_tau/Mass_n*(dens_neut)**R_2_3/two &
              + R_5_18*hbarc_square*c_tau/Mass_p*(dens_prot)**R_2_3/two &
              + R_10_9*c_tau*dummy_coeff*(dens_nucl/two)**R_5_3         &
              - coeff_b*dens_nucl - DOT_PRODUCT(dummy,dens_nucl**expt)

      dummy = dummy*coeff_delta

      L_sym = R_5_9*hbarc_square*c_tau/Mass_n*(dens_neut)**R_2_3/two &
            + R_5_9*hbarc_square*c_tau/Mass_p*(dens_prot)**R_2_3/two &
            + TEN*R_5_9*c_tau*dummy_coeff*(dens_nucl/two)**R_5_3     &
            - three*coeff_b*dens_nucl                                &
            - three*DOT_PRODUCT(dummy,dens_nucl**expt)

      dummy = dummy*(coeff_delta-one)

      K_sym = - R_5_9*hbarc_square*c_tau/Mass_n*(dens_neut)**R_2_3/two &
              - R_5_9*hbarc_square*c_tau/Mass_p*(dens_prot)**R_2_3/two &
              + (100.d0/9.d0)*c_tau*dummy_coeff*(dens_nucl/two)**R_5_3 &
              - 9.d0*DOT_PRODUCT(dummy,dens_nucl**expt)

      dummy = dummy*(coeff_delta-two)


      Q_sym = + (20.d0/9.d0)*hbarc_square*c_tau/Mass_n*(dens_neut)**R_2_3/two &
              + (20.d0/9.d0)*hbarc_square*c_tau/Mass_p*(dens_prot)**R_2_3/two &
              - (100.d0/9.d0)*c_tau*dummy_coeff*(dens_nucl/two)**R_5_3        &
              - 27.d0*DOT_PRODUCT(dummy,dens_nucl**expt)

      x0 = (dens_nucl - Nuc_Sat_Dens)/THREE

      write (11,"(2es20.12)") dble(i)/1.D3, &
             J_sym + L_sym*x0 + K_sym*x0*x0/2.D0 + Q_sym*x0*x0*x0/6.D0
    enddo
    close(11)

    dens_nucl = Nuc_Sat_Dens

    ! CALL SKYRME_BULK_PROPERTIES(log10_dens_neut,log10_dens_prot,Temperature,&
    ! dens_neut,dens_prot,dens_nucl,prot_frac,   &
    ! effective_Mass_n,effective_Mass_p,eta_neut,eta_prot, &
    ! tau_neut,tau_prot,v_neut,v_prot,mu_neut,mu_prot)

!   symmetry energy
    J_sym = + R_5_18*hbarc_square*c_tau/Mass_n*(dens_nucl/two)**R_2_3/two &
            + R_5_18*hbarc_square*c_tau/Mass_p*(dens_nucl/two)**R_2_3/two &
            + R_10_9*c_tau*dummy_coeff*(dens_nucl/two)**R_5_3 &
            - coeff_b*dens_nucl - DOT_PRODUCT(coeff_d,dens_nucl**(coeff_delta))

!   symmetry energy at half saturation density
    J_half= + R_5_18*hbarc_square*c_tau/Mass_n*(dens_nucl/four)**R_2_3/two &
            + R_5_18*hbarc_square*c_tau/Mass_p*(dens_nucl/four)**R_2_3/two &
            + R_10_9*c_tau*dummy_coeff*(dens_nucl/four)**R_5_3 &
            - coeff_b*dens_nucl/two &
            - DOT_PRODUCT(coeff_d,(dens_nucl/two)**(coeff_delta))

! slope of the symmetry energy of Symmetric Nuclear Matter
! L_V = 3n_sat dS(n)/dn|n=n_sat (often reffered to as J)
    L_sym = R_5_9*hbarc_square*c_tau/Mass_n*(dens_nucl/two)**R_2_3/two &
          + R_5_9*hbarc_square*c_tau/Mass_p*(dens_nucl/two)**R_2_3/two &
          + TEN*R_5_9*c_tau*dummy_coeff*(dens_nucl/two)**R_5_3 &
          - three*coeff_b*dens_nucl &
          - three*DOT_PRODUCT(coeff_d*(coeff_delta),dens_nucl**(coeff_delta))

! calculate the 2nd derivative w.r.t. n of the symmetry energy for symmetric Nuclear Matter
! K_S = 9n_sat² d²S(n)/dn²|n=n_sat (often reffered to as K_sym)

    K_sym = - R_5_9*hbarc_square*c_tau/Mass_n*(dens_nucl/two)**R_2_3/two &
            - R_5_9*hbarc_square*c_tau/Mass_p*(dens_nucl/two)**R_2_3/two &
            + (100.d0/9.d0)*c_tau*dummy_coeff*(dens_nucl/two)**R_5_3 &
            - 9.d0*DOT_PRODUCT(coeff_d*(coeff_delta)*(coeff_delta-one), &
            dens_nucl**(coeff_delta))

! calculate the 3rd derivative w.r.t. n of the symmetry energy for symmetric Nuclear Matter
! Q_S = 27n_sat³ d³S(n)/dn³|n=n_sat (often reffered to as Q_sym)

    Q_sym = + (20.d0/9.d0)*hbarc_square*c_tau/Mass_n*(dens_nucl/two)**R_2_3/two &
            + (20.d0/9.d0)*hbarc_square*c_tau/Mass_p*(dens_nucl/two)**R_2_3/two &
          - (100.d0/9.d0)*c_tau*dummy_coeff*(dens_nucl/two)**R_5_3 &
          - 27.d0*DOT_PRODUCT(coeff_d*(coeff_delta)*(coeff_delta-one)*(coeff_delta-two), &
          dens_nucl**(coeff_delta))

  END SUBROUTINE SKYRME_SYMMETRY_ENERGY

END MODULE Symmetry_Energy_Mod
