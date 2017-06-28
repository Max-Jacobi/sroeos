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
MODULE Determine_Nuclear_Properties_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  !USE Define_Operators_Mod
  USE Physical_Constants_Mod
  USE Skyrme_Coefficients_Mod
  USE Skyrme_Bulk_Mod
  USE Skyrme_Bulk_Observables_Mod
  USE Nuclear_Matter_Properties_Mod
  USE Symmetry_Energy_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE DETERMINE_NUCLEAR_PROPERTIES

    IMPLICIT NONE

    REAL(DP) :: log10_dens, log10_half_dens, low_T

    low_T = 0.01_DP

    CALL SATURATION_DENSITY(Nuc_Sat_Dens,Nuc_Bind_Ener)

!   Correct Nuc_Bind_Ener to include neutron proton mass difference
    Nuc_Bind_Ener = Nuc_Bind_Ener - Neut_Prot_Mass_Diff/2.D0

    log10_half_dens = log10(nuc_sat_dens/two)
    log10_dens = log10(nuc_sat_dens)

    PRESSURE_SNM = SKYRME_BULK_PRESSURE(log10_half_dens,log10_half_dens,low_T)
    PRESSURE_PNM = SKYRME_BULK_PRESSURE(log10_dens,-100.0_DP,low_T)

    Nuc_Ener_SNM_K0 = SKYRME_BULK_INCOMPRESSIBILITY(&
    log10_half_dens,log10_half_dens,low_T)

    Nuc_Ener_SNM_Q0 = SKYRME_BULK_SKEWNESS(&
    log10_half_dens,log10_half_dens,low_T)

    CALL SKYRME_SYMMETRY_ENERGY(Nuc_Sym_Param_J,Nuc_Sym_Param_J_half,&
    Nuc_Sym_Param_L,Nuc_Sym_Param_K,Nuc_Sym_Param_Q)

    CALL PRINT_NUCLEAR_PROPERTIES()

    ! WRITE(*,"(1ES20.12)") Nuc_Sat_Dens, Nuc_Bind_Ener,&
    ! PRESSURE_SNM, THREE*PRESSURE_PNM/Nuc_Sat_Dens/Nuc_Sym_Param_L, &
    ! Nuc_Ener_SNM_K0, Nuc_Ener_SNM_Q0, Nuc_Sym_Param_J,&
    ! Nuc_Sym_Param_J_half/Nuc_Sym_Param_J,&
    ! Nuc_Sym_Param_L,Nuc_Sym_Param_K,Nuc_Sym_Param_Q,&
    ! Mass_n_SNM_0/Mass_n,Mass_n_SNM_0-Mass_p_SNM_0

  END SUBROUTINE DETERMINE_NUCLEAR_PROPERTIES

! Simple subroutine to estimate nuclear saturation density
!  and binding energy of symmetric nuclear matter
  SUBROUTINE SATURATION_DENSITY(sat_dens,bind_e)

    IMPLICIT NONE

    REAL(DP), INTENT(OUT) :: sat_dens, bind_e
    REAL(DP) :: dens_nucl, log10_dens_neut, log10_dens_prot, low_T
    REAL(DP) :: Min_energy, Energy_density, Energy_per_nucleon
    INTEGER(I4B) :: i

    Min_energy = 1.d100
    low_T = 1.d-10

    DO i = -20000, 20000
      dens_nucl = 0.16d0*(one + dble(i)/2.d5)
      log10_dens_neut = log10(half*dens_nucl)
      log10_dens_prot = log10(half*dens_nucl)
      !CALL ENERGY_DENSITY()
      Energy_density = SKYRME_BULK_ENERGY(log10_dens_neut,log10_dens_neut,low_T)
      Energy_per_nucleon = Energy_density/dens_nucl
      Min_energy = min(Min_energy,Energy_per_nucleon)
      if (Min_energy == Energy_per_nucleon) sat_dens = dens_nucl
    ENDDO

    bind_e = Min_energy

    Mass_n_SNM_0 = Hbarc_square/TWO/Mass_n + &
                      (Coeff_alpha1 + Coeff_alpha2)*sat_dens/TWO
    Mass_n_SNM_0 = Hbarc_square/TWO/Mass_n_SNM_0

    Mass_p_SNM_0 = Hbarc_square/TWO/Mass_p + &
                      (Coeff_alpha1 + Coeff_alpha2)*sat_dens/TWO
    Mass_p_SNM_0 = Hbarc_square/TWO/Mass_p_SNM_0

  END SUBROUTINE SATURATION_DENSITY


  SUBROUTINE PRINT_NUCLEAR_PROPERTIES

    USE Make_Tables_Mod, ONLY : output_directory

    IMPLICIT NONE

    INTEGER(I4B) :: i
    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT   = "(1A10,1ES20.12,1A20)"
    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT_1 = "(1A2,1I1,1A7,1ES20.12,1A18)"
    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT_2 = "(1A6,1I1,1A3,1ES20.12,1A20)"
    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT_3 = "(1A3,1I1,1A6,1ES20.12,1A20)"
    CHARACTER(LEN=128) :: File_Name

    File_Name = ADJUSTL(TRIM(output_directory)) // &
                ADJUSTL(TRIM('/Nuclear_Properties.dat'))

    OPEN(10,file=File_Name)

!   Write nuclear matter properties
    WRITE(10,*)
    WRITE(10,*) 'Specific energy per baryon expansion terms. Eq.47 of SRO.'
    WRITE(10,*)
    WRITE(10,OUTPUT_FORMAT) 'n_0     = ', Nuc_Sat_Dens,    ' fm^-3'
    WRITE(10,OUTPUT_FORMAT) 'E_0     = ', Nuc_Bind_Ener,   ' MeV nucleon^-1'
    WRITE(10,OUTPUT_FORMAT) 'K_0     = ', Nuc_Ener_SNM_K0, ' MeV nucleon^-1'
    WRITE(10,OUTPUT_FORMAT) "K'      = ",-Nuc_Ener_SNM_Q0, ' MeV nucleon^-1'
    WRITE(10,*)
    WRITE(10,*) 'Symmetry energy expansion. Eq.50 of SRO.'
    WRITE(10,*)
    WRITE(10,OUTPUT_FORMAT) 'J       = ', Nuc_Sym_Param_J, ' MeV nucleon^-1'
    WRITE(10,OUTPUT_FORMAT) 'L       = ', Nuc_Sym_Param_L, ' MeV nucleon^-1'
    WRITE(10,OUTPUT_FORMAT) 'K_sym   = ', Nuc_Sym_Param_K, ' MeV nucleon^-1'
    WRITE(10,OUTPUT_FORMAT) "Q_sym   = ", Nuc_Sym_Param_Q, ' MeV nucleon^-1'
    WRITE(10,*)
    WRITE(10,*) ' Nucleon masses and effective masses'
    WRITE(10,*) '  at nuclear saturation density'
    WRITE(10,*)
    WRITE(10,OUTPUT_FORMAT) "m_n      = ", Mass_n
    WRITE(10,OUTPUT_FORMAT) "m_p      = ", Mass_p
    WRITE(10,OUTPUT_FORMAT) "m_n-m_p  = ", Mass_n - Mass_p
    WRITE(10,OUTPUT_FORMAT) "m*_n/m_n = ", Mass_n_SNM_0/Mass_n
    WRITE(10,OUTPUT_FORMAT) "m*_p/m_p = ", Mass_p_SNM_0/Mass_p
    WRITE(10,OUTPUT_FORMAT) "m*_n-m*_p= ", Mass_n_SNM_0-Mass_p_SNM_0
    WRITE(10,*)
    WRITE(10,*) 'Alpha particle properties'
    WRITE(10,*)
    WRITE(10,OUTPUT_FORMAT) "m_alpha  = ", m_alpha
    WRITE(10,OUTPUT_FORMAT) "b_alpha  = ", b_alpha
    WRITE(10,OUTPUT_FORMAT) "v_alpha  = ", v_alpha
    WRITE(10,*)

    CLOSE(10)

  END SUBROUTINE PRINT_NUCLEAR_PROPERTIES

END MODULE Determine_Nuclear_Properties_Mod
