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

MODULE Physical_Constants_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B

  IMPLICIT NONE

  SAVE

! some numbers
  REAL(DP), PARAMETER :: ZERO = 0.0_DP, ONE = 1.0_DP, TWO = 2.0_DP
  REAL(DP), PARAMETER :: THREE= 3.0_DP, FOUR= 4.0_DP, FIVE= 5.0_DP
  REAL(DP), PARAMETER :: SIX = 6.0_DP, TEN  =10.0_DP

  REAL(DP), PARAMETER :: PI = ACOS(-1.0_DP), TWO_PI_SQUARE = TWO*PI*PI
  REAL(DP), PARAMETER :: HALF = ONE/TWO, R_1_2 = HALF
  REAL(DP), PARAMETER :: R_1_3 = ONE/THREE, R_2_3 = TWO/THREE
  REAL(DP), PARAMETER :: R_3_2 = THREE/TWO, R_5_2 = FIVE/TWO
  REAL(DP), PARAMETER :: R_3_5 = THREE/FIVE, R_5_3 = FIVE/THREE
  REAL(DP), PARAMETER :: R_5_18 = FIVE/18.0_DP, R_5_9 = FIVE/9.0_DP
  REAL(DP), PARAMETER :: R_2_9 = TWO/9.0_DP, R_9_2 = ONE/R_2_9
  REAL(DP), PARAMETER :: R_10_9 = TWO*R_5_9

! electron, neutron and proton masses in (MeV c²)
  REAL(DP), PARAMETER :: Mass_e = 0.510999_DP
! neutron and proton mass set in the code.
  REAL(DP)            :: Mass_n
  REAL(DP)            :: Mass_p
  REAL(DP)            :: Neut_Prot_Mass_Diff

! volume (fm^3) and binding energy of alpha particles (MeV)
  REAL(DP)            :: v_alpha
  REAL(DP)            :: b_alpha
  REAL(DP)            :: M_alpha

! hbar in Mev fm
  REAL(DP), PARAMETER :: Hbarc = 197.3269788_DP,  Hbarc_Square = Hbarc*Hbarc
  REAL(DP), PARAMETER :: Alpha = 1.0_DP/137.036_DP
  REAL(DP), PARAMETER :: C_tau = R_3_5*(THREE*PI*PI)**R_2_3
  REAL(DP)            :: N_Q

! some constants used throughout the code
  REAL(DP), PARAMETER :: BETA_0 = ((THREE**FIVE)*(PI*ALPHA*Hbarc)/5.0_DP)**R_1_3

!   minimum/maximum densities used throughout to prevent underflow/overflow
  REAL(DP), PARAMETER :: log10_dens_min = -100.0_DP
  REAL(DP), PARAMETER :: log10_dens_max =  100.0_DP
  REAL(DP), PARAMETER :: dens_min = TEN**log10_dens_min
  REAL(DP), PARAMETER :: dens_max = TEN**log10_dens_max

!  for electron EoS subroutines
!  and conversion factor for tables
  !from NIST, also defined in lseos.2.7.f for EL_EOS
   REAL(DP), PARAMETER :: neutron_mass_MeV = 939.5654133d0 !neutron mass
   REAL(DP), PARAMETER :: proton_mass_MeV  = 938.2720813d0 ! proton mass

   REAL(DP), PARAMETER :: HBC = Hbarc
   REAL(DP), PARAMETER :: alpha_coul = 1.d0/137.035999139d0
   REAL(DP), PARAMETER :: clite = 29979245800.0d0

   REAL(DP), PARAMETER :: kb_erg = 1.38064852d-16
   REAL(DP), PARAMETER :: kb_mev = 8.6173303d-11
   REAL(DP), PARAMETER :: erg_to_mev = 6.241509099d5
   REAL(DP), PARAMETER :: mev_to_cgs = 1.6021766208d-6
   REAL(DP), PARAMETER :: temp_mev_to_kelvin = 1.16045221d10
   REAL(DP), PARAMETER :: press_EOS_to_cgs = 1.6021766208d33
   REAL(DP), PARAMETER :: energy_EOS_to_cgs = 1.6021766208d-6/1.674927471d-24
   REAL(DP), PARAMETER :: entropy_EOS_to_cgs = 1.38064852d-16/1.674927471d-24
   REAL(DP), PARAMETER :: rho_cgs_to_EOS = 5.970407778d-16 !neutron mass
   REAL(DP), PARAMETER :: press_cgs_to_EOS  = 6.24150964d-34
   REAL(DP), PARAMETER :: energy_cgs_to_EOS = 1.674927211d-24/1.60217649d-6
   REAL(DP), PARAMETER :: entropy_cgs_to_EOS = 1.674927211d-24/1.38064852d-16
   REAL(DP), PARAMETER :: planck_cgs = 6.626070040d-27
   REAL(DP), PARAMETER :: pi1 = 3.14159265358979d0
   REAL(DP), PARAMETER :: avo = 5.970407778d23 ! 1/m_n_cgs
   REAL(DP), PARAMETER :: sac_const = 2.0d0 * pi1 * (hbc)**2
  !   REAL(DP), PARAMETER :: avo = 6.0221367d23

   REAL(DP), PARAMETER :: m_n_cgs = neutron_mass_MeV*mev_to_cgs/clite**2
   REAL(DP), PARAMETER :: m_p_cgs = proton_mass_MeV*mev_to_cgs/clite**2

END MODULE Physical_Constants_Mod
