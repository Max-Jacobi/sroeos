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

MODULE Nuclear_Matter_Properties_Mod

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

  SAVE

  REAL(DP) :: Nuc_Bind_Ener
  REAL(DP) :: Nuc_Sat_Dens

  REAL(DP) :: Pressure_SNM, Pressure_PNM

  REAL(DP) :: Nuc_Ener_SNM_E0
  REAL(DP) :: Nuc_Ener_SNM_K0
  REAL(DP) :: Nuc_Ener_SNM_Q0

  REAL(DP) :: Nuc_Sym_Param_J
  REAL(DP) :: Nuc_Sym_Param_J_half
  REAL(DP) :: Nuc_Sym_Param_L
  REAL(DP) :: Nuc_Sym_Param_K
  REAL(DP) :: Nuc_Sym_Param_Q

  REAL(DP) :: Mass_n_SNM_0, Mass_p_SNM_0

  REAL(DP) :: Isospin_Incompressibility

  REAL(DP) :: Nucl_Sym_Half_Sat

  REAL(DP) :: Surface_Tension_Sigma_S
  REAL(DP) :: Surface_Sym_Ener
  REAL(DP) :: Surface_Level_Dens

END MODULE Nuclear_Matter_Properties_Mod
