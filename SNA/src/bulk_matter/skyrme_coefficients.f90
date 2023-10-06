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

MODULE Skyrme_Coefficients_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL

  IMPLICIT NONE

  SAVE

  INTEGER(I4B) :: Skyrme_parametrization_type
  INTEGER(I4B) :: non_local_terms
  INTEGER(I4B) :: Use_default_constants
  LOGICAL(LGCL) :: High_density_stiffening

! "Standard" Skyrme coefficients
  REAL(DP) :: Coeff_t0, Coeff_t1, Coeff_t2
  REAL(DP) :: Coeff_x0, Coeff_x1, Coeff_x2

  REAL(DP), DIMENSION(:), ALLOCATABLE :: Coeff_t3
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Coeff_x3
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Coeff_sigma

! "Simplified" Skyrme coefficients
  REAL(DP) :: Coeff_a, Coeff_b

  REAL(DP), DIMENSION(:), ALLOCATABLE :: Coeff_c, Coeff_d
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Coeff_delta

  REAL(DP) :: Coeff_alpha1, Coeff_alpha2, Coeff_alpha3
  REAL(DP) :: Coeff_beta1, Coeff_beta2, Coeff_beta3
  REAL(DP) :: Coeff_eps_n, Coeff_eps_p
  REAL(DP) :: Coeff_n_off

  REAL(DP) :: Coeff_qnn, Coeff_qpp
  REAL(DP) :: Coeff_qnp, Coeff_qpn

END MODULE Skyrme_Coefficients_Mod
