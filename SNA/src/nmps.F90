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

subroutine calc_nmps( &
    py_Coeff_a, &
    py_Coeff_b, &
    py_Coeff_c, &
    py_Coeff_d, &
    py_Coeff_delta, &
    py_Coeff_alpha1, &
    py_Coeff_alpha2, &
    py_non_local_terms, &
    py_Nuc_Bind_Ener, &
    py_Nuc_Sat_Dens, &
    py_Pressure_SNM, &
    py_Pressure_PNM, &
    py_Nuc_Ener_SNM_K0, &
    py_Nuc_Ener_SNM_Q0, &
    py_Nuc_Sym_Param_J, &
    py_Nuc_Sym_Param_L, &
    py_Nuc_Sym_Param_K, &
    py_Nuc_Sym_Param_Q, &
    py_Mass_n_SNM_0, &
    py_Mass_p_SNM_0)

  USE Determine_Nuclear_Properties_Mod
  USE Nuclear_Matter_Properties_Mod
  USE Skyrme_Coefficients_Mod
  USE Physical_Constants_Mod
  IMPLICIT NONE




  INTEGER, intent(in) :: py_non_local_terms
  REAL(8), intent(in) :: py_Coeff_a, py_Coeff_b
  REAL(8), intent(in), DIMENSION(*) :: py_Coeff_c, py_Coeff_d
  REAL(8), intent(in), DIMENSION(*) :: py_Coeff_delta
  REAL(8), intent(in) :: py_Coeff_alpha1, py_Coeff_alpha2

  REAL(8), intent(out) :: py_Nuc_Bind_Ener
  REAL(8), intent(out) :: py_Nuc_Sat_Dens
  REAL(8), intent(out) :: py_Pressure_SNM
  REAL(8), intent(out) :: py_Pressure_PNM
  REAL(8), intent(out) :: py_Nuc_Ener_SNM_K0
  REAL(8), intent(out) :: py_Nuc_Ener_SNM_Q0
  REAL(8), intent(out) :: py_Nuc_Sym_Param_J
  REAL(8), intent(out) :: py_Nuc_Sym_Param_L
  REAL(8), intent(out) :: py_Nuc_Sym_Param_K
  REAL(8), intent(out) :: py_Nuc_Sym_Param_Q
  REAL(8), intent(out) :: py_Mass_n_SNM_0
  REAL(8), intent(out) :: py_Mass_p_SNM_0

  ! set some parameters by hand
  Skyrme_parametrization_type = 1
  Use_default_constants = 0
  High_density_stiffening = .false.
  non_local_terms = py_non_local_terms

  Mass_n  = 939.56540D0
  Mass_p  = 938.27204D0
  Neut_Prot_Mass_Diff = Mass_n - Mass_p
  b_alpha = 28.29552D0 + TWO*Neut_Prot_Mass_Diff
  v_alpha = 24.0D0

  if (allocated(Coeff_c)) then
    DEALLOCATE(Coeff_c, Coeff_d, Coeff_delta)
  endif
  ALLOCATE(Coeff_c(non_local_terms))
  ALLOCATE(Coeff_d(non_local_terms))
  ALLOCATE(Coeff_delta(non_local_terms))
  coeff_qnn = 0
  coeff_qpp = 0
  coeff_qnp = 0
  coeff_qpn = 0

  Coeff_a = py_Coeff_a
  Coeff_b = py_Coeff_b
  Coeff_c = py_Coeff_c(1:non_local_terms)
  Coeff_d = py_Coeff_d(1:non_local_terms)
  Coeff_delta = py_Coeff_delta(1:non_local_terms)
  Coeff_alpha1 = py_Coeff_alpha1
  Coeff_alpha2 = py_Coeff_alpha2

  CALL DETERMINE_NUCLEAR_PROPERTIES

  py_Nuc_Bind_Ener                  = Nuc_Bind_Ener
  py_Nuc_Sat_Dens                   = Nuc_Sat_Dens
  py_Pressure_SNM                   = Pressure_SNM
  py_Pressure_PNM                   = Pressure_PNM
  py_Nuc_Ener_SNM_K0                = Nuc_Ener_SNM_K0
  py_Nuc_Ener_SNM_Q0                = Nuc_Ener_SNM_Q0
  py_Nuc_Sym_Param_J                = Nuc_Sym_Param_J
  py_Nuc_Sym_Param_L                = Nuc_Sym_Param_L
  py_Nuc_Sym_Param_K                = Nuc_Sym_Param_K
  py_Nuc_Sym_Param_Q                = Nuc_Sym_Param_Q
  py_Mass_n_SNM_0                   = Mass_n_SNM_0
  py_Mass_p_SNM_0                   = Mass_p_SNM_0

END subroutine calc_nmps
