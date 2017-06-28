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
module phase_space

  use Kind_Types_Mod

  implicit none

  save

! phase space limits
  REAL(DP) :: ye_min, ye_max, ye_step, ye_incr
  REAL(DP) :: n_min, n_max, n_step, n_incr
  REAL(DP) :: t_min, t_max, t_step, t_incr
  REAL(DP) :: ye_test, n_test, t_test
  INTEGER(I4B) ::  iye_ini, iye_fin, in_ini, in_fin, it_ini, it_fin

! transition values:
  REAL(DP) :: n_trans, n_delta, n_tol, nt_max, nt_min

end module phase_space
