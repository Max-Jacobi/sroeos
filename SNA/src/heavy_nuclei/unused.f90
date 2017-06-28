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
MODULE UNUSED_Mod

  USE Kind_Types_Mod, ONLY : DP

  IMPLICIT NONE

  SAVE

  REAL(DP) :: un_A, un_Z, un_rad, un_F, un_P, un_S, un_E
  REAL(DP) :: un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT
  REAL(DP) :: un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT
  REAL(DP) :: un_mu_no, un_mu_po, un_Meff_no, un_Meff_po
  REAL(DP) :: un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy
  REAL(DP) :: un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy
  REAL(DP) :: un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn
  REAL(DP) :: un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT

!$OMP THREADPRIVATE(un_A, un_Z, un_rad, un_F, un_P, un_S, un_E)
!$OMP THREADPRIVATE(un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT)
!$OMP THREADPRIVATE(un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT)
!$OMP THREADPRIVATE(un_mu_no, un_mu_po, un_Meff_no, un_Meff_po)
!$OMP THREADPRIVATE(un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy)
!$OMP THREADPRIVATE(un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy)
!$OMP THREADPRIVATE(un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn)
!$OMP THREADPRIVATE(un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT)

END MODULE UNUSED_Mod
