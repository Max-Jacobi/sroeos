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
MODULE Find_Non_Uniform_Matter_Derivatives_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, HALF, ONE, TWO, FOUR, v_alpha
  USE LU_decomposition

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Find_Non_Uniform_Matter_Derivatives ( n, Y, T, &
    u, n_ni, n_pi, n_i, y_i, n_no, n_po, n_o, n_alpha, &
    Dn_no_Deta_no,  Dn_po_Deta_no,  Dn_no_Deta_po,  Dn_po_Deta_po, &
    Dn_no_DT, Dn_po_DT, &
    Dmu_no_Deta_no, Dmu_po_Deta_no, Dmu_no_Deta_po, Dmu_po_Deta_po, &
    Dmu_no_DT, Dmu_po_DT, &
    DP_o_Deta_no, DP_o_Deta_po, DP_o_DT, &
    Dn_alpha_Deta_no, Dn_alpha_Deta_po, Dn_alpha_DT, &
    DP_alpha_Deta_no, DP_alpha_Deta_po, DP_alpha_DT, &
    Dmu_ni_Dn_ni, Dmu_pi_Dn_ni, Dmu_ni_Dn_pi, Dmu_pi_Dn_pi, &
    Dmu_ni_DT, Dmu_pi_DT, DP_i_Dn_ni, DP_i_Dn_pi, DP_i_DT, &
    DB1_Du, DB1_Dn_i, DB1_Dy_i, DB1_DT, &
    DB2_Du, DB2_Dn_i, DB2_Dy_i, DB2_DT, &
    DB3_Du, DB3_Dn_i, DB3_Dy_i, DB3_DT, &
    DU_DN, Dy_i_DN, Dn_i_DN, Deta_no_DN, Deta_po_DN, &
    DU_DT, Dy_i_DT, Dn_i_DT, Deta_no_DT, Deta_po_DT, &
    DU_DY, Dy_i_DY, Dn_i_DY, Deta_no_DY, Deta_po_DY )

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: n, y, T
    REAL(DP), INTENT(IN) :: u, n_ni, n_pi, n_i, y_i
    REAL(DP), INTENT(IN) :: n_no, n_po, n_o, n_alpha
    REAL(DP), INTENT(IN) :: Dn_no_Deta_no,  Dn_po_Deta_no
    REAL(DP), INTENT(IN) :: Dn_no_Deta_po,  Dn_po_Deta_po
    REAL(DP), INTENT(IN) :: Dn_no_DT, Dn_po_DT
    REAL(DP), INTENT(IN) :: Dmu_no_Deta_no, Dmu_po_Deta_no
    REAL(DP), INTENT(IN) :: Dmu_no_Deta_po, Dmu_po_Deta_po
    REAL(DP), INTENT(IN) :: Dmu_no_DT, Dmu_po_DT
    REAL(DP), INTENT(IN) :: DP_o_Deta_no, DP_o_Deta_po, DP_o_DT
    REAL(DP), INTENT(IN) :: Dmu_ni_Dn_ni, Dmu_pi_Dn_ni
    REAL(DP), INTENT(IN) :: Dmu_ni_Dn_pi, Dmu_pi_Dn_pi
    REAL(DP), INTENT(IN) :: Dmu_ni_DT, Dmu_pi_DT
    REAL(DP), INTENT(IN) :: DP_i_Dn_ni, DP_i_Dn_pi, DP_i_DT
    REAL(DP), INTENT(IN) :: Dn_alpha_Deta_no, Dn_alpha_Deta_po, Dn_alpha_DT
    REAL(DP), INTENT(IN) :: DP_alpha_Deta_no, DP_alpha_Deta_po, DP_alpha_DT
    REAL(DP), INTENT(IN) :: DB1_Du, DB1_Dn_i, DB1_Dy_i, DB1_DT
    REAL(DP), INTENT(IN) :: DB2_Du, DB2_Dn_i, DB2_Dy_i, DB2_DT
    REAL(DP), INTENT(IN) :: DB3_Du, DB3_Dn_i, DB3_Dy_i, DB3_DT

    REAL(DP), INTENT(OUT) :: DU_DN, Dy_i_DN, Dn_i_DN, Deta_no_DN, Deta_po_DN
    REAL(DP), INTENT(OUT) :: DU_DT, Dy_i_DT, Dn_i_DT, Deta_no_DT, Deta_po_DT
    REAL(DP), INTENT(OUT) :: DU_DY, Dy_i_DY, Dn_i_DY, Deta_no_DY, Deta_po_DY

    REAL(DP) :: exc_v_alpha
    REAL(DP) :: dAi_dxj(5,5), dAi_dyj(5,3), dAi_dxj_LU(5,5), dxi_dyj(5,3)

    INTEGER(I4B) :: ipvt(5)
!------------------------------------------------------------------!
! USE DERIVATIVES AND ITS CONSTRAINTS TO OBTAIN TOTAL DERIVATIVES
!  OF INTERNAL VARIABLES (u,x_i,n_i,eta_no,eta_po) w.r.t.
!   FIXED PARAMETERS (n, T, Y)
!------------------------------------------------------------------!
! dAi_dxj(I,J) = d A_I / d x_j where
!
!    A_1 = n - u*ni - (1-u)*(4*n_alpha+(n_no+n_po)*(1-n_alpha*v_alpha))
!    A_2 = n*Y - u*n*xi - (1-u)*(2*n_alpha+n_po*(1-n_alpha*v_alpha)
!    A_3 = mu_ni - B2 - mu_no
!    A_4 = mu_pi - B3 - mu_po
!    A_5 = P_i - B1 - P_o - P_alpha
!
!    x = (u, xi, ni, eta_po, eta_no)
!
! dAi_dyj(I,J) = - d A_I / d y_j where
!
!    y = (n, T, Y)
!
    exc_v_alpha = ONE - n_alpha*v_alpha
!
!    Equation 1 (Baryon conservation)
!
    dAi_dxj(1,1) = - n_i + (FOUR*n_alpha + n_o*exc_v_alpha)
    dAi_dxj(1,2) =   ZERO
    dAi_dxj(1,3) = - u
    dAi_dxj(1,4) = - (ONE-u)*exc_v_alpha*(Dn_po_Deta_po + Dn_no_Deta_po) &
               - (ONE-u)*(FOUR - v_alpha*n_o)*Dn_alpha_Deta_po
    dAi_dxj(1,5) = - (ONE-u)*exc_v_alpha*(Dn_no_Deta_no + Dn_po_Deta_no) &
               - (ONE-u)*(FOUR - v_alpha*n_o)*Dn_alpha_Deta_no
!
    dAi_dyj(1,1) = - ONE
    dAi_dyj(1,2) = (ONE-u)*( (FOUR - v_alpha*n_o )*Dn_alpha_DT &
                                 + exc_v_alpha*(Dn_no_DT + Dn_po_DT) )
    dAi_dyj(1,3) =   ZERO
!
!    Equation 2 (Charge conservation)
!
    dAi_dxj(2,1) = - n_i*y_i + (TWO*n_alpha+n_po*exc_v_alpha)
    dAi_dxj(2,2) = - u*n_i
    dAi_dxj(2,3) = - u*y_i
    dAi_dxj(2,4) = - (ONE-u)*exc_v_alpha*Dn_po_Deta_po &
               - (ONE-u)*(TWO-v_alpha*n_po)*Dn_alpha_Deta_po
    dAi_dxj(2,5) = - (ONE-u)*exc_v_alpha*Dn_po_Deta_no &
               - (ONE-u)*(TWO-v_alpha*n_po)*Dn_alpha_Deta_no

    dAi_dyj(2,1) = - Y
    dAi_dyj(2,2) = (ONE-u)*(TWO-v_alpha*n_po)*Dn_alpha_DT &
             + (ONE-u)*exc_v_alpha*Dn_po_DT
    dAi_dyj(2,3) = - n
!
!    Equation 3 (Proton chemical equilibrium)
!
    dAi_dxj(3,1) = - DB3_Du
    dAi_dxj(3,2) = - DB3_Dy_i + n_i*(Dmu_pi_Dn_pi - Dmu_pi_Dn_ni)
    dAi_dxj(3,3) = - DB3_Dn_i + (ONE-y_i)*Dmu_pi_Dn_ni + y_i*Dmu_pi_Dn_pi
    dAi_dxj(3,4) = - Dmu_po_Deta_po
    dAi_dxj(3,5) = - Dmu_po_Deta_no

    dAi_dyj(3,1) =   ZERO
    dAi_dyj(3,2) = - Dmu_pi_DT + DB3_DT + Dmu_po_DT
    dAi_dyj(3,3) =   ZERO
!
!    Equation 4 (Neutron chemical equilibrium)
!
    dAi_dxj(4,1) = - DB2_Du
    dAi_dxj(4,2) = - DB2_Dy_i + n_i*(Dmu_ni_Dn_pi - Dmu_ni_Dn_ni)
    dAi_dxj(4,3) = - DB2_Dn_i + (ONE-y_i)*Dmu_ni_Dn_ni + y_i*Dmu_ni_Dn_pi
    dAi_dxj(4,4) = - Dmu_no_Deta_po
    dAi_dxj(4,5) = - Dmu_no_Deta_no

    dAi_dyj(4,1) =   ZERO
    dAi_dyj(4,2) = - Dmu_ni_DT + DB2_DT + Dmu_no_DT
    dAi_dyj(4,3) =   ZERO
!
!    Equation 5 (Pressure equilibrium)
!
    dAi_dxj(5,1) = - DB1_Du
    dAi_dxj(5,2) = - DB1_Dy_i + n_i*(DP_i_Dn_pi-DP_i_Dn_ni)
    dAi_dxj(5,3) = - DB1_Dn_i + (ONE-y_i)*DP_i_Dn_ni + y_i*DP_i_Dn_pi
    dAi_dxj(5,4) = - DP_o_Deta_po - DP_alpha_Deta_po
    dAi_dxj(5,5) = - DP_o_Deta_no - DP_alpha_Deta_no

    dAi_dyj(5,1) =   ZERO
    dAi_dyj(5,2) = - DP_i_DT + DB1_DT + DP_o_DT + DP_alpha_DT
    dAi_dyj(5,3) =   ZERO
!
!    LU decomposition of dAi_dyj matrix
!
    CALL MATLUD(dAi_dxj,dAi_dxj_LU,5,ipvt)
    dxi_dyj = ZERO
!
!   Solve the LU decomposed linear system
!    to get derivatives w.r.t. density n
!
    CALL MLUSLV(dAi_dxj_LU,dxi_dyj(1:5,1),dAi_dyj(1:5,1),ipvt,5)
    WHERE (isnan(dxi_dyj)) dxi_dyj = 1.d-290

    DU_DN      = dxi_dyj(1,1)
    Dy_i_DN    = dxi_dyj(2,1)
    Dn_i_DN    = dxi_dyj(3,1)
    Deta_po_DN = dxi_dyj(4,1)
    Deta_no_DN = dxi_dyj(5,1)
!
!   Solve the LU decomposed linear system
!    to get derivatives w.r.t. temperature T
!
    CALL MLUSLV(dAi_dxj_LU,dxi_dyj(1:5,2),dAi_dyj(1:5,2),ipvt,5)
    WHERE (isnan(dxi_dyj)) dxi_dyj = 1.d-290

    DU_DT      = dxi_dyj(1,2)
    Dy_i_DT    = dxi_dyj(2,2)
    Dn_i_DT    = dxi_dyj(3,2)
    Deta_po_DT = dxi_dyj(4,2)
    Deta_no_DT = dxi_dyj(5,2)
!
!   solve the LU decomposed linear system
!    to get derivatives w.r.t. temperature Y
!
    CALL MLUSLV(dAi_dxj_LU,dxi_dyj(1:5,3),dAi_dyj(1:5,3),ipvt,5)

    WHERE (isnan(dxi_dyj)) dxi_dyj = 1.d-290

    DU_DY      = dxi_dyj(1,3)
    Dy_i_DY    = dxi_dyj(2,3)
    Dn_i_DY    = dxi_dyj(3,3)
    Deta_po_DY = dxi_dyj(4,3)
    Deta_no_DY = dxi_dyj(5,3)

    ! write (*,*)
    ! write (*,"(5ES20.12)") dAi_dyj
    ! write (*,*)
    ! write (*,"(5ES20.12)") dAi_dxj_LU
    ! write (*,*)
    ! write (*,"(5ES20.12)") dxi_dyj
    ! write (*,*)

  END SUBROUTINE Find_Non_Uniform_Matter_Derivatives

END MODULE Find_Non_Uniform_Matter_Derivatives_Mod
