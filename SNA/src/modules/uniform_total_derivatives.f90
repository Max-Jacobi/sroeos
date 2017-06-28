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
MODULE Find_Uniform_Matter_Derivatives_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, HALF, ONE, TWO, THREE, FOUR, TEN, &
                             v_alpha

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Find_Uniform_Matter_Derivatives ( n, T, Y, n_n, n_p, n_alpha, &
     Dn_n_Deta_n, Dn_p_Deta_n, Dn_n_Deta_p, Dn_p_Deta_p, Dn_n_DT, Dn_p_DT, &
     Dn_alpha_Deta_n, Dn_alpha_Deta_p, Dn_alpha_DT, &
     Deta_n_Dn, Deta_p_Dn, Deta_n_DT, Deta_p_DT, Deta_n_Dy, Deta_p_Dy )

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: n, y, T
    REAL(DP), INTENT(IN) :: n_n, n_p, n_alpha
    REAL(DP), INTENT(IN) :: Dn_n_Deta_n, Dn_n_Deta_p, Dn_p_Deta_n, Dn_p_Deta_p
    REAL(DP), INTENT(IN) :: Dn_n_DT, Dn_p_DT
    REAL(DP), INTENT(IN) :: Dn_alpha_Deta_n, Dn_alpha_Deta_p, Dn_alpha_DT

    REAL(DP), INTENT(OUT) :: Deta_n_Dn, Deta_p_Dn
    REAL(DP), INTENT(OUT) :: Deta_n_DT, Deta_p_DT
    REAL(DP), INTENT(OUT) :: Deta_n_Dy, Deta_p_Dy

    REAL(DP) :: exc_v_alpha
    REAL(DP) :: dAi_dxj(2,2), dAi_dyj(2,3), DET
!------------------------------------------------------------------------------!
! USE CONSTRAINTS TO OBTAIN TOTAL DERIVATIVES
!  OF INTERNAL VARIABLES (eta_no,eta_po) w.r.t.
!   FIXED PARAMETERS (n, T, Ye)
!------------------------------------------------------------------------------!
! DIJ(I,J) = d A_I / d x_j where
!
!    A_1 = n - u*ni - (1-u)*(4*n_alpha+(n_no+n_po)*(1-n_alpha*v_alpha))
!    A_2 = n*Ye - u*n*xi - (1-u)*(2*n_alpha+n_po*(1-n_alpha*v_alpha))
!
! if u-> 0:
!    A_1 = n - 4*n_alpha - (n_no+n_po)*(1-n_alpha*v_alpha)
!    A_2 = n*Ye - 2*n_alpha - n_po*(1-n_alpha*v_alpha)
!
!    x = (eta_po, eta_no)
!
! FIJ(I,J) = - d A_I / d y_j where
!
!    y = (n, T, Ye)
!
!    Equation 1 (Baryon conservation)
!
    exc_v_alpha = ONE - n_alpha*v_alpha

    dAi_dxj(1,1) = - exc_v_alpha*(Dn_p_Deta_p + Dn_n_Deta_p) &
                   - (FOUR-v_alpha*(n_n+n_p))*Dn_alpha_Deta_p
    dAi_dxj(1,2) = - exc_v_alpha*(Dn_p_Deta_n + Dn_n_Deta_n) &
                   - (FOUR-v_alpha*(n_n+n_p))*Dn_alpha_Deta_n
!
    dAi_dyj(1,1) = ONE
    dAi_dyj(1,2) = - ((FOUR-v_alpha*(n_n+n_p))*Dn_alpha_DT &
                   + exc_v_alpha*(Dn_n_DT+Dn_p_DT))
    dAi_dyj(1,3) = ZERO
!
!    Equation 2 (Charge conservation)
!
    dAi_dxj(2,1) = - exc_v_alpha*Dn_p_Deta_p - (TWO-v_alpha*n_p)*Dn_alpha_Deta_p
    dAi_dxj(2,2) = - exc_v_alpha*Dn_p_Deta_n - (TWO-v_alpha*n_p)*Dn_alpha_Deta_n
!
    dAi_dyj(2,1) = Y
    dAi_dyj(2,2) = - (TWO-v_alpha*n_p)*Dn_alpha_DT - exc_v_alpha*Dn_p_DT
    dAi_dyj(2,3) = n

!   DETERMINANT
    DET = dAi_dxj(1,1)*dAi_dxj(2,2) - dAi_dxj(1,2)*dAi_dxj(2,1)

!   Total derivatives w.r.t. n, T, Y
    Deta_p_Dn = (dAi_dxj(1,2)*dAi_dyj(2,1) - dAi_dxj(2,2)*dAi_dyj(1,1))/DET
    Deta_n_Dn = (dAi_dxj(2,1)*dAi_dyj(1,1) - dAi_dxj(1,1)*dAi_dyj(2,1))/DET

    Deta_p_DT = (dAi_dxj(1,2)*dAi_dyj(2,2) - dAi_dxj(2,2)*dAi_dyj(1,2))/DET
    Deta_n_DT = (dAi_dxj(2,1)*dAi_dyj(1,2) - dAi_dxj(1,1)*dAi_dyj(2,2))/DET

    Deta_p_Dy = (dAi_dxj(1,2)*dAi_dyj(2,3) - dAi_dxj(2,2)*dAi_dyj(1,3))/DET
    Deta_n_Dy = (dAi_dxj(2,1)*dAi_dyj(1,3) - dAi_dxj(1,1)*dAi_dyj(2,3))/DET

  END SUBROUTINE Find_Uniform_Matter_Derivatives

END MODULE Find_Uniform_Matter_Derivatives_Mod
