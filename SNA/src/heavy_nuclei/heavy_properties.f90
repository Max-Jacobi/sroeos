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
MODULE Heavy_Properties_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ONE, TWO, FOUR, R_3_2, &
      v_alpha, b_alpha, HC2 => Hbarc_Square, Neut_Prot_Mass_Diff
  USE Skyrme_Coefficients_Mod, ONLY : A => Coeff_a, B => Coeff_b, C => Coeff_c,&
                  D => Coeff_d, alpha1 => Coeff_alpha1, alpha2 => Coeff_alpha2,&
                  delta => Coeff_delta
  IMPLICIT NONE

CONTAINS

  SUBROUTINE HEAVY_PROPERTIES(u, T, n_ni, n_pi,           &
    mu_ni, mu_pi, tau_ni, tau_pi, Meff_ni, Meff_pi, P_i,  &
    Dmu_ni_Dn_i, Dmu_ni_Dy_i, Dmu_ni_DT, &
    Dmu_pi_Dn_i, Dmu_pi_Dy_i, Dmu_pi_DT, &
    Dtau_ni_Dn_i, Dtau_ni_Dy_i, Dtau_ni_DT, &
    Dtau_pi_Dn_i, Dtau_pi_Dy_i, Dtau_pi_DT, &
    DP_i_Dn_i, DP_i_Dy_i, DP_i_DT, Du_DT, Du_Dn, Du_Dy, &
    Dn_i_DT, Dn_i_Dn, Dn_i_Dy,     &
    Dy_i_DT, Dy_i_Dn, Dy_i_Dy,     &
    F_SC, DF_SC_Du, DF_SC_Dn_i, DF_SC_DT, DF_SC_Dy_i,                        &
    D2F_SC_Du2, D2F_SC_DuDn_i, D2F_SC_DuDT, D2F_SC_DuDy_i, D2F_SC_Dn_i2,     &
    D2F_SC_Dn_iDT, D2F_SC_Dn_iDy_i, D2F_SC_DT2, D2F_SC_DTDy_i, D2F_SC_Dy_i2, &
    F_TR, DF_TR_Du, DF_TR_Dn_i, DF_TR_DT, DF_TR_Dy_i,                        &
    D2F_TR_Du2, D2F_TR_DuDn_i, D2F_TR_DuDT, D2F_TR_DuDy_i, D2F_TR_Dn_i2,     &
    D2F_TR_Dn_iDT, D2F_TR_Dn_iDy_i, D2F_TR_DT2, D2F_TR_DTDy_i, D2F_TR_Dy_i2, &
    F_h, DF_h_DT, DF_h_Dn, DF_h_Dy, &
    S_h, DS_h_DT, DS_h_Dn, DS_h_Dy )

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: u, T
    REAL(DP), INTENT(IN) :: n_ni, n_pi, mu_ni, mu_pi
    REAL(DP), INTENT(IN) :: tau_ni, tau_pi, Meff_ni, Meff_pi
    REAL(DP), INTENT(IN) :: Dmu_ni_Dn_i, Dmu_ni_Dy_i, Dmu_ni_DT
    REAL(DP), INTENT(IN) :: Dmu_pi_Dn_i, Dmu_pi_Dy_i, Dmu_pi_DT
    REAL(DP), INTENT(IN) :: Dtau_ni_Dn_i, Dtau_ni_Dy_i, Dtau_ni_DT
    REAL(DP), INTENT(IN) :: Dtau_pi_Dn_i, Dtau_pi_Dy_i, Dtau_pi_DT
    REAL(DP), INTENT(IN) :: P_i, DP_i_Dn_i, DP_i_Dy_i, DP_i_DT
    REAL(DP), INTENT(IN) :: Du_DT, Du_Dn, Du_Dy
    REAL(DP), INTENT(IN) :: Dn_i_DT, Dn_i_Dn, Dn_i_Dy
    REAL(DP), INTENT(IN) :: Dy_i_DT, Dy_i_Dn, Dy_i_Dy


    REAL(DP), INTENT(IN) :: F_SC, DF_SC_Du, DF_SC_Dn_i, DF_SC_DT, DF_SC_Dy_i
    REAL(DP), INTENT(IN) :: D2F_SC_Du2, D2F_SC_DuDn_i, D2F_SC_DuDT
    REAL(DP), INTENT(IN) :: D2F_SC_DuDy_i, D2F_SC_Dn_i2
    REAL(DP), INTENT(IN) :: D2F_SC_Dn_iDT, D2F_SC_Dn_iDy_i, D2F_SC_DT2
    REAL(DP), INTENT(IN) :: D2F_SC_DTDy_i, D2F_SC_Dy_i2
    REAL(DP), INTENT(IN) :: F_TR, DF_TR_Du, DF_TR_Dn_i, DF_TR_DT, DF_TR_Dy_i
    REAL(DP), INTENT(IN) :: D2F_TR_Du2, D2F_TR_DuDn_i, D2F_TR_DuDT
    REAL(DP), INTENT(IN) :: D2F_TR_DuDy_i, D2F_TR_Dn_i2
    REAL(DP), INTENT(IN) :: D2F_TR_Dn_iDT, D2F_TR_Dn_iDy_i, D2F_TR_DT2
    REAL(DP), INTENT(IN) :: D2F_TR_DTDy_i, D2F_TR_Dy_i2

    REAL(DP), INTENT(OUT) :: F_h, DF_h_DT, DF_h_Dn, DF_h_Dy
    REAL(DP), INTENT(OUT) :: S_h, DS_h_DT, DS_h_Dn, DS_h_Dy

    REAL(DP) :: n_i, y_i
    REAL(DP) :: Fi, Ei, Si
    REAL(DP) :: DFi_Dy_i, DEi_Dy_i, DSi_Dy_i
    REAL(DP) :: DFi_Dn_i, DEi_Dn_i, DSi_Dn_i
    REAL(DP) :: DFi_DT,   DEi_DT,   DSi_DT

    REAL(DP) :: Ftr, Etr, Str
    REAL(DP) :: DFtr_Du,   DEtr_Du,   DStr_Du
    REAL(DP) :: DFtr_Dy_i, DEtr_Dy_i, DStr_Dy_i
    REAL(DP) :: DFtr_Dn_i, DEtr_Dn_i, DStr_Dn_i
    REAL(DP) :: DFtr_DT,   DEtr_DT,   DStr_DT

    REAL(DP) :: Fsc, Esc, Ssc
    REAL(DP) :: DFsc_Du,   DEsc_Du,   DSsc_Du
    REAL(DP) :: DFsc_Dy_i, DEsc_Dy_i, DSsc_Dy_i
    REAL(DP) :: DFsc_Dn_i, DEsc_Dn_i, DSsc_Dn_i
    REAL(DP) :: DFsc_DT,   DEsc_DT,   DSsc_DT


    n_i = n_ni + n_pi
    y_i = n_pi / n_i

!   Inside nucleons total free energy and its partial derivatives
    Fi          = n_ni*mu_ni+n_pi*mu_pi-P_i
    DFi_Dy_i    = n_i*(mu_pi-mu_ni) - DP_i_Dy_i &
                + n_i*(one-y_i)*Dmu_ni_Dy_i + n_i*y_i*Dmu_pi_Dy_i
    DFi_Dn_i    = (one-y_i)*mu_ni + y_i*mu_pi - DP_i_Dn_i &
                + n_i*(one-y_i)*Dmu_ni_Dn_i + n_i*y_i*Dmu_pi_Dn_i
    DFi_DT      = n_i*(one-y_i)*Dmu_ni_DT + n_i*y_i*Dmu_pi_DT - DP_i_DT
!   Inside nucleons total internal energy and its partial derivatives
    Ei          = HC2/TWO*(tau_ni/Meff_ni+tau_pi/Meff_pi) &
                + A*n_i**TWO + FOUR*B*n_ni*n_pi         &
                + DOT_PRODUCT(C,n_i**(ONE+delta)) &
                + DOT_PRODUCT(FOUR*D,n_i**(delta-ONE)*n_ni*n_pi) &
                - n_pi*Neut_Prot_Mass_diff

    DEi_Dy_i    = HC2/TWO*(Dtau_pi_Dy_i/Meff_pi + Dtau_ni_Dy_i/Meff_ni) &
                + n_i*(alpha1 - alpha2)*tau_pi + n_i*(alpha2 - alpha1)*tau_ni &
                + FOUR*B*(ONE-two*y_i)*n_i**TWO - n_i*Neut_Prot_Mass_diff &
                + FOUR*(ONE-TWO*y_i)*DOT_PRODUCT(D,n_i**(delta+ONE))

    DEi_Dn_i    = HC2/TWO*(Dtau_pi_Dn_i/Meff_pi + Dtau_ni_Dn_i/Meff_ni) &
                + (alpha1*y_i + alpha2*(one-y_i))*tau_pi &
                + (alpha2*y_i + alpha1*(one-y_i))*tau_ni &
                + TWO*(A + FOUR*B*y_i*(one-y_i))*n_i - y_i*Neut_Prot_Mass_diff &
                + DOT_PRODUCT(C*(delta+ONE),n_i**delta) &
                + FOUR*Y_I*(ONE-y_i)*DOT_PRODUCT(D*(delta+ONE),n_i**delta)

    DEi_DT      = HC2/TWO*(Dtau_ni_DT/Meff_ni+Dtau_pi_DT/Meff_pi)

!   Inside nucleons total entropy and its partial derivatives
    Si          = (Ei - Fi)/T
    DSi_Dy_i    = (DEi_Dy_i - DFi_Dy_i)/T
    DSi_Dn_i    = (DEi_Dn_i - DFi_Dn_i)/T
    DSi_DT      = (DEi_DT   - DFi_DT)  /T - (Ei - Fi)/T**TWO

!   Total Translational free energy and partial derivatives
    Ftr       = F_TR
    DFtr_Du   = DF_TR_Du
    DFtr_Dy_i = DF_TR_Dy_i
    DFtr_Dn_i = DF_TR_Dn_i
    DFtr_DT   = DF_TR_DT
!   Total Translational internal energy and partial derivatives
    Etr       = Ftr*(ONE  - T*DF_TR_DT/F_TR)
    DEtr_Du   = DFtr_DU   - T*D2F_TR_DuDT
    DEtr_Dy_i = DFtr_Dy_i - T*D2F_TR_DTDy_i
    DEtr_Dn_i = DFtr_Dn_i - T*D2F_TR_Dn_iDT
    DEtr_DT   =           - T*D2F_TR_DT2
!   Total Translational entropy and partial derivatives
    Str        = (Etr       - Ftr)      /T
    DStr_Du    = (DEtr_Du   - DFtr_Du)  /T
    DStr_Dy_i  = (DEtr_Dy_i - DFtr_Dy_i)/T
    DStr_Dn_i  = (DEtr_Dn_i - DFtr_Dn_i)/T
    DStr_DT    = (DEtr_DT   - DFtr_DT)  /T - (Etr - Ftr)/T**TWO

!   Total Surface+Coulomb free energy and partial derivatives
    Fsc       = F_SC
    DFsc_Du   = DF_SC_Du
    DFsc_Dy_i = DF_SC_Dy_i
    DFsc_Dn_i = DF_SC_Dn_i
    DFsc_DT   = DF_SC_DT
!   Total Surface+Coulomb internal energy and partial derivatives
    Esc       = Fsc*(ONE  - T*DF_SC_DT/F_SC)
    DEsc_Du   = DFsc_DU   - T*D2F_SC_DuDT
    DEsc_Dy_i = DFsc_Dy_i - T*D2F_SC_DTDy_i
    DEsc_Dn_i = DFsc_Dn_i - T*D2F_SC_Dn_iDT
    DEsc_DT   =           - T*D2F_SC_DT2
!   Total Surface+Coulomb entropy and partial derivatives
    Ssc        = (Esc       - Fsc)      /T
    DSsc_Du    = (DEsc_Du   - DFsc_Du)  /T
    DSsc_Dy_i  = (DEsc_Dy_i - DFsc_Dy_i)/T
    DSsc_Dn_i  = (DEsc_Dn_i - DFsc_Dn_i)/T
    DSsc_DT    = (DEsc_DT   - DFsc_DT)  /T - (Esc - Fsc)/T**TWO

!   total heavy nuclei free energy and entropy
    F_h = Ftr + Fsc + u*Fi
    S_h = Str + Ssc + u*Si

!   total free energy first derivatives w.r.t. (T,n,y)
    DF_h_DT = - u*Si - Str - Ssc

    DS_h_DT = Si*Du_DT + u*(DSi_DT + DSi_Dn_i*Dn_i_DT + DSi_Dy_i*Dy_i_DT)  &
            + DSsc_DT + DSsc_Dn_i*Dn_i_DT + DSsc_Dy_i*Dy_i_DT + DSsc_Du*Du_DT &
            + DStr_DT + DStr_Dn_i*Dn_i_DT + DStr_Dy_i*Dy_i_DT + DStr_Du*Du_DT

    DF_h_Dn = Fi*Du_Dn + u*(DFi_Dn_i*Dn_i_Dn + DFi_Dy_i*Dy_i_Dn)  &
            + DFsc_Dn_i*Dn_i_Dn + DFsc_Dy_i*Dy_i_Dn + DFsc_Du*Du_Dn &
            + DFtr_Dn_i*Dn_i_Dn + DFtr_Dy_i*Dy_i_Dn + DFtr_Du*Du_Dn

    DS_h_Dn = Si*Du_Dn + u*(DSi_Dn_i*Dn_i_Dn + DSi_Dy_i*Dy_i_Dn)  &
            + DSsc_Dn_i*Dn_i_Dn + DSsc_Dy_i*Dy_i_Dn + DSsc_Du*Du_Dn &
            + DStr_Dn_i*Dn_i_Dn + DStr_Dy_i*Dy_i_Dn + DStr_Du*Du_Dn
!
    DF_h_Dy = Fi*Du_Dy + u*(DFi_Dn_i*Dn_i_Dy + DFi_Dy_i*Dy_i_Dy)  &
            + DFsc_Dn_i*Dn_i_Dy + DFsc_Dy_i*Dy_i_Dy + DFsc_Du*Du_Dy &
            + DFtr_Dn_i*Dn_i_Dy + DFtr_Dy_i*Dy_i_Dy + DFtr_Du*Du_Dy

    DS_h_Dy = Si*Du_Dy + u*(DSi_Dn_i*Dn_i_Dy + DSi_Dy_i*Dy_i_Dy)  &
            + DSsc_Dn_i*Dn_i_Dy + DSsc_Dy_i*Dy_i_Dy + DSsc_Du*Du_Dy &
            + DStr_Dn_i*Dn_i_Dy + DStr_Dy_i*Dy_i_Dy + DStr_Du*Du_Dy

  END SUBROUTINE HEAVY_PROPERTIES

END MODULE Heavy_Properties_Mod
