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
MODULE Outside_Properties_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ONE, TWO, FOUR, R_3_2, &
      v_alpha, b_alpha, HC2 => Hbarc_Square, Neut_Prot_Mass_Diff
  USE Skyrme_Coefficients_Mod, ONLY : A => Coeff_a, B => Coeff_b, C => Coeff_c,&
                  D => Coeff_d, alpha1 => Coeff_alpha1, alpha2 => Coeff_alpha2,&
                  delta => Coeff_delta
  IMPLICIT NONE

CONTAINS

  SUBROUTINE OUTSIDE_PROPERTIES(u, T, n_no, n_po, mu_no, mu_po, tau_no, tau_po,&
    Meff_no, Meff_po, n_alpha, mu_alpha, P_o,  &
    Dn_no_Deta_no, Dn_no_Deta_po, Dn_no_DT,    &
    Dn_po_Deta_no, Dn_po_Deta_po, Dn_po_DT,    &
    Dmu_no_Deta_no, Dmu_no_Deta_po, Dmu_no_DT, &
    Dmu_po_Deta_no, Dmu_po_Deta_po, Dmu_po_DT, &
    Dtau_no_Deta_no, Dtau_no_Deta_po, Dtau_no_DT, &
    Dtau_po_Deta_no, Dtau_po_Deta_po, Dtau_po_DT, &
    Dn_alpha_Deta_no, Dn_alpha_Deta_po, Dn_alpha_DT,    &
    Dmu_alpha_Deta_no, Dmu_alpha_Deta_po, Dmu_alpha_DT, &
    DP_o_Deta_no, DP_o_Deta_po, DP_o_DT, Du_DT, Du_Dn, Du_Dy, &
    Deta_no_DT, Deta_no_Dn, Deta_no_Dy,     &
    Deta_po_DT, Deta_po_Dn, Deta_po_Dy,     &
    F_out, DF_out_DT, DF_out_Dn, DF_out_Dy, &
    S_out, DS_out_DT, DS_out_Dn, DS_out_Dy, &
    Dmu_nout_DT, Dmu_nout_Dn, Dmu_nout_Dy,  &
    Dmu_pout_DT, Dmu_pout_Dn, Dmu_pout_Dy )

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: u, T
    REAL(DP), INTENT(IN) :: n_no, n_po, mu_no, mu_po
    REAL(DP), INTENT(IN) :: tau_no, tau_po, Meff_no, Meff_po
    REAL(DP), INTENT(IN) :: n_alpha, mu_alpha, P_o
    REAL(DP), INTENT(IN) :: Dn_no_Deta_no, Dn_no_Deta_po, Dn_no_DT
    REAL(DP), INTENT(IN) :: Dn_po_Deta_no, Dn_po_Deta_po, Dn_po_DT
    REAL(DP), INTENT(IN) :: Dmu_no_Deta_no, Dmu_no_Deta_po, Dmu_no_DT
    REAL(DP), INTENT(IN) :: Dmu_po_Deta_no, Dmu_po_Deta_po, Dmu_po_DT
    REAL(DP), INTENT(IN) :: Dtau_no_Deta_no, Dtau_no_Deta_po, Dtau_no_DT
    REAL(DP), INTENT(IN) :: Dtau_po_Deta_no, Dtau_po_Deta_po, Dtau_po_DT
    REAL(DP), INTENT(IN) :: Dn_alpha_Deta_no, Dn_alpha_Deta_po, Dn_alpha_DT
    REAL(DP), INTENT(IN) :: Dmu_alpha_Deta_no, Dmu_alpha_Deta_po, Dmu_alpha_DT
    REAL(DP), INTENT(IN) :: DP_o_Deta_no, DP_o_Deta_po, DP_o_DT
    REAL(DP), INTENT(IN) :: Du_DT, Du_Dn, Du_Dy
    REAL(DP), INTENT(IN) :: Deta_no_DT, Deta_no_Dn, Deta_no_Dy
    REAL(DP), INTENT(IN) :: Deta_po_DT, Deta_po_Dn, Deta_po_Dy

    REAL(DP), INTENT(OUT) :: F_out, DF_out_DT, DF_out_Dn, DF_out_Dy
    REAL(DP), INTENT(OUT) :: S_out, DS_out_DT, DS_out_Dn, DS_out_Dy
    REAL(DP), INTENT(OUT) :: Dmu_nout_DT, Dmu_nout_Dn, Dmu_nout_Dy
    REAL(DP), INTENT(OUT) :: Dmu_pout_DT, Dmu_pout_Dn, Dmu_pout_Dy

    REAL(DP) :: dummy, n_o, odu, exc_v_alpha, u_out
    REAL(DP) :: Fa, Ea, Sa, Fo, Eo, So, E_out
    REAL(DP) :: DFa_Deta_po, DEa_Deta_po, DSa_Deta_po
    REAL(DP) :: DFa_Deta_no, DEa_Deta_no, DSa_Deta_no
    REAL(DP) :: DFa_DT,      DEa_DT,      DSa_DT
    REAL(DP) :: DFo_Deta_po, DEo_Deta_po, DSo_Deta_po
    REAL(DP) :: DFo_Deta_no, DEo_Deta_no, DSo_Deta_no
    REAL(DP) :: DFo_DT,      DEo_DT,      DSo_DT
    REAL(DP) :: Dn_a_DT, Dn_a_Dn, Dn_a_Dy

    n_o = n_no + n_po

    odu = one - u
    exc_v_alpha = ONE - v_alpha*n_alpha
    u_out = odu*exc_v_alpha


!   Outside nucleons total free energy and its partial derivatives
    Fo          = n_no*mu_no+n_po*mu_po-P_o
    DFo_Deta_no = n_no*Dmu_no_Deta_no + Dn_no_Deta_no*mu_no &
                + n_po*Dmu_po_Deta_no + Dn_po_Deta_no*mu_po - DP_o_Deta_no
    DFo_Deta_po = n_no*Dmu_no_Deta_po + Dn_no_Deta_po*mu_no &
                + n_po*Dmu_po_Deta_po + Dn_po_Deta_po*mu_po - DP_o_Deta_po
    DFo_DT      = n_no*Dmu_no_DT + Dn_no_DT*mu_no &
                + n_po*Dmu_po_DT + Dn_po_DT*mu_po - DP_o_DT

!   Outside nucleons total internal energy and its partial derivatives
    Eo          = HC2/TWO*(tau_no/Meff_no+tau_po/Meff_po) &
                + A*n_o**TWO + FOUR*B*n_no*n_po         &
                + DOT_PRODUCT(C,n_o**(ONE+delta)) &
                + DOT_PRODUCT(FOUR*D,n_o**(delta-ONE)*n_no*n_po) &
                - n_po*Neut_Prot_Mass_diff

    DEo_Deta_po = HC2/TWO*(Dtau_po_Deta_po/Meff_po + Dtau_no_Deta_po/Meff_no) &
                + (alpha1*Dn_po_Deta_po + alpha2*Dn_no_Deta_po)*tau_po &
                + (alpha1*Dn_no_Deta_po + alpha2*Dn_po_Deta_po)*tau_no &
                + (Dn_no_Deta_po + Dn_po_Deta_po)*( TWO*A*n_o        &
                + DOT_PRODUCT(C*(ONE+delta),n_o**delta)              &
                + DOT_PRODUCT(FOUR*D*(delta-ONE),n_o**(delta-TWO)*n_no*n_po)) &
                + (n_no*Dn_po_Deta_po+n_po*Dn_no_Deta_po)*(FOUR*B &
                + DOT_PRODUCT(FOUR*D,n_o**(delta-ONE))) &
                - Neut_Prot_Mass_diff*Dn_po_Deta_po

    DEo_Deta_no = HC2/TWO*(Dtau_po_Deta_no/Meff_po + Dtau_no_Deta_no/Meff_no) &
                + (alpha1*Dn_po_Deta_no + alpha2*Dn_no_Deta_no)*tau_po &
                + (alpha1*Dn_no_Deta_no + alpha2*Dn_po_Deta_no)*tau_no &
                + (Dn_no_Deta_no + Dn_po_Deta_no)*( TWO*A*n_o        &
                + DOT_PRODUCT(C*(ONE+delta),n_o**delta)              &
                + DOT_PRODUCT(FOUR*D*(delta-ONE),n_o**(delta-TWO)*n_no*n_po) ) &
                + (n_no*Dn_po_Deta_no+n_po*Dn_no_Deta_no)*( FOUR*B &
                + DOT_PRODUCT(FOUR*D,n_o**(delta-ONE)) ) &
                - Neut_Prot_Mass_diff*Dn_po_Deta_no

    DEo_DT      = HC2/TWO*(Dtau_no_DT/Meff_no+Dtau_po_DT/Meff_po) &
                + tau_po*(alpha1*Dn_po_DT + alpha2*Dn_no_DT) &
                + tau_no*(alpha1*Dn_no_DT + alpha2*Dn_po_DT) &
                + (Dn_no_DT + Dn_po_DT)*( TWO*A*n_o        &
                + DOT_PRODUCT(C*(ONE+delta),n_o**delta)              &
                + DOT_PRODUCT(FOUR*D*(delta-ONE),n_o**(delta-TWO)*n_no*n_po) ) &
                + (n_no*Dn_po_DT+n_po*Dn_no_DT)*( FOUR*B &
                + DOT_PRODUCT(FOUR*D,n_o**(delta-ONE)) ) &
                - Neut_Prot_Mass_diff*Dn_po_DT

!   Outside nucleons total entropy and its partial derivatives
    So          = (Eo - Fo)/T
    DSo_Deta_po = (DEo_Deta_po - DFo_Deta_po)/T
    DSo_Deta_no = (DEo_Deta_no - DFo_Deta_no)/T
    DSo_DT      = (DEo_DT      - DFo_DT     )/T - (Eo - Fo)/T**TWO

!   Alpha particle total free energy and partial derivatives
    dummy       = mu_alpha-b_alpha-T
    Fa          = n_alpha*dummy
    DFa_Deta_po = Dn_alpha_Deta_po*dummy + n_alpha*Dmu_alpha_Deta_po
    DFa_Deta_no = Dn_alpha_Deta_no*dummy + n_alpha*Dmu_alpha_Deta_no
    DFa_DT      = Dn_alpha_DT*dummy + n_alpha*(Dmu_alpha_DT-ONE)
!    Alpha particle internal energy and partial derivatives
    dummy       = R_3_2*T - b_alpha
    Ea          = n_alpha*dummy
    DEa_Deta_po = Dn_alpha_Deta_po*dummy
    DEa_Deta_no = Dn_alpha_Deta_no*dummy
    DEa_DT      = Dn_alpha_DT*dummy + R_3_2*n_alpha
!   Alpha particle entropy and partial derivatives
    Sa          = (Ea - Fa)/T
    DSa_Deta_po = (DEa_Deta_po - DFa_Deta_po)/T
    DSa_Deta_no = (DEa_Deta_no - DFa_Deta_no)/T
    DSa_DT      = (DEa_DT     - DFa_DT)      /T - (Ea - Fa)/T**TWO

!   total alpha particle free energy and entropy
    F_out = odu*(Fa + exc_v_alpha*Fo)
    S_out = odu*(Sa + exc_v_alpha*So)

!   total free energy first derivatives w.r.t. (T,n,y)
    Dn_a_DT = Dn_alpha_DT + Dn_alpha_Deta_no*Deta_no_DT &
                          + Dn_alpha_Deta_po*Deta_po_DT

    DF_out_DT = - S_out
    DS_out_DT = - Du_DT*(exc_v_alpha*So + Sa) - odu*v_alpha*Dn_a_DT*So &
            + odu * (DSa_DT + DSa_Deta_no*Deta_no_DT + DSa_Deta_po*Deta_po_DT) &
            + u_out*(DSo_DT + DSo_Deta_no*Deta_no_DT + DSo_Deta_po*Deta_po_DT)
!
    Dn_a_Dn = Dn_alpha_Deta_no*Deta_no_Dn + Dn_alpha_Deta_po*Deta_po_Dn

    DF_out_Dn = - Du_Dn*(exc_v_alpha*Fo + Fa) - odu*v_alpha*Dn_a_Dn*Fo &
              + odu * (DFa_Deta_no*Deta_no_Dn + DFa_Deta_po*Deta_po_Dn) &
              + u_out*(DFo_Deta_no*Deta_no_Dn + DFo_Deta_po*Deta_po_Dn)

    DS_out_Dn = - Du_Dn*(exc_v_alpha*So + Sa) - odu*v_alpha*Dn_a_Dn*So &
              + odu * (DSa_Deta_no*Deta_no_Dn + DSa_Deta_po*Deta_po_Dn) &
              + u_out*(DSo_Deta_no*Deta_no_Dn + DSo_Deta_po*Deta_po_Dn)
!
    Dn_a_Dy = Dn_alpha_Deta_no*Deta_no_Dy + Dn_alpha_Deta_po*Deta_po_Dy

    DF_out_Dy = - Du_Dy*(exc_v_alpha*Fo + Fa) - odu*v_alpha*Dn_a_Dy*Fo &
              + odu * (DFa_Deta_no*Deta_no_Dy + DFa_Deta_po*Deta_po_Dy) &
              + u_out*(DFo_Deta_no*Deta_no_Dy + DFo_Deta_po*Deta_po_Dy)

    DS_out_Dy = - Du_Dy*(exc_v_alpha*So + Sa) - odu*v_alpha*Dn_a_Dy*So &
              + odu * (DSa_Deta_no*Deta_no_Dy + DSa_Deta_po*Deta_po_Dy) &
              + u_out*(DSo_Deta_no*Deta_no_Dy + DSo_Deta_po*Deta_po_Dy)

!   chemical potential first derivatives w.r.t. (T,n,y)
    Dmu_nout_DN = Dmu_no_Deta_no*Deta_no_DN + Dmu_no_Deta_po*Deta_po_DN
    Dmu_pout_DN = Dmu_po_Deta_no*Deta_no_DN + Dmu_po_Deta_po*Deta_po_DN
    Dmu_nout_DY = Dmu_no_Deta_no*Deta_no_DY + Dmu_no_Deta_po*Deta_po_DY
    Dmu_pout_DY = Dmu_po_Deta_no*Deta_no_DY + Dmu_po_Deta_po*Deta_po_DY
    Dmu_nout_DT = Dmu_no_Deta_no*Deta_no_DT + Dmu_no_Deta_po*Deta_po_DT &
                + Dmu_no_DT
    Dmu_pout_DT = Dmu_po_Deta_no*Deta_no_DT + Dmu_po_Deta_po*Deta_po_DT &
                + Dmu_po_DT

  END SUBROUTINE OUTSIDE_PROPERTIES

END MODULE Outside_Properties_Mod
