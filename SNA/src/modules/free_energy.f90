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
MODULE Free_Energy_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL
  USE Physical_Constants_Mod, ONLY : ZERO, HALF, ONE, TWO, FOUR, TEN, R_3_2, &
      v_alpha, b_alpha, Hbarc_Square, PI, Mass_n
  USE Skyrme_Bulk_Mod
  USE Skyrme_Bulk_Density_Derivatives_Mod
  USE Skyrme_Bulk_Eta_Derivatives_Mod
  USE Skyrme_Bulk_Temperature_Derivatives_Mod
  USE Find_Uniform_Matter_Derivatives_Mod
  USE Find_Non_Uniform_Matter_Derivatives_Mod
  USE Alpha_Properties_Mod
  USE Free_Energy_Heavy_Mod
  USE Outside_Properties_Mod
  USE Heavy_Properties_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE FREE_ENERGY( log10_n_no, log10_n_po, log10_u, n, T, Yp, &
  F_o, F_i, F_alpha, F_TR, F_SC, &
  n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
  A, Z, rad, F, P, S, E, &
  DF_Dn, DF_Dy, DF_DT, DP_Dn, DP_Dy, DP_DT, &
  DS_Dn, DS_Dy, DS_DT, DE_Dn, DE_Dy, DE_DT, &
  mu_no, mu_po, Meff_no, Meff_po, &
  Dmu_no_DT, Dmu_no_Dn, Dmu_no_Dy, &
  Dmu_po_DT, Dmu_po_Dn, Dmu_po_Dy, &
  dLog10_n_no_dn, dLog10_n_po_dn, dLog10_u_dn, &
  dLog10_n_no_dT, dLog10_n_po_dT, dLog10_u_dT, &
  non_uniform_eq, non_uniform_jac, objective, error )

    IMPLICIT NONE
!
!   if objective = 0 only calculate free energies and densities
!                     used for solution of uniform system
!   if objective = 3 adds to calculation free energy
!                    first and second derivatives for uniform matter
!   if objective = 1 adds to calculation values for the
!                     non-uniform matter system of equations
!   if objective = 2 adds to calculation values for the
!                     analytical jacobian for non-uniform system of equations
!   if objective = 4 adds to calculation free energy
!                    first and second derivatives for non-uniform matter
!
!   error = 0 if all physical constranits are satisfied and u=0
!   error = 1 if all physical constranits are satisfied and u>0
!   error < 0 some physical constraint is not satisfied

    REAL(DP), INTENT(IN) :: log10_n_no, log10_n_po, log10_u
    REAL(DP), INTENT(IN) :: n, T, Yp
    REAL(DP), INTENT(OUT) :: F_o, F_i, F_alpha, F_TR, F_SC
    REAL(DP), INTENT(OUT) :: n_no, n_po, n_ni, n_pi, n_alpha, n_heavy
    REAL(DP), INTENT(OUT) :: dLog10_n_no_dn, dLog10_n_po_dn, dLog10_u_dn
    REAL(DP), INTENT(OUT) :: dLog10_n_no_dT, dLog10_n_po_dT, dLog10_u_dT

    REAL(DP), DIMENSION(3), INTENT(OUT) :: non_uniform_eq
    REAL(DP), DIMENSION(3,3), INTENT(OUT) :: non_uniform_jac

    INTEGER(I4B), INTENT(IN) :: objective
    INTEGER(I4B), INTENT(OUT) :: error

    REAL(DP), INTENT(OUT) :: A, Z, rad

    REAL(DP), INTENT(OUT) :: F, P, S, E
    REAL(DP), INTENT(OUT) :: DF_Dn, DF_Dy, DF_DT
    REAL(DP), INTENT(OUT) :: DP_Dn, DP_Dy, DP_DT
    REAL(DP), INTENT(OUT) :: DS_Dn, DS_Dy, DS_DT
    REAL(DP), INTENT(OUT) :: DE_Dn, DE_Dy, DE_DT

    REAL(DP), INTENT(OUT) :: Dmu_no_DT, Dmu_no_Dn, Dmu_no_Dy
    REAL(DP), INTENT(OUT) :: Dmu_po_DT, Dmu_po_Dn, Dmu_po_Dy

    REAL(DP) :: D2F_DT2, D2F_DTDn, D2F_DyDT, D2F_DnDy, D2F_Dn2, D2F_DY2

    REAL(DP) :: u, odu, B(3), JAC(3,3)
    ! inside region properties
    REAL(DP) :: log10_n_ni,log10_n_pi
    REAL(DP) :: n_i,y_i,P_i
    REAL(DP) :: Meff_ni,Meff_pi
    REAL(DP) :: eta_ni,eta_pi
    REAL(DP) :: tau_ni,tau_pi
    REAL(DP) :: v_ni,v_pi
    REAL(DP) :: mu_ni,mu_pi
    REAL(DP) :: G_ni,G_pi
    ! outside region properties
    REAL(DP) :: n_o,y_o,P_o
    REAL(DP) :: Meff_no,Meff_po
    REAL(DP) :: eta_no,eta_po
    REAL(DP) :: tau_no,tau_po
    REAL(DP) :: v_no,v_po
    REAL(DP) :: mu_no,mu_po
    REAL(DP) :: G_no,G_po
    ! derivatives of inside nucleon properties w.r.t. inside nucleon densities
    REAL(DP) :: Deta_ni_Dn_ni, Deta_ni_Dn_pi, Deta_pi_Dn_ni, Deta_pi_Dn_pi
    REAL(DP) :: Dtau_ni_Dn_ni, Dtau_ni_Dn_pi, Dtau_pi_Dn_ni, Dtau_pi_Dn_pi
    REAL(DP) :: Dmu_ni_Dn_ni,  Dmu_ni_Dn_pi,  Dmu_pi_Dn_ni,  Dmu_pi_Dn_pi
    REAL(DP) :: Dv_ni_Dn_ni,   Dv_ni_Dn_pi,   Dv_pi_Dn_ni,   Dv_pi_Dn_pi
    REAL(DP) :: DP_i_Dn_ni,    DP_i_Dn_pi
    REAL(DP) :: DP_i_Dn_i,     DP_i_Dy_i
    ! derivatives of outside nucleon properties w.r.t. outside nucleon densities
    REAL(DP) :: Deta_no_Dn_no, Deta_no_Dn_po, Deta_po_Dn_no, Deta_po_Dn_po
    REAL(DP) :: Dtau_no_Dn_no, Dtau_no_Dn_po, Dtau_po_Dn_no, Dtau_po_Dn_po
    REAL(DP) :: Dmu_no_Dn_no,  Dmu_no_Dn_po,  Dmu_po_Dn_no,  Dmu_po_Dn_po
    REAL(DP) :: Dv_no_Dn_no,   Dv_no_Dn_po,   Dv_po_Dn_no,   Dv_po_Dn_po
    REAL(DP) :: DP_o_Dn_no,      DP_o_Dn_po
    ! derivatives of outside nucleon properties w.r.t. nucleon degeneracy param.
    REAL(DP) :: Dn_no_Deta_no,Dn_no_Deta_po,Dn_po_Deta_no,Dn_po_Deta_po
    REAL(DP) :: Dtau_no_Deta_no,Dtau_no_Deta_po,Dtau_po_Deta_no,Dtau_po_Deta_po
    REAL(DP) :: Dmu_no_Deta_no, Dmu_no_Deta_po, Dmu_po_Deta_no, Dmu_po_Deta_po
    REAL(DP) :: DP_o_Deta_no,   DP_o_Deta_po
    ! derivatives of outside nucleon properties w.r.t. temperature
    REAL(DP) :: Dn_no_DT, Dn_po_DT, DP_o_DT
    REAL(DP) :: Dtau_no_DT, Dtau_po_DT
    ! derivatives of inside nucleon properties w.r.t. temperature
    REAL(DP) :: Dn_ni_DT, Dn_pi_DT, DP_i_DT
    REAL(DP) :: Dtau_ni_DT, Dtau_pi_DT, Dmu_ni_DT, Dmu_pi_DT
    ! alpha particle related quantities and properties
    REAL(DP) :: mass_Q, mu_alpha, P_alpha
    REAL(DP) :: exc_v_alpha
    ! alpha particle derivatives w.r.t. independent variable u
    REAL(DP) :: Dn_alpha_Du
    ! derivatives of alpha particle properties w.r.t. outside nucleon densities
    REAL(DP) :: Dn_alpha_Dn_no, Dn_alpha_Dn_po
    REAL(DP) :: DP_alpha_Dn_no, DP_alpha_Dn_po
    REAL(DP) :: Dmu_alpha_Dn_no, Dmu_alpha_Dn_po
    ! derivatives of alpha particle properties w.r.t. outside degen parameter
    REAL(DP) :: Dn_alpha_Deta_no, Dn_alpha_Deta_po
    REAL(DP) :: DP_alpha_Deta_no, DP_alpha_Deta_po
    REAL(DP) :: Dmu_alpha_Deta_no, Dmu_alpha_Deta_po
    ! derivatives of inside nucleon properties w.r.t. independent variable u
    REAL(DP) :: Dn_i_Du, Dy_i_Du
    ! derivatives of inside nucleon densities w.r.t. inside nucleon properties
    REAL(DP) :: Dmu_ni_Dn_i,  Dmu_ni_Dy_i,  Dmu_pi_Dn_i,  Dmu_pi_Dy_i
    REAL(DP) :: Dtau_ni_Dn_i, Dtau_ni_Dy_i, Dtau_pi_Dn_i, Dtau_pi_Dy_i
    ! derivatives of inside nucleon properties w.r.t. outside nucleon densities
    REAL(DP) :: Dn_i_Dn_no, Dn_i_Dn_po
    REAL(DP) :: Dy_i_Dn_no, Dy_i_Dn_po
    ! derivatives of pressure inside w.r.t. outside nucleon densities
    REAL(DP) :: DP_i_Dn_no, DP_i_Dn_po
    REAL(DP) :: DP_i_Dn_o, DP_i_Dy_o
    ! partial derivatives of inside variables w.r.t. independent variables
    REAL(DP) :: DA_in_Du(1:3), DA_in_Dn_no(1:3), DA_in_Dn_po(1:3)
    ! partial derivatives of outside variables w.r.t. independent variables
    REAL(DP) :: DA_out_Du(1:3), DA_out_Dn_no(1:3), DA_out_Dn_po(1:3)
    ! partial derivatives of equilibrium eqs. w.r.t. independent variables
    REAL(DP) :: DB_Dup(1:3), DB_Dn_i(1:3), DB_Dy_i(1:3), DB_DT(1:3)
    ! full derivatives of equilibrium eqs. w.r.t. independent variables
    REAL(DP) :: DB_Du(1:3), DB_Dn_no(1:3), DB_Dn_po(1:3)
    ! total density of nucleons and protons outside
    ! (includes nucleons in alphas)
    REAL(DP) :: tot_n_o, tot_n_po
    ! some auxiliary densities
    REAL(DP) :: u_out, n_1, n_2
    ! retivatives of w.r.t.
    REAL(DP) :: Dn_alpha_DT, Dmu_alpha_DT, DP_alpha_DT
    ! derivatives of translational and surface+coulomb free energies
    REAL(DP) :: DF_SC_Du, DF_SC_Dn_i, DF_SC_DT, DF_SC_Dy_i
    REAL(DP) :: D2F_SC_Du2, D2F_SC_DuDn_i, D2F_SC_DuDT, D2F_SC_DuDy_i
    REAL(DP) :: D2F_SC_Dn_i2, D2F_SC_Dn_iDT, D2F_SC_Dn_iDy_i
    REAL(DP) :: D2F_SC_DT2, D2F_SC_DTDy_i, D2F_SC_Dy_i2
    REAL(DP) :: DF_TR_Du, DF_TR_Dn_i, DF_TR_DT, DF_TR_Dy_i
    REAL(DP) :: D2F_TR_Du2, D2F_TR_DuDn_i, D2F_TR_DuDT, D2F_TR_DuDy_i
    REAL(DP) :: D2F_TR_Dn_i2, D2F_TR_Dn_iDT, D2F_TR_Dn_iDy_i
    REAL(DP) :: D2F_TR_DT2, D2F_TR_DTDy_i, D2F_TR_Dy_i2
    ! total derivatives of variables (u,n_i,y_i,eta_n,eta_p) w.r.t. (n,y,T)
    REAL(DP) :: DU_DN, Dy_i_DN, Dn_i_DN, Deta_no_DN, Deta_po_DN
    REAL(DP) :: DU_DT, Dy_i_DT, Dn_i_DT, Deta_no_DT, Deta_po_DT
    REAL(DP) :: DU_DY, Dy_i_DY, Dn_i_DY, Deta_no_DY, Deta_po_DY
    ! total outside derivatives of free energy, entropy and chemical potentials
    REAL(DP) :: F_out, DF_out_DT, DF_out_Dn, DF_out_Dy
    REAL(DP) :: S_out, DS_out_DT, DS_out_Dn, DS_out_Dy
    ! total heavy nuclei derivatives of free energy and entropy
    REAL(DP) :: F_h, DF_h_DT, DF_h_Dn, DF_h_Dy
    REAL(DP) :: S_h, DS_h_DT, DS_h_Dn, DS_h_Dy

    REAL(DP) ::  DF_o_DT, DF_o_Dn, DF_o_Dy
    REAL(DP) ::  DS_o_DT, DS_o_Dn, DS_o_Dy

!   set output to zero
    F_o = 1.d100 ; F_i = 1.d100; F_alpha = 1.d100
    F_TR = 1.d100; F_SC = 1.d100
    n_no = 1.d100; n_po = 1.d100
    n_ni = 1.d100; n_pi = 1.d100
    n_alpha = 1.d100; n_heavy = 1.d100
    dLog10_n_no_dn = 1.d100; dLog10_n_po_dn = 1.d100; dLog10_u_dn = 1.d100
    dLog10_n_no_dT = 1.d100; dLog10_n_po_dT = 1.d100; dLog10_u_dT = 1.d100

    F = 1.d100 ; P = 1.d100 ; S = 1.d100 ; E = 1.d100
    DF_Dn = 1.d100 ; DF_Dy = 1.d100 ; DF_DT = 1.d100
    DP_Dn = 1.d100 ; DP_Dy = 1.d100 ; DP_DT = 1.d100
    DS_Dn = 1.d100 ; DS_Dy = 1.d100 ; DS_DT = 1.d100
    DE_Dn = 1.d100 ; DE_Dy = 1.d100 ; DE_DT = 1.d100
    Dmu_no_DT = 1.d100 ; Dmu_no_Dn = 1.d100 ; Dmu_no_Dy = 1.d100
    Dmu_po_DT = 1.d100 ; Dmu_po_Dn = 1.d100 ; Dmu_po_Dy = 1.d100

    u = TEN**log10_u
    odu = ONE - u
    error = 0

!   check if volume fraction of heavy nuclei not in range 0<=u<=1.
    IF (u > ONE .OR. u < ZERO) THEN
      error = -1
      RETURN
    ENDIF

!   outside nucleons' properties
    CALL SKYRME_BULK_PROPERTIES(log10_n_no,log10_n_po,&
      T,n_no,n_po,n_o,y_o,Meff_no,Meff_po,eta_no,eta_po,tau_no,tau_po, &
      v_no,v_po,mu_no,mu_po,P_o,F_o)

!   check if outside density larger than total density
!    (add small tolerance to sum of n_no and n_po)
    IF (n_o > (one+5.d-2)*n) THEN
      error = -2
      RETURN
    ENDIF

    CALL ALPHA_PROPERTIES (objective, n_no, n_po, mu_no, mu_po, T, F_o, P_o, &
      n_alpha, mu_alpha, F_alpha, P_alpha)

!   chech if alpha particles occupy more than 100% of Wigner-Seiz cell volume
    exc_v_alpha = ONE - v_alpha*n_alpha
    IF (exc_v_alpha < ZERO) THEN
      error = -3
      RETURN
    ENDIF

!   if contribution due to heavy nuclei is very small ignore them
    IF (log10_u < -100.0D0) THEN
      error = 1
      F_i  = ZERO
      F_TR = ZERO
      F_SC = ZERO
      IF (objective == 0) RETURN
    ENDIF

    IF (objective == 2 .or. objective == 3 .or. objective ==4) THEN
      CALL SKYRME_BULK_DENSITY_DERIVATIVES( &
            log10_n_no, log10_n_po, T, G_no, G_po, &
            Deta_no_Dn_no, Deta_no_Dn_po, Deta_po_Dn_no, Deta_po_Dn_po, &
            Dtau_no_Dn_no, Dtau_no_Dn_po, Dtau_po_Dn_no, Dtau_po_Dn_po, &
            Dmu_no_Dn_no, Dmu_no_Dn_po,   Dmu_po_Dn_no, Dmu_po_Dn_po,   &
            Dv_no_Dn_no, Dv_no_Dn_po,     Dv_po_Dn_no, Dv_po_Dn_po,     &
            DP_o_Dn_no, DP_o_Dn_po)

      CALL SKYRME_BULK_ETA_DERIVATIVES(n, T, Yp, n_no, n_po, Meff_no, Meff_po, &
            eta_no, eta_po, tau_no, tau_po, v_no, v_po, mu_no, mu_po,   &
            Deta_no_Dn_no, Deta_no_Dn_po, Deta_po_Dn_no, Deta_po_Dn_po, &
            Dtau_no_Dn_no, Dtau_no_Dn_po, Dtau_po_Dn_no, Dtau_po_Dn_po, &
            Dmu_no_Dn_no,  Dmu_no_Dn_po,  Dmu_po_Dn_no,  Dmu_po_Dn_po,  &
            Dv_no_Dn_no,   Dv_no_Dn_po,   Dv_po_Dn_no,   Dv_po_Dn_po,   &
            DP_o_Dn_no,    DP_o_Dn_po,    DP_o_Deta_no,  DP_o_Deta_po,  &
            Dn_no_Deta_no,  Dn_no_Deta_po,  Dn_po_Deta_no,  Dn_po_Deta_po,  &
            Dmu_no_Deta_no, Dmu_no_Deta_po, Dmu_po_Deta_no, Dmu_po_Deta_po, &
            Dtau_no_Deta_no,Dtau_no_Deta_po,Dtau_po_Deta_no,Dtau_po_Deta_po )

      CALL SKYRME_BULK_TEMPEARTURE_DERIVATIVES(n_no, n_po, T, Meff_no, Meff_po,&
            eta_no, eta_po, tau_no, tau_po, v_no, v_po, mu_no, mu_po, &
            Dn_no_DT, Dn_po_DT, Dtau_no_DT, Dtau_po_DT, Dmu_no_DT, Dmu_po_DT, &
            DP_o_DT, 0)

      CALL ALPHA_DERIVATIVES  ( n_no, n_po, n_alpha, mu_alpha, T,            &
            Dmu_no_Dn_no,   Dmu_po_Dn_no,   Dmu_no_Dn_po,   Dmu_po_Dn_po,    &
            Dmu_no_Deta_no, Dmu_po_Deta_no, Dmu_no_Deta_po, Dmu_po_Deta_po,  &
            Dmu_no_DT, Dmu_po_DT, DP_o_DT,                                   &
            Dn_no_Deta_no,  Dn_po_Deta_no,  Dn_no_Deta_po,  Dn_po_Deta_po,   &
            DP_o_Dn_no, DP_o_Dn_po, DP_o_Deta_no, DP_o_Deta_po,              &
            Dn_alpha_Dn_no, Dn_alpha_Dn_po, Dmu_alpha_Dn_no, Dmu_alpha_Dn_po,&
            Dn_alpha_Deta_no, Dn_alpha_Deta_po,   &
            Dmu_alpha_Deta_no, Dmu_alpha_Deta_po, &
            DP_alpha_Deta_no, DP_alpha_Deta_po,   &
            Dn_alpha_DT, DP_alpha_DT, Dmu_alpha_DT )
    ENDIF

    IF (objective == 3) THEN

      CALL Find_Uniform_Matter_Derivatives ( n, T, Yp, n_no, n_po, n_alpha,  &
        Dn_no_Deta_no, Dn_po_Deta_no, Dn_no_Deta_po, Dn_po_Deta_po,          &
        Dn_no_DT, Dn_po_DT, Dn_alpha_Deta_no, Dn_alpha_Deta_po, Dn_alpha_DT,   &
        Deta_no_Dn, Deta_po_Dn, Deta_no_DT, Deta_po_DT, Deta_no_Dy, Deta_po_Dy )

      Du_DT = ZERO
      Du_Dn = ZERO
      Du_Dy = ZERO

      CALL OUTSIDE_PROPERTIES (u, T, n_no, n_po, mu_no, mu_po, tau_no, tau_po, &
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
          Dmu_no_DT, Dmu_no_Dn, Dmu_no_Dy,  &
          Dmu_po_DT, Dmu_po_Dn, Dmu_po_Dy )

      F     = F_out
      S     = S_out
      E     = F + T*S

      DF_Dn = DF_out_Dn
      DF_Dy = DF_out_Dy
      DF_DT = DF_out_DT

      DS_Dn = DS_out_Dn
      DS_Dy = DS_out_Dy
      DS_DT = DS_out_DT

      DE_Dn = DF_Dn + T*DS_Dn
      DE_Dy = DF_Dy + T*DS_Dy
      DE_DT = DF_DT + T*DS_DT + S

      D2F_DT2  = - DS_DT
      D2F_DTDn = (ONE-Yp)*Dmu_no_DT + Yp*Dmu_po_DT
      D2F_DyDT = - n*(Dmu_no_DT-Dmu_po_DT)
      D2F_DnDy = - (mu_no-mu_po) - n*(Dmu_no_Dn-Dmu_po_Dn)
      D2F_Dn2  = (ONE-Yp)*Dmu_no_Dn + Yp*Dmu_po_Dn
      D2F_DY2  = - n*(Dmu_no_Dy - Dmu_po_Dy)

      P     = n*DF_Dn - F
      DP_Dn = n*D2F_DN2
      DP_Dy = n*(mu_no - mu_po + D2F_DNDy)
      DP_DT = S + n*D2F_DTDn

      RETURN

    ENDIF

!   if heavy nuclei present obtain their properties
    u_out    = (one-u)*exc_v_alpha
    tot_n_o  =  n_o*exc_v_alpha + four*n_alpha
    tot_n_po = n_po*exc_v_alpha +  two*n_alpha

    n_i = (n-(ONE-u)*tot_n_o)/u
    n_1  = n_i - tot_n_o
    y_i = (n*Yp-(ONE-u)*tot_n_po)
    y_i = y_i/u/n_i
    n_2  = n_i*y_i-tot_n_po
    n_ni = n_i*(ONE-y_i)
    n_pi = n_i*(y_i)

!   check if inside neutron and proton densities smaller than outside densities
    IF (n_i < 1.05d0*n_o) THEN
      error = -3
      RETURN
    ENDIF

!   check if inside density too similar to total density
    IF (n_i < n*1.05D0) THEN
      error = -4
      RETURN
    ENDIF

!   check if proton fraction inside is outside (0,1) range
    IF (y_i <= 0.0D0 .OR. y_i >= ONE) THEN
      error = -5
      RETURN
    ENDIF

    log10_n_ni = log10(n_i*(one-y_i))
    log10_n_pi = log10(n_i*y_i)
!
!   inside nucleons' properties
    CALL SKYRME_BULK_PROPERTIES(log10_n_ni,log10_n_pi,&
      T,n_ni,n_pi,n_i,y_i,Meff_ni,Meff_pi,eta_ni,eta_pi,tau_ni,tau_pi, &
      v_ni,v_pi,mu_ni,mu_pi,P_i,F_i)

    CALL Free_Energy_Heavy( objective, n_ni, n_pi, u, T, A, Z, rad, &
      F_SC, DF_SC_Du, DF_SC_Dn_i, DF_SC_DT, DF_SC_Dy_i, &
      D2F_SC_Du2, D2F_SC_DuDn_i, D2F_SC_DuDT, D2F_SC_DuDy_i, D2F_SC_Dn_i2, &
      D2F_SC_Dn_iDT, D2F_SC_Dn_iDy_i, D2F_SC_DT2, D2F_SC_DTDy_i, D2F_SC_Dy_i2, &
      F_TR, DF_TR_Du, DF_TR_Dn_i, DF_TR_DT, DF_TR_Dy_i, &
      D2F_TR_Du2, D2F_TR_DuDn_i, D2F_TR_DuDT, D2F_TR_DuDy_i, D2F_TR_Dn_i2, &
      D2F_TR_Dn_iDT, D2F_TR_Dn_iDy_i, D2F_TR_DT2, D2F_TR_DTDy_i, D2F_TR_Dy_i2 )

!   if a heavy nucleus is too small, disregard it.
!   if F_TR or F_SC does not exist, disregard solution also
!    this happens whenever a heavy nuclei has a very small y_i (y_i ≲ 0.20)
!    and the surface tension \sigma(y_i,T) == NaN
    IF (A < 1.0d1 .OR. A > 1.0D3) THEN
      error = -6
      RETURN
    ENDIF

    IF (n_i > 0.5d0) THEN
      error = -8
      RETURN
    ENDIF

    IF (ieee_is_nan(F_TR) .OR. ieee_is_nan(F_SC)) THEN
      error = -7
      RETURN
    ENDIF

    n_heavy = (n_ni+n_pi)*TEN**log10_u/A

!   correct F_alpha and F_o due to exluded heavy nuclei volume "u"
    F_alpha = (one-u)*F_alpha
    F_o     = u_out*F_o
    F_i     = u*F_i

!   if OBJECTIVE = 1 then use this subroutine
!   to obtain equations to be minimimized for non_uniform matter
    IF (objective==0) RETURN

    B(1) = DF_TR_Du - n_i/u*DF_TR_Dn_i + DF_SC_Du - n_i/u*DF_SC_Dn_i

    B(2) = y_i/u/n_i*(DF_TR_Dy_i - n_i/y_i*DF_TR_Dn_i) &
         + y_i/u/n_i*(DF_SC_Dy_i - n_i/y_i*DF_SC_Dn_i)

    B(3) = -(ONE-y_i)/u/n_i*(DF_TR_Dy_i + DF_SC_Dy_i) &
                    - ONE/u*(DF_TR_Dn_i + DF_SC_Dn_i)


    IF (objective==1) THEN

      non_uniform_eq(1) = (P_i - P_o - P_alpha) - B(1)
      non_uniform_eq(2) = (mu_ni - mu_no) - B(2)
      non_uniform_eq(3) = (mu_pi - mu_po) - B(3)

      RETURN
    ENDIF

    Dy_i_Du = - (n_2-y_i*n_1)/u/n_i
    Dn_i_Du = - n_1/u

    CALL SKYRME_BULK_DENSITY_DERIVATIVES ( &
          log10_n_no, log10_n_po, T, G_no, G_po, &
          Deta_no_Dn_no, Deta_no_Dn_po, Deta_po_Dn_no, Deta_po_Dn_po, &
          Dtau_no_Dn_no, Dtau_no_Dn_po, Dtau_po_Dn_no, Dtau_po_Dn_po, &
          Dmu_no_Dn_no, Dmu_no_Dn_po,   Dmu_po_Dn_no, Dmu_po_Dn_po,   &
          Dv_no_Dn_no, Dv_no_Dn_po,     Dv_po_Dn_no, Dv_po_Dn_po,     &
          DP_o_Dn_no, DP_o_Dn_po)

    CALL ALPHA_DERIVATIVES  ( n_no, n_po, n_alpha, mu_alpha, T,            &
          Dmu_no_Dn_no,   Dmu_po_Dn_no,   Dmu_no_Dn_po,   Dmu_po_Dn_po,    &
          Dmu_no_Deta_no, Dmu_po_Deta_no, Dmu_no_Deta_po, Dmu_po_Deta_po,  &
          Dmu_no_DT, Dmu_po_DT, DP_o_DT,                                   &
          Dn_no_Deta_no,  Dn_po_Deta_no,  Dn_no_Deta_po,  Dn_po_Deta_po,   &
          DP_o_Dn_no, DP_o_Dn_po, DP_o_Deta_no, DP_o_Deta_po,              &
          Dn_alpha_Dn_no, Dn_alpha_Dn_po, Dmu_alpha_Dn_no, Dmu_alpha_Dn_po,&
          Dn_alpha_Deta_no, Dn_alpha_Deta_po,   &
          Dmu_alpha_Deta_no, Dmu_alpha_Deta_po, &
          DP_alpha_Deta_no, DP_alpha_Deta_po,   &
          Dn_alpha_DT, DP_alpha_DT, Dmu_alpha_DT )

    CALL SKYRME_BULK_DENSITY_DERIVATIVES ( &
          log10_n_ni, log10_n_pi, T, G_ni, G_pi, &
          Deta_ni_Dn_ni, Deta_ni_Dn_pi, Deta_pi_Dn_ni, Deta_pi_Dn_pi, &
          Dtau_ni_Dn_ni, Dtau_ni_Dn_pi, Dtau_pi_Dn_ni, Dtau_pi_Dn_pi, &
          Dmu_ni_Dn_ni, Dmu_ni_Dn_pi,   Dmu_pi_Dn_ni, Dmu_pi_Dn_pi,   &
          Dv_ni_Dn_ni, Dv_ni_Dn_pi,     Dv_pi_Dn_ni, Dv_pi_Dn_pi,     &
          DP_i_Dn_ni, DP_i_Dn_pi)

!   derivatives of total "inside" nucleon density w.r.t. "outside" densities
    Dn_i_Dn_no = -odu/u*(exc_v_alpha + (FOUR-v_alpha*n_o)*Dn_alpha_Dn_no)
    Dn_i_Dn_po = -odu/u*(exc_v_alpha + (FOUR-v_alpha*n_o)*Dn_alpha_Dn_po)

!   derivatives of "inside" nucleon chemical pot. w.r.t. "inside" properties
    Dmu_ni_Dn_i = (ONE-y_i)*Dmu_ni_Dn_ni + y_i*Dmu_ni_Dn_pi
    Dmu_pi_Dn_i = (ONE-y_i)*Dmu_pi_Dn_ni + y_i*Dmu_pi_Dn_pi

!   derivatives of "inside" nucleon chemical pot. w.r.t. "inside" properties
    Dmu_ni_Dy_i = n_i*(Dmu_ni_Dn_pi - Dmu_ni_Dn_ni)
    Dmu_pi_Dy_i = n_i*(Dmu_pi_Dn_pi - Dmu_pi_Dn_ni)

!   derivatives of "inside" nucleon tau w.r.t. "inside" properties
    Dtau_ni_Dn_i = (ONE-y_i)*Dtau_ni_Dn_ni + y_i*Dtau_ni_Dn_pi
    Dtau_pi_Dn_i = (ONE-y_i)*Dtau_pi_Dn_ni + y_i*Dtau_pi_Dn_pi

!   derivatives of "inside" nucleon tau w.r.t. "inside" properties
    Dtau_ni_Dy_i = n_i*(Dtau_ni_Dn_pi - Dtau_ni_Dn_ni)
    Dtau_pi_Dy_i = n_i*(Dtau_pi_Dn_pi - Dtau_pi_Dn_ni)

!   derivatives of total "inside" proton fraction w.r.t. "outside" densities
    Dy_i_Dn_no = (TWO-v_alpha*n_po)*Dn_alpha_Dn_no
    Dy_i_Dn_no = - ( u*y_i*Dn_i_Dn_no + odu*Dy_i_Dn_no)/n_i/u

    Dy_i_Dn_po = exc_v_alpha + (TWO-v_alpha*n_po)*Dn_alpha_Dn_po
    Dy_i_Dn_po = - ( u*y_i*Dn_i_Dn_po + odu*Dy_i_Dn_po)/n_i/u

!   derivative of inside pressure w.rt. inside properties
    DP_i_Dn_i = (ONE-y_i)*DP_i_Dn_ni + y_i*DP_i_Dn_pi
    DP_i_Dy_i = n_i*(DP_i_Dn_pi - DP_i_Dn_ni)

!   derivative of inside pressure w.r.t. outside nucleon densities
    DP_i_Dn_no = DP_i_Dy_i*Dy_i_Dn_no
    DP_i_Dn_po = DP_i_Dy_i*Dy_i_Dn_po
!   partial derivatives of Ai_in w.r.t. independent variables
!   A1_in = P_in ; A2_in = mu_ni ; A3_in = mu_pi
    DA_in_Du(1)    = DP_i_Dn_i*Dn_i_Du    + DP_i_Dy_i*Dy_i_Du
    DA_in_Dn_no(1) = DP_i_Dn_i*Dn_i_Dn_no + DP_i_Dy_i*Dy_i_Dn_no
    DA_in_Dn_po(1) = DP_i_Dn_i*Dn_i_Dn_po + DP_i_Dy_i*Dy_i_Dn_po

    DA_in_Du(2)    = Dmu_ni_Dn_i*Dn_i_Du    + Dmu_ni_Dy_i*Dy_i_Du
    DA_in_Dn_no(2) = Dmu_ni_Dn_i*Dn_i_Dn_no + Dmu_ni_Dy_i*Dy_i_Dn_no
    DA_in_Dn_po(2) = Dmu_ni_Dn_i*Dn_i_Dn_po + Dmu_ni_Dy_i*Dy_i_Dn_po

    DA_in_Du(3)    = Dmu_pi_Dn_i*Dn_i_Du    + Dmu_pi_Dy_i*Dy_i_Du
    DA_in_Dn_no(3) = Dmu_pi_Dn_i*Dn_i_Dn_no + Dmu_pi_Dy_i*Dy_i_Dn_no
    DA_in_Dn_po(3) = Dmu_pi_Dn_i*Dn_i_Dn_po + Dmu_pi_Dy_i*Dy_i_Dn_po
!   partial derivatives of Ai_out w.r.t. independent variables
!   A1_out = P_out ; A2_out = mu_no ; A3_out = mu_po
    DA_out_Du(1)    = ZERO
    DA_out_Dn_no(1) = DP_o_Dn_no +  DP_alpha_Dn_no
    DA_out_Dn_po(1) = DP_o_Dn_po +  DP_alpha_Dn_po
!
    DA_out_Du(2)    = ZERO
    DA_out_Dn_no(2) = Dmu_no_Dn_no
    DA_out_Dn_po(2) = Dmu_no_Dn_po
!
    DA_out_Du(3)    = ZERO
    DA_out_Dn_no(3) = Dmu_po_Dn_no
    DA_out_Dn_po(3) = Dmu_po_Dn_po
!   partial derivatives of B_i w.r.t. independent variables
    DB_Dup(1)     = D2F_TR_Du2 - n_i/u*D2F_TR_DuDn_i &
                  + D2F_SC_Du2 - n_i/u*D2F_SC_DuDn_i &
                  + n_i/u/u*(DF_TR_Dn_i + DF_SC_Dn_i)
    DB_Dn_i(1)    = D2F_TR_DuDn_i + D2F_SC_DuDn_i &
                  - n_i/u*( D2F_TR_Dn_i2 + D2F_SC_Dn_i2 ) &
                  - ( DF_TR_Dn_i + DF_SC_Dn_i )/u
    DB_Dy_i(1)    = D2F_TR_DuDy_i - n_i/u*D2F_TR_Dn_iDy_i &
                  + D2F_SC_DuDy_i - n_i/u*D2F_SC_Dn_iDy_i
    DB_DT(1)      = D2F_TR_DuDT + D2F_SC_DuDT &
                  - n_i/u*( D2F_TR_Dn_iDT + D2F_SC_Dn_iDT )
!
    DB_Dup(2)     = - B(2)/u + (y_i/n_i*D2F_TR_DuDy_i - D2F_TR_DuDn_i)/u &
                             + (y_i/n_i*D2F_SC_DuDy_i - D2F_SC_DuDn_i)/u
    DB_Dn_i(2)    = y_i/u/n_i*( ( D2F_TR_Dn_iDy_i - n_i/y_i*D2F_TR_Dn_i2   &
                                + D2F_SC_Dn_iDy_i - n_i/y_i*D2F_SC_Dn_i2 ) &
                                - ( DF_TR_Dy_i + DF_SC_Dy_i )/n_i )
    DB_Dy_i(2)    = y_i/u/n_i*( ( D2F_TR_Dy_i2 - n_i/y_i*D2F_TR_Dn_iDy_i   &
                                + D2F_SC_Dy_i2 - n_i/y_i*D2F_SC_Dn_iDy_i ) &
                                + ( DF_TR_Dy_i + DF_SC_Dy_i )/y_i )
    DB_DT(2)      = y_i/u/n_i*( ( D2F_TR_DTDy_i - n_i/y_i*D2F_TR_Dn_iDT   &
                                + D2F_SC_DTDy_i - n_i/y_i*D2F_SC_Dn_iDT ) )
!
    DB_Dup(3)     = - (D2F_SC_DuDn_i + (ONE-y_i)*D2F_SC_DuDy_i/n_i)/u &
                    - (D2F_TR_DuDn_i + (ONE-y_i)*D2F_TR_DuDy_i/n_i)/u &
                    - B(3)/u
    DB_Dn_i(3)    = - (D2F_TR_Dn_i2 + D2F_SC_Dn_i2)/u &
                    - (ONE-y_i)/u/n_i*(D2F_TR_Dn_iDy_i + D2F_SC_Dn_iDy_i) &
                    + (ONE-y_i)/u/n_i/n_i*(DF_TR_Dy_i + DF_SC_Dy_i)
    DB_Dy_i(3)    = - (ONE-y_i)/u/n_i*(D2F_TR_Dy_i2 + D2F_SC_Dy_i2) &
                    - ONE/u*(D2F_TR_Dn_iDy_i + D2F_SC_Dn_iDy_i) &
                    + ONE/u/n_i*(DF_TR_Dy_i + DF_SC_Dy_i)
    DB_DT(3)      = - (D2F_SC_Dn_iDT + (ONE-y_i)*D2F_SC_DTDy_i/n_i)/u &
                    - (D2F_TR_Dn_iDT + (ONE-y_i)*D2F_TR_DTDy_i/n_i)/u

    DB_Du       = DB_Dup   + DB_Dy_i*Dy_i_Du    + DB_Dn_i*Dn_i_Du
    DB_Dn_no    =            DB_Dy_i*Dy_i_Dn_no + DB_Dn_i*Dn_i_Dn_no
    DB_Dn_po    =            DB_Dy_i*Dy_i_Dn_po + DB_Dn_i*Dn_i_Dn_po


    IF (objective == 2 .OR. objective == 4) THEN
!   determine jacobian
      jac(1,1) = (DA_in_Dn_no(1) - DB_Dn_no(1) - DA_out_Dn_no(1))*n_no*Log(ten)
      jac(1,2) = (DA_in_Dn_po(1) - DB_Dn_po(1) - DA_out_Dn_po(1))*n_po*Log(ten)
      jac(1,3) = (DA_in_Du(1)    - DB_Du(1)    - DA_out_Du(1)   )*u   *Log(ten)

      jac(2,1) = (DA_in_Dn_no(2) - DB_Dn_no(2) - DA_out_Dn_no(2))*n_no*Log(ten)
      jac(2,2) = (DA_in_Dn_po(2) - DB_Dn_po(2) - DA_out_Dn_po(2))*n_po*Log(ten)
      jac(2,3) = (DA_in_Du(2)    - DB_Du(2)    - DA_out_Du(2)   )*u   *Log(ten)

      jac(3,1) = (DA_in_Dn_no(3) - DB_Dn_no(3) - DA_out_Dn_no(3))*n_no*Log(ten)
      jac(3,2) = (DA_in_Dn_po(3) - DB_Dn_po(3) - DA_out_Dn_po(3))*n_po*Log(ten)
      jac(3,3) = (DA_in_Du(3)    - DB_Du(3)    - DA_out_Du(3)   )*u   *Log(ten)

      non_uniform_jac = jac

      IF (objective == 2) RETURN
    ENDIF

    CALL SKYRME_BULK_TEMPEARTURE_DERIVATIVES( n_ni, n_pi, T, Meff_ni, Meff_pi, &
      eta_ni, eta_pi, tau_ni, tau_pi, v_ni, v_pi, mu_ni, mu_pi, &
      Dn_ni_DT, Dn_pi_DT, Dtau_ni_DT, Dtau_pi_DT, Dmu_ni_DT, Dmu_pi_DT, &
      DP_i_DT, 1)

    CALL Find_Non_Uniform_Matter_Derivatives ( n, Yp, T, &
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
      DB_Dup(1), DB_Dn_i(1), DB_Dy_i(1), DB_DT(1), &
      DB_Dup(2), DB_Dn_i(2), DB_Dy_i(2), DB_DT(2), &
      DB_Dup(3), DB_Dn_i(3), DB_Dy_i(3), DB_DT(3), &
      DU_DN, Dy_i_DN, Dn_i_DN, Deta_no_DN, Deta_po_DN, &
      DU_DT, Dy_i_DT, Dn_i_DT, Deta_no_DT, Deta_po_DT, &
      DU_DY, Dy_i_DY, Dn_i_DY, Deta_no_DY, Deta_po_DY )

    CALL OUTSIDE_PROPERTIES (u, T, n_no, n_po, mu_no, mu_po, tau_no, tau_po, &
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
        Dmu_no_DT, Dmu_no_Dn, Dmu_no_Dy,  &
        Dmu_po_DT, Dmu_po_Dn, Dmu_po_Dy )

    CALL HEAVY_PROPERTIES (u, T, n_ni, n_pi,           &
      mu_ni, mu_pi, tau_ni, tau_pi, Meff_ni, Meff_pi, P_i,  &
      Dmu_ni_Dn_i, Dmu_ni_Dy_i, Dmu_ni_DT, &
      Dmu_pi_Dn_i, Dmu_pi_Dy_i, Dmu_pi_DT, &
      Dtau_ni_Dn_i, Dtau_ni_Dy_i, Dtau_ni_DT, &
      Dtau_pi_Dn_i, Dtau_pi_Dy_i, Dtau_pi_DT, &
      DP_i_Dn_i, DP_i_Dy_i, DP_i_DT, &
      Du_DT, Du_Dn, Du_Dy,           &
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

      F     = F_h + F_out
      S     = S_h + S_out
      E     = F   + T*S

      DF_Dn = DF_h_Dn + DF_out_Dn
      DF_Dy = DF_h_Dy + DF_out_Dy
      DF_DT = DF_h_DT + DF_out_DT

      DS_Dn = DS_h_Dn + DS_out_Dn
      DS_Dy = DS_h_Dy + DS_out_Dy
      DS_DT = DS_h_DT + DS_out_DT

      DE_Dn = DF_Dn + T*DS_Dn
      DE_Dy = DF_Dy + T*DS_Dy
      DE_DT = DF_DT + T*DS_DT + S

      D2F_DT2  = - DS_DT
      D2F_DTDn = (ONE-Yp)*Dmu_no_DT + Yp*Dmu_po_DT
      D2F_DyDT = - n*(Dmu_no_DT-Dmu_po_DT)
      D2F_DnDy = - (mu_no-mu_po) - n*(Dmu_no_Dn-Dmu_po_Dn)
      D2F_Dn2  = (ONE-Yp)*Dmu_no_Dn + Yp*Dmu_po_Dn
      D2F_DY2  = - n*(Dmu_no_Dy - Dmu_po_Dy)

      P     = n*DF_Dn - F
      DP_Dn = n*D2F_Dn2
      DP_Dy = n*(mu_no - mu_po + D2F_DNDy)
      DP_DT = S + n*D2F_DTDn

      ! get x_nu derivatives w.r.t. n and T
      dLog10_n_no_dn = ONE/G_no*Deta_no_Dn/LOG(TEN)
      dLog10_n_po_dn = ONE/G_po*Deta_po_Dn/LOG(TEN)
      dLog10_u_dn    = DU_DN/u/LOG(TEN)
      dLog10_n_no_dT = (R_3_2/T + ONE/G_no*Deta_no_DT)/LOG(TEN)
      dLog10_n_po_dT = (R_3_2/T + ONE/G_po*Deta_po_DT)/LOG(TEN)
      dLog10_u_dT    = DU_DT/u/LOG(TEN)

  END SUBROUTINE FREE_ENERGY

END MODULE
