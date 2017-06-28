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
MODULE Free_Energy_Heavy_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL
  USE Physical_Constants_Mod, ONLY : ZERO, HALF, ONE, TWO, THREE, FOUR, &
      R_1_3, R_2_3, R_3_2, R_3_5, R_5_2, R_5_3, R_2_9, PI, N_Q, TEN
  USE Surface_Properties_Mod
  USE Delu_Mod, ONLY : DELU_SUB
  USE H_Mod
  USE Sigma_Mod
  USE Beta_Mod
  USE Nuclear_Radius_Mod
  USE Nuclear_Mass_Mod
  USE Mu_Mod
  USE Free_Energy_Heavy_Surface_Mod
  USE Free_Energy_Heavy_Translation_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Free_Energy_Heavy(objective, n_ni, n_pi, u, T, A, Z, Radius, &
    F_SC, DF_SC_DU, DF_SC_DN, DF_SC_DT, DF_SC_DY, &
    D2F_SC_DU2, D2F_SC_DUDN, D2F_SC_DUDT, D2F_SC_DUDY, &
    D2F_SC_DN2, D2F_SC_DNDT, D2F_SC_DNDY, D2F_SC_DT2, D2F_SC_DTDY, D2F_SC_DY2, &
    F_TR, DF_TR_DU, DF_TR_DN, DF_TR_DT, DF_TR_DY, &
    D2F_TR_DU2, D2F_TR_DUDN, D2F_TR_DUDT, D2F_TR_DUDY, &
    D2F_TR_DN2, D2F_TR_DNDT, D2F_TR_DNDY, D2F_TR_DT2, D2F_TR_DTDY, D2F_TR_DY2)

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: objective
    INTEGER(I4B) :: derivatives, derivatives_nuc_size

    REAL(DP), INTENT(IN) :: n_ni, n_pi, u, T
    REAL(DP) :: n, y

    REAL(DP), INTENT(OUT) :: A, Z, Radius
    REAL(DP), INTENT(OUT) :: F_SC, DF_SC_DU, DF_SC_DN, DF_SC_DT, DF_SC_DY
    REAL(DP), INTENT(OUT) :: D2F_SC_DU2, D2F_SC_DUDN, D2F_SC_DUDT, D2F_SC_DUDY
    REAL(DP), INTENT(OUT) :: D2F_SC_DN2, D2F_SC_DNDT, D2F_SC_DNDY
    REAL(DP), INTENT(OUT) :: D2F_SC_DT2, D2F_SC_DTDY, D2F_SC_DY2
    REAL(DP), INTENT(OUT) :: F_TR, DF_TR_DU, DF_TR_DN, DF_TR_DT, DF_TR_DY
    REAL(DP), INTENT(OUT) :: D2F_TR_DU2, D2F_TR_DUDN, D2F_TR_DUDT, D2F_TR_DUDY
    REAL(DP), INTENT(OUT) :: D2F_TR_DN2, D2F_TR_DNDT, D2F_TR_DNDY
    REAL(DP), INTENT(OUT) :: D2F_TR_DT2, D2F_TR_DTDY, D2F_TR_DY2

    REAL(DP) :: H, DH_DT, DH_DY, D2H_DT2, D2H_DY2, D2H_DTDY

    REAL(DP) :: SIGMA, DSIGMA_DT, DSIGMA_DY
    REAL(DP) :: D2SIGMA_DT2, D2SIGMA_DTDY, D2SIGMA_DY2

    REAL(DP) :: BETA, DBETA_DY, DBETA_DN, DBETA_DT
    REAL(DP) :: D2BETA_DY2, D2BETA_DN2, D2BETA_DT2
    REAL(DP) :: D2BETA_DNDY,D2BETA_DNDT,D2BETA_DTDY

    REAL(DP) :: DELU, DDELU_DU, D2DELU_DU2

    REAL(DP) :: R, DR_DU, DR_DN, DR_DT, DR_DY
    REAL(DP) :: D2R_DU2, D2R_DUDN, D2R_DUDT, D2R_DUDY
    REAL(DP) :: D2R_DN2, D2R_DNDT, D2R_DNDY
    REAL(DP) :: D2R_DT2, D2R_DTDY, D2R_DY2

    REAL(DP) :: DA_DU, DA_DN, DA_DT, DA_DY
    REAL(DP) :: D2A_DU2, D2A_DUDN, D2A_DUDT, D2A_DUDY
    REAL(DP) :: D2A_DN2, D2A_DNDT, D2A_DNDY
    REAL(DP) :: D2A_DT2, D2A_DTDY, D2A_DY2

    REAL(DP) :: MU, DMU_DU, DMU_DN, DMU_DT, DMU_DY
    REAL(DP) :: D2MU_DU2, D2MU_DUDN, D2MU_DUDT, D2MU_DUDY
    REAL(DP) :: D2MU_DN2, D2MU_DNDT, D2MU_DNDY
    REAL(DP) :: D2MU_DT2, D2MU_DTDY, D2MU_DY2

    n = n_ni + n_pi
    y = n_pi / n

    derivatives = objective

    IF (objective == 1) derivatives = 1
    IF (objective == 2) derivatives = 2
    IF (objective == 3) derivatives = 0
    IF (objective == 4) derivatives = 2

    derivatives_nuc_size = derivatives
!   calculate h(y,T) and its derivatives up to "derivatives" order
    CALL H_SUB (derivatives,y,T,H,DH_DT,DH_Dy,D2H_DT2,D2H_DTDy,D2H_Dy2)

    IF (H==ZERO) THEN
      A = ZERO
      RETURN
    ENDIF

!   calculates \sigma(y,T) and its derivatives up to "derivatives" order
!   depends on h(y,T) and its derivatives
    CALL SIGMA_SUB (derivatives,y,T,H,DH_DT,DH_Dy,D2H_DT2,D2H_Dy2,D2H_DTDy,  &
              SIGMA,DSIGMA_DT,DSIGMA_Dy,D2SIGMA_DT2,D2SIGMA_DTDy,D2SIGMA_Dy2)

!   calculates \beta(y,n,T) and its derivatives up to "derivatives" order
!   depends on h(y,T) and \sigma(y,T) and their derivatives
    CALL BETA_SUB (derivatives,y,n,H,DH_DT,DH_Dy,D2H_DT2,D2H_Dy2,D2H_DTDy,   &
          SIGMA,DSIGMA_DT,DSIGMA_Dy,D2SIGMA_DT2,D2SIGMA_DTDy,D2SIGMA_Dy2,    &
          BETA,DBETA_Dn,DBETA_DT,DBETA_Dy,D2BETA_Dn2,D2BETA_DnDT,D2BETA_DnDy,&
          D2BETA_DT2,D2BETA_DTDy,D2BETA_Dy2)

!   calculates Del(u) and its derivatives up to "derivatives" order
    CALL DELU_SUB (derivatives,u,DELU,DDELU_DU,D2DELU_DU2)
    !   calculate heavy nuclei generalized radius and its derivatives
    CALL NUCLEAR_RADIUS_SUB (derivatives,u,DELU,DDELU_DU,D2DELU_DU2,         &
            SIGMA,DSIGMA_DT,DSIGMA_Dy,D2SIGMA_DT2,D2SIGMA_Dy2,D2SIGMA_DTDy,  &
            BETA,DBETA_Dy,DBETA_Dn,DBETA_DT,D2BETA_Dy2,D2BETA_Dn2,D2BETA_DT2,&
            D2BETA_DnDy,D2BETA_DnDT,D2BETA_DTDy,R,DR_Du,DR_Dn,DR_DT,DR_Dy,   &
            D2R_Du2,D2R_DuDn,D2R_DuDT,D2R_DuDy,D2R_Dn2,D2R_DnDT,D2R_DnDy,    &
            D2R_DT2,D2R_DTDy,D2R_Dy2)

    !   calculate heavy nuclei nuclear mass number and its derivatives
    CALL NUCLEAR_MASS_SUB (derivatives,n,R,DR_DU,DR_DN,DR_DT,DR_DY,            &
          D2R_DU2,D2R_DUDN,D2R_DUDT,D2R_DUDY,D2R_DN2,D2R_DNDT,D2R_DNDY,        &
          D2R_DT2,D2R_DTDY,D2R_DY2,A,DA_DU,DA_DN,DA_DT,DA_DY,D2A_DU2,D2A_DUDN, &
          D2A_DUDT,D2A_DUDY,D2A_DN2,D2A_DNDT,D2A_DNDY,D2A_DT2,D2A_DTDY,D2A_DY2)

!   calculate heavy nuclei chamical potential and its derivatives
    CALL MU_SUB (derivatives,u,n,T,A,DA_DU,DA_DN,DA_DT,DA_DY,                 &
        D2A_DU2,D2A_DUDN,D2A_DUDT,D2A_DUDY,D2A_DN2,D2A_DNDT,D2A_DNDY,D2A_DT2, &
        D2A_DTDY,D2A_DY2,MU,DMU_DU,DMU_DN,DMU_DT,DMU_DY,D2MU_DU2,D2MU_DUDN,   &
        D2MU_DUDT,D2MU_DUDY,D2MU_DN2,D2MU_DNDT,D2MU_DNDY,D2MU_DT2,D2MU_DTDY,  &
        D2MU_DY2)

!   heavi nuclei nuclear charge and radius
    Z = y*A
    Radius = R

!   obtain free energy for surface+coulomb terms and its derivatives
    CALL FREE_ENERGY_SURFACE (derivatives,u,DELU,DDELU_DU,D2DELU_DU2,      &
        BETA,DBETA_Dn,DBETA_DT,DBETA_Dy,D2BETA_Dn2,D2BETA_DnDT,D2BETA_DnDy,&
        D2BETA_DT2,D2BETA_DTDy,D2BETA_Dy2,F_SC,DF_SC_DU,DF_SC_DN,DF_SC_DT, &
        DF_SC_DY,D2F_SC_DU2,D2F_SC_DUDN,D2F_SC_DUDT,D2F_SC_DUDY,D2F_SC_DN2,&
        D2F_SC_DNDT,D2F_SC_DNDY,D2F_SC_DT2,D2F_SC_DTDY,D2F_SC_DY2)

!   obtain free energy for translational terms and its derivatives
    CALL FREE_ENERGY_TRANSLATION (derivatives,n,T,u,DELU,DDELU_DU,D2DELU_DU2, &
        H,DH_DT,DH_Dy,D2H_DT2,D2H_Dy2,D2H_DTDy,A,DA_DU,DA_DN,DA_DT,DA_DY,     &
        D2A_DU2,D2A_DUDN,D2A_DUDT,D2A_DUDY,D2A_DN2,D2A_DNDT,D2A_DNDY,D2A_DT2, &
        D2A_DTDY,D2A_DY2,MU,DMU_DU,DMU_DN,DMU_DT,DMU_DY,D2MU_DU2,D2MU_DUDN,   &
        D2MU_DUDT,D2MU_DUDY,D2MU_DN2,D2MU_DNDT,D2MU_DNDY,D2MU_DT2,D2MU_DTDY,  &
        D2MU_DY2,F_TR,DF_TR_DU,DF_TR_DN,DF_TR_DT,DF_TR_DY,D2F_TR_DU2,        &
        D2F_TR_DUDN,D2F_TR_DUDT,D2F_TR_DUDY,D2F_TR_DN2,D2F_TR_DNDT,           &
        D2F_TR_DNDY,D2F_TR_DT2,D2F_TR_DTDY,D2F_TR_DY2)

  END SUBROUTINE Free_Energy_Heavy

END MODULE Free_Energy_Heavy_Mod
