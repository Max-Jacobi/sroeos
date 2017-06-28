!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    SRO_EOS is distributed in the hope that it will be USEful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with SRO_EOS.  If not, see <http://www.gnu.org/licenses/>.
!
MODULE FUNCTIONS_MOD

  USE Kind_types_Mod, ONLY : I4B, DP

  IMPLICIT NONE

CONTAINS

  FUNCTION calc_Ec(nb, T, Ye, A, Z, n0) RESULT(Ec)

    USE Physical_Constants_Mod, ONLY: alpha_coul, HBC, PI, R_1_3, R_3_5

    IMPLICIT NONE

    REAL(DP), INTENT(in)  :: nb,T,Ye,A,Z,n0
    REAL(DP) :: Ec, rr, uu

    IF (Z>2) THEN
      uu = Ye*nb/n0*A/Z
      rr = (3.d0*A/(4.d0*PI*n0))**R_1_3
      Ec = R_3_5*Z**2*alpha_coul/rr*(0.5d0*uu - 1.5d0*uu**R_1_3)*HBC
    ELSE
      Ec = 0.d0
    ENDIF
    RETURN

  END FUNCTION calc_Ec

  FUNCTION calc_dEcdln(nb, T, Ye, A, Z, n0, d2Ecdln2) RESULT(dEcdln)

    USE Physical_Constants_Mod, ONLY: alpha_coul, HBC, PI, R_3_5, R_1_3

    IMPLICIT NONE

    REAL(DP), INTENT(in)  :: nb,T,Ye,A,Z,n0
    REAL(DP), INTENT(out) :: d2Ecdln2
    REAL(DP) :: dEcdln
    REAL(DP) :: uu,rr

    IF (Z>2) THEN
      uu = Ye*nb/n0*A/Z
      rr = (3.d0*A/(4.d0*PI*n0))**(1.d0/3.d0)
      dEcdln   = R_3_5*Z**2*alpha_coul/rr*(0.5d0*uu-uu**R_1_3/2.d0)*HBC
      d2Ecdln2 = R_3_5*Z**2*alpha_coul/rr*(0.5d0*uu-uu**R_1_3/6.d0)*HBC
    ELSE
      dEcdln   = 0.d0
      d2Ecdln2 = 0.d0
    ENDIF
    RETURN

  END FUNCTION calc_dEcdln

  FUNCTION calc_mu(n, m, g, T, Ec) RESULT(mu)

    USE Physical_Constants_Mod, ONLY: sac_const, PI
    USE Fermi_Integrals_Mod, ONLY: inverse_fermi_one_half

    IMPLICIT NONE

    REAL(DP), INTENT(in) :: n ! number density in fm^{-3}
    REAL(DP), INTENT(in) :: m ! mass in MeV
    REAL(DP), INTENT(in) :: g ! partition FUNCTION
    REAL(DP), INTENT(in) :: T ! temperature in MeV
    REAL(DP), INTENT(in) :: Ec ! temperature in MeV
    REAL(DP) :: mu ! output chemical potential in MeV
    REAL(DP) :: delta, dmu
    REAL(DP) :: ifermi12

    mu = m + Ec + T * (log(n)-log(g)+1.5d0*(log(sac_const)-log(m)-log(T)))

    IF (m < 1000.0) THEN
      delta = SQRT(PI)*n/(2.d0*g*(m*T/sac_const)**1.5d0)
      mu = inverse_fermi_one_half(delta) * T  + m + Ec
    ENDIF

  END FUNCTION calc_mu

  FUNCTION calc_n(mu, m, g, T, Ec) RESULT(n)

    USE Physical_Constants_Mod, ONLY: sac_const, PI
    USE Fermi_Integrals_Mod, ONLY: fermi_one_half

    IMPLICIT NONE

    REAL(DP), INTENT(in) :: mu ! chemical potential in MeV
    REAL(DP), INTENT(in) :: m ! mass in MeV
    REAL(DP), INTENT(in) :: g ! partition FUNCTION
    REAL(DP), INTENT(in) :: T ! temperature in MeV
    REAL(DP), INTENT(in) :: Ec ! temperature in MeV
    REAL(DP) :: n ! output chemical potential in MeV
    REAL(DP) :: expt

    expt = (mu-m-Ec)/T

    IF (expt<-2.0D2) THEN
      n = 0.D0
    ELSEIF (expt>3.0D2) THEN
      n = g*(m*T/sac_const)**1.5D0*exp(3.0D2)
    ELSE
      n = g*(m*T/sac_const)**1.5D0*exp(expt)
    ENDIF

    IF (m < 1000.0) THEN
      n = g*(m*T/sac_const)**1.5d0 * 2.d0/SQRT(PI)*fermi_one_half((mu-m-Ec)/T)
    ENDIF

  END FUNCTION calc_n

END MODULE FUNCTIONS_MOD
