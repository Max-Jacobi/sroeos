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
MODULE Delu_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, HALF, ONE, TWO, THREE, &
      R_1_3, R_2_3, R_3_2, R_3_5, R_5_3, R_2_9, TEN

  IMPLICIT NONE

CONTAINS

  SUBROUTINE DELU_SUB(derivatives,u,DELU,DDELU_DU,D2DELU_DU2)

    INTEGER(I4B), INTENT(IN) :: derivatives
    REAL(DP), INTENT(IN) :: u
    REAL(DP), INTENT(OUT) :: DELU, DDELU_DU, D2DELU_DU2
    REAL(DP) :: Du, odu, Dodu, NUM, DEN
    REAL(DP) :: DDu_Du, DDodu_Du, DNUM_Du, DDEN_Du
    REAL(DP) :: D2Du_Du2, D2Dodu_Du2, D2NUM_Du2, D2DEN_Du2, AUXD1, AUXD2

    odu = one-u

    Du   = F_DU(u)
    Dodu = F_DU(odu)
    ! IF u<1.d-15 then odu -> 1 to within machine precision
    ! Thus, use Taylor series to evaluate of F_Du(odu)
    ! TODO: move this to function below
    IF (u<1.D-12) Dodu = u*u/6.0D0+5.0D0*u*u*u/54.d0

    NUM = odu*Du**R_1_3 + u*Dodu**R_1_3
    DEN = u*u + odu*odu + R_3_5*u*u*odu*odu

    DELU = u*odu*NUM/DEN

    DDELU_DU   = ZERO
    D2DELU_DU2 = ZERO

    IF ( derivatives == 0 ) RETURN

    DDu_Du   = F_DDu_Du(u)
    DDodu_Du = - F_DDu_Du(odu)
    ! IF u<1.d-15 then odu -> 1 to within machine precision
    ! Thus, use Taylor series to evaluate of F_DDu_Du(odu)
    ! TODO: move this to function below
    IF (u<1.D-12) DDodu_Du = -u/THREE

    DNUM_Du = - Du**R_1_3 + Dodu**R_1_3 &
            + R_1_3*(odu*DDu_Du/Du**R_2_3 + u*DDodu_Du/Dodu**R_2_3)
    DDEN_Du = TWO*(u - odu + R_3_5*u*odu*(odu-u))

    DDELU_DU = DELU*(ONE/u-ONE/odu-DDEN_Du/DEN+DNUM_Du/NUM)

    IF ( derivatives == 1 ) RETURN

    D2Du_Du2   = F_D2Du_Du2(u)
    D2Dodu_Du2 = F_D2Du_Du2(odu)

    D2DEN_Du2 = 4.0D0 + 1.2D0*odu*odu - 4.8D0*odu*u + 1.2D0*u*u
    D2NUM_Du2 = R_2_3*(DDodu_Du/Dodu**R_2_3 - DDu_Du/Du**R_2_3) &
            - R_2_9*(odu*DDu_Du**TWO/Du**R_5_3 + u*DDodu_Du**TWO/Dodu**R_5_3) &
            + R_1_3*(odu*D2Du_Du2/Du**R_2_3 + u*D2Dodu_Du2/Dodu**R_2_3)

    AUXD1 = -ONE/u**TWO-ONE/odu**TWO+(DDEN_Du/DEN)**TWO-D2DEN_Du2/DEN
    AUXD2 = D2NUM_Du2/NUM-(DNUM_Du/NUM)**TWO

    D2DELU_DU2 = DDELU_DU**TWO/DELU+DELU*(AUXD1+AUXD2)

    IF ( derivatives == 2 ) RETURN

  END SUBROUTINE DELU_SUB

! function to calculate D(u)
  FUNCTION F_DU(u) RESULT(D)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: u
    REAL(DP) :: odu, D

    IF (u>0.99999D0) THEN
      odu = ONE-u
      D = odu*odu/6.0D0+5.0D0*odu*odu*odu/54.d0
    ELSE
      D = ONE - R_3_2*u**R_1_3 + HALF*u
    ENDIF

    RETURN

  END FUNCTION F_DU

! function to calculate dD(u)/du
  FUNCTION F_DDu_Du(u) RESULT(DDuDu)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: u
    REAL(DP) :: odu, DDuDu

    DDuDu = HALF - HALF/u**R_2_3

    RETURN

  END FUNCTION F_DDu_Du

! function to calculate d²D(u)/du²
  FUNCTION F_D2Du_Du2(u) RESULT(D2DuDu2)

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: u
    REAL(DP) :: odu, D2DuDu2

    D2DuDu2 = R_1_3/u**R_5_3

    RETURN

  END FUNCTION F_D2Du_Du2

END MODULE Delu_Mod
