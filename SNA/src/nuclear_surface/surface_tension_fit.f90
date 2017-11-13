!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it AND/or modIFy
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
!    along with SRO_EOS.  IF not, see <http://www.gnu.org/licenses/>.
!
MODULE Fit_Surface_Tension_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, HALF

  IMPLICIT NONE

CONTAINS

  SUBROUTINE SURFACE_TENSION_FIT(m, n, x, fvec, iflag)

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN)  :: m, n
    REAL (DP), INTENT(IN)     :: x(:)
    REAL (DP), INTENT(IN OUT) :: fvec(:)
    INTEGER(I4B), INTENT(IN OUT) :: iflag
    INTEGER(I4B) :: i, j, k
    REAL(DP) :: DXI, DTEMP, TC, H, NUM, DEN, FUNC, P1, P2, P3
!      ----------
!      calculate the functions at x and return this vector in fvec.
!      ----------
    P1 = X(1)
    P2 = X(2)
    P3 = X(3)

    fvec = zero
    k = 0

    DO j = 1, 100
      DXI = dble(j)/2.d2
      DO i = 1, 200
        DTEMP = dble(i)/1.d1
        IF (lgttable(i,j)) THEN
          Tc = TCRIT*( TC0 + TC1*(ONE-TWO*DXI)**TWO &
             + TC2*(ONE-TWO*DXI)**FOUR + TC3*(ONE-TWO*DXI)**SIX)
          H = (ONE - (DTEMP/TC)**TWO)**P2
          IF (DTEMP>TC) H = ZERO
          NUM = TWO*TWO**P1+P3
          DEN = DXI**(-P1)+P3+(ONE-DXI)**(-P1)
          FUNC = H*NUM/DEN
          k = k + 1
          IF (FTABLE(I,J) == ZERO.or.ieee_is_nan(FTABLE(I,J))) THEN
            ! WRITE (*,*) I, J, H, NUM, DEN
            IF (FTABLE(I,J) == ZERO) THEN
              WRITE (*,*) I, J, ftable(i,j)
              WRITE (*,*) 'ERROR EVALUATING FVEC(K) AND/OR FTABLE(I,J)'
            ENDIF
            IF (ieee_is_nan(FTABLE(I,J)))  THEN
              WRITE (*,*)  'ERROR: FTABLE(I,J) IS NOT A NUMBER'
            ENDIF
            STOP
          ELSE
            fvec(k) = (FUNC-FTABLE(I,J))/FUNC
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE SURFACE_TENSION_FIT

END MODULE Fit_Surface_Tension_Mod
