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
MODULE LU_decomposition

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE

  IMPLICIT NONE

CONTAINS
!
!   VERSION OF LU DECOMPOSITION SUBROUTINES BY
!     DOUGLAS SWESTY 94
!
!   ADAPTED TO F95 BY A.S. SCHNEIDER IN 03/23/2017
!
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!***********************************************************************
!
!    MODULE:       MATLUD
!    TYPE:         SUBROUTINE
!    AUTHOR:       F. DOUGLAS SWESTY
!    DATE:         6/16/94
!
!    PURPOSE:      LU decomposes an NxN matrix
!
!    NOTE:         This subroutine employs Crout's algorithm with
!                  implicit pivoting as described in "Numerical
!                  Recipes", by Press, Flannery, Teukolsky, and
!                  Vetterling, (pub. Cambridge Univ. Press) first
!                  edition, whose authors desrve the credit for
!                  this algorithm.
!
!    CALL LINE:    CALL MATLUD(A,LU,N,IPVT)
!
!    INPUTS:       A = NxN array to be LU decomposed  (D)
!                  N = Size of A (I)
!
!    OUTPUTS:      LU = Array containing LU decomposition of A (D)
!                  IPVT = Vector of pivot indices (I)
!
!    CALLS :       None
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!

  SUBROUTINE MATLUD(A,LU,N,IPVT)

    USE, INTRINSIC :: IEEE_ARITHMETIC

    IMPLICIT NONE

!   Parameters
    INTEGER (I4B) ::  N
    INTEGER (I4B) ::  IPVT(N)
    REAL (DP) ::  A(N,N), LU(N,N), MIJ
!
!   Local variables
!
!   Loop variables
    INTEGER (I4B) ::  I, J, L
!   A small floating point value to prevent division by zero
    REAL (DP) ::  SMALLV
    PARAMETER(SMALLV=1.0D-20)
!   Value & original row index or largest L value
    REAL (DP) ::  E_MAX, TSTVAL
    INTEGER (I4B) ::  ROWPVT
!   A scratch variable to use when swapping rows
    REAL (DP) ::  SCRVAR
!   An array for the largest value of each row of A
    INTEGER (I4B), PARAMETER ::  NMAX = 100
    REAL (DP) ::  MAXVAL(NMAX)


!   Check if all elements on the matrix have real & finite values
    DO I = 1, N
      DO J = 1, N
        MIJ = A(I,J)
        IF (MIJ > 1.D100) MIJ =  1.D100
        IF (MIJ <-1.D100) MIJ = -1.D100
        IF (ieee_is_nan(MIJ)) MIJ = 0.D0
        A(I,J) = MIJ
      ENDDO
    ENDDO

!   Find the maximum absolute value in each row of A
    DO I=1,N,1
      MAXVAL(I) = - ONE
      DO J=1,N,1
        MAXVAL(I)  = MAX(SMALLV,ABS(A(I,J)),MAXVAL(I))
        LU(I,J) = A(I,J)
      ENDDO
    ENDDO
!   Now employ Crout's algorithm with implicit pivoting
    DO J=1,N,1
!     Calculate column J, U matrix elements
      DO I=1,J-1,1
        DO L=1,I-1,1
          LU(I,J) = LU(I,J)-LU(I,L)*LU(L,J)
        ENDDO
      ENDDO
!     Calculate column J, L matrix elements.  The element is
!     scaled by the largest element of the original matrix.
!     Also, the column from the diagonal down is searched
!     for the largest element.
      E_MAX = - ONE
      DO I=J,N,1
        DO L=1,J-1,1
          LU(I,J) = LU(I,J)-LU(I,L)*LU(L,J)
        ENDDO
        TSTVAL = ABS( LU(I,J) )/MAXVAL(I)
        IF(TSTVAL.GT.E_MAX) THEN
          E_MAX=TSTVAL
          ROWPVT = I
        ENDIF
      ENDDO
!     Keep track of which row was pivoted into row J
      IPVT(J) = ROWPVT
!     If the original diagonal element wasn't the
!     largest, then swap row J with row ROWPVT
      IF(ROWPVT.NE.J) THEN
        DO L=1,N,1
          SCRVAR = LU(ROWPVT,L)
          LU(ROWPVT,L) = LU(J,L)
          LU(J,L) = SCRVAR
        ENDDO
!       Set the MAXVAL for row ROWPVT to that of the one
!       that was just swapped in.
        MAXVAL(ROWPVT) = MAXVAL(J)
      ENDIF
!     Now divide the rest of the L column by the pivot element
      IF(ABS(LU(J,J)).GT.SMALLV )  THEN
        SCRVAR = LU(J,J)
      ELSE
        SCRVAR = SIGN(SMALLV,LU(J,J))
      ENDIF
      DO I=J+1,N,1
        LU(I,J) = LU(I,J)/SCRVAR
      ENDDO
    ENDDO

    RETURN

  END SUBROUTINE MATLUD
!***********************************************************************
!
!    MODULE:       MLUSLV
!    TYPE:         SUBROUTINE
!    AUTHOR:       F. DOUGLAS SWESTY
!    DATE:         6/16/94
!
!    PURPOSE:      LU decomposes an NxN matrix
!
!    NOTE:         This subroutine employs a forward & back substitution
!                  algorithm with implicit pivoting as described in
!                  "Numerical Recipes", by Press, Flannery, Teukolsky, &
!                  Vetterling, (pub. Cambridge Univ. Press) first
!                  edition, whose authors desrve the credit for
!                  this algorithm.
!
!    CALL LINE:    CALL MLUSLV(LU,X,B,IPVT,N)
!
!    INPUTS:       A = NxN LU decomposed array (D)
!                  B  = RHS vector (D)
!                  IPVT = Vector of pivots (I)
!                  N = Size of A (I)
!
!    OUTPUTS:      X = Solution vector(D)
!
!    CALLS :       None
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
  SUBROUTINE MLUSLV(LU,X,B,IPVT,N)

    IMPLICIT NONE
!   Parameters
    INTEGER (I4B) ::  N
    INTEGER (I4B) ::  IPVT(N)
    REAL (DP) ::  LU(N,N), X(N), B(N)
!
!   Local variables
!
!   Loop variables
    INTEGER (I4B) ::  I, J, L
!   A scratch variable
    REAL (DP) ::  SCRV
!   Copy the RHS into X so we don't destroy B by
!   unscrambling the pivots
    DO I=1,N,1
      X(I) = B(I)
    ENDDO
!   Do the forward substitution
    DO I=1,N,1
      L = IPVT(I)
      SCRV = X(L)
      X(L) = X(I)
      X(I) = SCRV
      DO J=1,I-1,1
        X(I) = X(I)-LU(I,J)*X(J)
      ENDDO
    ENDDO
!   Do the backward substitution
    DO I=N,1,-1
      DO J=I+1,N,1
        X(I) = X(I)-LU(I,J)*X(J)
      ENDDO
      X(I) = X(I)/LU(I,I)
    ENDDO

    RETURN

  END SUBROUTINE MLUSLV

END MODULE LU_decomposition
