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

MODULE Polynomial_Fit_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP

  IMPLICIT NONE

CONTAINS

  FUNCTION polyfit(vx, vy, d)
!
!    Adapted from the Rosetta Code
!     http://rosettacode.org/wiki/Polynomial_regression#Fortran
!
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN)              :: d
    REAL(DP), DIMENSION(d+1)              :: polyfit
    REAL(DP), DIMENSION(:), INTENT(IN)    :: vx, vy

    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: X
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: XT
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: XTX

    INTEGER(I4B) :: i, j

    INTEGER(I4B) :: n, lda, lwork
    INTEGER(I4B) :: info
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ipiv
    REAL(DP), DIMENSION(:), ALLOCATABLE :: work

    n = d+1
    lda = n
    lwork = n

    ALLOCATE(ipiv(n))
    ALLOCATE(work(lwork))
    ALLOCATE(XT(n, size(vx)))
    ALLOCATE(X(size(vx), n))
    ALLOCATE(XTX(n, n))

    ! prepare the matrix
    DO i = 0, d
       DO j = 1, size(vx)
          X(j, i+1) = vx(j)**i
       ENDDO
    ENDDO

    XT  = transpose(X)
    XTX = matmul(XT, X)

    ! calls to LAPACK subs DGETRF and DGETRI
    CALL DGETRF(n, n, XTX, lda, ipiv, info)
    IF ( info /= 0 ) THEN
       stop "problem in polynomial fit of T_crit 0"
       RETURN
    ENDIF

    CALL DGETRI(n, XTX, lda, ipiv, work, lwork, info)
    IF ( info /= 0 ) THEN
       stop "problem in polynomial fit of T_crit 1"
       RETURN
    ENDIF

    polyfit = matmul( matmul(XTX, XT), vy)

    DEALLOCATE(ipiv)
    DEALLOCATE(work)
    DEALLOCATE(X)
    DEALLOCATE(XT)
    DEALLOCATE(XTX)

  END FUNCTION polyfit

END MODULE Polynomial_Fit_Mod
