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

!
!  This file is originally part of the NLEQSLV of Berend Hasselman available at
!  https://cran.r-project.org/web/packages/nleqslv/index.html
!
!  It has been slightly modified by Andre da Silva Schneider
!   to fit SRO_EOS source code
!
MODULE limhpar_mod

USE Kind_Types_Mod
USE liqrev_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE limhpar(R, ldr, n, sdiag, qtf, dn, dnlen, &
                         glen, delta, mu, d, work)
      INTEGER(I4B) :: ldr, n
      REAL(DP) :: R(ldr,*), sdiag(*)
      REAL(DP) :: qtf(*),dn(*), dnlen, glen, d(*),work(*)
      REAL(DP) :: delta, mu

!-----------------------------------------------------------------------------
!
!     Arguments
!
!       Inout  R      Real(ldr,*)     upper triangular matrix R from QR (unaltered)
!                                     strict lower triangle CONTAINS
!                                        transposed strict upper triangle of the upper
!                                        triangular matrix S.
!
!       In     n      Integer         order of R.
!
!       In     ldr    Integer         leading dimension of the array R.
!
!       Out    sdiag  Real(*)         vector of size n, containing the
!                                     diagonal elements of the upper
!                                     triangular matrix S.
!
!       In     qtr    Real(*)         trans(Q)*f vector of size n
!       In     dn     Real(*)         Newton step
!       In     dnlen  Real            length Newton step
!       In     glen   Real            length gradient vector
!
!       Inout  mu     Real            Levenberg-Marquardt PARAMETER
!       In     delta  Real            size of trust region (euclidian norm)
!
!       Out    d      Real(*)         vector with step with norm very close to delta
!       Out    work   Real(*)         workspace of size n.
!
!     Description
!
!     determine Levenberg-Marquardt PARAMETER mu such that
!     norm[(R**T R + mu*I)**(-1) * qtf] - delta approximately 0
!     See description in liqrev.f for further details
!
!     Algorithm comes from More: The Levenberg-Marquardt algorithm, Implementation and Theory
!     Lecture Notes in Mathematics, 1978, no. 630.
!     USEs liqrev (in file liqrev.f) which is based on Nocedal's method (see comments in file)
!-----------------------------------------------------------------------------

      REAL(DP) :: phi, pnorm, qnorm, mulo, muhi,dmu, sqmu
      INTEGER(I4B) :: iter
      LOGICAL done
      REAL(DP) :: dnrm2
      REAL(DP), PARAMETER :: Rone = 1.0d0

      phi = dnlen - delta
      muhi = glen/delta

      CALL dcopy(n,dn,1,d,1)
      CALL dscal(n, Rone/dnlen, d, 1)

!     solve R**T * x = dn
      CALL dtrsv("U","T","N",n,R,ldr,d,1)
      qnorm = dnrm2(n,d,1)
      mulo = (phi/dnlen)/qnorm**2
      mu = mulo

      iter = 0
      done = .false.
      DO WHILE( .not. done )
          iter = iter + 1
          sqmu = sqrt(mu)
          CALL liqrev(n, R, ldr, sqmu, qtf, d, sdiag, work)
          pnorm = dnrm2(n,d,1)
          CALL dcopy(n,d,1,work,1)
          CALL dtrstt(R, ldr, n, sdiag, work)
          done = abs(pnorm-delta) .le. .1d0*delta .or. iter .gt. 5
          IF( .not. done ) THEN
              qnorm = dnrm2(n,work,1)
              IF( pnorm .gt. delta ) THEN
                  mulo = max(mulo,mu)
              ELSE IF( pnorm .lt. delta ) THEN
                  muhi = min(muhi,mu)
              ENDIF
              dmu = (pnorm-delta)/delta * (pnorm/qnorm)**2
              mu = max(mulo, mu + dmu)
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE

END MODULE
