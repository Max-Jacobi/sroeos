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
MODULE nwpure_mod

USE Kind_Types_Mod
USE nwutil_mod

IMPLICIT NONE

CONTAINS


      SUBROUTINE nwpure(n,xc,d,stepmx,scalex,fvec, &
                        xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)

      INTEGER(I4B) :: n,retcd,gcnt
      REAL(DP) :: stepmx,fpnorm
      REAL(DP) :: xc(*)
      REAL(DP) :: d(*),xp(*),fp(*),xw(*)
      REAL(DP) :: scalex(*)
      external fvec

      INTEGER(I4B) priter,iter

!-------------------------------------------------------------------------
!
!     Find a next iterate using geometric line search
!     along the newton direction
!
!     Arguments
!
!     In       n       Integer          dimension of problem
!     In       xc      Real(*)          current iterate
!     In       d       Real(*)          newton direction
!     In       stepmx  Real             maximum stepsize
!     In       scalex  Real(*)          diagonal scaling matrix for x()
!     In       fvec    Name             name of routine to calculate f()
!     In       xp      Real(*)          new x()
!     In       fp      Real(*)          new f(x)
!     In       fpnorm  Real             .5*||fp||**2
!     Out      xw      Real(*)          workspace for unscaling x
!
!     Out      retcd   Integer          RETURN code
!                                         0 new satisfactory x() found (!always)
!
!     Out      gcnt    Integer          number of steps taken
!     In       priter  Integer           >0 if intermediate steps to be printed
!                                        -1 if no printing
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i
      REAL(DP) :: oarg(3)
      REAL(DP) :: lambda
      REAL(DP) :: dnrm2
      REAL(DP) :: dlen

      INTEGER(I4B) :: idamax

      REAL(DP), PARAMETER :: Rone=1.0d0

!     safeguard initial step size

      dlen = dnrm2(n,d,1)
      IF( dlen .gt. stepmx ) THEN
          lambda = stepmx / dlen
      ELSE
          lambda = Rone
      ENDIF

      retcd  = 0
      gcnt   = 1

!     compute the next iterate xp

      DO i=1,n
         xp(i) = xc(i) + lambda*d(i)
      ENDDO

!     evaluate functions and the objective function at xp

      CALL nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)

!       IF( priter .gt. 0) THEN
!          oarg(1) = lambda
!          oarg(2) = fpnorm
!          oarg(3) = abs(fp(idamax(n,fp,1)))
!          CALL nwprot(iter,1,oarg)
!       ENDIF

      RETURN
      END SUBROUTINE

END MODULE nwpure_mod
