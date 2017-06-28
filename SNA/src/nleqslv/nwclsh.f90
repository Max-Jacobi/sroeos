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
MODULE nwclsh_mod

USE Kind_Types_Mod
USE nwutil_mod

IMPLICIT NONE

CONTAINS


      SUBROUTINE nwclsh(n,xc,fcnorm,d,g,stepmx,xtol,scalex,fvec, &
                        xp,fp,fpnorm,xw,retcd,gcnt,priter,iter)

      INTEGER(I4B) :: n,retcd,gcnt
      REAL(DP) :: stepmx,xtol,fcnorm,fpnorm
      REAL(DP) :: xc(*)
      REAL(DP) :: d(*),g(*),xp(*),fp(*),xw(*)
      REAL(DP) :: scalex(*)
      INTERFACE
      SUBROUTINE fvec ( x, r, n, flag )
        USE Kind_Types_Mod
        INTEGER(I4B), intent(in) :: n, flag
        REAL(DP), dimension(n), intent(in) :: x
        REAL(DP), dimension(n), intent(out) :: r
      END SUBROUTINE fvec
      END INTERFACE

      INTEGER(I4B) :: priter,iter

!-------------------------------------------------------------------------
!
!     Find a next acceptable iterate using a safeguarded cubic line search
!     along the newton direction
!
!     Arguments
!
!     In       n       Integer          dimension of problem
!     In       xc      Real(*)          current iterate
!     In       fcnorm  Real             0.5 * || f(xc) ||**2
!     In       d       Real(*)          newton direction
!     In       g       Real(*)          gradient at current iterate
!     In       stepmx  Real             maximum stepsize
!     In       xtol    Real             relative step size at which
!                                       successive iterates are considered
!                                       close enough to terminate algorithm
!     In       scalex  Real(*)          diagonal scaling matrix for x()
!     In       fvec    Name             name of routine to calculate f()
!     In       xp      Real(*)          new x()
!     In       fp      Real(*)          new f(x)
!     In       fpnorm  Real             .5*||fp||**2
!     Out      xw      Real(*)           workspace for unscaling x(*)
!
!     Out      retcd   Integer          RETURN code
!                                         0 new satisfactory x() found
!                                         1 no  satisfactory x() found
!                                           sufficiently distinct from xc()
!
!     Out      gcnt    Integer          number of steps taken
!     In       priter  Integer           >0 unit if intermediate steps to be printed
!                                        -1 if no printing
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i
      REAL(DP) :: slope,rsclen,oarg(4)
      REAL(DP) :: lambda,lamhi,lamlo,t
      REAL(DP) :: ddot,dnrm2, ftarg
      REAL(DP) :: dlen
      REAL(DP) :: a, b, disc, fpt, fpt0, fpnorm0, lambda0
      INTEGER(I4B) idamax
      LOGICAL firstback

      REAL(DP), PARAMETER :: alpha = 1.0d-4

      REAL(DP), PARAMETER :: Rzero=0.0d0, Rhalf=0.5d0, Rone=1.0d0, Rtwo=2.0d0, Rten=10.0d0, Rthree=3.0d0
      REAL(DP) :: t1,t2
      REAL(DP) :: dlamch

!     silence warnings issued by ftncheck

      lambda0 = Rzero
      fpnorm0 = Rzero

!     safeguard initial step size

      dlen = dnrm2(n,d,1)
      IF( dlen .gt. stepmx ) THEN
          lamhi  = stepmx / dlen
      ELSE
          lamhi  = Rone
      ENDIF

!     compute slope  =  g-trans * d

      slope = ddot(n,g,1,d,1)

!     compute the smallest value allowable for the damping
!     PARAMETER lambda ==> lamlo

      rsclen = nudnrm(n,d,xc)
      lamlo  = xtol / rsclen

!     initialization of retcd and lambda (linesearch length)

      retcd  = 2
      lambda = lamhi
      gcnt   = 0
      firstback = .true.

      DO WHILE( retcd .eq. 2 )

!        compute next x

         DO i=1,n
            xp(i) = xc(i) + lambda*d(i)
         ENDDO

!        evaluate functions and the objective function at xp

         CALL nwfvec(xp,n,scalex,fvec,fp,fpnorm,xw)
         gcnt = gcnt + 1
         ftarg = fcnorm + alpha * lambda * slope

!          IF( priter .gt. 0) THEN
!             oarg(1) = lambda
!             oarg(2) = ftarg
!             oarg(3) = fpnorm
!             oarg(4) = abs(fp(idamax(n,fp,1)))
!             CALL nwlsot(iter,1,oarg)
!          ENDIF

!        first is quadratic
!        test whether the standard step produces enough decrease
!        of the objective function.
!        If not update lambda and compute a new next iterate

         IF( fpnorm .le. ftarg ) THEN
             retcd = 0
         ELSE
            IF( fpnorm .gt. lamlo**2 * sqrt(dlamch('O')) ) THEN
!               safety against overflow in what follows (USE lamlo**2 for safety)
                lambda = lambda/Rten
                firstback = .true.
            ELSE
                IF( firstback ) THEN
                   t = ((-lambda**2)*slope/Rtwo)/ &
                              (fpnorm-fcnorm-lambda*slope)
                   firstback = .false.
                ELSE
                   fpt  = fpnorm -fcnorm - lambda*slope
                   fpt0 = fpnorm0-fcnorm - lambda0*slope
                   a = fpt/lambda**2 - fpt0/lambda0**2
                   b = -lambda0*fpt/lambda**2 + lambda*fpt0/lambda0**2
                   a = a /(lambda - lambda0)
                   b = b /(lambda - lambda0)
                   IF( abs(a) .le. dlamch('E') ) THEN
!                      not quadratic but linear
                       t = -slope/(2*b)
                   ELSE
!                      USE Higham procedure to compute roots acccurately
!                      Higham: Accuracy and Stability of Numerical Algorithms, second edition,2002, page 10.
!                      Actually solving 3*a*x^2+2*b*x+c=0 ==> (3/2)*a*x^2+b*x+(c/2)=0
                       disc = b**2 - Rthree * a * slope
                       t1 = -(b+sign(Rone,b)*sqrt(disc))/(Rthree*a)
                       t2 = slope/(Rthree*a)/t1
                       IF(a .gt. Rzero ) THEN
!                          upward opening parabola ==> rightmost is solution
                           t = max(t1,t2)
                       ELSE
!                          downward opening parabola ==> leftmost is solution
                           t = min(t1,t2)
                       ENDIF
                   ENDIF
                   t = min(t, Rhalf*lambda)
                ENDIF
                lambda0 = lambda
                fpnorm0 = fpnorm
                lambda = max(t,lambda/Rten)
                IF(lambda .lt. lamlo) THEN
                   retcd = 1
                ENDIF
            ENDIF
         ENDIF
      ENDDO

      RETURN
      END SUBROUTINE

END MODULE nwclsh_mod
