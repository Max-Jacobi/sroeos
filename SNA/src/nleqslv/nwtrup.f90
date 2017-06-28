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
MODULE nwtrup_mod

USE Kind_Types_Mod
USE nwutil_mod

IMPLICIT NONE

CONTAINS


      SUBROUTINE nwtrup(n,fcnorm,g,sc,nwtstep,stepmx,xtol,delta, &
                        fpred,retcd,xprev,fpnsav,fprev,xp,fp, &
                        fpnorm)

      INTEGER(I4B) :: n,retcd
      REAL(DP) :: fcnorm,stepmx,xtol,delta,fpred,fpnsav,fpnorm
      REAL(DP) :: xp(*),g(*)
      REAL(DP) :: sc(*),xprev(*),fprev(*),fp(*)
      LOGICAL nwtstep

!-------------------------------------------------------------------------
!
!     Decide whether to accept xp=xc+sc as the next iterate
!     and updates the trust region delta
!
!     Arguments
!
!     In       n       Integer         size of xc()
!     In       fcnorm  Real            .5*||f(xc)||**2
!     In       g       Real(*)         gradient at xc()
!     In       sc      Real(*)         current step
!     In       nwtstep Logical         true if sc is newton direction
!     In       stepmx  Real            maximum step size
!     In       xtol    Real            minimum step tolerance
!     Inout    delta   Real            trust region radius
!     In       fpred   Real            predicted value of .5*||f()||**2
!
!     Inout    retcd   Integer         RETURN code
!                                       0 xp accepted as next iterate;
!                                         delta trust region for next iteration.
!
!                                       1 xp unsatisfactory but
!                                         accepted as next iterate becaUSE
!                                         xp-xc .lt. smallest allowable
!                                         step length.
!
!                                       2 f(xp) too large.
!                                         continue current iteration with
!                                         new reduced delta.
!
!                                       3 f(xp) sufficiently small, but
!                                         quadratic model predicts f(xp)
!                                         sufficiently well to continue current
!                                         iteration with new doubled delta.
!
!                                      On first entry, retcd must be 4
!
!     Wk       xprev   Real(*)         (internal) work
!     Wk       fpnsav  Real            (internal) work
!     Wk       fprev   Real(*)         (internal) work
!     Inout    xp      Real(*)         new iterate x()
!     Inout    fp      Real(*)         new f(xp)
!     Inout    fpnorm  Real            new .5*||f(xp)||**2
!
!-------------------------------------------------------------------------

      REAL(DP) :: ared,pred,slope,sclen,rln,dltmp
      REAL(DP) :: dnrm2,ddot
      LOGICAL ret3ok

      REAL(DP), PARAMETER :: Rpten = 0.1d0, Rhalf=0.5d0, Rone=1.0d0, Rtwo=2.0d0, &
                             Rthree=3.0d0, Rfour=4.0d0, Rten=10.0d0

      REAL(DP), PARAMETER :: Rp99=Rone-Rten**(-2), Rp4=Rten**(-4), Rp75=Rthree/Rfour

      REAL(DP), PARAMETER :: alpha = Rp4

!     pred measures the predicted reduction in the function value

      ared  = fpnorm - fcnorm
      pred  = fpred  - fcnorm
      slope = ddot(n,g,1,sc,1)

      IF(retcd .ne. 3) THEN
         ret3ok = .false.
      ELSE
         ret3ok = fpnorm .ge. fpnsav .or. ared .gt. alpha * slope
      ENDIF

      IF(retcd .eq. 3 .and. ret3ok) THEN

!        reset xp,fp,fpnorm to saved values and terminate global step

         retcd = 0
         CALL dcopy(n,xprev,1,xp,1)
         CALL dcopy(n,fprev,1,fp,1)
         fpnorm = fpnsav
!        reset delta to initial value
!        but beware
!           if the trial step was a Newton step THEN delta is reset to
!           .5 * length(Newton step) which will be smaller
!           becaUSE delta is set to length(Newton step) ELSEwhere
!           see nwddlg.f and nwpdlg.f
         delta  = Rhalf*delta

      ELSEIF(ared .gt. alpha * slope) THEN

!        fpnorm too large (decrease not sufficient)

         rln = nudnrm(n,sc,xp)
         IF(rln .lt. xtol) THEN

!           cannot find satisfactory xp sufficiently distinct from xc

            retcd = 1

         ELSE

!           reduce trust region and continue current global step

            retcd = 2
            sclen = dnrm2(n,sc,1)
            dltmp = -slope*sclen/(Rtwo*(ared-slope))

            IF(dltmp .lt. Rpten*delta) THEN
               delta = Rpten*delta
            ELSE
               delta = min(Rhalf*delta, dltmp)
            ENDIF

         ENDIF

      ELSEIF(retcd .ne. 2 .and. (abs(pred-ared) .le. Rpten*abs(ared)) &
            .and. (.not. nwtstep) .and. (delta .le. Rp99*stepmx)) THEN

!        pred predicts ared very well
!        attempt a doubling of the trust region and continue global step
!        when not taking a newton step and trust region not at maximum

         CALL dcopy(n,xp,1,xprev,1)
         CALL dcopy(n,fp,1,fprev,1)
         fpnsav = fpnorm
         delta  = min(Rtwo*delta,stepmx)
         retcd  = 3

      ELSE

!        fpnorm sufficiently small to accept xp as next iterate.
!        Choose new trust region.

         retcd = 0
         IF(ared .ge. Rpten*pred) THEN

!           Not good enough. Decrease trust region for next iteration

            delta = Rhalf*delta
         ELSEIF( ared .le. Rp75*pred ) THEN

!           Wonderful. Increase trust region for next iteration

            delta = min(Rtwo*delta,stepmx)
         ENDIF

      ENDIF

      RETURN
      END SUBROUTINE

END MODULE nwtrup_mod
