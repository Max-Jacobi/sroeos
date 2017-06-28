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
MODULE nwddlg_mod

USE Kind_Types_Mod
USE nwtrup_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE nwddlg(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol, &
                        delta,qtf,scalex,fvec,d,xprev, &
                        ssd,v,wa,fprev,xp,fp,fpnorm,retcd,gcnt, &
                        priter,iter)

      INTEGER(I4B) :: ldr, n, retcd, gcnt, priter, iter
      REAL(DP) :: fcnorm, stepmx, xtol, fpnorm, delta
      REAL(DP) :: rjac(ldr,*), dn(*), g(*), xc(*), qtf(*)
      REAL(DP) :: scalex(*), d(*)
      REAL(DP) :: xprev(*), xp(*), fp(*)
      REAL(DP) :: ssd(*), v(*), wa(*), fprev(*)
      INTERFACE
      SUBROUTINE fvec ( x, r, n, flag )
        USE Kind_Types_Mod
        INTEGER(I4B), intent(in) :: n, flag
        REAL(DP), dimension(n), intent(in) :: x
        REAL(DP), dimension(n), intent(out) :: r
      END SUBROUTINE fvec
      END INTERFACE

!-------------------------------------------------------------------------
!
!     Find a next iterate xp by the double dogleg method
!
!     Arguments
!
!     In       n       Integer         size of problem: dimension x, f
!     In       Rjac    Real(ldr,*)     R of QR-factored jacobian
!     In       ldr     Integer         leading dimension of Rjac
!     Inout    dn      Real(*)         newton direction
!     Inout    g       Real(*)         gradient at current point
!                                      trans(jac)*f()
!     In       xc      Real(*)         current iterate
!     In       fcnorm  Real            .5*||f(xc)||**2
!     In       stepmx  Real            maximum stepsize
!     In       xtol    Real            x-tolerance (stepsize)
!     Inout    delta   Real            on input: initial trust region radius
!                                                if -1 THEN set to something
!                                                reasonable
!                                      on output: final value
!                                      ! Do not modify between CALLs WHILE
!                                        still iterating
!     In       qtf     Real(*)         trans(Q)*f(xc)
!     In       scalex  Real(*)         scaling factors for x()
!     In       fvec    Name            name of SUBROUTINE to evaluate f(x)
!                                      ! must be declared external in CALLer
!     Wk       d       Real(*)         work vector
!     Wk       xprev   Real(*)         work vector
!     Wk       ssd     Real(*)         work vector
!     Wk       v       Real(*)         work vector
!     Wk       wa      Real(*)         work vector
!     Wk       fprev   Real(*)         work vector
!     Out      xp      Real(*)         new x()
!     Out      fp      Real(*)         new f(xp)
!     Out      fpnorm  Real            new .5*||f(xp)||**2
!     Out      retcd   Integer         RETURN code
!                                       0  new satisfactory x() found
!                                       1  no  satisfactory x() found
!     Out      gcnt    Integer         number of steps taken
!     In       priter  Integer         print flag
!                                       -1 no intermediate printing
!                                       >0 yes for print of intermediate results
!     In       iter    Integer         current iteration (only USEd for above)
!
!     All vectors at least size n
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i
      REAL(DP) :: dnlen,ssdlen,alpha,beta,lambda,fpred
      REAL(DP) :: sqalpha,eta,gamma,fpnsav,oarg(7)
      REAL(DP) :: dnrm2,ddot
      LOGICAL nwtstep
      INTEGER(I4B) :: dtype

      INTEGER(I4B) :: idamax

      REAL(DP), PARAMETER :: Rhalf=0.5d0, Rone=1.0d0, Rtwo=2.0d0, &
      Rten=10.0d0, Rp2 = Rtwo/Rten, Rp8 = Rone - Rp2, Rzero=0.0d0

!     length newton direction

      dnlen = dnrm2(n, dn, 1)

!     steepest descent direction and length

      sqalpha = dnrm2(n,g,1)
      alpha   = sqalpha**2

      CALL dcopy(n, g, 1, d, 1)
      CALL dtrmv('U','N','N',n,rjac,ldr,d,1)
      beta = dnrm2(n,d,1)**2

      CALL dcopy(n, g, 1, ssd, 1)
      CALL dscal(n, -(alpha/beta), ssd, 1)

      ssdlen = alpha*sqalpha/beta

!     set trust radius to ssdlen or dnlen if required

      IF( delta .eq. -Rone ) THEN
         delta = min(ssdlen, stepmx)
      ELSEIF( delta .eq. -Rtwo ) THEN
         delta = min(dnlen, stepmx)
      ENDIF

!     calculate double dogleg PARAMETER

      gamma = alpha*alpha/(-beta*ddot(n,g,1,dn,1))
!      CALL dgdbg(gamma, alpha*alpha, -beta*ddot(n,g,1,dn,1))
!     precautionary (just in case)
      eta = max(Rzero, min(Rone,Rp2 + Rp8*gamma))

      retcd = 4
      gcnt  = 0

      DO WHILE( retcd .gt. 1 )
!        find new step by double dogleg algorithm

         CALL ddlgstp(n,dn,dnlen,delta,v, &
                     ssd,ssdlen,eta,d,dtype,lambda)
         nwtstep = dtype .eq. 4
!        compute the model prediction 0.5*||F + J*d||**2 (L2-norm)

         CALL dcopy(n,d,1,wa,1)
         CALL dtrmv('U','N','N',n,rjac,ldr,wa,1)
         CALL daxpy(n, Rone, qtf,1,wa,1)
         fpred = Rhalf * dnrm2(n,wa,1)**2

!        evaluate function at xp = xc + d

         DO i=1,n
            xp(i) = xc(i) + d(i)
         ENDDO

         CALL nwfvec(xp,n,scalex,fvec,fp,fpnorm,wa)
         gcnt = gcnt + 1

!        check whether the global step is acceptable

         oarg(2) = delta
         CALL nwtrup(n,fcnorm,g,d,nwtstep,stepmx,xtol,delta, &
                     fpred,retcd,xprev,fpnsav,fprev,xp,fp,fpnorm)

!          IF( priter .gt. 0 ) THEN
!             oarg(1) = lambda
!             oarg(3) = delta
!             oarg(4) = eta
!             oarg(5) = fpnorm
!             oarg(6) = abs(fp(idamax(n,fp,1)))
!             CALL nwdgot(iter,dtype,retcd,oarg)
!          ENDIF

      ENDDO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE ddlgstp(n,dn,dnlen,delta,v, &
                        ssd,ssdlen,eta,d,dtype,lambda)
      INTEGER(I4B) :: n
      REAL(DP) :: dn(*), ssd(*), v(*), d(*)
      REAL(DP) :: dnlen, delta, ssdlen, eta, lambda
      INTEGER(I4B) dtype

!-------------------------------------------------------------------------
!
!     Find a new step by the double dogleg algorithm
!     Internal routine for nwddlg
!
!     Arguments
!
!     In       n       Integer         size of problem
!     In       dn      Real(*)         current newton step
!     Out      dnlen   Real            length dn()
!     In       delta   Real            current trust region radius
!     Out      v       Real(*)         (internal) eta * dn() - ssd()
!     In       ssd     Real(*)         (internal) steepest descent direction
!     In       ssdlen  Real            (internal) length ssd
!     In       eta     Real            (internal) double dogleg PARAMETER
!     Out      d       Real(*)         new step for x()
!     Out      dtype   Integer         steptype
!                                       1 steepest descent
!                                       2 combination of dn and ssd
!                                       3 partial newton step
!                                       4 full newton direction
!     Out      lambda  Real            weight of eta*dn() in d()
!                                      closer to 1 ==> more of eta*dn()
!
!-----------------------------------------------------------------------

      INTEGER(I4B) :: i
      REAL(DP) :: vssd, vlen
      REAL(DP) :: dnrm2, ddot

      IF(dnlen .le. delta) THEN

!        Newton step smaller than trust radius ==> take it

         CALL dcopy(n, dn, 1, d, 1)
         delta = dnlen
         dtype = 4

      ELSEIF(eta*dnlen .le. delta) THEN

!        take partial step in newton direction

         CALL dcopy(n, dn, 1, d, 1)
         CALL dscal(n, delta / dnlen, d, 1)
         dtype = 3

      ELSEIF(ssdlen .ge. delta) THEN

!        take step in steepest descent direction

         CALL dcopy(n, ssd, 1, d, 1)
         CALL dscal(n, delta / ssdlen, d, 1)
         dtype = 1

      ELSE

!        calculate convex combination of ssd and eta*dn with length delta

         DO i=1,n
            v(i) = eta*dn(i) - ssd(i)
         ENDDO

         vssd = ddot(n,v,1,ssd,1)
         vlen = dnrm2(n,v,1)**2

         lambda =(-vssd+sqrt(vssd**2-vlen*(ssdlen**2-delta**2)))/vlen
         CALL dcopy(n, ssd, 1, d, 1)
         CALL daxpy(n, lambda, v, 1, d, 1)
         dtype = 2

      ENDIF

      RETURN
      END SUBROUTINE

END MODULE nwddlg_mod
