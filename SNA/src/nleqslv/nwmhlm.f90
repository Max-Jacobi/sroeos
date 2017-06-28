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
MODULE nwmhlm_mod

USE Kind_Types_Mod
USE limhpar_mod
USE nwutil_mod
USE nwtrup_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE nwmhlm(n,rjac,ldr,dn,g,xc,fcnorm,stepmx,xtol, &
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
!     Find a next iterate xp by the More-Hebden-Levenberg-Marquardt method
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
      REAL(DP) :: dnlen,glen,ssdlen,alpha,beta,mu,fpred
      REAL(DP) :: fpnsav,oarg(6)
      REAL(DP) :: dnrm2
      LOGICAL nwtstep
      INTEGER(I4B) :: dtype

      INTEGER(I4B) :: idamax

      REAL(DP), PARAMETER ::Rhalf=0.5d0, Rone=1.0d0, Rtwo=2.0d0

!     length newton direction

      dnlen = dnrm2(n, dn, 1)

!     gradient length and steepest descent direction and length

      glen  = dnrm2(n,g,1)
      alpha = glen**2

      CALL dcopy(n, g, 1, d, 1)
      CALL dtrmv('U','N','N',n,rjac,ldr,d,1)
      beta = dnrm2(n,d,1)**2

      CALL dcopy(n, g, 1, ssd, 1)
      CALL dscal(n, -(alpha/beta), ssd, 1)

      ssdlen = alpha*glen/beta

!     set trust radius to ssdlen or dnlen if required

      IF( delta .eq. -Rone ) THEN
         delta = min(ssdlen, stepmx)
      ELSEIF( delta .eq. -Rtwo ) THEN
         delta = min(dnlen, stepmx)
      ENDIF

      retcd = 4
      gcnt  = 0

      DO WHILE( retcd .gt. 1 )
!        find new step by More Hebden LM  algorithm
!        reUSE ssd as sdiag

         CALL nwmhstep(Rjac,ldr,n,ssd,qtf,dn,dnlen,glen,delta,mu, &
                       d, v, dtype)
         nwtstep = dtype .eq. 2
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
!             oarg(1) = mu
!             oarg(3) = delta
!             oarg(4) = dnrm2(n, d, 1)
!             oarg(5) = fpnorm
!             oarg(6) = abs(fp(idamax(n,fp,1)))
!             CALL nwmhot(iter,dtype,retcd,oarg)
!          ENDIF

      ENDDO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE nwmhstep(R,ldr,n,sdiag,qtf,dn,dnlen,glen,delta,mu, &
                        d, work, dtype)
      INTEGER(I4B) :: ldr, n
      REAL(DP) :: R(ldr,*)
      REAL(DP) :: sdiag(*), qtf(*), dn(*), d(*), work(*)
      REAL(DP) :: dnlen, glen, delta, mu
      INTEGER(I4B) :: dtype

!-------------------------------------------------------------------------
!
!     Find a new step by the More Hebden Levemberg Marquardt algorithm
!     Internal routine for nwmhlm
!
!     Arguments
!
!     In       R       Real(ldr,*)     R of QR-factored jacobian
!     In       ldr     Integer         leading dimension of R
!     In       n       Integer         size of problem
!     Out      sdiag   Real(*)         diagonal of LM lower triangular modified R
!     In       qtf     Real(*)         trans(Q)*f(xc)
!     In       dn      Real(*)         current newton step
!     Out      dnlen   Real            length dn()
!     In       glen    Real            length gradient
!     In       delta   Real            current trust region radius
!     Inout    mu      Real            Levenberg-Marquardt PARAMETER
!     Out      d       Real(*)         new step for x()
!     Work     work    Real(*)         work vector for limhpar
!     Out      dtype   Integer         steptype
!                                       1 LM step
!                                       2 full newton direction
!
!-----------------------------------------------------------------------

      REAL(DP), PARAMETER :: Rone=1.0D0

      IF(dnlen .le. delta) THEN

!        Newton step smaller than trust radius ==> take it

         CALL dcopy(n, dn, 1, d, 1)
         delta = dnlen
         dtype = 2

      ELSE

!        calculate LM step
         CALL limhpar(R, ldr, n, sdiag, qtf, dn, dnlen, glen, delta, &
                      mu, d, work)
!        change sign of step d (limhpar solves for trans(R)*R+mu*I)=qtf instead of -qtf)
         CALL dscal(n,-Rone,d,1)
         dtype = 1
      ENDIF

      RETURN
      END SUBROUTINE

END MODULE nwmhlm_mod
