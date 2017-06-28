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
MODULE nwutil_mod

USE Kind_Types_Mod
USE nwprnt_mod
USE lautil_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE nwtcvg(xplus,fplus,xc,xtol,retcd,ftol,iter, &
                        maxit,n,ierr,termcd)

      INTEGER(I4B) :: n,iter,maxit,ierr,termcd,retcd
      REAL(DP) :: xtol,ftol
      REAL(DP) :: xplus(*),fplus(*),xc(*)

!-------------------------------------------------------------------------
!
!     Decide whether to terminate the nonlinear algorithm
!
!     Arguments
!
!     In       xplus   Real(*)         new x values
!     In       fplus   Real(*)         new f values
!     In       xc      Real(*)         current x values
!     In       xtol    Real            stepsize tolerance
!     In       retcd   Integer         RETURN code from global search routines
!     In       ftol    Real            function tolerance
!     In       iter    Integer         iteration number
!     In       maxit   Integer         maximum number of iterations allowed
!     In       n       Integer         size of x and f
!     In       ierr    Integer         RETURN code of cndjac (condition estimation)
!
!     Out      termcd  Integer         termination code
!                                       0 no termination criterion satisfied
!                                         ==> continue iterating
!                                       1 norm of scaled function value is
!                                         less than ftol
!                                       2 scaled distance between last
!                                         two steps less than xtol
!                                       3 unsuccessful global strategy
!                                         ==> cannot find a better point
!                                       4 iteration limit exceeded
!                                       5 Jacobian too ill-conditioned
!                                       6 Jacobian singular
!                                       7 Jacobian not usable (all zero entries)
!-------------------------------------------------------------------------

      REAL(DP) :: fmax,rsx
      INTEGER(I4B) :: idamax

!     check whether function values are within tolerance

      termcd = 0

      IF( ierr .ne. 0 ) THEN
         termcd = 4 + ierr
         RETURN
      ENDIF

      fmax = abs(fplus(idamax(n,fplus,1)))
      IF( fmax .le. ftol) THEN
         termcd = 1
         RETURN
      ENDIF

!     initial check at start so there is no xplus
!     so only a check of function values is USEful
      IF(iter .eq. 0) RETURN

      IF(retcd .eq. 1) THEN
         termcd = 3
         RETURN
      ENDIF

!     check whether relative step length is within tolerance
!     Dennis Schnabel Algorithm A7.2.3

      rsx = nuxnrm(n, xplus, xc)
      IF(rsx .le. xtol) THEN
        termcd = 2
        RETURN
      ENDIF

!     check iteration limit

      IF(iter .ge. maxit) THEN
         termcd = 4
      ENDIF

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,iter)
      REAL(DP) :: xc(*),fc(*),fcnorm,xp(*),fp(*),fpnorm
      INTEGER(I4B) :: n, gcnt, priter, iter

!-------------------------------------------------------------------------
!
!     CALLing routine got an error in decomposition/update of Jacobian/Broyden
!     jacobian an singular or too ill-conditioned
!     prepare RETURN arguments
!
!     Arguments
!
!     In       n       Integer         size of x
!     In       xc      Real(*)         current (starting) x values
!     In       fc      Real(*)         function values f(xc)
!     In       fcnorm  Real            norm fc
!     Out      xp      Real(*)         final x values
!     Out      fp      Real(*)         function values f(xp)
!     Out      fpnorm  Real            final norm fp
!     Out      gcnt    Integer         # of backtracking steps (here set to 0)
!     In       priter  Integer         flag for type of output
!     In       iter    Integer         iteration counter
!
!-------------------------------------------------------------------------

      CALL dcopy(n,xc,1,xp,1)
      CALL dcopy(n,fc,1,fp,1)
      fpnorm = fcnorm
      gcnt   = 0
!       IF( priter .gt. 0 ) THEN
!          CALL nwjerr(iter)
!       ENDIF

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE chkjac1(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd)

      INTEGER(I4B) :: lda,n,termcd
      REAL(DP) :: A(lda,*),xc(*),fc(*)
      REAL(DP) :: epsm,scalex(*)
      REAL(DP) :: fz(*),wa(*),xw(*)
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
!     Check the USEr supplied jacobian against its finite difference approximation
!
!     Arguments
!
!     In       A       Real(lda,*)     USEr supplied jacobian
!     In       lda     Integer         leading dimension of ajanal
!     In       xc      Real(*)         vector of x values
!     In       fc      Real(*)         function values f(xc)
!     In       n       Integer         size of x
!     In       epsm    Real            machine precision
!     In       scalex  Real(*)         scaling vector for x()
!     Wk       fz      Real(*)         workspace
!     Wk       wa      Real(*)         workspace
!     Wk       xw      Real(*)         workspace
!     In       fvec    Name            name of routine to evaluate f(x)
!     Out      termcd  Integer         RETURN code
!                                        0  USEr supplied jacobian ok
!                                      -10  USEr supplied jacobian NOT ok
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i,j,errcnt
      REAL(DP) :: ndigit,p,h,xcj,dinf
      REAL(DP) :: tol
      integer idamax

      INTEGER(I4B), PARAMETER :: MAXERR=10

      REAL(DP), PARAMETER :: Rquart=0.25d0, Rten=10.0d0

      termcd = 0

!     compute the finite difference jacobian and check it against
!     the analytic one

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))
      tol    = epsm**Rquart

      errcnt = 0
      CALL dcopy(n,xc,1,xw,1)
      CALL vunsc(n,xw,scalex)

      DO j=1,n
         h = p + p * abs(xw(j))
         xcj   = xw(j)
         xw(j) = xcj + h

!        avoid (small) rounding errors
!        h = xc(j) - xcj but not here to avoid clever optimizers

         h = rnudIF(xw(j), xcj)

         CALL fvec(xw,fz,n,j)
         xw(j) = xcj

         DO i=1,n
            wa(i) = (fz(i)-fc(i))/h
         ENDDO

         dinf = abs(wa(idamax(n,wa,1)))

         DO i=1,n
            IF(abs(A(i,j)-wa(i)).gt.tol*dinf) THEN
               errcnt = errcnt + 1
               IF( errcnt .gt. MAXERR ) THEN
                  termcd = -10
                  RETURN
               ENDIF
               CALL nwckot(i,j,A(i,j),wa(i))
            ENDIF
         ENDDO
      ENDDO

!      CALL vscal(n,xc,scalex)

      IF( errcnt .gt. 0 ) THEN
         termcd = -10
      ENDIF
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE chkjac2(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd, &
                         dsub,dsuper)

      INTEGER(I4B) :: lda,n,termcd,dsub,dsuper
      REAL(DP) :: A(lda,*),xc(*),fc(*)
      REAL(DP) :: epsm,scalex(*)
      REAL(DP) :: fz(*),wa(*),xw(*)
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
!     Check the USEr supplied jacobian against its finite difference approximation
!
!     Arguments
!
!     In       A       Real(lda,*)     USEr supplied jacobian
!     In       lda     Integer         leading dimension of ajanal
!     In       xc      Real(*)         vector of x values
!     In       fc      Real(*)         function values f(xc)
!     In       n       Integer         size of x
!     In       epsm    Real            machine precision
!     In       scalex  Real(*)         scaling vector for x()
!     Wk       fz      Real(*)         workspace
!     Wk       wa      Real(*)         workspace
!     Wk       xw      Real(*)         workspace
!     In       fvec    Name            name of routine to evaluate f(x)
!     Out      termcd  Integer         RETURN code
!                                        0  USEr supplied jacobian ok
!                                      -10  USEr supplied jacobian NOT ok
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i,j,k,dsum,errcnt
      REAL(DP) ::  ndigit,p,h,dinf
      REAL(DP) ::  tol
      REAL(DP) :: w(n),xstep(n)

      INTEGER(I4B), PARAMETER :: MAXERR=10

      REAL(DP), PARAMETER :: Rquart=0.25d0, Rten=10.0d0, Rzero=0.0d0

      dsum = dsub + dsuper + 1

      termcd = 0

!     compute the finite difference jacobian and check it against
!     the USEr supplied one

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))
      tol    = epsm**Rquart

      errcnt = 0
      CALL dcopy(n,xc,1,xw,1)
      CALL vunsc(n,xw,scalex)

      DO j=1,n
          xstep(j) = p + p * abs(xw(j))
          w(j) = xw(j)
      ENDDO

      DO k=1,dsum
         DO j=k,n,dsum
            xw(j) = xw(j) + xstep(j)
         ENDDO

!        for non finite values error message will be wrong
         CALL fvec(xw,fz,n,n+k)

         DO j=k,n,dsum
             h = xstep(j)
             xw(j) = w(j)
             dinf = Rzero
             DO i=max(j-dsuper,1),min(j+dsub,n)
                wa(i) = (fz(i)-fc(i)) / h
                IF(abs(wa(i)).gt.dinf) dinf = abs(wa(i))
             ENDDO

             DO i=max(j-dsuper,1),min(j+dsub,n)
                IF(abs(A(i,j)-wa(i)).gt.tol*dinf) THEN
                   errcnt = errcnt + 1
                   IF( errcnt .gt. MAXERR ) THEN
                      termcd = -10
                      RETURN
                   ENDIF
                   CALL nwckot(i,j,A(i,j),wa(i))
                ENDIF
             ENDDO
         ENDDO
      ENDDO

!      CALL vscal(n,xc,scalex)

      IF( errcnt .gt. 0 ) THEN
         termcd = -10
      ENDIF
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE chkjac(A,lda,xc,fc,n,epsm,jacflg,scalex,fz,wa,xw,fvec, &
                        termcd)

      INTEGER(I4B) :: lda,n,termcd,jacflg(*)
      REAL(DP) :: A(lda,*),xc(*),fc(*)
      REAL(DP) :: epsm,scalex(*)
      REAL(DP) :: fz(*),wa(*),xw(*)
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
!     Check the USEr supplied jacobian against its finite difference approximation
!
!     Arguments
!
!     In       A       Real(lda,*)     USEr supplied jacobian
!     In       lda     Integer         leading dimension of ajanal
!     In       xc      Real(*)         vector of x values
!     In       fc      Real(*)         function values f(xc)
!     In       n       Integer         size of x
!     In       epsm    Real            machine precision
!     In       jacflg  Integer(*)      indicates how to compute jacobian
!                                      jacflg[1]:  0 numeric; 1 USEr supplied; 2 numerical banded
!                                                  3: USEr supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting jacobian allowed when
!                                                   singular or illconditioned
!     In       scalex  Real(*)         scaling vector for x()
!     Wk       fz      Real(*)         workspace
!     Wk       wa      Real(*)         workspace
!     Wk       xw      Real(*)         workspace
!     In       fvec    Name            name of routine to evaluate f(x)
!     Out      termcd  Integer         RETURN code
!                                        0  USEr supplied jacobian ok
!                                      -10  USEr supplied jacobian NOT ok
!
!-------------------------------------------------------------------------

      IF(jacflg(1) .eq. 3) THEN
!        USEr supplied and banded
         CALL chkjac2(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd, &
                      jacflg(2),jacflg(3))
      ELSE
         CALL chkjac1(A,lda,xc,fc,n,epsm,scalex,fz,wa,xw,fvec,termcd)
      ENDIF

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE fdjac0(xc,fc,n,epsm,fvec,fz,rjac,ldr)

      INTEGER(I4B) :: ldr,n
      REAL(DP) :: epsm
      REAL(DP) :: rjac(ldr,*),fz(*),xc(*),fc(*)
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
!     Compute the finite difference jacobian at the current point xc
!
!     Arguments
!
!     In       xc      Real(*)         current point
!     In       fc      Real(*)         function values at current point
!     In       n       Integer         size of x and f
!     In       epsm    Real            machine precision
!     In       fvec    Name            name of routine to evaluate f(x)
!     Wk       fz      Real(*)         workspace
!     Out      rjac    Real(ldr,*)     jacobian matrix at x
!                                        entry [i,j] is derivative of
!                                        f(i) wrt to x(j)
!     In       ldr     Integer         leading dimension of rjac
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i,j
      REAL(DP) :: ndigit,p,h,xcj

      REAL(DP), PARAMETER :: Rten=10d0

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))

      DO j=1,n
         h = p + p * abs(xc(j))

!        or as alternative h  = p * max(Rone, abs(xc(j)))

         xcj   = xc(j)
         xc(j) = xcj + h

!        avoid (small) rounding errors
!        h = xc(j) - xcj  but not here to avoid clever optimizers

         h = rnudIF(xc(j), xcj)
         CALL fvec(xc,fz,n,j)
         xc(j) = xcj
         DO i=1,n
            rjac(i,j) = (fz(i)-fc(i)) / h
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE fdjac2(xc,fc,n,epsm,fvec,fz,rjac,ldr,dsub,dsuper, &
                        w,xstep)

      INTEGER(I4B) :: ldr,n,dsub,dsuper
      REAL(DP) :: epsm
      REAL(DP) :: rjac(ldr,*),fz(*),xc(*),fc(*)
      REAL(DP) :: w(*), xstep(*)
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
!     Compute a banded finite difference jacobian at the current point xc
!
!     Arguments
!
!     In       xc      Real(*)         current point
!     In       fc      Real(*)         function values at current point
!     In       n       Integer         size of x and f
!     In       epsm    Real            machine precision
!     In       fvec    Name            name of routine to evaluate f(x)
!     Wk       fz      Real(*)         workspace
!     Out      rjac    Real(ldr,*)     jacobian matrix at x
!                                        entry [i,j] is derivative of
!                                        f(i) wrt to x(j)
!     In       ldr     Integer         leading dimension of rjac
!     In       dsub    Integer         number of subdiagonals
!     In       dsuper  Integer         number of superdiagonals
!     Internal w       Real(*)         for temporary saving of xc
!     Internal xstep   Real(*)         stepsizes
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i,j,k
      REAL(DP) ::  ndigit,p,h

      REAL(DP), PARAMETER :: Rten=10d0

      INTEGER(I4B) :: dsum

      dsum = dsub + dsuper + 1

      ndigit = -log10(epsm)
      p = sqrt(max(Rten**(-ndigit),epsm))

      DO k=1,n
         xstep(k) = p + p * abs(xc(k))
      ENDDO

      DO k=1,dsum
         DO j=k,n,dsum
            w(j) = xc(j)
            xc(j) = xc(j) + xstep(j)
         ENDDO

         CALL fvec(xc,fz,n,n+k)
         DO j=k,n,dsum
             CALL nuzero(n,rjac(1,j))
!            fdjac0 for why
!            doing this ensures that results for fdjac2 and fdjac0 will be identical
             h = rnudIF(xc(j),w(j))
             xc(j) = w(j)
             DO i=max(j-dsuper,1),min(j+dsub,n)
                rjac(i,j) = (fz(i)-fc(i)) / h
             ENDDO
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      function nudnrm(n, d, x)
      INTEGER(I4B) :: n
      REAL(DP) ::  d(*), x(*)
      REAL(DP) :: nudnrm

!-------------------------------------------------------------------------
!
!     calculate  max( abs(d[*]) / max(x[*],1) )
!
!     Arguments
!
!     In   n        Integer       number of elements in d() and x()
!     In   d        Real(*)       vector d
!     In   x        Real(*)       vector x
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i
      REAL(DP) :: t

      REAL(DP), PARAMETER :: Rzero=0.0d0, Rone=1.0d0

      t = Rzero
      DO i=1,n
         t = max(t, abs(d(i)) / max(abs(x(i)),Rone))
      ENDDO
      nudnrm = t

      RETURN
      END function

!-----------------------------------------------------------------------

      function nuxnrm(n, xn, xc)
      INTEGER(I4B) :: n
      REAL(DP) :: xn(*), xc(*)
      REAL(DP) ::nuxnrm

!-------------------------------------------------------------------------
!
!     calculate  max( abs(xn[*]-xc[*]) / max(xn[*],1) )
!
!     Arguments
!
!     In   n        Integer       number of elements in xn() and xc()
!     In   xn       Real(*)       vector xn
!     In   xc       Real(*)       vector xc
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i
      REAL(DP) :: t

      REAL(DP), PARAMETER :: Rzero=0.0d0, Rone=1.0d0

      t = Rzero
      DO i=1,n
         t = max(t, abs(xn(i)-xc(i)) / max(abs(xn(i)),Rone))
      ENDDO
      nuxnrm = t

      RETURN
      END function

!-----------------------------------------------------------------------

      function rnudIF(x, y)
      REAL(DP) :: x, y
      REAL(DP) :: rnudif

!-------------------------------------------------------------------------
!
!     Return difference of x and y (x - y)
!
!     Arguments
!
!     In   x  Real      argument 1
!     In   y  Real      argument 2
!
!-------------------------------------------------------------------------

      rnudif = x - y
      RETURN
      END function

!-----------------------------------------------------------------------

      SUBROUTINE compmu(r,ldr,n,mu,y,ierr)

      INTEGER(I4B) :: ldr,n,ierr
      REAL(DP) :: r(ldr,*),mu,y(*)

!-------------------------------------------------------------------------
!
!     Compute a small perturbation mu for the (almost) singular matrix R.
!     mu is USEd in the computation of the Levenberg-Marquardt step.
!
!     Arguments
!
!     In       R       Real(ldr,*)     upper triangular matrix from QR
!     In       ldr     Integer         leading dimension of R
!     In       n       Integer         column dimension of R
!     Out      mu      Real            sqrt(l1 norm of R * infinity norm of R
!                                      * n * epsm * 100) designed to make
!                                        trans(R)*R + mu * I not singular
!     Wk       y       Real(*)         workspace for dlange
!     Out      ierr    Integer         0 indicating mu ok
!                                      3 indicating mu much too small
!
!-------------------------------------------------------------------------

      REAL(DP) :: aifnrm,al1nrm,epsm
      REAL(DP) :: dlantr

      REAL(DP), PARAMETER :: Rhund=100d0

!     get the infinity norm of R
!     get the l1 norm of R
      ierr = 0
      aifnrm = dlantr('I','U','N',n,n,r,ldr,y)
      al1nrm = dlantr('1','U','N',n,n,r,ldr,y)
      epsm = epsmch()
      mu = sqrt(n*epsm*Rhund)*aifnrm*al1nrm
!     matrix consists of zero's or near zero's
!     LM correction in liqrev will not work
      IF( mu .le. Rhund*epsm ) THEN
         ierr = 3
      ENDIF
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE cndjac(n,r,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)
      INTEGER(I4B) :: n,ldr,icdwrk(*),ierr
      REAL(DP) :: cndtol,rcond,r(ldr,*),rcdwrk(*)

!---------------------------------------------------------------------
!
!     Check r for singularity and/or ill conditioning
!
!     Arguments
!
!     In       n       Integer         dimension of problem
!     In       r       Real(ldr,*)     upper triangular R from QR decomposition
!     In       ldr     Integer         leading dimension of rjac
!     In       cndtol  Real            tolerance of test for ill conditioning
!                                       when rcond <= cndtol THEN ierr is set to 1
!                                       cndtol should be >= machine precision
!     Out      rcond   Real            inverse condition  of r
!     Wk       rcdwrk  Real(*)         workspace (for dtrcon)
!     Wk       icdwrk  Integer(*)      workspace (fordtrcon)
!     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
!                                      1 indicating Jacobian too ill-conditioned
!                                      2 indicating Jacobian completely singular
!
!---------------------------------------------------------------------

      INTEGER(I4B) :: i,info
      LOGICAL rsing
      REAL(DP), PARAMETER :: Rzero=0.0d0

      ierr = 0

      rsing = .false.
      DO i=1,n
         IF( r(i,i) .eq. Rzero ) THEN
             rsing = .true.
         ENDIF
      ENDDO

      IF( rsing ) THEN
         ierr = 2
         rcond = Rzero
      ELSE
         CALL dtrcon('1','U','N',n,r,ldr,rcond,rcdwrk,icdwrk,info)
         IF( rcond .eq. Rzero ) THEN
             ierr = 2
         ELSEIF( rcond .le. cndtol ) THEN
             ierr = 1
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE nwfjac(x,scalex,f,fq,n,epsm,jacflg,fvec,mkjac,rjac, &
                        ldr,xw,w,xstep)

      INTEGER(I4B) :: ldr,n,jacflg(*)
      REAL(DP) :: epsm
      REAL(DP) :: x(*),f(*),scalex(*),xw(*),w(*),xstep(*)
      REAL(DP) :: rjac(ldr,*),fq(*)
      INTERFACE
      SUBROUTINE fvec ( x, r, n, flag )
        USE Kind_Types_Mod
        INTEGER(I4B), intent(in) :: n, flag
        REAL(DP), dimension(n), intent(in) :: x
        REAL(DP), dimension(n), intent(out) :: r
      END SUBROUTINE fvec
      SUBROUTINE mkjac ( jac, ldr, x, n )
        USE Kind_Types_Mod
        integer(i4b), intent(in) :: n, ldr
        REAL(DP), intent(in), dimension(n)  :: x(n)
        REAL(DP), intent(out) :: jac(ldr,*)
      END SUBROUTINE mkjac
      END INTERFACE

!-------------------------------------------------------------------------
!
!     Calculate the jacobian  matrix
!
!     Arguments
!
!     In       x       Real(*)         (scaled) current x values
!     In       scalex  Real(*)         scaling factors x
!     In       f       Real(*)         function values f(x)
!     Wk       fq      Real(*)         (internal) workspace
!     In       n       Integer         size of x and f
!     In       epsm    Real            machine precision
!     In       jacflg  Integer(*)      indicates how to compute jacobian
!                                      jacflg[1]:  0 numeric; 1 USEr supplied; 2 numerical banded
!                                                  3: USEr supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting jacobian allowed when
!                                                   singular or illconditioned
!     In       fvec    Name            name of routine to evaluate f()
!     In       mkjac   Name            name of routine to evaluate jacobian
!     Out      rjac    Real(ldr,*)     jacobian matrix (unscaled)
!     In       ldr     Integer         leading dimension of rjac
!     Internal xw      Real(*)         USEd for storing unscaled x
!     Internal w       Real(*)         workspace for banded jacobian
!     Internal xstep   Real(*)         workspace for banded jacobian
!
!-------------------------------------------------------------------------

!     compute the finite difference or analytic jacobian at x

      CALL dcopy(n,x,1,xw,1)
      CALL vunsc(n,xw,scalex)
      IF(jacflg(1) .eq. 0) THEN
         CALL fdjac0(xw,f,n,epsm,fvec,fq,rjac,ldr)
      ELSEIF(jacflg(1) .eq. 2) THEN
         CALL fdjac2(xw,f,n,epsm,fvec,fq,rjac,ldr,jacflg(2),jacflg(3), &
                     w,xstep)
      ELSE
         CALL mkjac(rjac,ldr,xw,n)
      ENDIF

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE nwscjac(n,rjac,ldr,scalex)
      INTEGER(I4B) :: n, ldr
      REAL(DP) :: rjac(ldr,*), scalex(*)

!-------------------------------------------------------------------------
!
!     Scale jacobian
!
!     Arguments
!
!     In       n       Integer         size of x and f
!     Inout    rjac    Real(ldr,*)     jacobian matrix
!     In       ldr     Integer         leading dimension of rjac
!     In       scalex  Real(*)         scaling factors for x
!
!-------------------------------------------------------------------------

      INTEGER(I4B) j
      REAL(DP) :: t
      REAL(DP), PARAMETER :: Rone=1.0d0

      DO j = 1,n
         t = Rone/scalex(j)
         CALL dscal(n,t,rjac(1,j),1)
      ENDDO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE nwunscjac(n,rjac,ldr,scalex)
      INTEGER(I4B) :: n, ldr
      REAL(DP) :: rjac(ldr,*), scalex(*)

!-------------------------------------------------------------------------
!
!     Unscale jacobian
!
!     Arguments
!
!     In       n       Integer         size of x and f
!     Inout    rjac    Real(ldr,*)     jacobian matrix
!     In       ldr     Integer         leading dimension of rjac
!     In       scalex  Real(*)         scaling factors for x
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: j
      REAL(DP) :: t

      DO j = 1,n
         t = scalex(j)
         CALL dscal(n,t,rjac(1,j),1)
      ENDDO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE nwcpsx(n,rjac,ldr,scalex,epsm, mode)

      INTEGER(I4B) :: ldr,n,mode
      REAL(DP) :: epsm
      REAL(DP) :: scalex(*)
      REAL(DP) :: rjac(ldr,*)

!-------------------------------------------------------------------------
!
!     Calculate scaling factors from the jacobian  matrix
!
!     Arguments
!
!     In       n       Integer         size of x and f
!     In       rjac    Real(ldr,*)     jacobian matrix
!     In       ldr     Integer         leading dimension of rjac
!     Out      scalex  Real(*)         scaling factors for x
!     In       epsm    Real            machine precision
!     In       mode    Integer         1: initialise, >1: adjust
!-------------------------------------------------------------------------

      INTEGER(I4B) :: k
      REAL(DP) :: dnrm2

      IF( mode .eq. 1 ) THEN
         DO k=1,n
            scalex(k) = dnrm2(n,rjac(1,k),1)
            IF( scalex(k) .le. epsm ) scalex(k) = 1
         ENDDO
      ELSE IF( mode .gt. 1 ) THEN
         DO k=1,n
            scalex(k) = max(scalex(k),dnrm2(n,rjac(1,k),1))
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE nwcpmt(n, x, scalex, factor, wrk, stepsiz)
      REAL(DP) :: x(*), scalex(*), wrk(*)
      REAL(DP) :: factor, stepsiz
      INTEGER(I4B) :: n

!-------------------------------------------------------------------------
!
!     Calculate maximum stepsize
!
!     Arguments
!
!     In       n       Integer     size of x
!     In       x       Real(*)     x-values
!     In       scalex  Real(*)     scaling factors for x
!     In       factor  Real        multiplier
!     Inout    wrk     Real(*)     workspace
!     Out      stepsiz Real        stepsize
!
!     Currently not USEd
!     Minpack USEs this to calculate initial trust region size
!     Not (yet) USEd in this code becaUSE it doesn't seem to help
!     Manually setting an initial trust region size works better
!
!-------------------------------------------------------------------------

      REAL(DP), PARAMETER :: Rzero=0.0d0

      REAL(DP) :: dnrm2

      CALL dcopy(n,x,1,wrk,1)
      CALL vscal(n,wrk,scalex)
      stepsiz = factor * dnrm2(n,wrk,1)
      IF( stepsiz .eq. Rzero ) stepsiz = factor
      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE vscal(n,x,sx)

      INTEGER(I4B) :: n
      REAL(DP) :: x(*),sx(*)

!-------------------------------------------------------------------------
!
!     Scale a vector x
!
!     Arguments
!
!     In       n       Integer         size of x
!     Inout    x       Real(*)         vector to scale
!     In       sx      Real(*)         scaling vector
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i

      DO i = 1,n
         x(i) = sx(i) * x(i)
      ENDDO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE vunsc(n,x,sx)

      INTEGER(I4B) :: n
      REAL(DP) :: x(*),sx(*)

!-------------------------------------------------------------------------
!
!     Unscale a vector x
!
!     Arguments
!
!     In       n       Integer         size of x
!     Inout    x       Real(*)         vector to unscale
!     In       sx      Real(*)         scaling vector
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i

      DO i = 1,n
         x(i) = x(i) / sx(i)
      ENDDO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE nwfvec(x,n,scalex,fvec,f,fnorm,xw)

      INTEGER(I4B) :: n
      REAL(DP) :: x(*),xw(*),scalex(*),f(*),fnorm
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
!     Evaluate the function at current iterate x and scale its value
!
!     Arguments
!
!     In       x       Real(*)         x
!     In       n       Integer         size of x
!     In       scalex  Real(*)         scaling vector for x
!     In       fvec    Name            name of routine to calculate f(x)
!     Out      f       Real(*)         f(x)
!     Out      fnorm   Real            .5*||f(x)||**2
!     Internal xw      Real(*)         USEd for storing unscaled xc
!
!-------------------------------------------------------------------------

      REAL(DP) :: dnrm2

      REAL(DP), PARAMETER :: Rhalf=0.5d0

      CALL dcopy(n,x,1,xw,1)
      CALL vunsc(n,xw,scalex)
      CALL fvec(xw,f,n,0)

      fnorm = Rhalf * dnrm2(n,f,1)**2

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      function epsmch()

!     Return machine precision
!     Use Lapack routine

      REAL(DP) :: epsmch
      REAL(DP) :: dlamch
      external dlamch

!     dlamch('e') RETURNs negeps (1-eps)
!     dlamch('p') RETURNs 1+eps

      epsmch = dlamch('p')

      RETURN
      END function

!-----------------------------------------------------------------------

      function dblhuge()

!     Return largest double precision number
!     Use Lapack routine

      REAL(DP) :: dblhuge
      REAL(DP) :: dlamch
      external dlamch

!     dlamch('o') RETURNs max double precision

      dblhuge = dlamch('o')

      RETURN
      END function

END MODULE nwutil_mod
