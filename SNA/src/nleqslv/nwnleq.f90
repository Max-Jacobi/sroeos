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
MODULE nwnleq_mod

USE Kind_Types_Mod
USE nwbrdn_mod
USE nwnwtn_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE nwnleq(x0,n,scalex,maxit, &
                        jacflg,xtol,ftol,btol,cndtol,method,global, &
                        xscalm,stepmx,delta,sigma,rjac,ldr, &
                        rwork,lrwork, &
                        rcdwrk,icdwrk,qrwork,qrwsiz,fjac,fvec,outopt,xp, &
                        fp,gp,njcnt,nfcnt,iter,termcd)

      INTEGER(I4B) :: n,jacflg(*),maxit,njcnt,nfcnt,iter,termcd,method
      INTEGER(I4B) :: global,xscalm,ldr,lrwork,qrwsiz
      INTEGER(I4B) :: outopt(*)
      REAL(DP) :: xtol,ftol,btol,cndtol,stepmx,delta,sigma
      REAL(DP) :: xp(*),fp(*),gp(*),x0(*)
      REAL(DP) :: rjac(ldr,*),rwork(*),rcdwrk(*),qrwork(*)
      REAL(DP) :: scalex(*)
      INTEGER(I4B) ::           icdwrk(*)
      INTERFACE
      SUBROUTINE fvec ( x, r, n, flag )
        USE Kind_Types_Mod
        INTEGER(I4B), intent(in) :: n, flag
        REAL(DP), dimension(n), intent(in) :: x
        REAL(DP), dimension(n), intent(out) :: r
      END SUBROUTINE fvec
      SUBROUTINE fjac ( jac, ldr, x, n )
        USE Kind_Types_Mod
        integer (i4b), intent(in) :: n, ldr
        REAL(DP), intent(in), dimension(n)  :: x(n)
        REAL(DP), intent(out) :: jac(ldr,*)
      END SUBROUTINE fjac
      END INTERFACE

!-------------------------------------------------------------------------
!
!     Solves systems of nonlinear equations using the Newton / Broyden
!     method with a global strategy either linesearch or double dogleg
!
!     In       x0      Real(*)         starting vector for x
!     In       n       Integer         dimension of problem
!     Inout    scalex  Real(*)         scaling factors x()
!     Inout    maxit   Integer         maximum number iterations
!     Inout    jacflg  Integer(*)      jacobian flag array
!                                      jacflg[1]:  0 numeric; 1 USEr supplied; 2 numerical banded
!                                                  3: USEr supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting step allowed when
!                                                   singular or illconditioned
!     Inout    xtol    Real            x tolerance
!     Inout    ftol    Real            f tolerance
!     Inout    btol    Real            x tolerance for backtracking
!     Inout    cndtol  Real            tolerance of test for ill conditioning
!     Inout    method  Integer         method to USE
!                                        0 Newton
!                                        1 Broyden
!     In       global  Integer         global strategy to USE
!                                        1 cubic linesearch
!                                        2 quadratic linesearch
!                                        3 geometric linesearch
!                                        4 double dogleg
!                                        5 powell dogleg
!                                        6 hookstep (More-Hebden Levenberg-Marquardt)
!     In       xscalm  Integer         scaling method
!                                        0 scale fixed and supplied by USEr
!                                        1 for scale from jac. columns a la Minpack
!     Inout    stepmx  Real            maximum stepsize
!     Inout    delta   Real            trust region radius
!                                        > 0.0 or special value for initial value
!                                        -1.0  ==> USE min(Cauchy length, stepmx)
!                                        -2.0  ==> USE min(Newton length, stepmx)
!     Inout    sigma   Real            reduction factor geometric linesearch
!     Inout    rjac    Real(ldr,*)     workspace jacobian
!                                         2*n*n for Broyden and n*n for Newton
!     In       ldr     Integer         leading dimension rjac
!     Out      rwork   Real(*)         real workspace (9n)
!     In       lrwork  Integer         size real workspace
!     In       rcdwrk  Real(*)         workspace for Dtrcon (3n)
!     In       icdwrk  Integer(*)      workspace for Dtrcon (n)
!     In       qrwork  Real(*)         workspace for Lapack QR routines (CALL liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!     In       fjac    Name            optional name of routine to calculate
!                                      USEr supplied jacobian
!     In       fvec    Name            name of routine to calculate f(x)
!     In       outopt  Integer(*)      output options
!                                       outopt(1)
!                                         0 no output
!                                         1 output an iteration report
!                                       outopt(2)
!                                         0 DO not check any USEr supplied jacobian
!                                         1 check USEr supplied jacobian if supplied
!     Out      xp      Real(*)         final values for x()
!     Out      fp      Real(*)         final values for f(x)
!     Out      gp      Real(*)         gradient of f() at xp()
!     Out      njcnt   Integer         number of jacobian evaluations
!     Out      nfcnt   Integer         number of function evaluations
!     Out      iter    Integer         number of (outer) iterations
!     Out      termcd  Integer         termination code
!                                       > 0 process terminated
!                                             1  function criterion near zero
!                                             2  no better point found
!                                             3  x-values within tolerance
!                                             4  iteration limit exceeded
!                                             5  singular/ill-conditioned jacobian
!                                             6  totally singular jacobian
!                                                (when allowSingular=TRUE)
!
!                                       < 0 invalid input PARAMETERs
!                                            -1  n not positive
!                                            -2  insufficient workspace rwork
!                                            -3  cannot check USEr supplied jacobian (not supplied)
!
!    The SUBROUTINE fvec must be declared as
!
!!        SUBROUTINE fvec(x,f,n,flag)
!         double precision x(*), f(*)
!         integer  n, flag
!
!         x() are the x values for which to calculate the function values f(*)
!         The dimension of these vectors is n
!         The flag argument is set to
!            0  for calculation of function values
!           >0  indicating that jacobian column <flag> is being computed
!               so that fvec can abort.
!
!    The SUBROUTINE fjac must be declared as
!
!!        SUBROUTINE mkjac(rjac,ldr,x,n)
!         integer ldr
!         double precision rjac(ldr,*), x(*)
!         integer  n
!
!         The routine calculates the jacobian in point x(*) of the
!         function. If any illegal values are encountered during
!         calculation of the jacobian it is the responsibility of
!         the routine to quit.

!-------------------------------------------------------------------------

      REAL(DP) :: epsm

!     check input PARAMETERs

      CALL nwpchk(n,lrwork,xtol,ftol,btol,cndtol,maxit, &
                  jacflg,method,global,stepmx,delta,sigma, &
                  epsm,outopt,scalex,xscalm,termcd)
      IF(termcd .lt. 0) THEN
         RETURN
      ENDIF

!     first argument of nwsolv/brsolv is leading dimension of rjac in those routines
!     should be at least n

      IF( method .eq. 0 ) THEN
         CALL nwsolv(ldr,x0,n,scalex,maxit,jacflg, &
                     xtol,ftol,btol,cndtol,global,xscalm, &
                     stepmx,delta,sigma, &
                     rjac, &
                     rwork(1    ),rwork(1+  n), &
                     rwork(1+2*n),rwork(1+3*n), &
                     rwork(1+4*n),rwork(1+5*n), &
                     rwork(1+6*n),rwork(1+7*n), &
                     rwork(1+8*n),rcdwrk,icdwrk,qrwork,qrwsiz,epsm, &
                     fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter,termcd)

      ELSEIF( method .eq. 1 ) THEN

         CALL brsolv(ldr,x0,n,scalex,maxit,jacflg, &
                     xtol,ftol,btol,cndtol,global,xscalm, &
                     stepmx,delta,sigma, &
                     rjac,rjac(1,n+1), &
                     rwork(1    ),rwork(1+  n), &
                     rwork(1+2*n),rwork(1+3*n), &
                     rwork(1+4*n),rwork(1+5*n), &
                     rwork(1+6*n),rwork(1+7*n), &
                     rwork(1+8*n),rcdwrk,icdwrk, &
                     qrwork,qrwsiz,epsm, &
                     fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter,termcd)

      ENDIF

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE nwpchk(n,lrwk, &
                        xtol,ftol,btol,cndtol,maxit,jacflg,method, &
                        global,stepmx,delta,sigma,epsm,outopt, &
                        scalex,xscalm,termcd)

      INTEGER(I4B), intent(in) :: n,lrwk
      INTEGER(I4B), intent(inout) :: jacflg(*)
      INTEGER(I4B) :: method,global,maxit,xscalm,termcd
      INTEGER(I4B) :: outopt(*)
      REAL(DP) :: xtol,ftol,btol,cndtol,stepmx,delta,sigma,epsm
      REAL(DP) :: scalex(*)

!-------------------------------------------------------------------------
!
!     Check input arguments for consistency and modify if needed/harmless
!
!     Arguments
!
!     In       n       Integer         dimension of problem
!     In       lrwk    Integer         size real workspace
!     Inout    xtol    Real            x tolerance
!     Inout    ftol    Real            f tolerance
!     Inout    btol    Real            x tolerance for backtracking
!     Inout    cndtol  Real            tolerance of test for ill conditioning
!     Inout    maxit   Integer         maximum number iterations
!     Inout    jacflg  Integer(*)      jacobian flag
!     Inout    method  Integer         method to USE (Newton/Broyden)
!     Inout    global  Integer         global strategy to USE
!     Inout    stepmx  Real            maximum stepsize
!     Inout    delta     Real            trust region radius
!     Inout    sigma   Real            reduction factor geometric linesearch
!     Out      epsm                    machine precision
!     Inout    scalex  Real(*)         scaling factors x()
!     Inout    xscalm  integer         0 for fixed scaling, 1 for automatic scaling
!     Out      termcd  Integer         termination code (<0 on errors)
!
!-------------------------------------------------------------------------

      INTEGER(I4B) :: i,len
      REAL(DP) :: Rhuge

      REAL(DP), PARAMETER :: Rzero=0.0d0, Rone=1.0d0, Rtwo=2.0d0, Rthree=3.0d0, Rhalf = 0.5d0
!     check that PARAMETERs only take on acceptable values
!     if not, set them to default values

!     initialize termcd to all ok

      termcd = 0

!     compute machine precision

      epsm = epsmch()

!     get largest double precision number
      Rhuge = dblhuge()

!     check dimensions of the problem

      IF(n .le. 0) THEN
         termcd = -1
         RETURN
      ENDIF

!     check dimensions of workspace arrays

      len = 9*n
!      +2*n*n
      IF(lrwk .lt. len) THEN
         termcd = -2
         RETURN
      ENDIF

!     check jacflg, method, and global

      IF(jacflg(1) .gt. 3 .or. jacflg(1) .lt. 0) jacflg(1) = 0

      IF(method .lt. 0 .or. method .gt. 1) method = 0

      IF(global .lt. 0 .or. global .gt. 6) global = 4

!     set outopt to correct values

      IF(outopt(1) .ne. 0 ) THEN
         outopt(1) = 1
      ENDIF

      IF(outopt(2) .ne. 0 ) THEN
         outopt(2) = 1
      ENDIF

!     check scaling scale matrices

      IF(xscalm .eq. 0) THEN
         DO i = 1,n
            IF(scalex(i) .lt. Rzero) scalex(i) = -scalex(i)
            IF(scalex(i) .eq. Rzero) scalex(i) = Rone
         ENDDO
      ELSE
         xscalm = 1
         DO i = 1,n
            scalex(i) = Rone
         ENDDO
      ENDIF
!     check step and function tolerances

      IF(xtol .lt. Rzero) THEN
         xtol = epsm**(Rtwo/Rthree)
      ENDIF

      IF(ftol .lt. Rzero) THEN
         ftol = epsm**(Rtwo/Rthree)
      ENDIF

      IF( btol .lt. xtol ) btol = xtol

      cndtol = max(cndtol, epsm)

!     check reduction in geometric linesearch

      IF( sigma .le. Rzero .or. sigma .ge. Rone ) THEN
         sigma = Rhalf
      ENDIF

!     check iteration limit

      IF(maxit .le. 0) THEN
         maxit = 150
      ENDIF

!     set stepmx

      IF(stepmx .le. Rzero) stepmx = Rhuge

!     check delta
      IF(delta .le. Rzero) THEN
         IF( delta .ne. -Rtwo ) THEN
            delta = -Rone
         ENDIF
      ELSEIF(delta .gt. stepmx) THEN
         delta = stepmx
      ENDIF

      RETURN
      END SUBROUTINE

END MODULE nwnleq_mod
