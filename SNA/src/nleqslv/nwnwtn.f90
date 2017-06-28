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
MODULE nwnwtn_mod

USE Kind_Types_Mod
USE nwutil_mod
USE nwbjac_mod
USE nwmhlm_mod
USE nwpdlg_mod
USE nwddlg_mod
USE nwglsh_mod
USE nwqlsh_mod
USE nwclsh_mod
USE nwpure_mod
USE liqrup_mod
USE nwnjac_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE nwsolv(ldr,xc,n,scalex,maxit, &
                        jacflg,xtol,ftol,btol,cndtol,global,xscalm, &
                        stepmx,delta,sigma, &
                        rjac,wrk1,wrk2,wrk3,wrk4,fc,fq,dn,d,qtf, &
                        rcdwrk,icdwrk,qrwork,qrwsiz,epsm, &
                        fjac,fvec,outopt,xp,fp,gp,njcnt,nfcnt,iter, &
                        termcd)

      INTEGER(I4B) :: ldr,n,termcd,njcnt,nfcnt,iter
      INTEGER(I4B) :: maxit,jacflg(*),global,xscalm,qrwsiz
      INTEGER(I4B) :: outopt(*)
      REAL(DP) :: xtol,ftol,btol,cndtol
      REAL(DP) :: stepmx,delta,sigma,fpnorm,epsm
      REAL(DP) :: rjac(ldr,*)
      REAL(DP) :: xc(*),fc(*),xp(*),fp(*),dn(*),d(*)
      REAL(DP) :: wrk1(*),wrk2(*),wrk3(*),wrk4(*)
      REAL(DP) :: qtf(*),gp(*),fq(*)
      REAL(DP) :: scalex(*)
      REAL(DP) :: rcdwrk(*),qrwork(*)
      INTEGER(I4B) :: icdwrk(*)
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

!-----------------------------------------------------------------------
!
!     Solve system of nonlinear equations with Newton and global strategy
!
!
!     Arguments
!
!     In       ldr     Integer         leading dimension of rjac
!     In       xc      Real(*)         initial estimate of solution
!     In       n       Integer         dimensions of problem
!     Inout    scalex  Real(*)         scaling factors x(*)
!     In       maxit   Integer         maximum number of allowable iterations
!     In       jacflg  Integer(*)      jacobian flag array
!                                      jacflg[1]:  0 numeric; 1 USEr supplied; 2 numerical banded
!                                                  3: USEr supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting step allowed when
!                                                   singular or illconditioned
!     In       xtol    Real            tolerance at which successive iterates x()
!                                      are considered close enough to
!                                      terminate algorithm
!     In       ftol    Real            tolerance at which function values f()
!                                      are considered close enough to zero
!     Inout    btol    Real            x tolerance for backtracking
!     Inout    cndtol  Real            tolerance of test for ill conditioning
!     In       global  Integer         global strategy to USE
!                                        1 cubic linesearch
!                                        2 quadratic linesearch
!                                        3 geometric linesearch
!                                        4 double dogleg
!                                        5 powell dogleg
!                                        6 hookstep (More-Hebden Levenberg-Marquardt)
!     In       xscalm  Integer         x scaling method
!                                        1 from column norms of first jacobian
!                                          increased if needed after first iteration
!                                        0 scaling USEr supplied
!     In       stepmx  Real            maximum allowable step size
!     In       delta   Real            trust region radius
!     In       sigma   Real            reduction factor geometric linesearch
!     Inout    rjac    Real(ldr,*)     jacobian (n columns)
!     Wk       wrk1    Real(*)         workspace
!     Wk       wrk2    Real(*)         workspace
!     Wk       wrk3    Real(*)         workspace
!     Wk       wrk4    Real(*)         workspace
!     Inout    fc      Real(*)         function values f(xc)
!     Wk       fq      Real(*)         workspace
!     Wk       dn      Real(*)         workspace
!     Wk       d       Real(*)         workspace
!     Wk       qtf     Real(*)         workspace
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!     In       qrwork  Real(*)         workspace for Lapack QR routines (CALL liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!     In       epsm    Real            machine precision
!     In       fjac    Name            name of routine to calculate jacobian
!                                      (optional)
!     In       fvec    Name            name of routine to calculate f()
!     In       outopt  Integer(*)      output options
!     Out      xp      Real(*)         final x()
!     Out      fp      Real(*)         final f(xp)
!     Out      gp      Real(*)         gradient at xp()
!     Out      njcnt   Integer         number of jacobian evaluations
!     Out      nfcnt   Integer         number of function evaluations
!     Out      iter    Integer         number of (outer) iterations
!     Out      termcd  Integer         termination code
!
!-----------------------------------------------------------------------

      INTEGER(I4B) :: gcnt,retcd,ierr
      REAL(DP) ::  dum(2),fcnorm,rcond
      LOGICAL fstjac
      INTEGER(I4B) :: priter

      INTEGER(I4B) :: idamax

!     initialization

      retcd = 0
      iter  = 0
      njcnt = 0
      nfcnt = 0
      ierr  = 0

      dum(1) = 0

      IF( outopt(1) .eq. 1 ) THEN
         priter = 1
      ELSE
         priter = -1
      ENDIF

!     evaluate function

      CALL vscal(n,xc,scalex)
      CALL nwfvec(xc,n,scalex,fvec,fc,fcnorm,wrk1)

!     evaluate USEr supplied or finite difference jacobian and check USEr supplied
!     jacobian, if requested

      fstjac = .false.
      IF(mod(jacflg(1),2) .eq. 1) THEN

        IF( outopt(2) .eq. 1 ) THEN
           fstjac = .true.
           njcnt = njcnt + 1
           CALL nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,fvec,fjac,rjac, &
                       ldr,wrk1,wrk2,wrk3)
           CALL chkjac(rjac,ldr,xc,fc,n,epsm,jacflg,scalex, &
                       fq,wrk1,wrk2,fvec,termcd)
           IF(termcd .lt. 0) THEN
!              copy initial values
               CALL dcopy(n,xc,1,xp,1)
               CALL dcopy(n,fc,1,fp,1)
               CALL vunsc(n,xp,scalex)
               fpnorm = fcnorm
               RETURN
           ENDIF
        ENDIF

      ENDIF

!     check stopping criteria for input xc

      CALL nwtcvg(xc,fc,xc,xtol,retcd,ftol,iter,maxit,n,ierr,termcd)

      IF(termcd .gt. 0) THEN
          CALL dcopy(n,xc,1,xp,1)
          CALL dcopy(n,fc,1,fp,1)
          fpnorm = fcnorm
          IF( outopt(3) .eq. 1 .and. .not. fstjac ) THEN
             njcnt = njcnt + 1
             CALL nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,fvec,fjac,rjac, &
                         ldr,wrk1,wrk2,wrk3)
          ENDIF
          RETURN
      ENDIF

!       IF( priter .gt. 0 ) THEN

!          dum(1) = fcnorm
!          dum(2) = abs(fc(idamax(n,fc,1)))

!          IF( global .eq. 0 ) THEN
!             CALL nwprot(iter, -1, dum)
!          ELSEIF( global .le. 3 ) THEN
!             CALL nwlsot(iter,-1,dum)
!          ELSEIF( global .eq. 4 ) THEN
!             CALL nwdgot(iter,-1,0,dum)
!          ELSEIF( global .eq. 5 ) THEN
!             CALL nwpwot(iter,-1,0,dum)
!          ELSEIF( global .eq. 6 ) THEN
!             CALL nwmhot(iter,-1,0,dum)
!          ENDIF

!       ENDIF
      DO WHILE( termcd .eq. 0 )
         iter = iter + 1
         CALL nwnjac(rjac,ldr,n,xc,fc,fq,fvec,fjac,epsm,jacflg,wrk1, &
                     wrk2,wrk3, &
                     xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn, &
                     qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)
!        - choose the next iterate xp by a global strategy
         IF( ierr .gt. 0 ) THEN
!           jacobian singular or too ill-conditioned
            CALL nweset(n,xc,fc,fcnorm,xp,fp,fpnorm,gcnt,priter,iter)
         ELSEIF(global .eq. 0) THEN
            CALL nwpure(n,xc,dn,stepmx,scalex, &
                        fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
                        priter,iter)
         ELSEIF(global .eq. 1) THEN
            CALL nwclsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex, &
                        fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
                        priter,iter)
         ELSEIF(global .eq. 2) THEN
            CALL nwqlsh(n,xc,fcnorm,dn,gp,stepmx,btol,scalex, &
                        fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
                        priter,iter)
         ELSEIF(global .eq. 3) THEN
            CALL nwglsh(n,xc,fcnorm,dn,gp,sigma,stepmx,btol,scalex, &
                        fvec,xp,fp,fpnorm,wrk1,retcd,gcnt, &
                        priter,iter)
         ELSEIF(global .eq. 4) THEN
            CALL nwddlg(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx, &
                        btol,delta,qtf,scalex, &
                        fvec,d,fq,wrk1,wrk2,wrk3,wrk4, &
                        xp,fp,fpnorm,retcd,gcnt,priter,iter)
         ELSEIF(global .eq. 5) THEN
            CALL nwpdlg(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx, &
                        btol,delta,qtf,scalex, &
                        fvec,d,fq,wrk1,wrk2,wrk3,wrk4, &
                        xp,fp,fpnorm,retcd,gcnt,priter,iter)
         ELSEIF(global .eq. 6) THEN
            CALL nwmhlm(n,rjac,ldr,dn,gp,xc,fcnorm,stepmx, &
                        btol,delta,qtf,scalex, &
                        fvec,d,fq,wrk1,wrk2,wrk3,wrk4, &
                        xp,fp,fpnorm,retcd,gcnt,priter,iter)
         ENDIF

         nfcnt = nfcnt + gcnt

!        - check stopping criteria for the new iterate xp

         CALL nwtcvg(xp,fp,xc,xtol,retcd,ftol,iter,maxit,n,ierr,termcd)

         IF(termcd .eq. 0) THEN
!           update xc, fc, and fcnorm
            CALL dcopy(n,xp,1,xc,1)
            CALL dcopy(n,fp,1,fc,1)
            fcnorm = fpnorm
         ENDIF

      ENDDO

      IF( outopt(3) .eq. 1 ) THEN
         CALL nwfjac(xp,scalex,fp,fq,n,epsm,jacflg,fvec,fjac,rjac, &
                     ldr,wrk1,wrk2,wrk3)
      ENDIF

      CALL vunsc(n,xp,scalex)

      RETURN
      END SUBROUTINE

END MODULE nwnwtn_mod
