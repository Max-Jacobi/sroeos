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
MODULE nwbjac_mod

USE Kind_Types_Mod
USE nwutil_mod
USE lautil_mod
USE liqrev_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE nwbjac(rjac,r,ldr,n,xc,fc,fq,fvec,fjac,epsm,jacflg, &
                        wrk1,wrk2,wrk3, &
                        xscalm,scalex,gp,cndtol,rcdwrk,icdwrk,dn, &
                        qtf,rcond,qrwork,qrwsiz,njcnt,iter,fstjac,ierr)

!-----------------------------------------------------------------------
!
!     Compute Jacobian matrix in xc, fc
!     scale it, compute gradient in xc and generate QR decomposition
!     calculate Newton step
!
!     Arguments
!
!     Out      rjac    Real(ldr,*)     jacobian (n columns) and USEd for storing full Q from Q
!     Out      r       Real(ldr,*)     USEd for storing R from QR factorization
!     In       ldr     Integer         leading dimension of rjac
!     In       n       Integer         dimensions of problem
!     In       xc      Real(*)         initial estimate of solution
!     Inout    fc      Real(*)         function values f(xc)
!     Wk       fq      Real(*)         workspace
!     In       fjac    Name            name of routine to calculate jacobian
!                                      (optional)
!     In       fvec    Name            name of routine to calculate f()
!     In       epsm    Real            machine precision
!     In       jacflg  Integer(*)      jacobian flag array
!                                      jacflg[1]:  0 numeric; 1 USEr supplied; 2 numerical banded
!                                                  3: USEr supplied banded
!                                      jacflg[2]: number of sub diagonals or -1 if not banded
!                                      jacflg[3]: number of super diagonals or -1 if not banded
!                                      jacflg[4]: 1 if adjusting step allowed when
!                                                   singular or illconditioned
!     Wk       wrk1    Real(*)         workspace
!     Wk       wrk2    Real(*)         workspace
!     Wk       wrk3    Real(*)         workspace
!     In       xscalm  Integer         x scaling method
!                                        1 from column norms of first jacobian
!                                          increased if needed after first iteration
!                                        0 scaling USEr supplied
!     Inout    scalex  Real(*)         scaling factors x(*)
!     Out      gp      Real(*)         gradient at xp()
!     In       cndtol  Real            tolerance of test for ill conditioning
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!     Out      dn      Real(*)         Newton step
!     Out      qtf     Real(*)         workspace for nwnstp
!     Out      rcond   Real            estimated inverse condition of R from QR
!     In       qrwork  Real(*)         workspace for Lapack QR routines (CALL liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!     Out      njcnt   Integer         number of jacobian evaluations
!     In       iter    Integer         iteration counter (USEd in scaling)
!     Inout    fstjac  LOGICAL         .true. if initial jacobian is available
!                                      on exit set to .false.
!     Out      ierr    Integer         error code
!                                        0 no error
!                                       >0 error in nwnstp (singular ...)
!
!-----------------------------------------------------------------------

      INTEGER(I4B) :: ldr,n,iter, njcnt, ierr
      INTEGER(I4B) :: jacflg(*),xscalm,qrwsiz
      LOGICAL fstjac
      REAL(DP) :: epsm, cndtol, rcond
      REAL(DP) :: rjac(ldr,*),r(ldr,*)
      REAL(DP) :: xc(*),fc(*),dn(*)
      REAL(DP) :: wrk1(*),wrk2(*),wrk3(*)
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

      LOGICAL stepadj
      REAL(DP), PARAMETER :: Rzero = 0.0d0, Rone = 1.0d0

!     evaluate the jacobian at the current iterate xc

      IF( .not. fstjac ) THEN
         CALL nwfjac(xc,scalex,fc,fq,n,epsm,jacflg,fvec,fjac,rjac, &
                     ldr,wrk1,wrk2,wrk3)
         njcnt = njcnt + 1
      ELSE
         fstjac = .false.
      ENDIF

!     if requested calculate x scale from jacobian column norms a la Minpack

      IF( xscalm .eq. 1 ) THEN
         CALL vunsc(n,xc,scalex)
         CALL nwcpsx(n,rjac,ldr,scalex,epsm,iter)
         CALL vscal(n,xc,scalex)
      ENDIF

      CALL nwscjac(n,rjac,ldr,scalex)

!     evaluate the gradient at the current iterate xc
!     gp = trans(Rjac) * fc
      CALL dgemv('T',n,n,Rone,rjac,ldr,fc,1,Rzero,gp,1)

!     get broyden (newton) step
      stepadj = jacflg(4) .eq. 1
      CALL dcopy(n,fc,1,fq,1)
      CALL brdstp(rjac,r,ldr,fq,n,cndtol, stepadj, &
                  wrk1,dn,qtf,ierr,rcond, &
                  rcdwrk,icdwrk,qrwork,qrwsiz)

!     save some data about jacobian for later output
!      CALL nwsnot(0,ierr,rcond)

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------

      SUBROUTINE brdstp(rjac,r,ldr,fn,n,cndtol, stepadj, &
                        qraux,dn,qtf,ierr,rcond, &
                        rcdwrk,icdwrk,qrwork,qrwsiz)

      INTEGER(I4B) :: ldr,n,ierr,qrwsiz
      REAL(DP) ::  cndtol,rjac(ldr,*),r(ldr,*),qraux(*),fn(*)
      REAL(DP) ::dn(*),qtf(*)
      REAL(DP) :: rcdwrk(*),qrwork(*)
      INTEGER(I4B) :: icdwrk(*)
      REAL(DP) :: rcond
      LOGICAL stepadj

!-----------------------------------------------------------------------
!
!     Calculate the newton step
!
!     Arguments
!
!     Inout    rjac    Real(ldr,*)     jacobian matrix at current iterate; on RETURN full Q
!     Inout    r       Real(ldr,*)     jacobian matrix at current iterate; on RETURN R fom QR
!     In       ldr     Integer         leading dimension of rjac
!     In       fn      Real(*)         function values at current iterate
!     In       n       Integer         dimension of problem
!     In       cndtol  Real            tolerance of test for ill conditioning
!     In       stepadj Logical         allow adjusting step for singular/illconditioned jacobian
!     Inout    qraux   Real(*)         QR info from liqrfa (CALLing Lapack dgeqrf)
!     Out      dn      Real(*)         Newton direction
!     Out      qtf     Real(*)         trans(Q)*f()
!     Out      ierr    Integer         0 indicating Jacobian not ill-conditioned or singular
!                                      1 indicating Jacobian ill-conditioned
!                                      2 indicating Jacobian completely singular
!                                      3 indicating almost zero LM correction
!     Out      rcond   Real            inverse condition of upper triangular R of QR
!     Wk       rcdwrk  Real(*)         workspace
!     Wk       icdwrk  Integer(*)      workspace
!     In       qrwork  Real(*)         workspace for Lapack QR routines (CALL liqsiz)
!     In       qrwsiz  Integer         size of qrwork
!
!-----------------------------------------------------------------------

      INTEGER(I4B) info,k

      REAL(DP), PARAMETER :: Rone=1.0d0
      REAL(DP) :: mu

!     perform a QR factorization of rjac (simple Lapack routine)
!     check for singularity or ill conditioning
!     form qtf = trans(Q) * fn

      CALL liqrfa(rjac,ldr,n,qraux,qrwork,qrwsiz,ierr)

!     check for singularity or ill conditioning

      CALL cndjac(n,rjac,ldr,cndtol,rcond,rcdwrk,icdwrk,ierr)

!     compute qtf = trans(Q)*fn

      CALL dcopy(n,fn,1,qtf,1)
      CALL liqrqt(rjac, ldr, n, qraux, qtf, qrwork, qrwsiz, info)

!     copy the upper triangular part of the QR decomposition
!     contained in Rjac into R[*, 1..n].
!     form Q from the QR decomposition (taur/qraux in wrk1)

      CALL dlacpy('U',n,n,rjac,ldr,r,ldr)
      CALL liqrqq(rjac,ldr,qraux,n,qrwork,qrwsiz,info)

!     now Rjac[* ,1..n] holds expanded Q
!     now R[* ,1..n] holds full upper triangle R

      IF( ierr .eq. 0 ) THEN
!         Normal Newton step
!         solve Jacobian*dn  =  -fn
!         ==> R*dn = - qtf

          CALL dcopy(n,qtf,1,dn,1)
          CALL dtrsv('U','N','N',n,r,ldr,dn,1)
          CALL dscal(n, -Rone, dn, 1)

      ELSEIF( stepadj ) THEN
!         Adjusted Newton step
!         approximately from pseudoinverse(Jac+)
!         USE mu to solve (trans(R)*R + mu*I*mu*I) * x = - trans(R) * fn
!         directly from the QR decomposition of R stacked with mu*I
!         a la Levenberg-Marquardt
          CALL compmu(r,ldr,n,mu,rcdwrk,ierr)
          IF( ierr .eq. 0 ) THEN
             CALL liqrev(n,r,ldr,mu,qtf,dn, &
                         rcdwrk(1+n),rcdwrk(2*n+1))
             CALL dscal(n, -Rone, dn, 1)

!            copy lower triangular R to upper triangular
             DO k=1,n
                CALL dcopy (n-k+1,r(k,k),1,r(k,k),ldr)
                r(k,k) = rcdwrk(1+n+k-1)
             ENDDO
          ENDIF
      ENDIF

      RETURN
      END SUBROUTINE

END MODULE nwbjac_mod
