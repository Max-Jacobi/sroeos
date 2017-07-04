!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it AND/or modIFy
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
!    along with SRO_EOS.  IF not, see <http://www.gnu.org/licenses/>.
!
MODULE Critical_Point_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO
  USE Nuclear_Matter_Properties_Mod, ONLY : Nuc_Sat_Dens
  USE Skyrme_Bulk_Mod, ONLY : SKYRME_BULK_PRESSURE
  USE nwnleq_mod
  USE lautil_mod

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE FIND_CRITICAL_POINT ( x_sol, prot_frac, lgt )
!   x_sol(1) = critical point density
!   x_sol(2) = critical point temperature
!   lgt = .true. IF solution is actually a critical point
    IMPLICIT NONE

    REAL(DP), DIMENSION(2), INTENT(INOUT) :: x_sol
    REAL(DP) :: prot_frac
    REAL(DP), DIMENSION(2) :: x1, x2, r
    REAL(DP), DIMENSION(2,2) :: g
    LOGICAL, INTENT(OUT) :: lgt
    INTEGER(I4B) :: FLAG, retry
    REAL(DP) :: RES1, RES2

!   parameters for non-linear equation solver "nleqslv"
    INTEGER(I4B) :: N, MAXIT, JACFLG(1:4),  SIZE
    INTEGER(I4B) :: METHOD, GLOBAL, XSCALM, LDR, LRWORK
    INTEGER(I4B) :: NJCNT, NFCNT, ITER, TERMCD, QRWSIZ, OUTOPT(2)
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ICDWRK
    REAL(DP) :: XTOL, FTOL, BTOL, CNDTOL
    REAL(DP) :: STEPMX, DELTA, SIGMA, TRACE, DSUB
    REAL(DP), DIMENSION(:), ALLOCATABLE :: RJAC, RWORK, RCDWRK
    REAL(DP), DIMENSION(:), ALLOCATABLE :: QRWORK, SCALEX

    retry = 0

!   set parameters for non-linear equation solver
    N = 2 ; MAXIT = 250 ; JACFLG(1:4) = (/0,-1,-1, 1/) ; OUTOPT(1:2) = (/0,1/)
    METHOD = 0; GLOBAL = 4; XSCALM = 1; LDR = N; LRWORK = 9*N
    XTOL = 1.D-15; FTOL = 1.D-12; BTOL = 1.D-4; CNDTOL = 1.D-7
    STEPMX = -1.D0; DELTA = -1.D0; SIGMA = 0.5D0; TRACE = 1.D0; DSUB = 0.D0
    SIZE = (METHOD+1)*N*N

10  CONTINUE

    ALLOCATE(RJAC(SIZE),RWORK(9*N),RCDWRK(3*N),ICDWRK(N),SCALEX(N))
    CALL liqsiz(n,qrwsiz)
    ALLOCATE(qrwork(qrwsiz))

    RJAC = ZERO ; RWORK = ZERO; RCDWRK = ZERO ; QRWORK = ZERO ; ICDWRK = 0
    SCALEX = ONE

!   Set initial guesses for
!   critical density
    x2(1) = 0.5d0
!   critical temperature
    x2(2) = 10.d0

    x1 = x2

!   solve equation for equilibrium inside AND outside nuclei
    CALL nwnleq(x2,n,scalex,maxit,jacflg,xtol,ftol,btol,cndtol,method,global, &
                xscalm,stepmx,delta,sigma,rjac,ldr,rwork,lrwork, &
                rcdwrk,icdwrk,qrwork,qrwsiz,&
                jacob_critical_point,critical_point,outopt,x1, &
                r,g,njcnt,nfcnt,iter,termcd)
    deallocate(RJAC, RWORK, RCDWRK, ICDWRK, QRWORK, SCALEX)
!   chech whether solution x1 found is actually a solution to eq being solved.
!   sotemimes output for x1 is not a solution, but a point where nlwleq stalled.
    flag = 0
    CALL critical_point( x1, r, n, flag )

    lgt = .FALSE.
!   IF x1 is a solution THEN set as so
    IF (DOT_PRODUCT(r,r)<1.d-12) THEN
      x_sol = x1
      IF (x_sol(1)>zero .AND. x_sol(1)<one .AND. x_sol(2) > zero) lgt = .TRUE.
      IF (lgt) RETURN
    ENDIF

!   sometimes solution is not found with one method of nleqslv
!   but a different method can find the solution
    IF (retry<3) THEN
      retry = retry+1
!        IF (retry==1) nlglobal = 0
!        IF (retry==2) nlglobal = 2
!        IF (retry==2) nlglobal = 4
      IF (retry==1) global = 5
      IF (retry==2) global = 6
      GOTO 10
    ENDIF

    RETURN

CONTAINS

  SUBROUTINE critical_point ( x, s, m, flag )
!
!    Given the value for n_ni, n_pi, n_no, n_po,
!    this routine solves the equilibrium equations
!    where pressure and chemical potentials are
!    the same inside and outside nuclei.
!
    IMPLICIT NONE

    INTEGER (I4B), INTENT(IN) :: m, flag
    REAL(DP), INTENT(IN)  :: x(m)
    REAL(DP), INTENT(OUT) :: s(m)
    REAL(DP) :: XI, DNNO, DNPO
    REAL(DP) :: P_up, P_0, P_dn, EPS, XRHO, XR
    REAL(DP) :: log10_dens_neut, log10_dens_prot, Temperature

    EPS  = 1.D-2

!   tempearture of the system
    TEMPERATURE = X(2)

!   high density phase proton fraction
    XI   = prot_frac

!   low density phase density, neutron density and proton density
    XRHO = X(1)*Nuc_Sat_Dens
    XR   = XRHO
    DNNO = (ONE-XI)*XRHO
    DNPO = XI*XRHO

    log10_dens_neut = log10(DNNO)
    log10_dens_prot = log10(DNPO)

    !   calculate nuclear pressure
    P_0 = SKYRME_BULK_PRESSURE(log10_dens_neut, log10_dens_prot, Temperature)

!   low density phase density, neutron density and proton density
    XRHO = XR*(ONE-EPS)
    DNNO = (ONE-XI)*XRHO
    DNPO = XI*XRHO

    log10_dens_neut = log10(DNNO)
    log10_dens_prot = log10(DNPO)

!   calculate nuclear pressure
    P_dn = SKYRME_BULK_PRESSURE(log10_dens_neut, log10_dens_prot, Temperature)

    XRHO = XR*(ONE+EPS)
    DNNO = (ONE-XI)*XRHO
    DNPO = XI*XRHO

    log10_dens_neut = log10(DNNO)
    log10_dens_prot = log10(DNPO)

!   calculate chemical potentials of neutrons and protons and their pressure
    P_up = SKYRME_BULK_PRESSURE(log10_dens_neut, log10_dens_prot, Temperature)

    S(1) = (P_up-P_dn)/(TWO*XR*EPS)
    S(2) = (P_up-TWO*P_0+P_dn)/(XR*EPS)**TWO

    RETURN
  END SUBROUTINE critical_point

  SUBROUTINE jacob_critical_point ( jac, ldr, x1, m )
!
!   jacobian set to zero
!   since we obtain the jacobian numerically
!
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: m, ldr
    REAL(DP), INTENT(IN)  :: x1(m)
    REAL(DP), INTENT(OUT) :: jac(ldr,*)

    jac(1:ldr,1:m) = zero

    RETURN
  END SUBROUTINE jacob_critical_point

  END SUBROUTINE FIND_CRITICAL_POINT

END MODULE Critical_Point_Mod
