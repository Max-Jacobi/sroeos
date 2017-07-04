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
MODULE Surface_Equilibrium_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, TEN
  USE Skyrme_Bulk_Mod, ONLY : SKYRME_BULK_PROPERTIES
  USE nwnleq_mod
  USE lautil_mod

  IMPLICIT NONE

  REAL(DP) :: prot_frac_inside, temp

CONTAINS

  SUBROUTINE FIND_SURFACE_EQUILIBRIUM ( x_sol, prot_frac,  &
                                        temperature, residue, lgt )
!   x_sol(1) = inside neutron density
!   x_sol(2) = inside proton density
!   x_sol(3) = outside neutron density
!   x_sol(4) = outside proton density
!   lgt = .true. IF inside and outside proton densities different enough
    IMPLICIT NONE

    REAL(DP), DIMENSION(4), INTENT(INOUT) :: x_sol, residue
    REAL(DP), INTENT(IN) :: prot_frac, temperature
    REAL(DP), DIMENSION(4) :: x1, x2, r
    REAL(DP), DIMENSION(4,4) :: g
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

    prot_frac_inside = prot_frac
    temp = temperature

    retry = 0

!   set parameters for non-linear equation solver
    N = 4 ; MAXIT = 250 ; JACFLG(1:4) = (/0,-1,-1, 1/) ; OUTOPT(1:2) = (/0,1/)
    METHOD = 0; GLOBAL = 4; XSCALM = 1; LDR = N; LRWORK = 9*N
    XTOL = 1.D-12; FTOL = 1.D-8; BTOL = 1.D-4; CNDTOL = 1.D-7
    STEPMX = -1.D0; DELTA = -1.D0; SIGMA = 0.5D0; TRACE = 1.D0; DSUB = 0.D0
    SIZE = (METHOD+1)*N*N

10  CONTINUE

    ALLOCATE(RJAC(SIZE),RWORK(9*N),RCDWRK(3*N),ICDWRK(N),SCALEX(N))
    CALL liqsiz(n,qrwsiz)
    ALLOCATE(qrwork(qrwsiz))

    RJAC = ZERO ; RWORK = ZERO; RCDWRK = ZERO ; QRWORK = ZERO ; ICDWRK = 0
    SCALEX = ONE

!   Set initial guesses for
!    inside neutron fraction
    x2(1) = - 1.1d0
!   inside proton fraction
    x2(2) = x2(1) + log10(prot_frac)
!   outside neutron fraction
    x2(3) = -10.d0
!   outside proton fraction
    x2(4) = x2(3)*(two*prot_frac)

!   if initial value not zero
!    then use them as initial guesses
    if (sum(x_sol(1:2))/=zero) x2(1:2) = x_sol(1:2)
    if (sum(x_sol(3:4))/=zero) x2(3:4) = x_sol(3:4)

    x1 = ZERO

!   solve equation for equilibrium inside AND outside nuclei
    CALL nwnleq(x2,n,scalex,maxit,jacflg,xtol,ftol,btol,cndtol,method,global, &
                xscalm,stepmx,delta,sigma,rjac,ldr,rwork,lrwork, &
                rcdwrk,icdwrk,qrwork,qrwsiz,&
                jacob_surface_equilibrium,surface_equilibrium,outopt,x1, &
                r,g,njcnt,nfcnt,iter,termcd)
    deallocate(RJAC, RWORK, RCDWRK, ICDWRK, QRWORK, SCALEX)
!   chech whether solution x1 found is actually a solution to eq being solved.
!   sotemimes output for x1 is not a solution, but a point where nlwleq stalled.
    flag = 0
    CALL surface_equilibrium( x1, r, n, flag )
    residue = r

    lgt = .FALSE.
!   IF x1 is a solution THEN set as so
    IF (DOT_PRODUCT(r,r)<1.d-12) THEN
      x_sol = x1
      ! only set as true solution with inside/outside matter
      ! if inside/outside densities are different
      RES1 = ABS((x_sol(1)-x_sol(3))/x_sol(1))
      RES2 = ABS((x_sol(2)-x_sol(4))/x_sol(2))
      if (res1>1.d-6 .and. res2>1.d-6) lgt = .true.
      if (lgt) return
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

  SUBROUTINE surface_equilibrium ( x, s, m, surface_flag )
!
    IMPLICIT NONE

    INTEGER (I4B), INTENT(IN) :: m, surface_flag
    REAL(DP), INTENT(IN)  :: x(m)
    REAL(DP), INTENT(OUT) :: s(m)
    ! inside region properties
    REAL(DP) :: log10_dens_neut_in,log10_dens_prot_in
    REAL(DP) :: dens_neut_in,dens_prot_in,dens_nucl_in,prot_frac_in
    REAL(DP) :: effective_mass_neut_in,effective_mass_prot_in
    REAL(DP) :: eta_neut_in,eta_prot_in
    REAL(DP) :: tau_neut_in,tau_prot_in
    REAL(DP) :: v_neut_in,v_prot_in
    REAL(DP) :: mu_neut_in,mu_prot_in
    REAL(DP) :: pressure_in, free_energy_in
    ! outside region properties
    REAL(DP) :: log10_dens_neut_out,log10_dens_prot_out
    REAL(DP) :: dens_neut_out,dens_prot_out,dens_nucl_out,prot_frac_out
    REAL(DP) :: effective_mass_neut_out,effective_mass_prot_out
    REAL(DP) :: eta_neut_out,eta_prot_out
    REAL(DP) :: tau_neut_out,tau_prot_out
    REAL(DP) :: v_neut_out,v_prot_out
    REAL(DP) :: mu_neut_out,mu_prot_out
    REAL(DP) :: pressure_out, free_energy_out

    log10_dens_neut_in = X(1)
    log10_dens_prot_in = X(2)

    ! pressure_in = SKYRME_BULK_PRESSURE(log10_dens_neut_in,&
    !                                    log10_dens_prot_in,Temp)
    CALL SKYRME_BULK_PROPERTIES(log10_dens_neut_in,log10_dens_prot_in, &
    Temp,dens_neut_in,dens_prot_in,dens_nucl_in,prot_frac_in,   &
    effective_mass_neut_in,effective_mass_prot_in,   &
    eta_neut_in,eta_prot_in,tau_neut_in,tau_prot_in, &
    v_neut_in,v_prot_in,mu_neut_in,mu_prot_in,       &
    pressure_in,free_energy_in)

    log10_dens_neut_out = X(3)
    log10_dens_prot_out = X(4)

    dens_neut_out = TEN**log10_dens_neut_out
    dens_prot_out = TEN**log10_dens_prot_out

    ! pressure_out = SKYRME_BULK_PRESSURE(log10_dens_neut_out,&
    !                                     log10_dens_prot_out,Temp)
    CALL SKYRME_BULK_PROPERTIES(log10_dens_neut_out,log10_dens_prot_out, &
    Temp,dens_neut_out,dens_prot_out,dens_nucl_out,prot_frac_out, &
    effective_mass_neut_out,effective_mass_prot_out,     &
    eta_neut_out,eta_prot_out,tau_neut_out,tau_prot_out, &
    v_neut_out,v_prot_out,mu_neut_out,mu_prot_out,       &
    pressure_out,free_energy_out)

    S(1) = (pressure_in  - pressure_out)*1.d1
    S(2) = (mu_neut_in   - mu_neut_out)
    S(3) = (mu_prot_in   - mu_prot_out)
    S(4) = (prot_frac_in - prot_frac_inside)

    RETURN
  END SUBROUTINE surface_equilibrium

  SUBROUTINE jacob_surface_equilibrium ( jac, ldr, x1, m )
!
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: m, ldr
    REAL(DP), INTENT(IN)  :: x1(m)
    REAL(DP), INTENT(OUT) :: jac(ldr,*)

    jac(1:ldr,1:m) = zero

    RETURN
  END SUBROUTINE jacob_surface_equilibrium

  END SUBROUTINE FIND_SURFACE_EQUILIBRIUM

END MODULE Surface_Equilibrium_Mod
