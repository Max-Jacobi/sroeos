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
MODULE Find_Uniform_Solution_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL
  USE Global_Variables_Mod, ONLY : IS_TEST
  USE Physical_Constants_Mod, ONLY : ZERO, HALF, ONE, TWO, FOUR, TEN, &
                             v_alpha
  USE Free_Energy_Mod
  USE Skyrme_Bulk_Density_Derivatives_Mod
  USE nwnleq_mod
  USE lautil_mod
  USE unused_mod
  USE phase_space_point_mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Find_Uniform_Solution( n_in, T_in, Yp_in, x_guess, x_sol, &
                                    Ftotal, residue, retry, x_lgcl)

    IMPLICIT NONE
!
    REAL(DP), INTENT(IN) :: n_in, T_in, Yp_in, x_guess(1)
    REAL(DP), INTENT(OUT) :: x_sol(1), Ftotal, residue
    LOGICAL(LGCL), INTENT(OUT) :: x_lgcl
    LOGICAL(LGCL), INTENT(IN)  :: retry
    REAL(DP), DIMENSION(1,1) :: g
    REAL(DP) :: x_grid(1), x1(1), x2(1), r(1), X_p
    REAL(DP) :: log10_n_no, log10_n_po, log10_u
    REAL(DP) :: dLog10_n_no_dn, dLog10_n_po_dn, dLog10_u_dn
    REAL(DP) :: dLog10_n_no_dT, dLog10_n_po_dT, dLog10_u_dT
    REAL(DP) :: F_out, F_in, F_alpha, F_heavy, F_surf_coul
    REAL(DP) :: n_no, n_po
    REAL(DP) :: n_ni, n_pi
    REAL(DP) :: n_alpha, n_heavy
    REAL(DP) :: non_uniform_eq(3), non_uniform_jac(3,3)
    INTEGER(I4B) :: error

!   parameters for non-linear equation solver "nleqslv"
    INTEGER(I4B) :: N_eq, MAXIT, JACFLG(1:4),  SIZE
    INTEGER(I4B) :: METHOD, GLOBAL, XSCALM, LDR, LRWORK
    INTEGER(I4B) :: NJCNT, NFCNT, ITER, TERMCD, QRWSIZ, OUTOPT(2)
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ICDWRK
    REAL(DP) :: XTOL, FTOL, BTOL, CNDTOL
    REAL(DP) :: STEPMX, DELTA, SIGMA, TRACE, DSUB
    REAL(DP), DIMENSION(:), ALLOCATABLE :: RJAC, RWORK, RCDWRK
    REAL(DP), DIMENSION(:), ALLOCATABLE :: QRWORK, SCALEX

    n = n_in
    T = T_in
    Yp = Yp_in
!
!   set parameters for non-linear equation solver
    N_eq = 1; MAXIT = 500 ; JACFLG(1:4) = (/0,-1,-1, 1/) ; OUTOPT(1:2) = (/0,1/)
    METHOD = 0; GLOBAL = 4; XSCALM = 1; LDR = N_eq; LRWORK = 9*N_eq
    XTOL = 1.D-18; FTOL = 1.D-8; BTOL = 1.D-6; CNDTOL = 1.D-9
    STEPMX = -1.D0; DELTA = -1.D0; SIGMA = 0.5D0; TRACE = 1.D0; DSUB = 0.D0
    SIZE = (METHOD+1)*N_eq*N_eq

10  CONTINUE

    ALLOCATE(RJAC(SIZE),RWORK(9*N_eq),RCDWRK(3*N_eq),ICDWRK(N_eq),SCALEX(N_eq))
    CALL liqsiz(N_eq,qrwsiz)
    ALLOCATE(qrwork(qrwsiz))

    IF (retry) THEN !.AND.n<1.d-5) THEN
      JACFLG(1:4) = (/1,-1,-1,1/)
      OUTOPT(1:2) = (/1,1/)
    endif

    RJAC = ZERO ; RWORK = ZERO; RCDWRK = ZERO ; QRWORK = ZERO ; ICDWRK = 0
    SCALEX = ONE

!   if no initial guess given look for one in the "grid"
    IF (sum(x_guess)==zero) then
      ! if retry is true grid will have twice as many points.
      CALL GRID (x_grid,retry)
      ! replace guess with point from grid
      x2 = x_grid
    ELSE
      x2 = x_guess
    ENDIF

    x1 = zero

!   solve equation for equilibrium for uniform matter
!   if matter is neutron rich
    IF (Yp<=half) &
    CALL nwnleq(x2,N_eq,scalex,maxit,jacflg,xtol,ftol,btol,cndtol,method,global, &
                xscalm,stepmx,delta,sigma,rjac,ldr,rwork,lrwork, &
                rcdwrk,icdwrk,qrwork,qrwsiz,&
                JACOB_SOLVE_UNIFORM_NEUTRON_RICH,SOLVE_UNIFORM_NEUTRON_RICH,&
                outopt,x1,r,g,njcnt,nfcnt,iter,termcd)
!   if matter is proton rich
    IF (Yp> half) &
    CALL nwnleq(x2,N_eq,scalex,maxit,jacflg,xtol,ftol,btol,cndtol,method,global, &
                xscalm,stepmx,delta,sigma,rjac,ldr,rwork,lrwork, &
                rcdwrk,icdwrk,qrwork,qrwsiz,&
                JACOB_SOLVE_UNIFORM_PROTON_RICH,SOLVE_UNIFORM_PROTON_RICH,&
                outopt,x1,r,g,njcnt,nfcnt,iter,termcd)

    deallocate(RJAC, RWORK, RCDWRK, ICDWRK, QRWORK, SCALEX)
!   chech whether solution x1 found is actually a solution to eq being solved.
!   sotemimes output for x1 is not a solution, but a point where nlwleq stalled.

    x_lgcl = .FALSE.
    residue = DOT_PRODUCT(r,r)
    IF (residue<1.d-12) x_lgcl = .TRUE.

    X_p = TEN**x1(1)
    IF (Yp>HALF) THEN
      log10_n_po = log10(n*(-X_p*(n*v_alpha)*(Yp) &
                      - two*((one-two*Yp)-X_p))/(two-n*v_alpha*(one-Yp)))
      log10_n_no = log10(n*X_p)
      log10_u    = -290.D0
    ELSE
      log10_n_no = log10(n*(-X_p*(n*v_alpha)*(one-Yp) &
                          + two*((one-two*Yp)+X_p))/(two-n*v_alpha*Yp))
      log10_n_po = log10(n*X_p)
      log10_u    = -290.0D0
    ENDIF

    IF (x_lgcl) x_sol = x1

    CALL FREE_ENERGY(log10_n_no, log10_n_po, log10_u, n, T, Yp, &
          F_out, F_in, F_alpha, F_heavy, F_surf_coul, &
          n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
          un_A, un_Z, un_rad, un_F, un_P, un_S, un_E, &
          un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT, &
          un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT, &
          un_mu_no, un_mu_po, un_Meff_no, un_Meff_po, &
          un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy, &
          un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy, &
          un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn, &
          un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT, &
          non_uniform_eq, non_uniform_jac, 0, error)

    IF (IS_TEST) THEN
      WRITE (*,*) ' SOLVING UNIFORM SYSTEM '
      WRITE (*,*)
      WRITE (*,*) '  INITIAL GUESS:   ', x_guess, ' (0.0 if none)'
      WRITE (*,*) '  SOLUTION:        ', x1 
      WRITE (*,*) '  RESIDUE:         ', r
      WRITE (*,*) '  FREE ENERGY:     ', (F_out+F_alpha)/n
      WRITE (*,*) '  JACOBIAN:        ', g
      WRITE (*,*) '  TERMINATION CODE:', termcd
      WRITE (*,*) '  ERROR:           ', error
      WRITE (*,*)
      IF (termcd < 0 ) WRITE (*,*) '   INVALID VALUES FOR SOLUTION!  '
      IF (termcd == 1 .AND. error == 1) WRITE (*,*) '   SOLUTION FOUND!  '
      IF (termcd == 1 .AND. error /= 1) WRITE (*,*) '   SOLUTION NOT-PHYSICAL!  '
      IF (termcd == 2) WRITE (*,*) '   SOLUTION MAY BE INACCURATE!'
      IF (termcd == 3) WRITE (*,*) '   SOLUTION MAY BE INACCURATE!'
      IF (termcd == 4) WRITE (*,*) '   ITERATION LIMIT EXCEEDED AND SOLUTION NOT FOUND!'
      IF (termcd == 5) WRITE (*,*) '   PROBLEM WITH JACOBIAN'
      IF (termcd == 6) WRITE (*,*) '   PROBLEM WITH JACOBIAN'
      WRITE (*,*)
    ENDIF

    Ftotal = (F_out + F_in + F_alpha + F_heavy + F_surf_coul)/n

    RETURN

  END SUBROUTINE Find_Uniform_Solution

! checks on a grid which point for independent variable Xp
!  gets closer to solving system of equations for uniform matter
  SUBROUTINE GRID ( x_inout, lgcl_retry )

    IMPLICIT NONE

    REAL(DP), DIMENSION(1), INTENT(INOUT) :: x_inout
    LOGICAL(LGCL), INTENT(IN) :: lgcl_retry
    INTEGER(I4B) :: i, imin, imax, loc, location(1)
    REAL(DP) :: Xp, denom, residue(1), array_dummy(2)
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temporary
    LOGICAL(LGCL) :: sol_lgcl
    INTEGER(I4B) :: error, flag
    INTEGER(I4B), PARAMETER :: n_dim = 1

    IF (x_inout(1) /= zero) THEN
      ! if initial guess for solution not zero
      ! then try points near initial guess
      ! Yp/10 <= Xp <= Yp*10 with dXp = 10^(1/20)
      imin = -12 ; imax = +12
      denom = 240.0d0
    ELSE
      ! if initial guess for solution not given (xp = zero)
      ! then try points in range below. (should be enough)
      ! Yp*10^-50 <= Xp <= Yp*10 with dXp = 10^(1/200)
      imin = -10 ; imax = 20*int(TEN**(-log10(T)+2))
      IF (T>1.D0) imax = 2000
      IF (T>5.D0) imax = 1000
      IF (T>1.D1) imax = 200
      IF (T>2.D1) imax = 100
      denom = 200.0d0
      ! for low temperatures no need to
      ! try very low neutron/proton fraction
      ! (check later if actually true)
      ! IF (T<0.1_DP) IMAX = 1
    ENDIF

!   if this is a retry then double
!   number of grid points to look for solution
    IF (lgcl_retry) THEN
      denom = two*denom ; imax = 2*imax
    ENDIF

!   set a temporary array to store values of
!    function we're trying to minimize
    ALLOCATE(Temporary(imin:imax,1:2))

    IF (IS_TEST) WRITE (*,*) log10(Yp)-dble(imin)/denom, ' < X_p < ', &
                             log10(Yp)-dble(imax)/denom, imin, imax

    IF (x_inout(1) == zero) THEN
      DO i = imin, imax
!       set initial guess for nucleon densities
!       based on whether system is neutron or proton rich
        IF (Yp<HALF) THEN
          Xp = log10(Yp)-dble(i)/denom
          CALL SOLVE_UNIFORM_NEUTRON_RICH( (/Xp/), residue, n_dim, flag )
        ELSE
          Xp = log10(one-Yp)-dble(i)/denom
          CALL SOLVE_UNIFORM_PROTON_RICH( (/Xp/), residue, n_dim, flag )
        ENDIF
        Temporary(i,1:2) = (/ Xp, abs(residue) /)
      ENDDO
      location = MINLOC(Temporary(:,2))
      loc = location(1) - 1 + imin
      IF (Temporary(loc,1)>=1.d100) RETURN
      x_inout = Temporary(loc,1)
    ELSE
      CONTINUE
      DO i = imin, imax
        if (i<0)  Xp = x_inout(1)*(one+dble(i)/denom)
        if (i>0)  Xp = x_inout(1)*(one-dble(-i)/denom)
        if (i==0) Xp = x_inout(1)
        IF (Yp<HALF) THEN
          CALL SOLVE_UNIFORM_NEUTRON_RICH( (/Xp/), residue, n_dim, flag )
        ELSE
          CALL SOLVE_UNIFORM_PROTON_RICH( (/Xp/), residue, n_dim, flag )
        ENDIF
      ENDDO
    ENDIF

    DEALLOCATE(Temporary)

    RETURN

  END SUBROUTINE GRID

  SUBROUTINE SOLVE_UNIFORM_NEUTRON_RICH ( x, s, m, uniform_flag )
!   given initial input for independent parameter 'x'
!   determine residue 's' for Eq that solves EoS of
!   uniform matter for neutron rich matter
    IMPLICIT NONE

    INTEGER (I4B), INTENT(IN) :: m, uniform_flag
    REAL(DP), INTENT(IN)  :: x(m)
    REAL(DP), INTENT(OUT) :: s(m)

    REAL(DP) :: X_p, log10_n_no, log10_n_po, log10_u
    REAL(DP) :: F_out, F_in, F_alpha, F_heavy, F_surf_coul
    REAL(DP) :: n_no, n_po, n_o, n_ni, n_pi, n_alpha, n_heavy
    REAL(DP) :: exc_v_alpha
    REAL(DP) :: non_uniform_eq(3), non_uniform_jac(3,3)
    INTEGER(I4B) :: error

    !  REAL(DP) :: y_o,Meff_no,Meff_po,eta_no,eta_po,tau_no,tau_po, &
    !  v_no,v_po,mu_no,mu_po,P_o,F_o,P_alpha,mu_alpha

    X_p = TEN**x(1)

    log10_n_no = log10(n*(-X_p*(n*v_alpha)*(one-Yp) &
                        + two*((one-two*Yp)+X_p))/(two-n*v_alpha*Yp))
    log10_n_po = log10(n*X_p)
    log10_u = -290.0D0

    IF (ISNAN(log10_n_no) .OR. ISNAN(log10_n_po)) THEN
      S(1) = 1.D100
      RETURN
    ENDIF

    CALL FREE_ENERGY(log10_n_no, log10_n_po, log10_u, n, T, Yp, &
          F_out, F_in, F_alpha, F_heavy, F_surf_coul, &
          n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
          un_A, un_Z, un_rad, un_F, un_P, un_S, un_E, &
          un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT, &
          un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT, &
          un_mu_no, un_mu_po, un_Meff_no, un_Meff_po, &
          un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy, &
          un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy, &
          un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn, &
          un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT, &
          non_uniform_eq, non_uniform_jac, 0, error)

    !  CALL SKYRME_BULK_PROPERTIES(log10_n_no,log10_n_po,&
    !    T,n_no,n_po,n_o,y_o,Meff_no,Meff_po,eta_no,eta_po,tau_no,tau_po, &
    !    v_no,v_po,mu_no,mu_po,P_o,F_o)
     !
    !  CALL ALPHA_PROPERTIES (0, n_no, n_po, mu_no, mu_po, T, F_o, P_o, &
    !          n_alpha, mu_alpha, F_alpha, P_alpha)

    IF (error<0) THEN
      ! IF (IS_TEST) WRITE (*,*) log10_n_no,log10_n_po,log10_u,n_no, n_po, n_alpha,error
      S(1) = 1.D100
      RETURN
    ENDIF

    n_o = n_no + n_po
    exc_v_alpha = one - n_alpha*v_alpha

    S(1) = ONE - (n_o*exc_v_alpha + FOUR*n_alpha)/n
    ! IF (IS_TEST) WRITE (*,*) log10_n_no,log10_n_po,log10_u,un_F,S,error

  END SUBROUTINE SOLVE_UNIFORM_NEUTRON_RICH

  SUBROUTINE SOLVE_UNIFORM_PROTON_RICH ( x, s, m, uniform_flag )
!   given initial input for independent parameter 'x'
!   determine residue 's' for Eq that solves EoS of
!   uniform matter for neutron rich matter
    IMPLICIT NONE

    INTEGER (I4B), INTENT(IN) :: m, uniform_flag
    REAL(DP), INTENT(IN)  :: x(m)
    REAL(DP), INTENT(OUT) :: s(m)

    REAL(DP) :: X_p, log10_n_no, log10_n_po, log10_u
    REAL(DP) :: F_out, F_in, F_alpha, F_heavy, F_surf_coul
    REAL(DP) :: n_no, n_po, n_o, n_ni, n_pi, n_alpha, n_heavy
    REAL(DP) :: exc_v_alpha
    REAL(DP) :: non_uniform_eq(3), non_uniform_jac(3,3)
    INTEGER(I4B) :: error

    X_p = TEN**x(1)

    log10_n_po = log10(n*(-X_p*(n*v_alpha)*(Yp) &
                    - two*((one-two*Yp)-X_p))/(two-n*v_alpha*(one-Yp)))
    log10_n_no = log10(n*X_p)
    log10_u = -290.0D0

    IF (ISNAN(log10_n_no) .OR. ISNAN(log10_n_po)) THEN
      S(1) = 1.D100
      RETURN
    ENDIF

    CALL FREE_ENERGY(log10_n_no, log10_n_po, log10_u, n, T, Yp, &
          F_out, F_in, F_alpha, F_heavy, F_surf_coul, &
          n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
          un_A, un_Z, un_rad, un_F, un_P, un_S, un_E, &
          un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT, &
          un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT, &
          un_mu_no, un_mu_po, un_Meff_no, un_Meff_po, &
          un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy, &
          un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy, &
          un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn, &
          un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT, &
          non_uniform_eq, non_uniform_jac, 0, error)

    IF (error<0) THEN
      S(1) = 1.D100
      RETURN
    ENDIF

    n_o = n_no + n_po
    exc_v_alpha = one - n_alpha*v_alpha

    S(1) = ONE - (n_o*exc_v_alpha + FOUR*n_alpha)/n

  END SUBROUTINE SOLVE_UNIFORM_PROTON_RICH

  SUBROUTINE JACOB_SOLVE_UNIFORM_NEUTRON_RICH ( jac, ldr, x1, m )
!   given initial input for independent parameter 'x'
!   determine residue 's' for Eq that solves EoS of
!   uniform matter for neutron rich matter
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: m, ldr
    REAL(DP), INTENT(IN)  :: x1(m)
    REAL(DP), INTENT(OUT) :: jac(ldr,*)

    REAL(DP) :: X_p, log10_n_no, log10_n_po, log10_u
    REAL(DP) :: G_n, G_p
    REAL(DP) :: log10_n_n,log10_n_p
    REAL(DP) :: Deta_n_Dn_n, Deta_n_Dn_p, Deta_p_Dn_n, Deta_p_Dn_p
    REAL(DP) :: Dtau_n_Dn_n, Dtau_n_Dn_p, Dtau_p_Dn_n, Dtau_p_Dn_p
    REAL(DP) :: Dmu_n_Dn_n, Dmu_n_Dn_p, Dmu_p_Dn_n, Dmu_p_Dn_p
    REAL(DP) :: Dv_n_Dn_n, Dv_n_Dn_p, Dv_p_Dn_n, Dv_p_Dn_p
    REAL(DP) :: DP_Dn_n, DP_Dn_p
    REAL(DP) :: F_out, F_in, F_alpha, F_heavy, F_surf_coul
    REAL(DP) :: n_no, n_po, n_o, n_ni, n_pi, n_alpha, n_heavy
    REAL(DP) :: exc_v_alpha
    REAL(DP) :: DP_alpha_Dn_n, DP_alpha_Dn_p
    REAL(DP) :: Dn_alpha_Dn_n, Dn_alpha_Dn_p
    REAL(DP) :: Dn_n_Dxp, Dn_p_Dxp, Dn_alpha_Dxp
    REAL(DP) :: non_uniform_eq(3), non_uniform_jac(3,3)
    INTEGER(I4B) :: error

    X_p = TEN**x1(1)

    log10_n_no = log10(n*(-X_p*(n*v_alpha)*(one-Yp) &
                        + two*((one-two*Yp)+X_p))/(two-n*v_alpha*Yp))
    log10_n_po = log10(n*X_p)
    log10_u = -290.0D0

    IF (ISNAN(log10_n_no) .OR. ISNAN(log10_n_po)) THEN
      JAC(1,1) = ZERO
      RETURN
    ENDIF

    CALL FREE_ENERGY(log10_n_no, log10_n_po, log10_u, n, T, Yp, &
          F_out, F_in, F_alpha, F_heavy, F_surf_coul, &
          n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
          un_A, un_Z, un_rad, un_F, un_P, un_S, un_E, &
          un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT, &
          un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT, &
          un_mu_no, un_mu_po, un_Meff_no, un_Meff_po, &
          un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy, &
          un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy, &
          un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn, &
          un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT, &
          non_uniform_eq, non_uniform_jac, 0, error)

    IF (error<0) THEN
      JAC(1,1) = ZERO
      RETURN
    ENDIF

    n_o = n_no + n_po
    exc_v_alpha = one - n_alpha*v_alpha

    CALL SKYRME_BULK_DENSITY_DERIVATIVES( &
          log10_n_no, log10_n_po, T, G_n, G_p, &
          Deta_n_Dn_n, Deta_n_Dn_p, Deta_p_Dn_n, Deta_p_Dn_p, &
          Dtau_n_Dn_n, Dtau_n_Dn_p, Dtau_p_Dn_n, Dtau_p_Dn_p, &
          Dmu_n_Dn_n, Dmu_n_Dn_p,   Dmu_p_Dn_n, Dmu_p_Dn_p,   &
          Dv_n_Dn_n, Dv_n_Dn_p,     Dv_p_Dn_n, Dv_p_Dn_p,     DP_Dn_n, DP_Dn_p)

    Dn_p_Dxp = n*X_p*LOG(TEN)

    Dn_n_Dxp = n*(two-n*v_alpha*(one-Yp))/&
                            (two-n*v_alpha*Yp)*X_p*LOG(TEN)
!
    DP_alpha_Dn_n = n_alpha*( (two-n_no*v_alpha)*Dmu_n_Dn_n + &
                              (two-n_po*v_alpha)*Dmu_p_Dn_n )

    DP_alpha_Dn_p = n_alpha*( (two-n_no*v_alpha)*Dmu_n_Dn_p + &
                              (two-n_po*v_alpha)*Dmu_p_Dn_p )

    Dn_alpha_Dn_n = DP_alpha_Dn_n/T
    Dn_alpha_Dn_p = DP_alpha_Dn_p/T

    Dn_alpha_Dxp = Dn_alpha_Dn_n*Dn_n_Dxp + Dn_alpha_Dn_p*Dn_p_Dxp

    JAC(1,1) = -ONE/n*( (Dn_p_Dxp+Dn_n_Dxp)*exc_v_alpha - &
                            (v_alpha*n_o-four)*Dn_alpha_Dxp )

  END SUBROUTINE JACOB_SOLVE_UNIFORM_NEUTRON_RICH

  SUBROUTINE JACOB_SOLVE_UNIFORM_PROTON_RICH ( jac, ldr, x1, m )
!   given initial input for independent parameter 'x'
!   determine residue 's' for Eq that solves EoS of
!   uniform matter for neutron rich matter
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: m, ldr
    REAL(DP), INTENT(IN)  :: x1(m)
    REAL(DP), INTENT(OUT) :: jac(ldr,*)

    REAL(DP) :: X_p, log10_n_no, log10_n_po, log10_u
    REAL(DP) :: G_n, G_p
    REAL(DP) :: log10_n_n,log10_n_p
    REAL(DP) :: Deta_n_Dn_n, Deta_n_Dn_p, Deta_p_Dn_n, Deta_p_Dn_p
    REAL(DP) :: Dtau_n_Dn_n, Dtau_n_Dn_p, Dtau_p_Dn_n, Dtau_p_Dn_p
    REAL(DP) :: Dmu_n_Dn_n, Dmu_n_Dn_p, Dmu_p_Dn_n, Dmu_p_Dn_p
    REAL(DP) :: Dv_n_Dn_n, Dv_n_Dn_p, Dv_p_Dn_n, Dv_p_Dn_p
    REAL(DP) :: DP_Dn_n, DP_Dn_p
    REAL(DP) :: F_out, F_in, F_alpha, F_heavy, F_surf_coul
    REAL(DP) :: n_no, n_po, n_o, n_ni, n_pi, n_alpha, n_heavy
    REAL(DP) :: exc_v_alpha
    REAL(DP) :: DP_alpha_Dn_n, DP_alpha_Dn_p
    REAL(DP) :: Dn_alpha_Dn_n, Dn_alpha_Dn_p
    REAL(DP) :: Dn_n_Dxp, Dn_p_Dxp, Dn_alpha_Dxp
    REAL(DP) :: non_uniform_eq(3), non_uniform_jac(3,3)
    INTEGER(I4B) :: error

    X_p = TEN**x1(1)

    log10_n_po = log10(n*(-X_p*(n*v_alpha)*(Yp) &
                    - two*((one-two*Yp)-X_p))/(two-n*v_alpha*(one-Yp)))
    log10_n_no = log10(n*X_p)
    log10_u = -290.0D0

    IF (ISNAN(log10_n_no) .OR. ISNAN(log10_n_po)) THEN
      JAC(1,1) = ZERO
      RETURN
    ENDIF

    CALL FREE_ENERGY(log10_n_no, log10_n_po, log10_u, n, T, Yp, &
          F_out, F_in, F_alpha, F_heavy, F_surf_coul, &
          n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
          un_A, un_Z, un_rad, un_F, un_P, un_S, un_E, &
          un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT, &
          un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT, &
          un_mu_no, un_mu_po, un_Meff_no, un_Meff_po, &
          un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy, &
          un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy, &
          un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn, &
          un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT, &
          non_uniform_eq, non_uniform_jac, 0, error)

    IF (error<0) THEN
      JAC(1,1) = ZERO
      RETURN
    ENDIF

    n_o = n_no + n_po
    exc_v_alpha = one - n_alpha*v_alpha

    CALL SKYRME_BULK_DENSITY_DERIVATIVES( &
          log10_n_no, log10_n_po, T, G_n, G_p, &
          Deta_n_Dn_n, Deta_n_Dn_p, Deta_p_Dn_n, Deta_p_Dn_p, &
          Dtau_n_Dn_n, Dtau_n_Dn_p, Dtau_p_Dn_n, Dtau_p_Dn_p, &
          Dmu_n_Dn_n, Dmu_n_Dn_p,   Dmu_p_Dn_n, Dmu_p_Dn_p,   &
          Dv_n_Dn_n, Dv_n_Dn_p,     Dv_p_Dn_n, Dv_p_Dn_p,     DP_Dn_n, DP_Dn_p)

    Dn_p_Dxp = n*(two-n*v_alpha*(Yp))/&
                            (two-n*v_alpha*(one-Yp))*X_p*LOG(TEN)

    Dn_n_Dxp = n*X_p*LOG(TEN)

    DP_alpha_Dn_n = n_alpha*( (two-n_no*v_alpha)*Dmu_n_Dn_n + &
                              (two-n_po*v_alpha)*Dmu_p_Dn_n )

    DP_alpha_Dn_p = n_alpha*( (two-n_no*v_alpha)*Dmu_n_Dn_p + &
                              (two-n_po*v_alpha)*Dmu_p_Dn_p )

    Dn_alpha_Dn_n = DP_alpha_Dn_n/T
    Dn_alpha_Dn_p = DP_alpha_Dn_p/T

    Dn_alpha_Dxp = Dn_alpha_Dn_n*Dn_n_Dxp + Dn_alpha_Dn_p*Dn_p_Dxp

    JAC(1,1) = -ONE/n*( (Dn_p_Dxp+Dn_n_Dxp)*exc_v_alpha - &
                            (v_alpha*n_o-four)*Dn_alpha_Dxp )

  END SUBROUTINE JACOB_SOLVE_UNIFORM_PROTON_RICH

END MODULE Find_Uniform_Solution_Mod
