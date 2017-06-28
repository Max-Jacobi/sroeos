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
MODULE Find_Non_Uniform_Solution_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL
  USE Global_Variables_Mod, ONLY : IS_TEST
  USE Physical_Constants_Mod, ONLY : ZERO, HALF, ONE, TWO, THREE, FOUR, TEN, &
                             V_Alpha, energy_cgs_to_EOS, temp_mev_to_kelvin, &
                             press_cgs_to_EOS, rho_cgs_to_EOS
  USE Nuclear_Matter_Properties_Mod, ONLY : Nuc_Sat_Dens
  USE Free_Energy_Mod
  USE Skyrme_Bulk_Density_Derivatives_Mod
  USE nwnleq_mod
  USE lautil_mod
  USE unused_mod
  USE wrap
  USE phase_space_point_mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Find_Non_Uniform_Solution( n_in, T_in, Yp_in, x_guess, x_sol, &
                                        Ftotal, residue, retry, x_lgcl)

    IMPLICIT NONE
!
    REAL(DP), INTENT(IN) :: n_in, T_in, Yp_in
    REAL(DP), INTENT(OUT) :: x_sol(3), Ftotal, residue
    REAL(DP), INTENT(INOUT) :: x_guess(3)
    LOGICAL(LGCL), INTENT(OUT) :: x_lgcl
    LOGICAL(LGCL), INTENT(IN)  :: retry
    LOGICAL(LGCL) :: x_residue, x_gamma, x_A, x_Z
    REAL(DP), DIMENSION(3,3) :: g
    REAL(DP) :: x_grid(500,4), x1(3), x2(3)
    REAL(DP) :: r(3), s(3), dAdu(3,3), jacobian(3,3), x_step(3)
    INTEGER(I4B) :: i, objective, error, count_max, eval_max
    REAL(DP) :: Fmin, F_out, F_in, F_alpha, F_trans, F_surf_coul
    REAL(DP) :: dummy, solution(3), solution_jac(3,3)
    REAL(DP) :: n_no, n_po
    REAL(DP) :: n_ni, n_pi
    REAL(DP) :: n_alpha, n_heavy
    REAL(DP) :: x_no, x_po, x_alpha, x_heavy, Log10_u
    REAL(DP) :: DF_Dn, DF_Dy, DF_DT, DP_Dn, DP_Dy, DP_DT
    REAL(DP) :: DS_Dn, DS_Dy, DS_DT, DE_Dn, DE_Dy, DE_DT
    REAL(DP) :: A, Z, P, P_tot, Abar, Zbar, exc_v_alpha, u
    REAL(DP) :: DE_DT_tot, DP_DT_tot, DP_Dn_tot
    REAL(DP) :: rho,T_K,etot,ptot,stot
    REAL(DP) :: dedT,dpdT,dsdT,dedr,dpdr,dsdr,dedy,dpdy,dsdy,gamma,etaele,sound

!   parameters for non-linear equation solver "nleqslv"
    INTEGER(I4B) :: N_eq, MAXIT, JACFLG(1:4),  SIZE
    INTEGER(I4B) :: METHOD, GLOBAL, XSCALM, LDR, LRWORK
    INTEGER(I4B) :: NJCNT, NFCNT, ITER, TERMCD, QRWSIZ, OUTOPT(2)
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ICDWRK
    REAL(DP) :: XTOL, FTOL, BTOL, CNDTOL
    REAL(DP) :: STEPMX, DELTA, SIGMA, TRACE, DSUB
    REAL(DP), DIMENSION(:), ALLOCATABLE :: RJAC, RWORK, RCDWRK
    REAL(DP), DIMENSION(:), ALLOCATABLE :: QRWORK, SCALEX

!   set parameters for non-linear equation solver
    N_eq = 3; MAXIT = 150 ; JACFLG(1:4) = (/1,-1,-1, 1/) ; OUTOPT(1:2) = (/1,1/)
    METHOD = 0; GLOBAL = 4; XSCALM = 1; LDR = N_eq; LRWORK = 9*N_eq
    XTOL = 1.D-15; FTOL = 1.D-8; BTOL = 1.D-6; CNDTOL = 1.D-8
    STEPMX = -1.D0; DELTA = -1.D0; SIGMA = 0.5D0; TRACE = 1.D0; DSUB = 0.D0
    SIZE = (METHOD+1)*N_eq*N_eq

    n = n_in
    T = T_in
    Yp = Yp_in

    ALLOCATE(RJAC(SIZE),RWORK(9*N_eq),RCDWRK(3*N_eq),ICDWRK(N_eq),SCALEX(N_eq))
    CALL liqsiz(N_eq,qrwsiz)
    ALLOCATE(qrwork(qrwsiz))

    IF (retry.AND.n<1.d-5) THEN
      JACFLG(1:4) = (/1,-1,-1,1/)
      OUTOPT(1:2) = (/1,1/)
      objective = 0
    endif

    RJAC = ZERO ; RWORK = ZERO; RCDWRK = ZERO ; QRWORK = ZERO ; ICDWRK = 0
    SCALEX = ONE

!   get an array of points with increasing free energy
!    in a grid near the initial guess
!    or in some reasonable large grid where the solution may be
!   (if retry is true grid will have twice as many points.)
    CALL GRID (x_guess,x_grid,retry)

    x2 = x_guess
    x1 = zero

    i = 0
!   solve equation for equilibrium for non uniform matter
    x_lgcl    = .FALSE.
    x_residue = .FALSE.
    x_gamma   = .FALSE.
    x_A       = .FALSE.
    x_Z       = .FALSE.

    eval_max = 0
    count_max = 500
    IF (n<0.01d0) count_max = 200

    DO WHILE (i < count_max .and. eval_max < 25 .and. .not. x_lgcl)
      i = i + 1
      x2(1:3) = x_grid(i,1:3)
      IF (x_grid(i,4)>1.d99) CYCLE
      CALL nwnleq(x2,N_eq,scalex,maxit,jacflg,xtol,ftol,btol,cndtol,method,&
                  global,xscalm,stepmx,delta,sigma,rjac,ldr,rwork,lrwork, &
                  rcdwrk,icdwrk,qrwork,qrwsiz,&
                  jacob_SOLVE_NON_UNIFORM,SOLVE_NON_UNIFORM,&
                  outopt,x1,r,g,njcnt,nfcnt,iter,termcd)

      x_step(1:3) = x_grid(i,1:3)-x1(1:3)

      residue = DOT_PRODUCT(r,r)
      IF (residue<1.d-12) THEN
        eval_max = eval_max + 1
        x_residue = .TRUE.
      ELSE
        CYCLE
      ENDIF

      objective = 4

      CALL FREE_ENERGY(x1(1), x1(2), x1(3), n, T, Yp,&
          F_out, F_in, F_alpha, F_trans, F_surf_coul,&
          n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
          A, Z, un_rad, un_F, P, un_S, un_E, &
          DF_Dn, DF_Dy, DF_DT, DP_Dn, DP_Dy, DP_DT, &
          DS_Dn, DS_Dy, DS_DT, DE_Dn, DE_Dy, DE_DT, &
          un_mu_no, un_mu_po, un_Meff_no, un_Meff_po, &
          un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy, &
          un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy, &
          un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn, &
          un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT, &
          solution, solution_jac, objective, error)

!     heavy nuclei must have A>=10.
      IF (A>10.0D0 .AND. A<1.0D3) THEN
        x_A = .TRUE.
      ELSE
        CYCLE
      ENDIF

!     heavy nuclei must have Z >= 6. 
!     Uless Yp < 1%, then Z > 2.
      IF (Yp_in >= 0.01D0 .AND. Z> 6.0D0 .AND. Z<1.0D3) THEN
        x_Z = .TRUE.
      ELSEIF (Yp_in <= 0.01D0 .AND. Z > 2.0D0 .AND. Z<1.0D3) THEN
        x_Z = .TRUE.
      ELSE
        CYCLE
      ENDIF

!     proton fraction of heavy nuclei must be in (0,1) range
      IF (Z/A <= 0.0D0 .OR. Z/A >= ONE) THEN
        CYCLE
      ENDIF

      Fmin = F_out + F_in + F_alpha + F_trans + F_surf_coul
!     get adiabatic index
!     TODO: make subroutine
!     convert temperature to K and density to g/cm^3 before call to wrap_timmes
      T_K = T*temp_mev_to_kelvin
      rho = n/rho_cgs_to_EOS
      !   occupied volume
      Log10_u = x1(3)
      u = TEN**Log10_u
      !   set u to zero if uniform atter is the solution
      IF (u < 1.d-100) u = zero
      !   volume excluded by alpha particles
      exc_v_alpha = one - n_alpha*v_alpha
      !   get number fractions
      x_no = (ONE-u)*exc_v_alpha*n_no / n
      x_po = (ONE-u)*exc_v_alpha*n_po / n
      x_alpha = (ONE-u)*FOUR*n_alpha / n
      x_heavy = (n_ni+n_pi)*ten**log10_u/n
      IF (n_heavy>zero .AND. A>zero) THEN
        Abar = one/(x_no+x_po+x_alpha/four+x_heavy/A)
      ELSE
        Abar = one/(x_no+x_po+x_alpha/four)
      ENDIF
      Zbar = Abar*Yp
      CALL WRAP_TIMMES(Abar,Zbar,rho,T_K,etot,ptot,stot,&
          dedT,dpdT,dsdT,dedr,dpdr,dsdr,dedy,dpdy,dsdy,gamma,etaele,sound)
      P_tot = P + ptot*press_cgs_to_EOS
      DP_Dn_tot = DP_Dn + dpdr*press_cgs_to_EOS/rho_cgs_to_EOS
      DP_DT_tot = DP_DT + dpdT*press_cgs_to_EOS*temp_mev_to_kelvin
      De_DT_tot = DE_DT/n + dedt*energy_cgs_to_EOS*temp_mev_to_kelvin
      Gamma = n/P_tot*DP_Dn_tot + T*(DP_DT_tot*DP_DT_tot)/n/P_tot/DE_DT_tot

      IF (gamma > ZERO) THEN
        x_gamma = .TRUE.
      ELSEIF (gamma > -0.1D0) THEN
        gamma = ZERO
        x_gamma = .TRUE.
!     accept gamma < 0 if Yp < 0.01.
!     helps with the LS parametrization
      ELSEIF (Yp<0.01d0 .AND. gamma > -1.D0) THEN
        gamma = ZERO
        x_gamma = .TRUE.
      ELSE
        CYCLE
      ENDIF

      IF (x_residue .AND. x_gamma .AND. x_A .AND. x_Z) x_lgcl = .TRUE.
    ENDDO

    DEALLOCATE(RJAC, RWORK, RCDWRK, ICDWRK, QRWORK, SCALEX)

    IF (IS_TEST) THEN
      WRITE (*,*) A, Z, GAMMA, x_gamma, x_A, x_Z, x_residue
      WRITE (*,*) ' SOLVING NON-UNIFORM SYSTEM '
      WRITE (*,*)
      WRITE (*,*) '   INITIAL GUESS:   ', x_guess
      WRITE (*,*) '   SOLUTION:        ', x1
      WRITE (*,*) '   RESIDUE:         ', r
      WRITE (*,*) '   FREE ENERGY:     ', Fmin/n
      WRITE (*,*) '   JACOBIAN(1,1:3): ', solution_jac(1,:)
      WRITE (*,*) '   JACOBIAN(2,1:3): ', solution_jac(2,:)
      WRITE (*,*) '   JACOBIAN(3,1:3): ', solution_jac(3,:)
      WRITE (*,*) '   TERMINATION CODE:', termcd
      WRITE (*,*) '   ERROR:           ', error
      WRITE (*,*)
      IF (termcd < 0 ) WRITE (*,*) '   INVALID VALUES FOR SOLUTION!  '
      IF (termcd == 1 .AND. error == 0) WRITE (*,*) '   SOLUTION FOUND!  '
      IF (termcd == 1 .AND. error /= 0) WRITE (*,*) '   SOLUTION NOT-PHYSICAL!  '
      IF (termcd == 2) WRITE (*,*) '   SOLUTION MAY BE INACCURATE!'
      IF (termcd == 3) WRITE (*,*) '   SOLUTION MAY BE INACCURATE!'
      IF (termcd == 4) WRITE (*,*) '   ITERATION LIMIT EXCEEDED AND SOLUTION NOT FOUND!'
      IF (termcd == 5) WRITE (*,*) '   NO SOLUTION: ILL CONDITIONED JACOBIAN!'
      IF (termcd == 6) WRITE (*,*) '   NO SOLUTION: SINGULAR JACOBIAN!'
      WRITE (*,*)
    ENDIF

    ! make restriction on gamma less strict

    IF (x_lgcl .AND. Gamma < zero) THEN
      WRITE (*,"(A36,3ES16.7,A10,ES16.7)") &
        'WARNING! GAMMA < 0 FOR (Yp,T,n) = ', Yp, T, n, ' GAMMA = ', GAMMA
      IF (GAMMA < -ONE .OR. n > 0.05d0 .OR. T > ONE) THEN
        x_lgcl = .FALSE.
        IF (.NOT. retry) THEN
          WRITE (*,"(A36)") ' WILL TRY TO FIND OTHER SOLUTION!'
        ELSE
          WRITE (*,"(A36)") ' WILL SET SOLUTION TO UNIFORM!'
        ENDIF
      ELSE
        GAMMA = 1.D-10
        WRITE (*,"(A36)") ' KEPT SOLUTION, BUT SET GAMMA TO 1.D-10!'
      ENDIF
    ENDIF

    IF (x_lgcl .AND. A<10.D0) THEN
      WRITE (*,"(A36,3ES16.7,A10,ES16.7)") &
        'WARNING! A < 10 FOR (Yp,T,n) = ', Yp, T, n, ' A = ', A
      WRITE (*,"(A66)") ' IGNORED SOLUTION BECAUSE &
                          HEAVY NUCLEI MASS NUMBER IS TOO SMALL!'
      x_lgcl = .FALSE.
    ENDIF

    IF (x_lgcl) THEN
      x_sol   = x1
      Ftotal  = Fmin/n
    ENDIF
!   chech whether solution x1 found is actually a solution to eq being solved.
!   sotemimes output for x1 is not a solution, but a point where nlwleq stalled.

    ! IF (IS_TEST) write (*,"(10ES15.6,L4)") x_grid(i,1:3), x1, r, fmin, x_lgcl
!   TODO: MAKE FOLLOWING TEST OUTPUT A SUBROUTINE
    IF (IS_TEST .AND. x_lgcl) THEN
      WRITE (*,*) 'NO OTHER ERRORS FOUND!'
      WRITE (*,*)
      CALL  SOLVE_NON_UNIFORM ( x_sol, R, n_eq, n_eq )
      x_sol(3) = x_sol(3)+1.d-8
      CALL  SOLVE_NON_UNIFORM ( x_sol, S, n_eq, n_eq )
      dAdu(1:3,3) = (s(1:3)-r(1:3))/(1.d-8)
      x_sol(3) = x_sol(3)-1.d-8
      x_sol(2) = x_sol(2)+1.d-8
      CALL  SOLVE_NON_UNIFORM ( x_sol, S, n_eq, n_eq )
      dAdu(1:3,2) = (s(1:3)-r(1:3))/(1.d-8)
      x_sol(2) = x_sol(2)-1.d-8
      x_sol(1) = x_sol(1)+1.d-8
      CALL  SOLVE_NON_UNIFORM ( x_sol, S, n_eq, n_eq )
      dAdu(1:3,1) = (s(1:3)-r(1:3))/(1.d-8)
      x_sol(1) = x_sol(1)-1.d-8
      CALL JACOB_SOLVE_NON_UNIFORM ( jacobian, ldr, x_sol, n_eq )
      write (*,*)
      write (*,*) 'Comparison of numerical and analytical derivatives'
      write (*,*)
      write (*,*) '(x1,x2,x3) = (log10(nno),log10(npo),log10(u))'
      write (*,*)
      write (*,*) x_sol
      write (*,*)
      write (*,*) '(A1,A2,A3)'
      write (*,*)
      write (*,*) r(1:3)
      write (*,*)
      write (*,*) 'dA_i/dx_j analytical'
      write (*,*)
      write (*,*) jacobian(1:3,1)
      write (*,*) jacobian(1:3,2)
      write (*,*) jacobian(1:3,3)
      write (*,*)
      write (*,*) 'dA_i/dx_j numerical'
      write (*,*)
      write (*,*) dadu(1:3,1)
      write (*,*) dadu(1:3,2)
      write (*,*) dadu(1:3,3)
      write (*,*)
      write (*,*) 'dA_i/dx_j (analytical - numerical)'
      write (*,*)
      write (*,*) jacobian(1:3,1) - dadu(1:3,1)
      write (*,*) jacobian(1:3,2) - dadu(1:3,2)
      write (*,*) jacobian(1:3,3) - dadu(1:3,3)
      write (*,*)
    endif

    RETURN

  END SUBROUTINE Find_Non_Uniform_Solution

! checks on a grid which point for independent variables
!  Xp(1:3) = (log10(n_no),log10(n_po),log10(u))
!  gets closer to solving system of equations for non-uniform matter
  SUBROUTINE GRID ( x_inout, x_grid, lgcl_retry )

    IMPLICIT NONE

    REAL(DP), DIMENSION(3), INTENT(INOUT) :: x_inout
    LOGICAL(LGCL), INTENT(IN) :: lgcl_retry
    INTEGER(I4B) :: Xp_min(3), Xp_max(3), Xp_int
    INTEGER(I4B) :: nn, np, u, i, location(1:3)
    REAL(DP) :: denomn, denomp
    REAL(DP) :: Xp_real, Xp_fac(3), X_input(3)
    REAL(DP) :: residue(3), array_dummy(2)
    REAL(DP) :: log10_n_no, log10_n_po, log10_u
    ! used G instead of F below to differentiate from main subroutine
    REAL(DP) :: G_total, G_min
    REAL(DP) :: G_out, G_in, G_alpha, G_heavy, G_surf_coul
    ! used m instead of n below to differentiate from main subroutine
    REAL(DP) :: m_no, m_po, m_ni, m_pi, m_alpha, m_heavy
    REAL(DP) :: non_uniform_eq(3), non_uniform_jac(3,3)
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: temporary
    REAL(DP), DIMENSION(1:500,4), INTENT(OUT) :: x_grid
    LOGICAL(LGCL) :: sol_lgcl
    INTEGER(I4B) :: error, objective, flag
    INTEGER(I4B), PARAMETER :: n_dim = 3

!   set mimimum in free energy to a very high value
    G_min = 1.d100
    x_grid = ZERO
!   set flag to only compute free energies
    objective = 0
!   set some helpful variables
    denomp = 12.d0
    denomn = 12.d0 - 4.d0*24.d0*(min(Yp,HALF)-HALF)/11.d0
    IF (Yp<ONE/24.d0) THEN
      denomn  = MIN(FOUR/THREE/Yp/TWO,30.d0)
    ENDIF

    IF ( x_inout(3) /= ZERO ) THEN
      ! if initial guess for solution not zero
      ! then try points near initial guess
      Xp_min(1) =  998
      Xp_max(1) = 1002
      Xp_fac(1) = x_inout(1)/1.0d3
      Xp_min(2) =  998
      Xp_max(2) = 1002
      Xp_fac(2) = x_inout(2)/1.0d3
      Xp_min(3) =  998
      Xp_max(3) = 1002
      Xp_fac(3) = x_inout(3)/1.0d3
    ELSE
      ! if initial guess for solution not given
      ! then try points in range below. (should be enough)
      ! functions mostly determined by trial and error
      Xp_int = INT(60.D0 - 30.D0*LOG10(T) + 6.D0*LOG10(T)**TWO)
      Xp_min(1) = - MIN(Xp_int,int(2.d2*(MIN(Yp,HALF)/(HALF)))+1)
      Xp_max(1) =    0
      Xp_fac(1) = LOG10(n*(ONE-Yp))
      Xp_min(2) = - Xp_int
      Xp_max(2) = MIN(0,Xp_min(2)+80)
      Xp_fac(2) = LOG10(n*Yp)
      Xp_min(3) = - 100
      Xp_max(3) =     0
      Xp_fac(3) = 1.0d2
    ENDIF

!   if this is a retry then double number of grid points
!   for each independent variable to look for solution
    IF (lgcl_retry) THEN
      Xp_min = 2*Xp_min ; Xp_max = 2*Xp_max ; Xp_fac(3) = TWO*Xp_fac(3)
      DENOMN = TWO*DENOMN ; DENOMP = TWO*DENOMP
!     for very low proton fraction try even lower proton densities
      IF (Yp<0.05d0) Xp_min(2) = 2*Xp_min(2)
    ENDIF

!   set a temporary array to store values of
!    function we're trying to minimize
    ALLOCATE(Temporary(Xp_min(1):Xp_max(1),Xp_min(2):Xp_max(2),&
                       Xp_min(3):Xp_max(3),1:4))
    temporary = zero

!   from all initial guesses of 'x' in the grid
!   find the one with the lowest free energy
    ! IF (IS_TEST) THEN
    !   write (*,*)
    !   WRITE (*,"(A12,3es20.12)") 'x_inout = ', x_inout
    !   WRITE (*,"(3I5)") xp_min, xp_max, Xp_int
    !   WRITE (*,*) 'limits on log10(nn), log10(np) and log10(u)'
    !   WRITE (*,*) Xp_fac(1)-two**(-dble(Xp_min(1))/DENOMN)+one,&
    !               Xp_fac(1)-two**(-dble(Xp_max(1))/DENOMN)+one
    !   WRITE (*,*) Xp_fac(2)-two**(-dble(Xp_min(2))/DENOMP)+one,&
    !               Xp_fac(2)-two**(-dble(Xp_max(2))/DENOMP)+one
    !   Xp_real = -log10(Nuc_Sat_Dens)+log10(n)-2.2d0
    !   WRITE (*,*) Xp_real-2.5d0*dble(Xp_min(3))/Xp_fac(3), &
    !               Xp_real-2.5d0*dble(Xp_max(3))/Xp_fac(3), Nuc_Sat_Dens
    !   WRITE (*,*)
    ! ENDIF

    Xp_real = -log10(Nuc_Sat_Dens)+log10(n)-2.2d0

    IF ( x_inout(3) /= ZERO ) THEN
      DO u = Xp_min(3), Xp_max(3)
        log10_u = Xp_fac(3)*DBLE(u)
        DO np = Xp_min(2), Xp_max(2)
          log10_n_po = Xp_fac(2)*DBLE(np)
          DO nn = Xp_min(1), Xp_max(1)
            log10_n_no = Xp_fac(1)*DBLE(nn)
            CALL FREE_ENERGY(log10_n_no, log10_n_po, log10_u, n, T, Yp, &
              G_out, G_in, G_alpha, G_heavy, G_surf_coul, &
              m_no, m_po, m_ni, m_pi, m_alpha, m_heavy, &
              un_A, un_Z, un_rad, un_F, un_P, un_S, un_E, &
              un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT, &
              un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT, &
              un_mu_no, un_mu_po, un_Meff_no, un_Meff_po, &
              un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy, &
              un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy, &
              un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn, &
              un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT, &
              non_uniform_eq, non_uniform_jac, objective, error)
            IF (error==0) THEN
              G_total = (G_out + G_in + G_alpha + G_heavy + G_surf_coul)/n
            ELSE
              G_total = 1.d100
            ENDIF
            G_min = MIN(G_min,G_total)
            Temporary(nn,np,u,1:4) = (/log10_n_no,log10_n_po,log10_u,G_total/)
          ENDDO
        ENDDO
      ENDDO
    ELSE
      DO u = Xp_min(3), Xp_max(3)
        log10_u = Xp_real-2.5d0*dble(u)/Xp_fac(3)
        DO np = Xp_max(2), Xp_min(2), -1
          log10_n_po = Xp_fac(2)-TWO**(-DBLE(np)/DENOMP)+ONE
          DO nn = Xp_max(1), Xp_min(1), -1
            log10_n_no =  Xp_fac(1)-TWO**(-DBLE(nn)/DENOMN)+ONE
            CALL FREE_ENERGY(log10_n_no, log10_n_po, log10_u, n, T, Yp, &
              G_out, G_in, G_alpha, G_heavy, G_surf_coul, &
              m_no, m_po, m_ni, m_pi, m_alpha, m_heavy, &
              un_A, un_Z, un_rad, un_F, un_P, un_S, un_E, &
              un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT, &
              un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT, &
              un_mu_no, un_mu_po, un_Meff_no, un_Meff_po, &
              un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy, &
              un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy, &
              un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn, &
              un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT, &
              non_uniform_eq, non_uniform_jac, objective, error)
            IF (error==0) THEN
              G_total = (G_out + G_in + G_alpha + G_heavy + G_surf_coul)/n
            ELSE
              G_total = 1.d100
            ENDIF
            if (G_total<G_min) then
              G_min = G_total
            endif
            Temporary(nn,np,u,1:4) = (/log10_n_no,log10_n_po,log10_u,G_total/)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    objective = 0

!   order count_max points with lowest free energy
    DO i = 1, 500
      location(1:3) = MINLOC(Temporary(:,:,:,4))
      nn = location(1) - 1 + Xp_min(1)
      np = location(2) - 1 + Xp_min(2)
      u  = location(3) - 1 + Xp_min(3)
      x_grid(i,1:4) = Temporary(nn,np,u,1:4)
      log10_n_no = x_grid(i,1)
      log10_n_po = x_grid(i,2)
      log10_u    = x_grid(i,3)
!      CALL FREE_ENERGY(log10_n_no, log10_n_po, log10_u, n, T, Yp, &
!        G_out, G_in, G_alpha, G_heavy, G_surf_coul, &
!        m_no, m_po, m_ni, m_pi, m_alpha, m_heavy, &
!        un_A, un_Z, un_rad, un_F, un_P, un_S, un_E, &
!        un_DF_Dn, un_DF_Dy, un_DF_DT, un_DP_Dn, un_DP_Dy, un_DP_DT, &
!        un_DS_Dn, un_DS_Dy, un_DS_DT, un_DE_Dn, un_DE_Dy, un_DE_DT, &
!        un_mu_no, un_mu_po, un_Meff_no, un_Meff_po, &
!        un_Dmu_no_DT, un_Dmu_no_Dn, un_Dmu_no_Dy, &
!        un_Dmu_po_DT, un_Dmu_po_Dn, un_Dmu_po_Dy, &
!        un_dLog10_n_no_dn, un_dLog10_n_po_dn, un_dLog10_u_dn, &
!        un_dLog10_n_no_dT, un_dLog10_n_po_dT, un_dLog10_u_dT, &
!        non_uniform_eq, non_uniform_jac, objective, error)
        Temporary(nn,np,u,4) = 1.d100
    ENDDO

    DEALLOCATE(Temporary)
801 format (4es15.5,i4,5L2)
    RETURN

  END SUBROUTINE GRID

  SUBROUTINE SOLVE_NON_UNIFORM ( x, s_eq, m, uniform_flag )
!   given initial input for independent parameter 'x'
!   determine residue 's' for Eq that solves EoS of
!   uniform matter for neutron rich matter
    IMPLICIT NONE

    INTEGER (I4B), INTENT(IN) :: m, uniform_flag
    REAL(DP), INTENT(IN)  :: x(m)
    REAL(DP), INTENT(OUT) :: S_eq(m)

    REAL(DP) :: log10_n_no, log10_n_po
    REAL(DP) :: log10_u
    REAL(DP) :: F_out, F_in, F_alpha, F_heavy, F_surf_coul
    REAL(DP) :: n_no, n_po, dens_nucl_out
    REAL(DP) :: n_ni, n_pi
    REAL(DP) :: n_alpha, n_heavy
    REAL(DP) :: excluded_alpha_volume, B(m), S_jac(m,m)
    INTEGER(I4B) :: error

    log10_n_no = X(1)
    log10_n_po = X(2)
    log10_u = X(3)

    IF (ISNAN(log10_n_no) .OR. ISNAN(log10_n_po) .OR. ISNAN(log10_u)) THEN
      S_eq = 1.D100
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
          S_eq, S_jac, 1, error)

    IF (error<0) THEN
      S_eq = 1.D100
      RETURN
    ENDIF

    RETURN

  END SUBROUTINE SOLVE_NON_UNIFORM

  SUBROUTINE JACOB_SOLVE_NON_UNIFORM ( jac, ldr, x1, m )
!   given initial input for independent parameter 'x'
!   determine residue 's' for Eq that solves EoS of
!   uniform matter for neutron rich matter
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: m, ldr
    REAL(DP), INTENT(IN)  :: x1(m)
    REAL(DP), INTENT(OUT) :: jac(ldr,*)

    REAL(DP) :: X_p, log10_n_no, log10_n_po, log10_u
    REAL(DP) :: log10_n_n,log10_n_p,Temperature
    REAL(DP) :: F_out, F_in, F_alpha, F_heavy, F_surf_coul
    REAL(DP) :: n_no, n_po, dens_nucl_out
    REAL(DP) :: n_ni, n_pi
    REAL(DP) :: n_alpha, n_heavy
    REAL(DP) :: S_eq(3), S_jac(3,3)
    INTEGER(I4B) :: objective, error

    log10_n_no = x1(1)
    log10_n_po = x1(2)
    log10_u = x1(3)

    IF (ISNAN(log10_n_no) .OR. ISNAN(log10_n_po) .OR. ISNAN(log10_u)) THEN
      JAC(1:ldr,1:3) = ZERO
      RETURN
    ENDIF

    objective = 2

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
          S_eq, S_jac, objective, error)

    IF (error<0) THEN
      JAC(ldr,1:3) = ZERO
      RETURN
    ENDIF

    JAC(1:ldr,1:3) = S_JAC(1:ldr,1:3)

  END SUBROUTINE JACOB_SOLVE_NON_UNIFORM

END MODULE Find_NoN_Uniform_Solution_Mod
