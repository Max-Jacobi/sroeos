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
MODULE Output_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL
  USE Physical_Constants_Mod, ONLY : ONE, TWO, TEN, HALF, v_alpha
  USE Global_Variables_Mod, only : IS_TEST
  USE Free_Energy_Mod
  USE Print_Test_Mod
  USE Main_output_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE GET_OUTPUT(n, T, Yp, x_u, x_nu, f_u, f_nu, &
                        sol_u, sol_nu, nu_is_true_sol  )

    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: n, Yp, T, x_u(1), x_nu(3), f_u, f_nu
    LOGICAL(LGCL), INTENT(IN) :: sol_u, sol_nu
    LOGICAL(LGCL), INTENT(OUT) :: nu_is_true_sol

    REAL(DP)      :: X_p, log10_n_no, log10_n_po, log10_u
    REAL(DP)      :: non_uniform_eq(3), non_uniform_jac(3,3)
    REAL(DP)      :: exc_v_alpha

    INTEGER(I4B)  :: objective, error

    nu_is_true_sol = .FALSE.
    objective = 0

    IF (sol_u .AND. sol_nu) THEN
      IF (IS_TEST) write (*,*) 'BOTH UNIFORM AND NON UNIFORM SOLUTION FOUND'
      IF (f_nu > f_u) THEN
        IF (IS_TEST) write (*,*) 'UNIFORM MATTER HAS LOWEST FREE ENERGY'
        X_p = TEN**x_u(1)
        objective = 3
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
      ELSE
        IF (IS_TEST) write (*,*) 'NON UNIFORM MATTER HAS LOWEST FREE ENERGY'
        nu_is_true_sol = .TRUE.
        objective = 4
        log10_n_no = x_nu(1)
        log10_n_po = x_nu(2)
        log10_u    = x_nu(3)
      ENDIF
      !CALL COMPARE_FREE_ENERGIES()
    ELSEIF (sol_u .AND. .not. sol_nu) THEN
      IF (IS_TEST) write (*,*) 'ONLY UNIFORM SOLUTION FOUND!'
      objective = 3
      X_p = TEN**x_u(1)
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
    ELSEIF (.not. sol_u .AND. sol_nu) THEN
      IF (IS_TEST) write (*,*) 'ONLY NON UNIFORM SOLUTION FOUND!'
      nu_is_true_sol = .TRUE.
      objective = 4
      log10_n_no = x_nu(1)
      log10_n_po = x_nu(2)
      log10_u    = x_nu(3)
    ELSEIF (.not. sol_u .AND. .not. sol_nu) THEN
      ! WRITE (filenumber2,"(3ES16.8,6ES20.12,2L2)") n, T, Yp, &
      !                      x_u, x_nu, f_u, f_nu, sol_u, sol_nu
      CONTINUE
    ENDIF

    CALL FREE_ENERGY( log10_n_no, log10_n_po, log10_u, n, T, Yp, &
                      F_o, F_i, F_alpha, F_TR, F_SC, &
                      n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
                      A_heavy, Z_heavy, rad, F, P, S, E, &
                      DF_Dn, DF_Dy, DF_DT, DP_Dn, DP_Dy, DP_DT, &
                      DS_Dn, DS_Dy, DS_DT, DE_Dn, DE_Dy, DE_DT, &
                      mu_no, mu_po, Meff_no, Meff_po, &
                      Dmu_no_DT, Dmu_no_Dn, Dmu_no_Dy, &
                      Dmu_po_DT, Dmu_po_Dn, Dmu_po_Dy, &
                      dLog10_n_no_dn, dLog10_n_po_dn, dLog10_u_dn, &
                      dLog10_n_no_dT, dLog10_n_po_dT, dLog10_u_dT, &
                      non_uniform_eq, non_uniform_jac, objective, error )

    ! WRITE (filenumber1,"(3ES16.8,6ES20.12,2L2)") n, T, Yp, &
    !                      x_u, x_nu, f_u, f_nu, sol_u, sol_nu

!   occupied volume
    u = TEN**Log10_u
!   heavy nuclei generalized radius
    r = rad
!   set u to zero if uniform atter is the solution
    IF (u < 1.d-100) u = zero
!   volume excluded by alpha particles
    exc_v_alpha = one - n_alpha*v_alpha
!   get number fractions (should add to one)
    x_no = (ONE-u)*exc_v_alpha*n_no / n
    x_po = (ONE-u)*exc_v_alpha*n_po / n
    x_alpha = (ONE-u)*FOUR*n_alpha / n
    x_heavy = (n_ni+n_pi)*ten**log10_u/n
    IF (x_heavy < 1.d-50) x_heavy = zero

  IF (IS_TEST) CALL PRINT_TEST( n, T, Yp, F_o, F_i, F_alpha, F_TR, F_SC, &
                   x_no, x_po, n_ni, n_pi, x_alpha, x_heavy, &
                   A_heavy, Z_heavy, rad, F, P, S, E, &
                   DF_Dn, DF_Dy, DF_DT, DP_Dn, DP_Dy, DP_DT, &
                   DS_Dn, DS_Dy, DS_DT, DE_Dn, DE_Dy, DE_DT, &
                   mu_no, mu_po, Meff_no, Meff_po, &
                   Dmu_no_DT, Dmu_no_Dn, Dmu_no_Dy, &
                   Dmu_po_DT, Dmu_po_Dn, Dmu_po_Dy)

  END SUBROUTINE GET_OUTPUT

END MODULE Output_Mod
