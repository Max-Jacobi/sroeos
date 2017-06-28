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
MODULE Nuclear_Surface_Tension_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP
  USE Make_Tables_Mod, ONLY : output_directory
  USE Polynomial_Fit_Mod
  USE Surface_Properties_Mod
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, R_5_3, THREE, TEN
  USE Minimize_surface_tension_Mod
  USE Surface_Equilibrium_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Nuclear_Surface_Tension(data_points,Surface_Tension_array)

    USE Skyrme_Bulk_Mod, ONLY : SKYRME_BULK_PRESSURE

    IMPLICIT NONE

    INTEGER(I4B) :: i, j, imax, flag, prot_frac_step, max_temperature
    INTEGER(I4B), INTENT(OUT) :: data_points
    REAL(DP) :: prot_frac, Temperature, Tcrit
    REAL(DP) :: x1(1), x2(2), x3(3), x4(4), residue(4)
    REAL(DP) :: tension_residue1(1), tension_residue3(3)
    REAL(DP), DIMENSION(200,100,4) :: xn, xt
    REAL(DP), DIMENSION(10000,3), INTENT(OUT) :: Surface_Tension_Array
    REAL(DP), DIMENSION(100) :: neutron_excess_array
    REAL(DP), DIMENSION(100) :: critical_temperature_array
    REAL(DP), DIMENSION(100) :: critical_density_array
    REAL(DP) :: log10_dens_neut_in,  log10_dens_prot_in
    REAL(DP) :: log10_dens_neut_out, log10_dens_prot_out
    REAL(DP) :: pressure_in, pressure_out
    REAL(DP) :: tol, neut_tol, prot_tol
    REAL(DP) :: Surf_Tension
    REAL(DP) :: r_neut, r_prot, a_neut, a_prot
    LOGICAL(LGCL) :: surf_eq, surf_ten

!   given temperature T and internal proton fraction xi
!   obtain inside and outside densities n_ni, n_pi, n_no, n_po
!   for which matter is in equilibrium. We chose:
!   T in the range  (0.10 to T_crit(xi)) in 0.01 MeV intervals and
!   xi in the range (0.10 to 0.50)       in 0.005 intervals.
    WRITE (*,*)
    WRITE (*,"('Fitting parameters (a,q,p) for surface tension s(x,T)')")
    WRITE (*,"(' for given nuclei proton fraction x and temperature T:')")
    WRITE (*,"('s(x,T) = sigma_s*h(T/Tc(x))*(2*2^a+q)/(x^-a+q+(1-x)^-a)')")
    WRITE (*,"('where h(T/Tc(x)) = (1+(T/Tc(x))²)^p')")
    WRITE (*,*)
    WRITE (*,"('THIS MAY TAKE A COUPLE OF MINUTES!!')")
    WRITE (*,*)
    WRITE (*,"('Values will be printed in OUTPUT_DIR/surface_properties.dat')")
    WRITE (*,*)

    data_points = 0  ! saves array size for surface tension fitting procedure
    x4 = zero ! array with neutron and proton inside and outside densities
    xn = zero ! array to save x4 solution for a given prot frac and temperature

!   tabulate surface tension for later fitting
!   Surface_Tension_Array(:,1:3) = (prot_frac,temperature,surface_tension)
    Surface_Tension_Array(:,:) = zero

!   open files where surface properties are output
    open (20,file=trim(adjustl(output_directory))//"/surface_density.dat")
    open (21,file=trim(adjustl(output_directory))//"/surface_tension.dat")
!   open (22,file=trim(adjustl(output_directory))//"/surface_properties.dat")

!   set proton fraction from 0.50 down to minimum to be determined every 0.02
    prot_frac_step = 2

    DO j = 100, 1, - prot_frac_step
      prot_frac = DBLE(j)/200.d0
      max_temperature = INT(ten*T_table(j))
      write (*,*) j, prot_frac,T_table(j)
      ! set temperature (T_min = 0.1 MeV, T_max T_crit)
      DO i = 1, max_temperature-1
        Temperature = DBLE(i)/10.D0
        x4 = ZERO
        ! if a solution has been found for lower T or higher Yp
        !  use it to find solution for new T, Yp
        IF (j < 100) x4(1:4) = xn(i,j+prot_frac_step,1:4)
        IF (i > 1) THEN
          IF (sum(xn(i-1,j,1:4))/=zero) x4(1:4) = xn(i-1,j,1:4)
        ENDIF
        ! Determine inside/outside densities for
        ! pressure and chemical potential equilibria
        !write (*,"(12ES14.6)") prot_frac, Temperature
        CALL FIND_SURFACE_EQUILIBRIUM (x4,prot_frac,Temperature,residue,surf_eq)
        !write (*,"(12ES14.6)") x4, residue
        xn(i,j,1:4) = x4(1:4)

        ! if solution has been found calculate skin thickness and
        ! haf density radii assuming a woods-saxon like distribution.
        log10_dens_neut_in  = X4(1)
        log10_dens_prot_in  = X4(2)
        log10_dens_neut_out = X4(3)
        log10_dens_prot_out = X4(4)
        ! check if high and low density phases have
        ! different proton and neutron densities
        neut_tol = ABS(log10_dens_neut_in-log10_dens_neut_out) &
                  /ABS(log10_dens_neut_in)
        IF (neut_tol<tol) surf_eq = .false.
        prot_tol = ABS(log10_dens_prot_in-log10_dens_prot_out) &
                  /ABS(log10_dens_prot_in)
        IF (prot_tol<tol) surf_eq = .false.
        pressure_in = SKYRME_BULK_PRESSURE &
                      (log10_dens_neut_in,log10_dens_prot_in,Temperature)
        pressure_out = SKYRME_BULK_PRESSURE &
                      (log10_dens_neut_out,log10_dens_prot_out,Temperature)
        IF (surf_eq) T_table(j) = Temperature
        IF (surf_eq) WRITE (20,"(12ES15.6)") prot_frac, temperature, x4, &
                residue, pressure_in, pressure_out
        FLUSH(20)
        ! if densities are too similar assume matter
        ! uniform (no surface tension)
        IF (.NOT. surf_eq) CYCLE
!       otherwise compute surface tension
        IF (surf_eq) THEN
          ! given the solution above minimize
          ! surface thermodynamic potential
          ! for symmetric nuclear matter (j=100)
          IF (j==100) THEN
            IF (i==1) x1 = zero
            CALL MINIMIZE_SURFACE_TENSION (1, x1, x4, &
                      prot_frac, Temperature, tension_residue1, surf_ten)
            x3(1:3) = (/zero,x1(1),x1(1)/)
            tension_residue3 = (/zero,tension_residue1(1),tension_residue1(1)/)
            xt(i,j,1:3) = x3(1:3)
          ELSE
          ! for neutron rich matter (j<100)
          ! initial guess for T_min
            x3(1:3) = zero
          ! if T = T_min use solution for proton fraction prot_frac(j)
          !  to find solution for proton fraction prot_frac(j-prot_frac_step)
          IF (i==1) THEN
            IF (j < 100) x3(1:3) = xt(i,j+prot_frac_step,1:3)
            IF (j < 100-prot_frac_step) x3(1:3) = xt(i,j+prot_frac_step,1:3) &
                  - (xt(i,j+2*prot_frac_step,1:3) - xt(i,j+prot_frac_step,1:3))
          ENDIF
          ! if T > T_min try to use solution for lower temperature T(i-1)
          ! as initial guess to find solution for T(i)
          IF (i/=1) THEN
            IF (SUM(xt(i-1,j,1:3))/=zero) x3(1:3) = xt(i-1,j,1:3)
            ! If solution for both T(i-2) and T(i-1) exist
            ! make simple extrapolation to use as initial guess
            ! to find solutin at T(i).
            IF (i>2) THEN
              IF (SUM(xt(i-1,j,1:3))/=zero .AND. sum(xt(i-2,j,1:3))/=zero) &
                x3(1:3) = xt(i-1,j,1:3) + (xt(i-1,j,1:3) - xt(i-2,j,1:3))
              ENDIF
            ENDIF
            CALL MINIMIZE_SURFACE_TENSION (3, x3, x4, &
                      prot_frac, Temperature, tension_residue3, surf_ten)
            xt(i,j,1:3) = x3(1:3)
          ENDIF
          ! compute surface tension
          r_neut = x3(1) ; r_prot = zero ; a_neut = x3(2) ; a_prot = x3(3)
          Surf_Tension = Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                log10_dens_neut_in, log10_dens_prot_in,  &
                log10_dens_neut_out, log10_dens_prot_out, temperature)
          WRITE(21,"( 9ES15.6,L2)") prot_frac, temperature, x3, &
                                    tension_residue3, Surf_Tension, surf_ten
          FLUSH(21)
          ! log10_dens_neut_in  = X4(1)
          ! log10_dens_prot_in  = X4(2)
          ! log10_dens_neut_out = X4(3)
          ! log10_dens_prot_out = X4(4)
          !WRITE(22,"(16ES15.6,L2)") prot_frac, temperature, x4, r_neut, r_prot,&
          !                        a_neut, a_prot, Surf_Tension, tension_residue3
          !FLUSH(22)
          IF (surf_ten) THEN
            data_points = data_points + 1
            Surface_Tension_Array(data_points,1:3) = &
                                      (/prot_frac,temperature,Surf_Tension/)
          ENDIF
        ENDIF
      ENDDO
      WRITE (20,*) ; WRITE (21,*) !; WRITE (22,*) 
    ENDDO
    CLOSE(20) ; CLOSE(21)

  END SUBROUTINE Nuclear_Surface_Tension

END MODULE Nuclear_Surface_Tension_Mod
