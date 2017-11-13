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
  USE Surface_Properties_Mod, ONLY : T_table
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, R_5_3, THREE, TEN
  USE Minimize_surface_tension_Mod
  USE Surface_Equilibrium_Mod
  USE Skyrme_Bulk_Mod, ONLY : Skyrme_Bulk_Pressure
  USE OMP_LIB

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Nuclear_Surface_Tension(data_points,Surface_Tension_array)

!    USE Skyrme_Bulk_Mod, ONLY : SKYRME_BULK_PRESSURE

    IMPLICIT NONE

    INTEGER(I4B) :: i, j, k, thread, num_threads
    INTEGER(I4B) :: prot_frac_step, max_temperature
    INTEGER(I4B), INTENT(OUT) :: data_points
    REAL(DP) :: prot_frac, Temperature, Tcrit
    REAL(DP) :: x1(1), x2(2), x3(3), x4(4), residue(4)
    REAL(DP) :: tension_residue1(1), tension_residue3(3)
    REAL(DP), DIMENSION(200,4) :: xn, xt
    REAL(DP), DIMENSION(200,100,3) :: tension_array
    REAL(DP), DIMENSION(10000,3), INTENT(OUT) :: Surface_Tension_Array
    REAL(DP), DIMENSION(100) :: neutron_excess_array
    REAL(DP), DIMENSION(100) :: critical_temperature_array
    REAL(DP), DIMENSION(100) :: critical_density_array
    REAL(DP), DIMENSION(200,100,10) :: Surface_Properties
    REAL(DP) :: log10_dens_neut_in,  log10_dens_prot_in
    REAL(DP) :: log10_dens_neut_out, log10_dens_prot_out
    REAL(DP) :: pressure_in, pressure_out
    REAL(DP) :: tol, neut_tol, prot_tol
    REAL(DP) :: Surf_Tension
    REAL(DP) :: r_neut, r_prot, a_neut, a_prot
    LOGICAL(LGCL) :: surf_eq, surf_ten, error
    LOGICAL(LGCL) :: logical_array(200,100)

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

!   tabulate surface tension for later fitting
!   Surface_Tension_Array(:,1:3) = (prot_frac,temperature,surface_tension)
    Surface_Tension_Array(:,:) = zero
    Tension_Array(:,:,:) = zero
    Logical_Array(:,:) = .false.

!   open files where surface properties are output
!    open (20,file=trim(adjustl(output_directory))//"/surface_density.dat")
!    open (21,file=trim(adjustl(output_directory))//"/surface_tension.dat")
    open (20,file=trim(adjustl(output_directory))//"/surface_properties.dat")

!   set proton fraction from 0.50 down to minimum to be determined every 0.02
    prot_frac_step = 2
    tol = 1.d-2
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) DEFAULT(NONE) &
!$OMP FIRSTPRIVATE(prot_frac_step,T_table,tol) &
!$OMP PRIVATE(prot_frac,max_temperature,temperature) &
!$OMP PRIVATE(i,thread,surf_eq,surf_ten,x1,x3,x4,xn,xt) &
!$OMP PRIVATE(residue,tension_residue1,tension_residue3) &
!$OMP PRIVATE(log10_dens_neut_in, log10_dens_prot_in)  &
!$OMP PRIVATE(log10_dens_neut_out,log10_dens_prot_out) &
!$OMP PRIVATE(neut_tol,prot_tol,r_neut,r_prot,a_neut,a_prot) &
!$OMP PRIVATE(pressure_in,pressure_out,surf_tension,error) &
!$OMP SHARED(Tension_Array,Logical_Array,Surface_Properties)
    LOOP_Y: DO j = 100, 20, -prot_frac_step
      thread = omp_get_thread_num()
      prot_frac = DBLE(j)/200.d0
      max_temperature = INT(ten*T_table(j))
      write (*,*) thread, j, prot_frac, T_table(j)
      error = .false. 
      i = 0
      ! set temperature (T_min = 0.1 MeV, T_max = T_crit)
      LOOP_T: DO i = 1, max_temperature-1
        if (error) CYCLE
        Temperature = DBLE(i)/10.D0
        x4 = ZERO
        IF (i/=1) THEN
          IF (SUM(xn(i-1,1:4))/=zero) x4(1:4) = xn(i-1,1:4)
          ! If solution for both T(i-2) and T(i-1) exist
          ! make simple extrapolation to use as initial guess
          ! to find solutin at T(i).
          IF (i>2) THEN
            IF (SUM(xn(i-1,1:4))/=zero .AND. sum(xn(i-2,1:4))/=zero) &
              x4(1:4) = xn(i-1,1:4)+(xn(i-1,1:4)-xn(i-2,1:4))
          ENDIF
        ENDIF
        ! Determine inside/outside densities for
        ! pressure and chemical potential equilibria
        CALL FIND_SURFACE_EQUILIBRIUM (x4,prot_frac,Temperature,residue,surf_eq)
        xn(i,1:4) = x4(1:4)
        ! if solution has been found calculate skin thickness and
        ! haf density radii assuming a woods-saxon like distribution.
        log10_dens_neut_in  = X4(1)
        log10_dens_prot_in  = X4(2)
        log10_dens_neut_out = X4(3)
        log10_dens_prot_out = X4(4)
        write (2000+thread,"(2I4,4ES15.6,L4)") i, j, x4, surf_eq
        Surface_Properties(i,j,1:2) = X4(1:2)
        Surface_Properties(i,j,4:5) = X4(3:4)
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
        Surface_Properties(i,j,3:6:3) = (/pressure_in,pressure_out/)
!        IF (surf_eq) T_table(j) = Temperature
!        IF (surf_eq) WRITE (20,"(12ES15.6)") prot_frac, temperature, x4, &
!                residue, pressure_in, pressure_out
!        FLUSH(20)
        ! if densities are too similar assume matter
        ! uniform (no surface tension)
        IF (.NOT. surf_eq) error = .true.
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
            xt(i,1:3) = x3(1:3)
            IF (.not.surf_ten) error = .true.
            IF (error) CYCLE
          ELSE
          ! for neutron rich matter (j<100)
          ! initial guess for T_min
            x3(1:3) = zero
          ! if T > T_min try to use solution for lower temperature T(i-1)
          ! as initial guess to find solution for T(i)
            IF (i/=1) THEN
              IF (SUM(xt(i-1,1:3))/=zero) x3(1:3) = xt(i-1,1:3)
              ! If solution for both T(i-2) and T(i-1) exist
              ! make simple extrapolation to use as initial guess
              ! to find solutin at T(i).
              IF (i>2) THEN
                IF (SUM(xt(i-1,1:3))/=zero .AND. sum(xt(i-2,1:3))/=zero) &
                  x3(1:3) = xt(i-1,1:3)+(xt(i-1,1:3)-xt(i-2,1:3))
              ENDIF
            ENDIF
            CALL MINIMIZE_SURFACE_TENSION (3, x3, x4, &
                    prot_frac, Temperature, tension_residue3, surf_ten)
            IF (.not.surf_ten) error = .true.
            IF (error) CYCLE
            xt(i,1:3) = x3(1:3)
          ENDIF
          ! compute surface tension
          r_neut = x3(1) ; r_prot = zero ; a_neut = x3(2) ; a_prot = x3(3)
          Surf_Tension = Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                log10_dens_neut_in, log10_dens_prot_in,  &
                log10_dens_neut_out, log10_dens_prot_out, temperature)
          Surface_Properties(i,j,7:9) = x3(1:3)
          Surface_Properties(i,j,10)  = surf_tension
          write (3000+thread,"(2I12,3ES15.6,L4)") i, j, x3, surf_ten
!          WRITE(21,"( 9ES15.6,L2)") prot_frac, temperature, x3, &
!                                    tension_residue3, Surf_Tension, surf_ten
!          FLUSH(21)
!$OMP CRITICAL
          Tension_Array(i,j,1:3) = (/prot_frac,temperature,Surf_Tension/)
          Logical_Array(i,j)     = surf_ten
!$OMP END CRITICAL
        ENDIF
      ENDDO LOOP_T
      write (*,*) thread, j, prot_frac, T_table(j), 'end'
    ENDDO LOOP_Y
!$OMP END PARALLEL DO


    data_points = 0
    DO j = 100, 1, -prot_frac_step
      prot_frac = dble(j)/200.d0
      DO i = 1, 200
        temperature = dble(i)/10.d0
        IF (logical_array(i,j) .and. tension_array(i,j,3)>zero) THEN
          data_points = data_points + 1
          Surface_Tension_Array(data_points,1:3) = Tension_Array(i,j,1:3)
          ! write (*,*) Tension_Array(i,j,1:3)
          write (20,"(2ES9.2,10ES15.7)") prot_frac, temperature, &
                                       Surface_properties(i,j,1:10)
        ENDIF
      ENDDO
      WRITE (20,*)
    ENDDO

    CLOSE(20) !; CLOSE(21)


  END SUBROUTINE Nuclear_Surface_Tension

END MODULE Nuclear_Surface_Tension_Mod
