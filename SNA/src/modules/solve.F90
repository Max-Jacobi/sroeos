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
MODULE Solve_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL
  USE Physical_Constants_Mod, ONLY : ZERO, TEN
  USE Phase_Space_Input_Mod
  USE Find_Uniform_Solution_Mod
  USE Find_Non_Uniform_Solution_Mod
  USE Output_Mod
  USE Main_output_Mod
  USE Make_Tables_Mod
  USE Save_to_Table_Mod
#ifdef _OPENMP
  USE OMP_LIB
#endif
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE SOLVE

    IMPLICIT NONE

    REAL(DP) :: UNIFORM_SOL(1), NON_UNIFORM_SOL(3)
    REAL(DP) :: UNIFORM_SOL_GUESS(1), NON_UNIFORM_SOL_GUESS(3)
    LOGICAL(LGCL) :: UNIFORM_GUESS, NON_UNIFORM_GUESS
    REAL(DP) :: n, Yp, T
    REAL(DP) :: F_uniform, F_non_uniform, uniform_residue, non_uniform_residue
    REAL(DP) :: dble_n, increment_n, Delta_n, Delta_T
    REAL(DP) :: t_transition, n_transition, n_transition_min, n_transition_max

    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: x_nu_guess_n, x_nu_guess_t
    LOGICAL(LGCL), DIMENSION(:,:), ALLOCATABLE :: u_solution, nu_solution

    INTEGER(I4B) :: i_Yp, i_n, i_T, count_u, count_nu
    INTEGER(I4B) :: flag, thread

    LOGICAL(LGCL) :: uniform_solution, non_uniform_solution, retry
    LOGICAL(LGCL) :: check_uniform, check_non_uniform, check_non_uniform_flag
    LOGICAL(LGCL) :: nu_flag, nu_true
    LOGICAL(LGCL) :: stop_check_uniform
    LOGICAL(LGCL) :: T_LARGER_THAN_Tcrit, n_LARGER_THAN_nmax, n_SMALLER_THAN_nmin
    LOGICAL(LGCL) :: non_uniform_sol_found_for_T

    CHARACTER(LEN=256) :: filename1, filename2, string_Yp, outdir
    INTEGER(I4B) :: filenumber1, filenumber2, nthreads

    REAL(DP) :: wtime
    REAL(SP) :: xtime
    INTEGER(I4B) :: hours, minutes, seconds

    ! open file where output solutions will be written
    !  n, T, Ye, log10(n_no), log10(n_po), log10(u)
    outdir = output_directory

#ifdef _OPENMP    
    nthreads = omp_get_max_threads()
#else
    nthreads = 1
#endif

    write(6,*) "OpenMP parallization: density-temperature slices at fixed proton fraction."
    write(6,"(A11,I4,A7)") " There are ",Yp_fin," slices."
    if (nthreads > Yp_fin) then 
       write(6,*) "Requesting less proton fraction slices than OpenMP threads. This will fail!"
       write(6,*) "Aborting!"
       stop
    endif
 
    open (9,file=trim(adjustl(output_directory))//"/timer.dat")

!$OMP PARALLEL DO SCHEDULE (DYNAMIC,1) DEFAULT(NONE) &
!$OMP FIRSTPRIVATE(Yp_ini,Yp_fin,Yp_step,Yp_min,outdir,write_solutions_to_file) &
!$OMP FIRSTPRIVATE(n_ini,n_fin,steps_per_decade_in_n,Log10n_min) &
!$OMP FIRSTPRIVATE(T_ini,T_fin,steps_per_decade_in_T,Log10T_min) &
!$OMP PRIVATE(thread,i_T,i_n,i_Yp,dble_n,increment_n,count_u,count_nu) &
!$OMP PRIVATE(n,Yp,T,Delta_n,Delta_T,t_transition,non_uniform_sol_found_for_T) &
!$OMP PRIVATE(n_transition,n_transition_max,n_transition_min) &
!$OMP PRIVATE(UNIFORM_SOL,NON_UNIFORM_SOL) &
!$OMP PRIVATE(F_uniform,F_non_uniform,uniform_residue,non_uniform_residue) &
!$OMP PRIVATE(UNIFORM_SOL_GUESS,NON_UNIFORM_SOL_GUESS) &
!$OMP PRIVATE(UNIFORM_GUESS,NON_UNIFORM_GUESS) &
!$OMP PRIVATE(check_uniform,check_non_uniform) &
!$OMP PRIVATE(check_non_uniform_flag,stop_check_uniform) &
!$OMP PRIVATE(T_LARGER_THAN_Tcrit, n_LARGER_THAN_nmax, n_SMALLER_THAN_nmin) &
!$OMP PRIVATE(uniform_solution,non_uniform_solution,retry,nu_flag,nu_true) &
!$OMP PRIVATE(x_nu_guess_n,x_nu_guess_T,u_solution,nu_solution) &
!$OMP PRIVATE(filename1,filename2,filenumber1,filenumber2,string_Yp) &
!$OMP PRIVATE(wtime,hours,minutes,seconds)
! loop over proton fractions
    DO i_Yp = Yp_fin, Yp_ini, -1
      ! start timer
#ifdef _OPENMP 
      wtime = omp_get_wtime ( )
#else
      call cpu_time ( xtime )
#endif

      ! set proton fraction
      Yp = Yp_min + dble(i_Yp-1)*Yp_step
      ! n_transition is the maximum density where we expect a transition
      ! from uniform matter to non-uniform matter to occur.
      ! this value is updated every i_T loop and used at a lower temperature
      n_transition_max = 1.d-1
      n_transition_min = 2.d-2
!       once temperature is very low and density lower than transition to
!       non-uniform matter stop looking for uniform solution
      stop_check_uniform = .FALSE.
      ! t_transition is an approximate temperature for two-phase coexistence
      ! slightly overestimated value is calculated by Critical_Temperature(Yp)
      ! TODO: CHECK FUNCTION SAFETY
      t_transition = Critical_Temperature(Yp)
      ! allocate arrays used as initial guesses for EoS
      allocate(x_nu_guess_T(n_ini-1:n_fin,t_ini-1:t_fin,1:3))
      allocate(x_nu_guess_n(n_ini-1:n_fin,t_ini-1:t_fin,1:3))
      allocate( u_solution(n_ini:n_fin+1,t_ini:t_fin+1))
      allocate(nu_solution(n_ini:n_fin+1,t_ini:t_fin+1))
      ! values stored to guess solution based on higher temperature/density
      x_nu_guess_t = zero
      x_nu_guess_n = zero
      u_solution = .FALSE.
      nu_solution = .FALSE.
!     Set file names for output.
!      TODO: Adjust filename for y < 0.001
      WRITE (string_Yp,"(1F5.3)") Yp
#ifdef _OPENMP
      thread = omp_get_thread_num()
#else
      thread = 1
#endif
      IF (write_solutions_to_file) THEN
        filename1 = TRIM(adjustl(outdir))//"/SOL/"//trim(adjustl(string_Yp))
        filenumber1 = 2000*(thread+1)+11
        filename2 = TRIM(adjustl(outdir))//"/NO_SOL/"//trim(adjustl(string_Yp))
        filenumber2 = 2000*(thread+1)+12
!     Open ASCII files to output solutions and for when no solution is found
!     TODO: Add something to check for NaNs.
!     TODO: Add flag whether want to print solutions as it may use to much space.
        OPEN(filenumber1,file=filename1,status='replace')
        OPEN(filenumber2,file=filename2,status='replace')
      ENDIF
!$OMP CRITICAL
      WRITE (*,"('Thread',1i4,' computing loop ',1i4,' of ',1i4, ' for y = ',1F8.6)") &
          thread, i_Yp, Yp_fin, Yp
!$OMP END CRITICAL
!     loop in Temperature from high to low
      DO i_T = T_fin, T_ini, -1
        ! TODO: make this a function
        T = TEN**(Log10T_min + dble(i_T-1)/dble(steps_per_decade_in_T))
        Delta_T = T*(ONE - TEN**(ONE/dble(steps_per_decade_in_T)))
!       set flag to check uniform solution to true
        check_uniform = .TRUE.
!       set initial uniform solution guess to false/zero
        UNIFORM_GUESS = .FALSE.
        UNIFORM_SOL_GUESS = ZERO
        non_uniform_solution = .FALSE.
!       set initial non-uniform solution guess to false/zero
        NON_UNIFORM_GUESS = .FALSE.
        NON_UNIFORM_SOL_GUESS = ZERO
        uniform_solution = .FALSE.
!       check if T > T_crit estimated for each Yp above
        IF (T > t_transition) THEN
          T_LARGER_THAN_Tcrit = .TRUE.
        ELSE
          T_LARGER_THAN_Tcrit = .FALSE.
        ENDIF
!       set flag of non-uniform solution found for Yp, T set to false
        non_uniform_sol_found_for_T = .FALSE.
!       loop in density from high to low
        DO i_n = n_fin, n_ini, -1
          ! set debug counters to zero
          count_nu = 0 ; count_u = 0
          ! TODO: make this a function
          n = TEN**(Log10n_min+dble(i_n-1)/dble(steps_per_decade_in_n))
          Delta_n = n*(ONE-TEN**(ONE/dble(steps_per_decade_in_n)))
!------------------------------------------------------------------------------!
!         look for uniform matter solution
!         TODO: make this a subroutine
!------------------------------------------------------------------------------!
          uniform_solution = .FALSE.
          uniform_residue = 1.d10
          F_uniform = 1.d100
!         check if temperature and density are low enough to skip
!         looking for uniform solution, which may be time consuming at low T
          IF (non_uniform_sol_found_for_T .AND. stop_check_uniform) &
            check_uniform=.FALSE.
          ! write (*,*) Yp, T, n, check_uniform
!         if initial guess provided, use it
          IF (check_uniform .AND. UNIFORM_GUESS) THEN
            count_u = count_u + 1
            RETRY = .FALSE.
            CALL Find_Uniform_Solution(n,T,Yp,UNIFORM_SOL_GUESS,UNIFORM_SOL,&
                              F_uniform,uniform_residue,retry,uniform_solution)
            ! if solution found use as guess for next density
            IF (uniform_solution) THEN
              UNIFORM_GUESS = .TRUE.
              UNIFORM_SOL_GUESS = UNIFORM_SOL
            ENDIF
          ENDIF
!         if no initial guess or solution not found
          IF (check_uniform .AND. .NOT. uniform_solution) THEN
            count_u = count_u + 1
            UNIFORM_GUESS = .FALSE.
            UNIFORM_SOL_GUESS = ZERO
            RETRY = .FALSE.
            CALL Find_Uniform_Solution(n,T,Yp,UNIFORM_SOL_GUESS,UNIFORM_SOL,&
                              F_uniform,uniform_residue,retry,uniform_solution)
            ! if no solution found double grid points
            ! TODO: add flag to only do this for lower densities
            IF (.NOT. uniform_solution) THEN
              count_u = count_u + 1
              RETRY = .TRUE.
              CALL Find_Uniform_Solution(n,T,Yp,UNIFORM_SOL_GUESS,UNIFORM_SOL,&
                              F_uniform,uniform_residue,retry,uniform_solution)
            ENDIF
            ! if solution found use as guess for next density
            IF (uniform_solution) THEN
              UNIFORM_GUESS = .TRUE.
              UNIFORM_SOL_GUESS = UNIFORM_SOL
            ENDIF
          ENDIF
          IF (uniform_solution) u_solution(i_n,i_T) = .TRUE.
!------------------------------------------------------------------------------!
!       look for non-uniform matter solution
!       TODO: make this a subroutine
!------------------------------------------------------------------------------!
          non_uniform_solution = .FALSE.
          non_uniform_sol = ZERO
          non_uniform_residue = 1.d10
          F_non_uniform = 1.d99
!         check if n is between critical densities
!         TODO: make this subroutine
          IF (n > n_transition_max) THEN
            n_LARGER_THAN_nmax = .TRUE.
          ELSE
            n_LARGER_THAN_nmax = .FALSE.
          ENDIF
          IF (n < n_transition_min) THEN
            n_SMALLER_THAN_nmin = .TRUE.
          ELSE
            n_SMALLER_THAN_nmin = .FALSE.
          ENDIF
!         check whether or not to look for uniform solution
          check_non_uniform = .TRUE.
          IF (T_LARGER_THAN_Tcrit) check_non_uniform = .FALSE.
          IF (n_LARGER_THAN_nmax ) check_non_uniform = .FALSE.
          IF (n_SMALLER_THAN_nmin) check_non_uniform = .FALSE.
          check_non_uniform_flag = check_non_uniform
!         if a solution was found at a higher density use it as initial guess
          IF (check_non_uniform .AND. abs(x_nu_guess_n(i_n,i_t,3))>zero) THEN
            ! write (11,*) Yp, T, n, x_nu_guess_n(i_n,i_T,1:3), 'density'
            count_nu = count_nu + 1
            RETRY = .FALSE.
            NON_UNIFORM_GUESS = .TRUE.
            NON_UNIFORM_SOL_GUESS(1:3) = x_nu_guess_n(i_n,i_t,1:3)
            CALL Find_Non_Uniform_Solution(n,T,Yp,NON_UNIFORM_SOL_GUESS,&
                                           NON_UNIFORM_SOL,F_non_uniform,&
                                           non_uniform_residue,retry,&
                                           non_uniform_solution)
          ENDIF
!         if solution found at a higher temperature use it as initial guess
          IF (non_uniform_solution) check_non_uniform = .FALSE.
          IF (check_non_uniform .AND. abs(x_nu_guess_T(i_n,i_t,3))>zero) THEN
            ! write (11,*) Yp, T, n, x_nu_guess_t(i_n,i_T,1:3), 'Temperature'
            count_nu = count_nu + 2
            RETRY = .FALSE.
            NON_UNIFORM_GUESS = .TRUE.
            NON_UNIFORM_SOL_GUESS(1:3) = x_nu_guess_T(i_n,i_t,1:3)
            CALL Find_Non_Uniform_Solution(n,T,Yp,NON_UNIFORM_SOL_GUESS,&
                                           NON_UNIFORM_SOL,F_non_uniform,&
                                           non_uniform_residue,retry,&
                                           non_uniform_solution)
          ENDIF
!         if no guesses above then try to find solution from scratch
          IF (non_uniform_solution) check_non_uniform = .FALSE.
          IF (check_non_uniform .AND. .NOT. non_uniform_solution) THEN
            count_nu = count_nu + 4
            NON_UNIFORM_GUESS = .FALSE.
            NON_UNIFORM_SOL_GUESS = ZERO
            RETRY = .FALSE.
            CALL Find_Non_Uniform_Solution(n,T,Yp,NON_UNIFORM_SOL_GUESS,&
                                           NON_UNIFORM_SOL,F_non_uniform,&
                                           non_uniform_residue,retry,&
                                           non_uniform_solution)
            ! If no solution found and one is expected, try again with twice the
            ! number of points for each independent variable
            nu_flag = .FALSE.
            IF (nu_solution(i_n  ,i_T+1)) nu_flag = .TRUE.
            IF (nu_solution(i_n+1,i_T+1)) nu_flag = .TRUE.
            IF (.NOT.u_solution(i_n,i_T)) nu_flag = .TRUE.
            IF (.NOT. non_uniform_solution .AND. nu_flag) THEN
              count_nu = count_nu + 8
              RETRY = .TRUE.
              CALL Find_Non_Uniform_Solution(n,T,Yp,NON_UNIFORM_SOL_GUESS,&
                                             NON_UNIFORM_SOL,F_non_uniform,&
                                             non_uniform_residue,retry,&
                                             non_uniform_solution)
            ENDIF
          ENDIF
          IF (non_uniform_solution) nu_solution(i_n,i_T) = .TRUE.
!         get output to save to table
!         last six values are used to udate guess to find solutions.
          CALL GET_OUTPUT(n, T, Yp, uniform_sol, non_uniform_sol, F_uniform,   &
                        F_non_uniform, uniform_solution, non_uniform_solution, &
                        nu_true )

!         If a non-uniform solution is found update values for n_transition.
!         n_transition_* determines range to look for a solution at a given n
          IF (non_uniform_solution) THEN
!           n_transition_max only affects the next (lower) T to be calculated.
!           n_transition_max is increased slightly if T > 1 MeV up to n_max
!           and kept constant otherwise
            n_transition = n
            dble_n = dble(i_n)
            increment_n = dble(steps_per_decade_in_n)
            IF (.NOT. non_uniform_sol_found_for_T) THEN
              n_transition = n
              IF (T>=1.0d0) n_transition = n*TEN**(1.D0/increment_n)
              IF (T> 3.0d0) n_transition = n*TEN**(2.D0/increment_n)
              IF (T> 5.0d0) n_transition = n*TEN**(3.D0/increment_n)
              n_transition_max = min(n_transition,0.12d0)
            ENDIF
!           n_transition_min determines whether to try to check for solution at
!           lower density. The choice made here ensures the code tried to find a
!           solution down to 0.01 fm^-3, if no solution found at higher density
!           or down to another two points in density phase space if a solution
!           is found at n > 0.01 fm^-3.
            n_transition_min = n/TEN**(TWO/increment_n)
            non_uniform_sol_found_for_T = .TRUE.
            IF (n_transition_min<TEN**Log10n_min) stop_check_uniform = .TRUE.
            ! write (*,*) Yp, T, n, n_transition_min, n_transition_max, &
            ! uniform_solution, non_uniform_solution, F_uniform, F_non_uniform
          ENDIF
!------------------------------------------------------------------------------!
!         write solutions to output file
          IF (write_solutions_to_file) &
                WRITE (filenumber1,"(3ES15.8,8ES16.8,I2,I3,6L2)") n, T, Yp, &
                uniform_sol, non_uniform_sol, F_uniform, F_non_uniform,  &
                uniform_residue, non_uniform_residue, count_u, count_nu, &
                uniform_solution, non_uniform_solution, &
                check_uniform, check_non_uniform_flag 

!         check for points where no solution was found
!         TODO: save these points and interpolate a solution later
          IF ( .NOT. uniform_solution .AND. .NOT. non_uniform_solution) THEN
            IF (write_solutions_to_file) WRITE (filenumber2,"(5ES20.12,8L2)") &
             n, T, Yp, uniform_residue, non_uniform_residue, uniform_solution,&
             non_uniform_solution, check_uniform, check_non_uniform_flag
          ENDIF

          IF (nu_true) THEN
!           save to use as a guess for non-uniform solution at a lower T
            x_nu_guess_T(i_n,i_t-1,1)=NON_UNIFORM_SOL(1)+dLog10_n_no_dT*delta_T
            x_nu_guess_T(i_n,i_t-1,2)=NON_UNIFORM_SOL(2)+dLog10_n_po_dT*delta_T
            x_nu_guess_T(i_n,i_t-1,3)=NON_UNIFORM_SOL(3)+dLog10_u_dT*delta_T
!           save to use as a guess for non-uniform solution at a lower n
            x_nu_guess_n(i_n-1,i_t,1)=NON_UNIFORM_SOL(1)+dLog10_n_no_dn*delta_n
            x_nu_guess_n(i_n-1,i_t,2)=NON_UNIFORM_SOL(2)+dLog10_n_po_dn*delta_n
            x_nu_guess_n(i_n-1,i_t,3)=NON_UNIFORM_SOL(3)+dLog10_u_dn*delta_n
          ELSE
            x_nu_guess_t(i_n,i_T-1,:) = zero
            x_nu_guess_n(i_n-1,i_T,:) = zero
          ENDIF
          CALL SAVE_TO_TABLE (i_n, i_T, i_Yp, n, Yp, T)

          IF (write_solutions_to_file) THEN 
            FLUSH(filenumber1) ; FLUSH(filenumber2)
          ENDIF
        ENDDO
      ENDDO
      IF (write_solutions_to_file) THEN 
        CLOSE(filenumber1) ; CLOSE(filenumber2)
      ENDIF
      deallocate(x_nu_guess_T,x_nu_guess_n,u_solution,nu_solution)
      ! start timer
#ifdef _OPENMP 
      wtime = omp_get_wtime ( ) - wtime
      hours   = int(wtime)/3600
      minutes = mod(int(wtime),3600)/60
      seconds = mod(mod(int(wtime),3600),60)
!$OMP CRITICAL
      WRITE (9,"('Thread',1i4,' computed EOS for y = ',1F8.6, ' in ',1I4,':',1I2.2,':',1I2.2)") &
          thread, Yp, hours, minutes, seconds
!$OMP END CRITICAL
#else
      call cpu_time ( xtime )
      hours   = int(xtime)/3600
      minutes = mod(int(xtime),3600)/60
      seconds = mod(mod(int(xtime),3600),60)
      WRITE (9,"('Thread',1i4,' computed EOS for y = ',1F8.6, ' in',1I3,':',1I2.2,':',1I2.2)") 
          thread, Yp, hours, minutes, seconds
#endif
    ENDDO
!$OMP end parallel do
    CLOSE(9)

  END SUBROUTINE SOLVE

END MODULE Solve_Mod
