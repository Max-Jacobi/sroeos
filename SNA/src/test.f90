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

PROGRAM Test
  ! Use modules
  USE Kind_Types_Mod,       only : I4B, DP
  USE Global_Variables_Mod, only : IS_TEST
  USE Read_Input_Mod,       only : READ_TEST_INPUT
  USE Test_Input_Mod
  USE Read_Skyrme_Coefficients_Mod
  USE Determine_Nuclear_Properties_Mod
  USE Determine_Surface_Properties_Mod
  USE Surface_Observable_Mod
  USE Find_Uniform_Solution_Mod
  USE Find_Non_Uniform_Solution_Mod
  USE Output_Mod
  USE Main_output_Mod
  USE phase_space_point_mod
  ! No implicit typing
  IMPLICIT NONE
  ! Declare variables
  REAL(DP)      :: UNIFORM_SOL(1), NON_UNIFORM_SOL(3)
  REAL(DP)      :: F_uniform, F_non_uniform
  REAL(DP)      :: uniform_residue, non_uniform_residue
  ! REAL(DP)      :: n, Yp, T

  LOGICAL(LGCL) :: uniform_solution, non_uniform_solution, retry
  LOGICAL(LGCL) :: non_uniform_is_true_solution

  IS_TEST = .TRUE.

  WRITE (*,*)
  WRITE (*,*) 'Step 1: Read in Skyrme Coefficients.'
  WRITE (*,*)
  CALL READ_SKYRME_COEFFICIENTS

  WRITE (*,*) 'Step 2: Determine nuclear matter properties.'
  WRITE (*,*)
  CALL DETERMINE_NUCLEAR_PROPERTIES

  WRITE (*,*) 'Step 3: Determine Surface properties.'
  WRITE (*,*) '        Could take ~5 minutes if surface_fit = .true.!'
  WRITE (*,*)
  CALL DETERMINE_SURFACE_PROPERTIES
  CALL SURFACE_OBSERVABLES

  WRITE (*,*) 'Step 4: Read test input.'
  WRITE (*,*)
  CALL READ_TEST_INPUT

  n  = density
  T  = Temperature
  Yp = Proton_Fraction

  IF (.NOT.GUESS) THEN
    UNIFORM_SOL_GUESS = ZERO
    NON_UNIFORM_SOL_GUESS = ZERO
  ENDIF

! NO SOLUTION FOUND YET
  uniform_solution = .FALSE.
  ! SOLVE EOS FOR UNIFORM MATTER IF INITIAL GUESS GIVEN
  IF (GUESS) THEN
    RETRY = .FALSE.
    CALL Find_Uniform_Solution(n,T,Yp,UNIFORM_SOL_GUESS,UNIFORM_SOL,&
                              F_uniform,uniform_residue,retry,uniform_solution)
  ENDIF
! IF GUESS NOT GIVEN OR NOT HELPFUL,
!  TRY TO FIND SOLUTION WITHOUT INITIAL GUESS
  IF (.NOT. uniform_solution) THEN
    UNIFORM_SOL_GUESS = ZERO
    RETRY = .FALSE.
    CALL Find_Uniform_Solution(n,T,Yp,UNIFORM_SOL_GUESS,UNIFORM_SOL,&
                              F_uniform,uniform_residue,retry,uniform_solution)
    ! IF STILL NO SOLUTION FOUND DOUBLE ARRAY RANGE
    IF (.NOT. uniform_solution) THEN
      RETRY = .TRUE.
      CALL Find_Uniform_Solution(n,T,Yp,UNIFORM_SOL_GUESS,UNIFORM_SOL,&
                                F_uniform,uniform_residue,retry,uniform_solution)
    ENDIF
  ENDIF

  ! write (*,*) n,T,Yp,UNIFORM_SOL_GUESS,UNIFORM_SOL,&
  ! F_uniform,uniform_residue,retry,uniform_solution

! NO SOLUTION FOUND YET
  non_uniform_solution = .FALSE.
  F_non_uniform = 1.d100
  ! SOLVE EOS FOR UNIFORM MATTER IF INITIAL GUESS GIVEN
  IF (GUESS) THEN
    RETRY = .FALSE.
    CALL Find_Non_Uniform_Solution(n,T,Yp,NON_UNIFORM_SOL_GUESS,&
                                   NON_UNIFORM_SOL,F_non_uniform,&
                                   non_uniform_residue,retry,&
                                   non_uniform_solution)
  ENDIF
! IF GUESS NOT GIVEN OR NOT HELPFUL,
!  TRY TO FIND SOLUTION WITHOUT INITIAL GUESS
  IF (.NOT. non_uniform_solution) THEN
    NON_UNIFORM_SOL_GUESS = ZERO
    RETRY = .FALSE.
    CALL Find_Non_Uniform_Solution(n,T,Yp,NON_UNIFORM_SOL_GUESS,&
                                   NON_UNIFORM_SOL,F_non_uniform,&
                                   non_uniform_residue,retry,&
                                   non_uniform_solution)
    ! IF STILL NO SOLUTION FOUND DOUBLE ARRAY RANGE
    IF (.NOT. non_uniform_solution) THEN
      RETRY = .TRUE.
      CALL Find_Non_Uniform_Solution(n,T,Yp,NON_UNIFORM_SOL_GUESS,&
                                     NON_UNIFORM_SOL,F_non_uniform,&
                                     non_uniform_residue,retry,&
                                     non_uniform_solution)
    ENDIF
  ENDIF

  CALL GET_OUTPUT (n,T,Yp,UNIFORM_SOL,NON_UNIFORM_SOL,F_uniform,F_non_uniform,&
                   uniform_solution,non_uniform_solution,       &
                   non_uniform_is_true_solution )

END PROGRAM Test
