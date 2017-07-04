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
MODULE Read_Input_Mod

  ! Use modules
  USE Kind_Types_Mod, ONLY : I4B, DP
  USE Physical_Constants_Mod, ONLY : ZERO, HALF

  ! No implicit typing
  IMPLICIT NONE

CONTAINS

  SUBROUTINE READ_TEST_INPUT

    USE Test_Input_Mod

    IMPLICIT NONE

    REAL(DP), DIMENSION(1) :: X1
    REAL(DP), DIMENSION(3) :: X3

    NAMELIST /INPUT_TEST_LIST/ &
                 DENSITY,TEMPERATURE,PROTON_FRACTION,&
                 UNIFORM_SOL_GUESS,NON_UNIFORM_SOL_GUESS,GUESS

    GUESS = .FALSE.
    X1 = 0.D0
    X3 = 0.D0

    OPEN(10,FILE='input/test.in',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
    READ(10,NML=INPUT_TEST_LIST)
    CLOSE(10)

    IF (ABS(PROTON_FRACTION-HALF)<1.D-7) THEN
      WRITE (*,*)
      WRITE (*,*) 'WARNING!! Code often fails for Yp = 1/2.'
      WRITE (*,*) 'Given your choices for proton_fraction '
      WRITE (*,*) 'code will try to calculate EOS within 10^-7 of 0.5.'
      WRITE (*,*) 'Example: Yp = 0.50001D0 or Yp = 0.49999D0.'
      CALL CONTINUE_EXCUTION()
    ENDIF

    IF (DENSITY < 1.D-14) THEN
      WRITE (*,*) 'Code is unreliable for densities n < 1.d-14 fm^-3'
      CALL CONTINUE_EXCUTION()
    ENDIF

    IF (DENSITY > 10.D0) THEN
      WRITE (*,*) 'Code is unreliable for densities n > 10.d0 fm^-3'
      CALL CONTINUE_EXCUTION()
    ENDIF

    IF (TEMPERATURE < 1.D-4) THEN
      WRITE (*,*) 'Code is unreliable for temperatures n < 1.d-4 MeV'
      CALL CONTINUE_EXCUTION()
    ENDIF

    IF (TEMPERATURE > 350.D0) THEN
      WRITE (*,*) 'Code is unreliable for temperatures T > 300.0 MeV'
      CALL CONTINUE_EXCUTION()
    ENDIF

    X1   = UNIFORM_SOL_GUESS
    X3   = NON_UNIFORM_SOL_GUESS

    CALL PRINT_TEST_INPUT

  END SUBROUTINE READ_TEST_INPUT


  SUBROUTINE READ_MAIN_INPUT

    USE Phase_Space_Input_Mod
    USE Global_Variables_Mod

    IMPLICIT NONE

    REAL(DP) :: n_old, t_old, Yp_chk, n_chk, t_chk, Yp
    INTEGER(I4B) :: i_Yp

    NAMELIST /SPACE_LIST/ Yp_min, Yp_max, steps_in_Yp, Yp_spacing, &
                Log10n_min, Log10n_max, steps_per_decade_in_n, log10n_spacing, &
                Log10T_min, Log10T_max, steps_per_decade_in_T, log10T_spacing

!   set all spacings to zero
    steps_in_Yp = 0
    Yp_spacing = ZERO
    steps_per_decade_in_n = ZERO
    log10n_spacing = ZERO
    steps_per_decade_in_T = ZERO
    log10T_spacing = ZERO

!   read input file
    OPEN(10,FILE=trim(adjustl(input_space_filename)), &
         FORM='FORMATTED',STATUS='OLD',ACTION='READ')
    READ(10,NML=SPACE_LIST)
    CLOSE(10)

!   check range of proton fraction Yp.
    CALL CHECK_RANGE_DP('Yp_min',Yp_min,1.D-3,1.D0)
    CALL CHECK_RANGE_DP('Yp_max',Yp_max,1.D-3,1.D0)
    CALL CHECK_ORDER(Yp_min,'Yp_min',Yp_max,'Yp_max')

    IF (Yp_spacing .NE. ZERO) THEN
      CALL CHECK_RANGE_DP('Yp_spacing',Yp_spacing,1.D-3,0.02D0)
      Yp_step = Yp_spacing
    ELSEIF (steps_in_Yp .NE. 0) THEN
      CALL CHECK_RANGE_INT('steps_in_Yp',steps_in_Yp,2,200)
      Yp_step = (Yp_max - Yp_min)/DBLE(steps_in_Yp-1)
    ELSE
       STOP "Need at least one of Yp_spacing or steps_in_Yp to be given!"
    ENDIF

    Yp_ini = 1
    Yp_chk = (Yp_max - Yp_min)/Yp_step
    Yp_fin = idnint(Yp_chk) + 1

!   check for Yp = 1/2.
    DO i_Yp = Yp_ini, Yp_fin
      ! set proton fraction
      Yp = Yp_min + dble(i_Yp-1)*Yp_step
      IF (ABS(Yp-HALF)<1.D-7) THEN
        WRITE (*,*)
        WRITE (*,*) 'WARNING!! Code often fails for Yp = 1/2.'
        WRITE (*,*) 'Given your choices for Yp_min, Yp_max, and Yp_spacing/steps_in_Yp '
        WRITE (*,*) 'code will try to calculate EOS within 10^-7 of 0.5.'
        WRITE (*,*) 'To avoid problems add/subtract ~ 10^-5 to Yp_min and Yp_max.'
        WRITE (*,"(A20,1F9.5)") 'Example: Yp_min = ', Yp_min - 1.D-5
        WRITE (*,"(A20,1F9.5)") '         Yp_max = ', Yp_max - 1.D-5
        CALL CONTINUE_EXCUTION()
      ENDIF
    ENDDO

!   check range of density n
    CALL CHECK_RANGE_DP('Log10n_min',Log10n_min,-14.D0,1.D0)
    CALL CHECK_RANGE_DP('Log10n_max',Log10n_max,-14.D0,1.D0)
    CALL CHECK_ORDER(Log10n_min,'Log10n_min',Log10n_max,'Log10n_max')

    IF (log10n_spacing .NE. ZERO) THEN
      CALL CHECK_RANGE_DP('log10n_spacing',log10n_spacing,1.D-3,1.D-1)
      steps_per_decade_in_n = 1.D0/log10n_spacing
    ELSEIF (steps_per_decade_in_n .NE. 0) THEN
      CALL CHECK_RANGE_DP('steps_per_decade_in_n',&
                           steps_per_decade_in_n,5.D0,100.D0)
    ELSE
      STOP &
         "Need at least one of log10n_spacing or steps_per_decade_in_n to be given!"
    ENDIF

    n_ini = 1
    n_chk = (Log10n_max-Log10n_min)*(steps_per_decade_in_n)
! check if n_chk is integer.
! only necessary to not mess-up ASCII output files
! modify Log10n_min slightly to adjust for that if necessary
    IF (abs(n_chk-dble(idnint(n_chk)))>1.d-8) THEN
      WRITE (*,*)
      WRITE (*,*) 'Steps in density computed to be', n_chk
      n_old = Log10n_min
      Log10n_min = Log10n_max - dble(int(n_chk)+1) / (steps_per_decade_in_n)
      WRITE (*,*) 'Adjusting Log10n_min so steps in density is ', int(n_chk)
      WRITE (*,*) 'Log10n_min changed from', n_old,' to ', Log10n_min
    ENDIF
    n_fin = idnint((Log10n_max-Log10n_min)*steps_per_decade_in_n) + 1

!   check range of temperature T
    CALL CHECK_RANGE_DP('Log10T_min',Log10T_min,-4.D0,2.5D0)
    CALL CHECK_RANGE_DP('Log10T_max',Log10T_max,-4.D0,2.5D0)
    CALL CHECK_ORDER(Log10T_min,'Log10T_min',Log10T_max,'Log10T_max')

    IF (log10n_spacing .NE. ZERO) THEN
      CALL CHECK_RANGE_DP('log10T_spacing',log10T_spacing,1.D-3,1.D-1)
      steps_per_decade_in_T = 1.D0/log10T_spacing
    ELSEIF (steps_per_decade_in_T .NE. 0) THEN
      CALL CHECK_RANGE_DP('steps_per_decade_in_T',&
                           steps_per_decade_in_T,5.D0,100.D0)
    ELSE
       STOP "Need at least one of 'log10T_spacing' or 'steps_per_decade_in_T' to be given!"
    ENDIF

    t_ini = 1
    t_chk = (Log10T_max-Log10T_min)*(steps_per_decade_in_T)
! check if t_chk is integer.
! only necessary to not mess-up ASCII output files
! modify Log10T_min slightly to adjust for that if necessary
    IF (abs(t_chk-dble(idnint(t_chk)))>1.d-8) THEN
      WRITE (*,*)
      WRITE (*,*) 'Steps in temperature computed to be', t_chk
      t_old = Log10T_min
      Log10T_min = Log10T_max - dble(int(t_chk)+1) / (steps_per_decade_in_T)
      WRITE (*,*) 'Adjusting Log10n_min so steps in temperature is ', int(t_chk)
      WRITE (*,*) 'Log10T_min changed from', t_old,' to ', Log10T_min
    ENDIF
    T_fin = idnint((Log10T_max-Log10T_min)*steps_per_decade_in_T) + 1

    CALL PRINT_MAIN_INPUT

  END SUBROUTINE READ_MAIN_INPUT

  SUBROUTINE PRINT_TEST_INPUT

    USE Test_Input_Mod

    WRITE (*,*)
    WRITE (*,*) 'INPUT FOR TEST:'
    WRITE (*,*)
    WRITE (*,"(A20,1ES18.10,A8)") 'TEMPERATURE      : ', TEMPERATURE, ' MeV   '
    WRITE (*,"(A20,1ES18.10,A8)") 'DENSITY          : ', DENSITY,     ' fm^-3 '
    WRITE (*,"(A20,1ES18.10,A8)") 'PROTON FRACTION  : ', PROTON_FRACTION
    WRITE (*,*)
    WRITE (*,"(A20,L1)"         ) 'USE INITIAL GUESS? ', GUESS
    IF (GUESS) THEN
      WRITE (*,"(A30,1ES18.10,A8)") 'UNIFORM SOLUTION GUESS    : ', &
                                     UNIFORM_SOL_GUESS
      WRITE (*,"(A30,3ES18.10,A8)") 'NON-UNIFORM SOLUTION GUESS: ', &
                                     NON_UNIFORM_SOL_GUESS
    ENDIF
    WRITE (*,*)

  END SUBROUTINE PRINT_TEST_INPUT


  SUBROUTINE PRINT_MAIN_INPUT

    USE Phase_Space_Input_Mod

    USE Physical_Constants_Mod, ONLY : TEN

    WRITE (*,*)
    WRITE (*,*) 'PHASE SPACE FOR EOS:'
    WRITE (*,*)
    WRITE (*,"(A20,1ES18.10,A4,1ES18.10,2A8,1ES18.10,A20)") &
      'TEMPERATURE        : ', TEN**Log10T_min, &
      ' to ',                  TEN**Log10T_max, ' MeV   ', &
      ' with ', steps_per_decade_in_T, ' steps per decade. '

    WRITE (*,"(A20,1ES18.10,A4,1ES18.10,2A8,1ES18.10,A20)") &
      'DENSITY            : ', TEN**Log10n_min, &
      ' to ',                  TEN**Log10n_max, ' fm^-3 ', &
      ' with ', steps_per_decade_in_n, ' steps per decade. '

    WRITE (*,"(A20,1ES18.10,A4,1ES18.10,2A8,1ES18.10,A20)") &
      'PROTON FRACTION    : ', Yp_min, &
      ' to ',                  Yp_max, '        ', &
      ' with ', Yp_step, ' steps.'
    WRITE (*,*)

    WRITE (*,*)

  END SUBROUTINE PRINT_MAIN_INPUT

  SUBROUTINE CONTINUE_EXCUTION()

    CHARACTER(LEN=1) :: continue_flag

    WRITE (*,*) 'CONTINUE ANYWAY? (Y/N)'
    READ  (*,*) continue_flag
    continue_flag = ADJUSTL(TRIM(continue_flag))
    IF     (TRIM(continue_flag) == 'Y' .OR. TRIM(continue_flag) == 'y') THEN
      WRITE (*,*) 'User chose to continue code execution despite warning.'
    ELSEIF (TRIM(continue_flag) == 'N' .OR. TRIM(continue_flag) == 'n') THEN
      STOP "Stopped code execution."
    ELSE
      WRITE (*,*) 'Invalid option ', continue_flag
      WRITE (*,*) 'Will stop code execution.'
      STOP "Stopped code execution."
    ENDIF

  END SUBROUTINE CONTINUE_EXCUTION

  SUBROUTINE CHECK_RANGE_DP(char_in,input_value,lower_limit,upper_limit)

    CHARACTER(LEN=*), INTENT(IN) :: char_in
    REAL(DP), INTENT(IN) :: input_value, lower_limit, upper_limit

    IF (input_value < lower_limit) THEN
      WRITE (*,*)
      WRITE (*,*) 'Input variable ',  char_in, ' set to ', input_value, &
                  ' is lower than reliable limit ', lower_limit
      CALL CONTINUE_EXCUTION()
    ENDIF

    IF (input_value > upper_limit) THEN
      WRITE (*,*)
      WRITE (*,*) 'Input variable ',  char_in, ' set to ', input_value, &
                  ' is larger than reliable limit ', upper_limit
      CALL CONTINUE_EXCUTION()
    ENDIF

    RETURN

  END SUBROUTINE CHECK_RANGE_DP

  SUBROUTINE CHECK_RANGE_INT(char_in,input_value,lower_limit,upper_limit)

    CHARACTER(LEN=*), INTENT(IN) :: char_in
    INTEGER(I4B), INTENT(IN) :: input_value, lower_limit, upper_limit

    IF (input_value < lower_limit) THEN
      WRITE (*,*)
      WRITE (*,*) 'Input variable ', char_in, ' set to ', input_value, &
                  ' is lower than reliable limit ', lower_limit
      CALL CONTINUE_EXCUTION()
    ENDIF

    IF (input_value > upper_limit) THEN
      WRITE (*,*)
      WRITE (*,*) 'Input variable ',  char_in, ' set to ', input_value, &
                  ' is larger than reliable limit ', upper_limit
      CALL CONTINUE_EXCUTION()
    ENDIF

    RETURN

  END SUBROUTINE CHECK_RANGE_INT

  SUBROUTINE CHECK_ORDER(lower_limit,char_lower,upper_limit,char_upper)

    CHARACTER(LEN=*), INTENT(IN) :: char_lower, char_upper
    REAL(DP), INTENT(IN) :: lower_limit, upper_limit

    IF (upper_limit <= lower_limit) THEN
      WRITE (*,*)
      WRITE (*,*) 'Input variable ', char_lower, &
                  'should not be larger than input variable limit ', char_upper
      STOP
    ENDIF

    RETURN

  END SUBROUTINE CHECK_ORDER

END MODULE Read_Input_Mod
