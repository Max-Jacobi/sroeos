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

  include "../version.inc"

CONTAINS

  SUBROUTINE READ_OUTPUT_INPUT

    USE Tables_Input_Mod
    USE Phase_Space_Input_Mod, ONLY : Yp_fin, n_fin, T_fin
    USE Input_Files_Mod, ONLY : input_tables

    CHARACTER(LEN=20) :: string_rho, string_temp, string_yp
    CHARACTER(LEN=100) :: gitinfo

    INTEGER(I4B) :: values(8)
    CHARACTER(LEN=8) :: date

    NAMELIST /TABLES/   sna_eos, nse_eos, &
                        make_hdf5_table, merge_base, &
                        only_SNA, only_NSE

    make_hdf5_table  = .TRUE.
    merge_base = ' '
    sna_eos = ' '
    nse_eos = ' '
    only_SNA = .FALSE.
    only_NSE = .FALSE.

    ! nchanged is in ../version.inc, included in
    ! Make_Tables_Mod
    if(nchanged > 0) then
       gitinfo = "M"//trim(adjustl(git_version))
    else
       gitinfo = git_version
    endif

!   read input file
    OPEN(10,FILE=trim(adjustl(input_tables)), &
         FORM='FORMATTED',STATUS='OLD',ACTION='READ')
    READ(10,NML=TABLES)
    CLOSE(10)

    IF (only_SNA .AND. only_NSE) THEN
      STOP 'Parameters only_SNA and only_NSE cannot be both true. '
    ENDIF

!   check date and time and write timestamp and filename
    CALL date_and_time(DATE=date,VALUES=values)
    WRITE (string_rho ,*) n_fin
    WRITE (string_temp,*) t_fin
    WRITE (string_yp,  *) Yp_fin

    h5_filename = &
        trim(adjustl(merge_base))//"_rho"//trim(adjustl(string_rho))// &
        "_temp"//trim(adjustl(string_temp))//"_ye"//trim(adjustl(string_yp))// &
        "_git"//trim(adjustl(gitinfo))//"_"//trim(adjustl(date))//".h5"

!   TODO: CHECK IF SNA AND NSE TABLES EXIST

  END SUBROUTINE READ_OUTPUT_INPUT

  SUBROUTINE READ_TRANSITION_INPUT

    USE Transition_Input_Mod
    USE Tables_Input_Mod, ONLY : only_SNA, only_NSE
    USE Input_Files_Mod, ONLY : input_transition
    USE Physical_Constants_Mod, ONLY : ONE, TWO, TEN

    NAMELIST /TRANSITION/   n_transition, n_delta, n_tolerance

    n_transition = zero
    n_delta      = zero
    n_tolerance  = zero

!   read input file
    OPEN(10,FILE=trim(adjustl(input_transition)), &
         FORM='FORMATTED',STATUS='OLD',ACTION='READ')
    READ(10,NML=TRANSITION)
    CLOSE(10)

    IF (.not.only_sna .AND. .not. only_nse) THEN
      IF (n_transition==zero) THEN
        WRITE (*,*) 'n_transition in transition.in file not set.'
        WRITE (*,*) 'Using default value n_transition = -4.D0'
        n_transition = -4.00d0
      ENDIF
      IF (n_delta==zero) THEN
        WRITE (*,*) 'n_delta in transition.in file not set.'
        WRITE (*,*) 'Using default value n_delta = 0.33'
        n_delta      =  0.33d0
      ENDIF
      IF (n_tolerance==zero) THEN
        WRITE (*,*) 'n_tolerance in transition.in file not set.'
        WRITE (*,*) 'Using default value n_delta = 1.D-4'
      ENDIF
    ENDIF

    IF (only_SNA) THEN
      WRITE (*,*)
      WRITE (*,*) 'Paramater only_SNA was set to TRUE.'
      WRITE (*,*) ' Setting value of n_transition to -100 to avoid issues.'
      WRITE (*,*) ' Setting value of n_delta      to  0.5 to avoid issues.'
      WRITE (*,*)
      n_transition = -100.0d0
      n_delta      =  0.5d0
      n_tolerance  =  1.0d-6
    ENDIF

    IF (only_NSE) THEN
      WRITE (*,*)
      WRITE (*,*) 'Paramater only_NSE was set to TRUE.'
      WRITE (*,*) ' Setting value of n_transition to +100 to avoid issues.'
      WRITE (*,*) ' Setting value of n_delta      to  0.5 to avoid issues.'
      WRITE (*,*)
      n_transition = +100.0d0
      n_delta      =  0.5d0
      n_tolerance  =  1.0d-6
    ENDIF

    Log10nt_min = n_transition + n_delta*atanh(two*n_tolerance-one)
    Log10nt_max = n_transition - n_delta*atanh(two*n_tolerance-one)

    write (*,"(A22,1ES12.4)") 'Transition density:   ', n_transition
    write (*,"(A22,1ES12.4)") 'Transition thickness: ', n_delta
    write (*,"(A22,1ES12.4)") 'Transition tolerance: ', n_tolerance
    write (*,"(A30,2ES12.4)") 'Log10  merge region (fm^-3): ', &
                              Log10nt_min, Log10nt_max
    write (*,"(A30,2ES12.4)") 'Linear merge region (fm^-3): ', &
                              ten**Log10nt_min, ten**Log10nt_max

    ! TODO: CHECK IF N_TRANSITION, N_DELTA, AND N_DELTA
    !       ARE WITHIN ACCEPTABLE RANGES.

  END SUBROUTINE READ_TRANSITION_INPUT


  SUBROUTINE READ_SPACE_INPUT

    USE Phase_Space_Input_Mod
    USE Input_Files_Mod, ONLY : input_space
    USE Table_Sizes_Mod, ONLY : final_tab, nvars

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
    OPEN(10,FILE=trim(adjustl(input_space)), &
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
      STOP "Need at least one of 'Yp_spacing' or 'steps_in_Yp' to be given!"
    ENDIF

    Yp_ini = 1
    Yp_chk = (Yp_max - Yp_min)/Yp_step
    Yp_fin = INT(Yp_chk) + 1

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
      STOP "Need at least one of 'log10n_spacing' or 'steps_per_decade_in_n' &
              to be given!"
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
      STOP "Need at least one of 'log10T_spacing' or 'steps_per_decade_in_T' &
              to be given!"
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

    ALLOCATE(final_tab(n_fin,t_fin,yp_fin,nvars))

    CALL PRINT_MAIN_INPUT

  END SUBROUTINE READ_SPACE_INPUT

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
