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
MODULE Determine_Surface_Properties_Mod

  USE Global_Variables_Mod, ONLY: input_skyrme_filename
  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL
  USE Surface_Properties_Mod
  USE Fit_Critical_Temperature_Mod
  USE Nuclear_Surface_Tension_Mod
  USE Nuclear_Surface_Fit_Mod
  USE Surface_Observable_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE DETERMINE_SURFACE_PROPERTIES

    IMPLICIT NONE

    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT   = "(1A10,1ES20.12,1A10)"
    CHARACTER(LEN=256) :: File_Name

    INTEGER(I4B) :: data_points, io_stat
    REAL(DP), DIMENSION(10000,1:3) :: TAB_SURF_TENS_VAL
    REAL(DP) :: S_S, A_S

    NAMELIST /SURFACE_PROPERTIES/ fit_surface, surface_sigma, &
                                  surface_p, surface_q, surface_lambda, &
                                  fit_temperature, T_crit, T_crit_coeffs

    NAMELIST /HEAVY_NUCLEI_SIZE/  fix_heavy_nuclei_size, A0, A00

!   set default for procedures below
    fit_temperature = .true.
    fit_surface = .true.
    fix_heavy_nuclei_size = .false.

    OPEN(10,FILE=input_skyrme_filename,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
    READ(10,NML=SURFACE_PROPERTIES,IOSTAT=io_stat)
    READ(10,NML=HEAVY_NUCLEI_SIZE,IOSTAT=io_stat)
    CLOSE(10)

    IF (fit_temperature) THEN
      CALL FIT_CRITICAL_TEMPERATURE
      WRITE (*,*) '        Fitted Critical temperature based on Skyrme parametrization.'
    ELSE
      WRITE (*,*) '        Used user input for critical temperature fitting.'
    ENDIF

    IF (fit_surface) THEN
      WRITE (*,*) '        Fitted Surface tension based on Skyrme parametrization.'
      Surface_sigma  = zero
      Surface_q      = zero
      surface_p      = zero
      surface_lambda = zero
      CALL NUCLEAR_SURFACE_TENSION(data_points,TAB_SURF_TENS_VAL)
      WRITE (*,*) '        Fitted Surface tension based on Skyrme parametrization.'
      CALL NUCLEAR_SURFACE_FIT(data_points,TAB_SURF_TENS_VAL)
    ELSE
      WRITE (*,*) '        Used user input for Surface fitting.'
    ENDIF

    CALL SURFACE_OBSERVABLES

  END SUBROUTINE DETERMINE_SURFACE_PROPERTIES

END MODULE Determine_Surface_Properties_Mod
