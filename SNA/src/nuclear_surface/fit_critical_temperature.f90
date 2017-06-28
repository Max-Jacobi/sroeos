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
MODULE Fit_Critical_Temperature_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP
  USE Polynomial_Fit_Mod
  USE Surface_Properties_Mod
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO
  USE Critical_Point_Mod
  USE Make_Tables_Mod, ONLY : output_directory

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Fit_Critical_Temperature()

    IMPLICIT NONE

    INTEGER(I4B) :: i, j, imax

    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT   = "(1A10,1ES20.12,1A10)"
    CHARACTER(LEN=64) :: File_Name

    REAL(DP) :: prot_frac, temp, Tcrit
    REAL(DP), DIMENSION(2) :: x2
    REAL(DP), DIMENSION(100) :: neutron_excess_array
    REAL(DP), DIMENSION(100) :: critical_temperature_array
    REAL(DP), DIMENSION(100) :: critical_density_array
    LOGICAL(LGCL) :: critical_solution

    neutron_excess_array = zero
    critical_temperature_array = zero
    T_table = zero
    i = 0
    DO j = 100, 1, -1
      prot_frac = dble(j)/2.d2
!     For a given proton fraction find the maximum temperature and density
!     where two phases of different density and proton fractions can coexist
      CALL FIND_CRITICAL_POINT (x2, prot_frac, critical_solution)
      ! WRITE (*,*) prot_frac, x2, critical_solution
!     if a solution is found set values to array
      IF (critical_solution) THEN
        i = i + 1
        Temp = x2(2)
        T_table(j) = Temp
        IF (j==100) Tcrit = T_table(100)
        ! array of neutron excess in high density region
        neutron_excess_array(i) = (one-two*prot_frac)**two
        ! array of critical densities
        critical_density_array(i) = x2(1)
        ! array of (normalized) critical temperatures
        critical_temperature_array(i) = Temp/Tcrit
        ! use lin interpolation to determine whether T < 0 in next iteration
        ! if so, leave loop as no critical point exists for next iteration
        ! of the proton fraction
        IF (j<99) THEN
          IF (T_table(j+1)==zero .and. T_table(j+2)/=zero) &
            T_table(j+1) = (T_table(j+2) + T_table(j))/two
          IF (TWO*Temp-T_table(j+1)<ZERO) EXIT
        ENDIF
      ENDIF
    ENDDO
!    DO j = 100, 1, -1
!      IF (T_table(j) /= zero) WRITE (10,*) j, T_table(j)
!    ENDDO
    imax = i
    CLOSE(10)
! !   Check if all intermediate values for critical_temperature_array and
! !   critical_density_array were found
! !    If not interpolate from ones found
!     DO i = 2, imax
!       IF (critical_density_array(i)<tolerance .and. critical_density_array(i-1) > tolerance .and. critical_density_array(i+1) > tolerance ) then
!         xarray(i) = (xarray(i-1) + xarray(i+1)) / two
!         yarray(i) = (yarray(i-1) + yarray(i+1)) / two
!       endif
!     enddo

!   Fit critical tempearture with function
!   Tc(x) = Tc(0.5)*(a+b(1-2*x)^2+c(1-2*x)^4+d(1-2*x)^6)
!   where x is proton fraction
!   Use temporary arrays (x,T) -> ((1-2x)²,T) created above
    T_crit_coeffs(1:4) = polyfit(neutron_excess_array(1:imax), &
                                 critical_temperature_array(1:imax), 3)
    T_crit = Tcrit

!   output this in some decent format somewhere
    File_Name = ADJUSTL(TRIM(output_directory)) // &
                ADJUSTL(TRIM('/Critical_temperature_fit.dat'))
    OPEN(10,file=File_Name)

    WRITE (10,*)
    WRITE (10,*) 'Critical temperature fit. See Eq. 22 of SRO.'
    WRITE (10,*)
    WRITE (10,OUTPUT_FORMAT) 'T_c     =  ', T_crit, ' MeV '
    WRITE (10,OUTPUT_FORMAT) 'a       =  ', T_crit_coeffs(1)
    WRITE (10,OUTPUT_FORMAT) 'b       =  ', T_crit_coeffs(2)
    WRITE (10,OUTPUT_FORMAT) 'c       =  ', T_crit_coeffs(3)
    WRITE (10,OUTPUT_FORMAT) 'd       =  ', T_crit_coeffs(4)
    WRITE (10,*)
    WRITE (10,*)
    CLOSE(10)

  END SUBROUTINE Fit_Critical_Temperature

END MODULE Fit_Critical_Temperature_Mod
