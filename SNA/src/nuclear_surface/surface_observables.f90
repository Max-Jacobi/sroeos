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
MODULE Surface_Observable_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, THREE, FOUR, PI, R_2_3
  USE H_Mod
  USE Sigma_Mod
  USE Nuclear_Matter_Properties_Mod, ONLY : Nuc_Sat_Dens
  USE Nuclear_Surface_Fit_Mod
  USE Make_Tables_Mod

CONTAINS

  SUBROUTINE SURFACE_OBSERVABLES()

    IMPLICIT NONE

    CHARACTER(LEN=256) :: File_Name
    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT   = "(1A20,1ES23.15,1A20)"

    INTEGER(I4B) :: derivatives
    REAL(DP) :: y,T,H,DH_DT,DH_DY,D2H_DT2,D2H_DTDY,D2H_DY2
    REAL(DP) :: SIGMA,DSIGMA_DT,DSIGMA_DY,D2SIGMA_DT2,D2SIGMA_DTDY,D2SIGMA_DY2
    REAL(DP) :: A_S, S_S

!
!   output this in some decent format somewhere
    File_Name = ADJUSTL(TRIM(output_directory)) // &
                ADJUSTL(TRIM('/Surface_fit.dat'))
    OPEN(10,file=File_Name)

    WRITE (10,*)
    WRITE (10,*) 'Surface tenstion fit. See Eqs. 19 and 20 of SRO.'
    WRITE (10,*)

    WRITE (10,OUTPUT_FORMAT) 'surface_sigma   = ', surface_sigma, '! in  MeV fm^-2'
    WRITE (10,OUTPUT_FORMAT) 'surface_q       = ', surface_q
    WRITE (10,OUTPUT_FORMAT) 'surface_lambda  = ', surface_lambda
    WRITE (10,OUTPUT_FORMAT) 'surface_p       = ', surface_p
    WRITE (10,*)

    derivatives = 2

    y = 0.5D0
    T = 0.001D0

    CALL H_SUB (derivatives,y,T,H,DH_DT,DH_DY,D2H_DT2,D2H_DTDY,D2H_DY2)
    CALL SIGMA_SUB (derivatives,y,T,H,DH_DT,DH_DY,D2H_DT2,D2H_DY2,D2H_DTDY,&
              SIGMA,DSIGMA_DT,DSIGMA_DY,D2SIGMA_DT2,D2SIGMA_DTDY,D2SIGMA_DY2)

    S_S = - PI/TWO*(THREE/FOUR/PI/Nuc_Sat_Dens)**(R_2_3)*D2SIGMA_DY2
    A_S = - PI*TWO*(THREE/FOUR/PI/Nuc_Sat_Dens)**(R_2_3)*D2SIGMA_DT2

    WRITE (10,*)
    WRITE (10,*) 'Surface tension properties. See Eqs. 53 of SRO.'
    WRITE (10,*)

    WRITE (10,OUTPUT_FORMAT) 'A_S     =  ', A_S
    WRITE (10,OUTPUT_FORMAT) 'S_S     =  ', S_S
    WRITE (10,*)

    CLOSE(10)

  END SUBROUTINE SURFACE_OBSERVABLES


END MODULE Surface_Observable_Mod
