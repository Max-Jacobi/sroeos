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

!
!  This file is originally part of the NLEQSLV of Berend Hasselman available at
!  https://cran.r-project.org/web/packages/nleqslv/index.html
!
!  It has been slightly modified by Andre da Silva Schneider
!   to fit SRO_EOS source code
!
MODULE nwprnt_mod

  USE Kind_Types_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE nwckot(i,j,aij,wi)

    IMPLICIT NONE
!
!   output for a single incorrect jacobian entry
!
    INTEGER(I4B) :: i, j, k
    REAL(DP) ::  aij, wi

    k = 0
    ! removed print statements
!    write (17,801) i, j, aij
!    write (17,802) i, j, wi

801   format ('Chkjac possible error in jacobian[',I2,',',I2,']=', 1ES20.13)
802   format ('                        Estimated[',I2,',',I2,']=', 1ES20.13)

  END SUBROUTINE nwckot

END MODULE nwprnt_mod
