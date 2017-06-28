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
MODULE Kind_Types_Mod

  ! No implicit typing
  IMPLICIT NONE

  ! Explicit visibility declaration
  PRIVATE
  !PUBLIC :: Int_Byte, Int_Short, Int_Long
  !PUBLIC :: Real_Single, Real_Double
  PUBLIC :: SP, DP, I4B, LGCL

  ! Integer kinds
  INTEGER, PARAMETER :: Int_Byte  = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: Int_Short = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: Int_Long  = SELECTED_INT_KIND(9)

  ! Real kinds
  INTEGER, PARAMETER :: Real_Single = SELECTED_REAL_KIND(p=6)
  INTEGER, PARAMETER :: Real_Double = SELECTED_REAL_KIND(p=15)

  ! Logical kind
  INTEGER, PARAMETER :: LOGICAL = KIND(.true.)

  ! Generic kind used throughout
  INTEGER, PARAMETER :: SP   = Real_Single
  INTEGER, PARAMETER :: DP   = Real_Double
  INTEGER, PARAMETER :: I4B  = Int_Long
  INTEGER, PARAMETER :: LGCL = LOGICAL

!  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
!  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
!  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
!  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
!  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
!  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
!  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
!  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp

END MODULE Kind_Types_Mod
