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

MODULE Global_Variables_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL

  IMPLICIT NONE

  SAVE

  LOGICAL(LGCL) :: IS_TEST
  character(len=256) :: input_space_filename
  character(len=256) :: input_skyrme_filename

END MODULE Global_Variables_Mod
