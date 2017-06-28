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

MODULE Input_Files_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL

  IMPLICIT NONE

  SAVE

  character(len=256), parameter :: input_space = "input/space.in"
  character(len=256), parameter :: input_partition = "input/partition.in"
  character(len=256), parameter :: input_iso_properties = "input/isotopes.in"
  character(len=256), parameter :: input_iso_list = "input/list.in"
  character(len=256), parameter :: input_out_list = "input/output.in"

END MODULE Input_Files_Mod
