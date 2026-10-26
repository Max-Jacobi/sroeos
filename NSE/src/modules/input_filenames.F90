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
!    Input file names. The stock files under input/ are the defaults;
!    READ_COMMANDLINE overrides all five when exactly five arguments are
!    given (same convention as SNA's sro_eos and MERGE's merge_eos).
!
MODULE Input_Files_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL

  IMPLICIT NONE

  SAVE

  character(len=256) :: input_space = "input/space.in"
  character(len=256) :: input_partition = "input/partition.in"
  character(len=256) :: input_iso_properties = "input/isotopes.in"
  character(len=256) :: input_iso_list = "input/list.in"
  character(len=256) :: input_out_list = "input/output.in"

CONTAINS

  SUBROUTINE READ_COMMANDLINE
    integer :: num_args

    num_args = command_argument_count()

    if (num_args == 0) return

    if (num_args /= 5) then
      write(*,*) 'Usage: main [<space.in> <partition.in> <isotopes.in>', &
                 ' <list.in> <output.in>]'
      write(*,*) 'With no arguments the files are read from input/.'
      stop
    end if

    call get_command_argument(1, input_space)
    call get_command_argument(2, input_partition)
    call get_command_argument(3, input_iso_properties)
    call get_command_argument(4, input_iso_list)
    call get_command_argument(5, input_out_list)
  END SUBROUTINE READ_COMMANDLINE

END MODULE Input_Files_Mod
