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

MODULE Read_Commandline_Mod

  USE GLOBAL_VARIABLES_MOD

  IMPLICIT NONE

  contains

  SUBROUTINE READ_COMMANDLINE
    integer :: num_args

    num_args = command_argument_count()

    if (num_args .ne. 2) then
      write(*,*) 'Usage: sro_eos input_space_filename input_skyrme_filename'
      stop
    end if

    call get_command_argument(1, input_space_filename)
    call get_command_argument(2, input_skyrme_filename)
  END SUBROUTINE READ_COMMANDLINE

END MODULE Read_Commandline_Mod
