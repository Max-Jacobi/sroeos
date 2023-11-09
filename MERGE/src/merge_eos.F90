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
PROGRAM MERGE_EOS

    USE Read_Commandline_Mod
    USE Read_Input_Mod
    use Read_Table_Mod
    use MERGE_TABLES_MOD
    use write_table
#ifdef _OPENMP
    use omp_lib
#endif
!*****************************************************************************80
!
!  Ths program merges three EoS tables
!    - SNA EOS table
!    - NSE EOS
!    - lepton/photon EOS


  IMPLICIT NONE

  CHARACTER(LEN=4) :: flag
  INTEGER(I4B) :: i
  integer :: nthreads
#ifdef _OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif

  ! output banner
  write(6,*) "******************************************************************************"
  write(6,*) "*  Schneider-Roberts-Ott EOS Code: MERGE                                     *"
  write(6,*) "******************************************************************************"
  write(6,*)
  write(6,*) "git revision: ",adjustl(trim(git_version))
  write(6,*) "Version date: ",adjustl(trim(version_date))
  if (nchanged > 0) then
     write(6,*) "Warning: uncommitted changes."
  else
     write(6,*) "No uncommitted changes."
  endif
  write(6,*) "Source code and inputs will be included in the HDF5 table."
  write(6,*)
  write(6,"(A19,I3,A10)") " Using OpenMP with ",nthreads," thread(s)"
  write(6,*)
  write(6,*) "******************************************************************************"
  call sleep(5)

!  call subroutine that reads input parameters
!   that define the range of validity of each table
!   and the names of the tables to be merged

  CALL READ_COMMANDLINE

  WRITE (*,*)
  WRITE (*,*)  'Read file space input  file'
  CALL READ_SPACE_INPUT


  WRITE (*,*)
  WRITE (*,*)  'Read file tables input file'
  CALL READ_OUTPUT_INPUT

  WRITE (*,*)
  WRITE (*,*)  'Read file transition input file'
  CALL READ_TRANSITION_INPUT

!  call subroutine that reads NSE EOS
  IF (only_NSE .OR. .not.only_SNA) THEN
    WRITE (*,*)
    WRITE (*,*)  'Reading NSE table ', nse_eos
    flag = "NSE"
    CALL READTABLE(nse_eos,flag)
  ENDIF

  IF (only_SNA .OR. .not.only_NSE) THEN
    WRITE (*,*)
    WRITE (*,*)  'Reading SNA table ', sna_eos
    flag = "SNA"
    CALL READTABLE(sna_eos,flag)
  ENDIF

! this is for testing
!  call copy_src_input(trim(adjustl(nse_eos))//C_NULL_CHAR,&
!       trim(adjustl(sna_eos))//C_NULL_CHAR, &
!       trim(adjustl("bogus"))//C_NULL_CHAR,i)

!  call subroutine that merges the tables
  WRITE (*,*)
  WRITE (*,*) 'Merging SNA and NSE tables.'
  CALL MERGE

!  call subroutine to output table
  WRITE (*,*)
  WRITE (*,*) 'Write output files'
  CALL OUTPUT

END PROGRAM MERGE_EOS
