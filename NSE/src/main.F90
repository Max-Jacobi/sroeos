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

PROGRAM Main
  ! Use modules
  USE Kind_Types_Mod, ONLY : I4B, DP, LGCL
  USE Read_Input_Mod, ONLY : READ_OUTPUT_INPUT, READ_SPACE_INPUT
  USE READ_NUCLEAR_DATA_TABLE_MOD
  USE Make_Tables_Mod
  USE Allocate_Mod
  USE Solve_Mod
  USE Write_to_table_Mod
#ifdef _OPENMP
  USE omp_lib
#endif
  
  ! No implicit typing
  IMPLICIT NONE

  integer :: nthreads

#ifdef _OPENMP
  nthreads = omp_get_max_threads()
#else
  nthreads = 1
#endif
  
  ! output banner
  write(6,*) "******************************************************************************"
  write(6,*) "*  Schneider-Roberts-Ott EOS Code: NSE MAIN                                  *"
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

  WRITE (*,*)
  WRITE (*,*) 'Step 1: Reading input.'
  CALL READ_OUTPUT_INPUT

  WRITE (*,*)
  WRITE (*,*) 'Step 2: Reading table range and resolution.'
  CALL READ_SPACE_INPUT

  WRITE (*,*)
  WRITE (*,*) 'Step 3: Reading list of isotopes.'
  CALL READ_NUCLEAR_DATA

  WRITE (*,*)
  WRITE (*,*) 'Step 4: Obtaining Equation of State.'

! allocate array sizes for HDF5 table
  WRITE (*,*)
  WRITE (*,*) ' Step 4a:  Allocate arrays for HDF5 output.'
  IF (make_hdf5_table)  CALL ALLOCATE_OUTPUT

! get NSE EOS solution
  WRITE (*,*)
  WRITE (*,*) ' Step 4b:  Solving equation of state.'
  CALL SOLVE
  !
  CALL WRITE_TO_TABLE

END PROGRAM Main
