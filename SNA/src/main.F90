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
  USE Read_Input_Mod, only : READ_MAIN_INPUT
  USE Read_Skyrme_Coefficients_Mod, ONLY : READ_SKYRME_COEFFICIENTS
  USE Determine_Nuclear_Properties_Mod
  USE Determine_Surface_Properties_Mod
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
  
  IS_TEST = .FALSE.

  ! output banner
  write(6,*) "******************************************************************************"
  write(6,*) "*  Schneider-Roberts-Ott EOS Code: SNA MAIN                                  *"
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
  WRITE (*,*) 'Step 1: Read in Skyrme Coefficients.'
  CALL READ_SKYRME_COEFFICIENTS

  WRITE (*,*)
  WRITE (*,*) 'Step 2: Determine nuclear matter properties.'
  CALL DETERMINE_NUCLEAR_PROPERTIES

  WRITE (*,*)
  WRITE (*,*) 'Step 3: Determine Surface properties.'
  WRITE (*,*) '        Could take ~ 5 minutes.'
  CALL DETERMINE_SURFACE_PROPERTIES

  WRITE (*,*)
  WRITE (*,*) 'Step 4: Read table range and resolution.'
  CALL READ_MAIN_INPUT

! allocate array sizes for HDF5 table
  WRITE (*,*)
  WRITE (*,*) 'Step 5: Allocate arrays for output.'
  IF (make_hdf5_table) CALL ALLOCATE_OUTPUT

  WRITE (*,*)
  WRITE (*,*) 'Step 6: Obtain EoS.'
  CALL SOLVE

  WRITE (*,*)
  WRITE (*,*) 'Step 7: Write output table(s).'
  CALL WRITE_TO_TABLE

END PROGRAM Main
