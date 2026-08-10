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

MODULE Make_Tables_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL

  IMPLICIT NONE

  SAVE

  LOGICAL(LGCL) ::  make_hdf5_table
  LOGICAL(LGCL) ::  write_solutions_to_file
  LOGICAL(LGCL) ::  redefine_print_parameters

  CHARACTER(LEN=256)  :: filename_base, output_directory

! optional per-cell seeds for the non-uniform solve
! (log10 n_no, log10 n_po, log10 u), read from seed_filename
  CHARACTER(LEN=256)  :: seed_filename
  LOGICAL(LGCL) :: use_seeds
  REAL(DP),      DIMENSION(:,:,:,:), ALLOCATABLE :: cell_seed
  LOGICAL(LGCL), DIMENSION(:,:,:),   ALLOCATABLE :: cell_seed_mask

  include "../version.inc"

END MODULE Make_Tables_Mod
