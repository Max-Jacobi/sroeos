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
MODULE HDF5_Writer_Mod

  USE HDF5

  IMPLICIT NONE

CONTAINS

  FUNCTION h5_output_double(lfid,var_name,var_array,ldims,lrank) result(lierr)

    IMPLICIT NONE

    integer :: lierr
    integer(HID_T), intent(in) :: lfid
    character(LEN=*), intent(in)          :: var_name
    real(8), dimension(*), intent(in) :: var_array
    integer(HSIZE_T), dimension(:),     intent(in) :: ldims
    integer, intent(in) :: lrank

    integer(HID_T) :: ldspace_id,ldset_id,mem_type
    integer :: lerror
    mem_type = H5T_NATIVE_DOUBLE
    call h5screate_simple_f(lrank, ldims, ldspace_id, lerror)
    lierr = lerror
    call h5dcreate_f(lfid, TRIM(var_name), mem_type, &
         ldspace_id, ldset_id, lerror)
    lierr = lerror + lierr
    call h5dwrite_f(ldset_id, mem_type, var_array, &
         ldims, lerror)
    lierr = lerror + lierr
    call h5dclose_f(ldset_id, lerror)
    lierr = lerror + lierr
    call h5sclose_f(ldspace_id, lerror)
    lierr = lerror + lierr
    if (lierr .ne. 0) then
      write(6,'(a,a,a)') 'Error writing variable ',TRIM(var_name),' to file.'
    else
      write(6,'(a,a,a)') TRIM(var_name),' written to file.'
    endif

    RETURN

  END FUNCTION h5_output_double

  FUNCTION h5_output_int(lfid,var_name,var_array,ldims,lrank) result(lierr)

    IMPLICIT NONE

    integer :: lierr
    integer(HID_T), intent(in) :: lfid
    character(LEN=*), intent(in)          :: var_name
    integer, dimension(*), intent(in) :: var_array
    integer(HSIZE_T), dimension(:),     intent(in) :: ldims
    integer, intent(in) :: lrank

    integer(HID_T) :: ldspace_id,ldset_id,mem_type
    integer :: lerror

    mem_type = H5T_NATIVE_INTEGER
    call h5screate_simple_f(lrank, ldims, ldspace_id, lerror)
    lierr = lerror
    call h5dcreate_f(lfid, TRIM(var_name), mem_type, &
         ldspace_id, ldset_id, lerror)
    lierr = lerror + lierr
    call h5dwrite_f(ldset_id, mem_type, var_array, &
         ldims, lerror)
    lierr = lerror + lierr
    call h5dclose_f(ldset_id, lerror)
    lierr = lerror + lierr
    call h5sclose_f(ldspace_id, lerror)
    lierr = lerror + lierr
    if (lierr .ne. 0) then
      write(6,'(a,a,a)') 'Error writing variable ',TRIM(var_name),' to file.'
    else
      write(6,'(a,a,a)') TRIM(var_name),' written to file.'
    endif

    RETURN

  END FUNCTION h5_output_int

END MODULE HDF5_Writer_Mod
