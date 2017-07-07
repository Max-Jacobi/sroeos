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
MODULE Write_to_table_Mod

  USE Kind_Types_Mod,        ONLY : I4B, DP
  USE Physical_Constants_Mod
  USE HDF5
  USE HDF5_Writer_Mod
  USE Make_Tables_Mod
  USE Allocatable_Tables_Mod
  USE Input_Files_Mod
  USE iso_c_binding

  IMPLICIT NONE

CONTAINS

  SUBROUTINE write_to_table

!   Subroutine that writes results of EoS to an
!    HDF5 table to be read by GR1D code
!   Almost all of this subroutine was copied from
!    EOSmaker v 1.0 by Evan O'Connor and Christian Ott

    IMPLICIT NONE

    INTEGER(I4B) :: values(8)
    CHARACTER(LEN=8) :: date
    REAL(DP) :: timestamp

    CHARACTER(LEN=128) h5_filename
    CHARACTER(LEN=100) :: string_n, string_T, string_Ye
    CHARACTER(LEN=100) :: gitinfo

  ! Variables needed for output
    INTEGER(I4B) :: error, rank, ierr
    INTEGER(HID_T) :: file_id
    INTEGER(HSIZE_T) dims1(1), dims2(2), dims3(3)

!   check date and time and write timestamp and filename
    call date_and_time(DATE=date,VALUES=values)

!   number of points n logrho, logT and Ye.
    write(string_n, *) pointsrho
    write(string_T, *) pointsTemp
    write(string_Ye,*) pointsYe

    if(nchanged > 0) then
       gitinfo = "M"//trim(adjustl(git_version))
    else
       gitinfo = git_version
    endif

    timestamp =   dble(values(1))*10000.0d0 &
                + dble(values(2))*100.0&
                + dble(values(3)) &
                + (dble(values(5)) + dble(values(6))/60.0d0 &
                + dble(values(7))/3600.0d0)/24.0

    h5_filename = trim(adjustl(output_directory))//"/"// & 
                  trim(adjustl(filename_base)) &
                  //"_rho"//trim(adjustl(string_n))  &
                  //"_temp"//trim(adjustl(string_T)) &
                  //"_ye"//trim(adjustl(string_Ye))  &
                  //"_git"//trim(adjustl(gitinfo)) &
                  //"_"//trim(adjustl(date))//".h5"


    write( *,"(A)")  "Writing HDF5 EOS Table: "
    write( *,"(A)")  trim(adjustl(h5_filename))
    write( *,"(A9,F16.7)")  "at time: ", timestamp

! Write to hdf5 file
    call h5open_f(error)
    call h5fcreate_f(h5_filename,H5F_ACC_TRUNC_F,file_id,error)

    rank = 1

    dims1(1) = 1

    ierr = ierr + h5_output_int(file_id,"pointsrho",(/pointsrho/),dims1,rank)
    ierr = ierr + h5_output_int(file_id,"pointstemp",(/pointstemp/),dims1,rank)
    ierr = ierr + h5_output_int(file_id,"pointsye",(/pointsYe/),dims1,rank)

    rank = 1

    dims1(1) = pointsrho
    ierr = ierr + h5_output_double(file_id,"logn",log10n_tab,dims1,rank)

    dims1(1) = pointstemp
    ierr = ierr + h5_output_double(file_id,"logtemp",log10t_tab,dims1,rank)

    dims1(1) = pointsYe
    ierr = ierr + h5_output_double(file_id,"ye",yp_tab,dims1,rank)

    rank=3

    dims3(1) = pointsrho
    dims3(2) = pointstemp
    dims3(3) = pointsYe

    ierr = ierr + h5_output_double(file_id,"p",p_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"s",s_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"e",e_tab,dims3,rank)

    mun_tab = mun_tab - neutron_mass_MeV
    mup_tab = mup_tab - neutron_mass_MeV

    ierr = ierr + h5_output_double(file_id,"muh",muh_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"mun",mun_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"mup",mup_tab,dims3,rank)

    ierr = ierr + h5_output_double(file_id,"dpdn",dpdn_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"dsdn",dsdn_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"dmudn",dmudn_tab,dims3,rank)

    ierr = ierr + h5_output_double(file_id,"dpdt",dpdt_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"dsdt",dsdt_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"dmudt",dmudt_tab,dims3,rank)

    ierr = ierr + h5_output_double(file_id,"dpdy",dpdy_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"dsdy",dsdy_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"dmudy",dmudy_tab,dims3,rank)

    ierr = ierr + h5_output_double(file_id,"zh",zh_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"ah",ah_tab,dims3,rank)

    ierr = ierr + h5_output_double(file_id,"zl",zl_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"al",al_tab,dims3,rank)

    ierr = ierr + h5_output_double(file_id,"xn",xn_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"xp",xp_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"xa",xa_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"xh",xh_tab,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"xl",xl_tab,dims3,rank)

  !  ierr = ierr + h5_output_int(file_id,"error",ierror,dims3,rank)

    call h5fclose_f(file_id,error)
    call h5close_f(error)

    ! now add source code and input file to table:
    call include_src_input(trim(adjustl(h5_filename))//C_NULL_CHAR, &
         trim(adjustl(input_space))//C_NULL_CHAR,      &
         trim(adjustl(input_partition))//C_NULL_CHAR,  &
         trim(adjustl(input_iso_properties))//C_NULL_CHAR,  &
         trim(adjustl(input_iso_list))//C_NULL_CHAR,  &
         trim(adjustl(input_out_list))//C_NULL_CHAR)

    write(*,*) "Done! :-)"

    RETURN

  END SUBROUTINE write_to_table

END MODULE Write_to_table_Mod
