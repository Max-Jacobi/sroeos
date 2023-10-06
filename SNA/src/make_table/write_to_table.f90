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
  USE Global_Variables_Mod
  USE Print_Parameters_Mod
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

    CHARACTER(LEN=256) h5_filename
    CHARACTER(LEN=100) :: gitinfo, string_n, string_T, string_Ye

  ! Variables needed for output
    INTEGER(I4B) :: error, rank, ierr
    INTEGER(HID_T) :: file_id
    INTEGER(HSIZE_T) dims1(1), dims2(2), dims3(3)

    ! nchanged is in ../version.inc, included in
    ! Make_Tables_Mod
    if(nchanged > 0) then
       gitinfo = "M"//trim(adjustl(git_version))
    else
       gitinfo = git_version
    endif

!   check date and time and write timestamp and filename
    call date_and_time(DATE=date,VALUES=values)

!   number of points n logrho, logT and Ye.
    write(string_n, *) pointsrho
    write(string_T, *) pointsTemp
    write(string_Ye,*) pointsYe

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

    if (print_p) ierr = ierr + h5_output_double(file_id,"p",p_tab,dims3,rank)
    if (print_s) ierr = ierr + h5_output_double(file_id,"s",s_tab,dims3,rank)
    if (print_e) ierr = ierr + h5_output_double(file_id,"e",e_tab,dims3,rank)

    if (print_muh) ierr = ierr + h5_output_double(file_id,"muh",muh_tab,dims3,rank)
    if (print_mun) ierr = ierr + h5_output_double(file_id,"mun",mun_tab,dims3,rank)
    if (print_mup) ierr = ierr + h5_output_double(file_id,"mup",mup_tab,dims3,rank)

    if (print_dpdn)  ierr = ierr + h5_output_double(file_id,"dpdn",dpdn_tab,dims3,rank)
    if (print_dsdn)  ierr = ierr + h5_output_double(file_id,"dsdn",dsdn_tab,dims3,rank)
    if (print_dmuhdn) ierr = ierr + h5_output_double(file_id,"dmudn",dmudn_tab,dims3,rank)

    if (print_dpdt)  ierr = ierr + h5_output_double(file_id,"dpdt",dpdt_tab,dims3,rank)
    if (print_dsdt)  ierr = ierr + h5_output_double(file_id,"dsdt",dsdt_tab,dims3,rank)
    if (print_dmuhdt) ierr = ierr + h5_output_double(file_id,"dmudt",dmudt_tab,dims3,rank)

    if (print_dpdy)  ierr = ierr + h5_output_double(file_id,"dpdy",dpdy_tab,dims3,rank)
    if (print_dsdy)  ierr = ierr + h5_output_double(file_id,"dsdy",dsdy_tab,dims3,rank)
    if (print_dmuhdy) ierr = ierr + h5_output_double(file_id,"dmudy",dmudy_tab,dims3,rank)

    if (print_zbar) ierr = ierr + h5_output_double(file_id,"zbar",zbar_tab,dims3,rank)
    if (print_abar) ierr = ierr + h5_output_double(file_id,"abar",abar_tab,dims3,rank)

    if (print_u) ierr = ierr + h5_output_double(file_id,"u",u_tab,dims3,rank)
    if (print_r) ierr = ierr + h5_output_double(file_id,"r",r_tab,dims3,rank)

    if (print_meff_n) ierr = ierr + h5_output_double(file_id,"meff_n",Meff_n_tab,dims3,rank)
    if (print_meff_p) ierr = ierr + h5_output_double(file_id,"meff_p",Meff_p_tab,dims3,rank)

    if (print_xn) ierr = ierr + h5_output_double(file_id,"xn",xn_tab,dims3,rank)
    if (print_xp) ierr = ierr + h5_output_double(file_id,"xp",xp_tab,dims3,rank)
    if (print_xa) ierr = ierr + h5_output_double(file_id,"xa",xa_tab,dims3,rank)
    if (print_xh) ierr = ierr + h5_output_double(file_id,"xh",xh_tab,dims3,rank)

  !  ierr = ierr + h5_output_int(file_id,"error",ierror,dims3,rank)

    call h5fclose_f(file_id,error)
    call h5close_f(error)

    ! now add source code and input file to table:
    call include_src_input(trim(adjustl(h5_filename))//C_NULL_CHAR, &
         trim(adjustl(input_space_filename))//C_NULL_CHAR, & 
         trim(adjustl(input_skyrme_filename))//C_NULL_CHAR)

    write(*,*) "Done! :-)"

    RETURN

  END SUBROUTINE write_to_table

END MODULE Write_to_table_Mod
