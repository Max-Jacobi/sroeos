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
module write_table

  USE Kind_Types_Mod
  USE Tables_Input_Mod
  USE Table_Sizes_Mod
  USE hdf5
  USE h5_writer
  USE physical_constants_mod, ONLY : energy_EOS_to_cgs
  use iso_c_binding  ! this is a built-in F2003 module

  implicit none

  contains

  subroutine output

    implicit none

   ! Variables needed for output
    integer :: error, rank, ierr, i, j, k, itemp(1)
    integer(HID_T) :: file_id
    integer(HSIZE_T) :: dims1(1), dims2(2), dims3(3)

    integer :: inse ! this is needed for passing logical to C++
    integer :: isna ! this is needed for passing logical to C++

    ! Write to hdf5
    call h5open_f(error)
    call h5fcreate_f(h5_filename,H5F_ACC_TRUNC_F,file_id,error)

    ! number of points across density, temperature and proton fraction dimensions
    rank = 1

    dims1(1) = 1

    ierr = ierr + h5_output_int(file_id,"pointsrho", (/nrho/), dims1,rank)
    ierr = ierr + h5_output_int(file_id,"pointstemp",(/ntemp/),dims1,rank)
    ierr = ierr + h5_output_int(file_id,"pointsye",  (/nyp/),  dims1,rank)

    ! Here we are outputting an integer flag to indicate to the
    ! interpolation code that the speed of sound in this table
    ! is already divided by the specific enthalpy. This is to distinguish
    ! SRO tables from the legacy tables provided by O'Connor & Ott on
    ! stellarcollapse.org
    itemp = 1
    ierr = ierr + h5_output_int(file_id,"have_rel_cs2", itemp, dims1,rank)

    ierr = ierr + h5_output_double(file_id,"energy_shift", &
         (/energy_shift*energy_EOS_to_cgs/),dims1,rank)

    ! points accros density, temperature and proton fraction dimensions
    rank = 1

    dims1(1) = nrho
    ierr = ierr + h5_output_double(file_id,"logrho",logrho,dims1,rank)

    dims1(1) = ntemp
    ierr = ierr + h5_output_double(file_id,"logtemp",logtemp,dims1,rank)

    dims1(1) = nyp
    ierr = ierr + h5_output_double(file_id,"ye",yp,dims1,rank)

    ! thermo quantities
    rank=3

    dims3(1) = nrho
    dims3(2) = ntemp
    dims3(3) = nyp

! index variable mapping for final table (final_tab):
  !  1 -> press
  !  2 -> energy
  !  3 -> entropy
  !  4 -> muhat
  !  5 -> mun
  !  6 -> mup
  !  7 -> mue
  !  8 -> munu
  !  9 -> dedT
  ! 10 -> dpdrho|e
  ! 11 -> dpde|rho
  ! 12 -> gamma
  ! 13 -> cs²
  ! 14 ->
  ! 15 ->
  ! 16 -> xn
  ! 17 -> xp
  ! 18 -> xa
  ! 19 -> xh
  ! 20 -> abar
  ! 21 -> zbar
  ! 22 -> xl
  ! 23 -> albar
  ! 24 -> zlbar
  ! 25 -> meffn
  ! 26 -> meffp

    ierr = ierr + h5_output_double(file_id,"logpress" ,final_tab(:,:,:,1),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"logenergy",final_tab(:,:,:,2),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"entropy"  ,final_tab(:,:,:,3),dims3,rank)

!    do i = 1, nrho
!      do j = 1, ntemp
!        do k = 1, nyp
!          write (1234,"(3I4,12ES14.6)") i, j, k, final_tab(i,j,k,1)
!        enddo
!      enddo
!    enddo

    ierr = ierr + h5_output_double(file_id,"muhat",final_tab(:,:,:,4),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"mu_n" ,final_tab(:,:,:,5),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"mu_p" ,final_tab(:,:,:,6),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"mu_e" ,final_tab(:,:,:,7),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"munu" ,final_tab(:,:,:,8),dims3,rank)

    ierr = ierr + h5_output_double(file_id,"dedt"   ,final_tab(:,:,:,9) ,dims3,rank)
    ierr = ierr + h5_output_double(file_id,"dpdrhoe",final_tab(:,:,:,10),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"dpderho",final_tab(:,:,:,11),dims3,rank)

    ierr = ierr + h5_output_double(file_id,"gamma",final_tab(:,:,:,12),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"cs2"  ,final_tab(:,:,:,13),dims3,rank)

    ierr = ierr + h5_output_double(file_id,"r",final_tab(:,:,:,14),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"u",final_tab(:,:,:,15),dims3,rank)

    ierr = ierr + h5_output_double(file_id,"Xn",final_tab(:,:,:,16),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"Xp",final_tab(:,:,:,17),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"Xa",final_tab(:,:,:,18),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"Xh",final_tab(:,:,:,19),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"Xl",final_tab(:,:,:,22),dims3,rank)

    ierr = ierr + h5_output_double(file_id,"Abar",final_tab(:,:,:,20),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"Zbar",final_tab(:,:,:,21),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"Albar",final_tab(:,:,:,23),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"Zlbar",final_tab(:,:,:,24),dims3,rank)

    ierr = ierr + h5_output_double(file_id,"meffn",final_tab(:,:,:,25),dims3,rank)
    ierr = ierr + h5_output_double(file_id,"meffp",final_tab(:,:,:,26),dims3,rank)

    call h5fclose_f(file_id,error)
    call h5close_f(error)

    ! copy the SNA and NSE source/input file datasets
    if(.not.only_SNA) then
       inse = 1
    else
       inse = 0
    endif

    if(.not.only_NSE) then
       isna = 1
    else
       isna = 0
    endif

    call copy_src_input(trim(adjustl(nse_eos))//C_NULL_CHAR,&
         trim(adjustl(sna_eos))//C_NULL_CHAR,&
         trim(adjustl(h5_filename))//C_NULL_CHAR,inse,isna)
    
    ! now include MERGE code and input files
    call include_src_input(trim(adjustl(h5_filename))//C_NULL_CHAR,&
         trim(adjustl("src.tar.gz"))//C_NULL_CHAR,&
         trim(adjustl("MERGE-src.tar.gz"))//C_NULL_CHAR)

    call include_src_input(trim(adjustl(h5_filename))//C_NULL_CHAR,&
         trim(adjustl("input/space.in"))//C_NULL_CHAR,&
         trim(adjustl("MERGE-space.in"))//C_NULL_CHAR)

    call include_src_input(trim(adjustl(h5_filename))//C_NULL_CHAR,&
         trim(adjustl("input/tables.in"))//C_NULL_CHAR,&
         trim(adjustl("MERGE-tables.in"))//C_NULL_CHAR)

    call include_src_input(trim(adjustl(h5_filename))//C_NULL_CHAR,&
         trim(adjustl("input/transition.in"))//C_NULL_CHAR,&
         trim(adjustl("MERGE-transition.in"))//C_NULL_CHAR)

  end subroutine output


end module write_table
