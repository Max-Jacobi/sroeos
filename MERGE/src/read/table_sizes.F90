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
MODULE Table_Sizes_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP

  IMPLICIT NONE

  SAVE

  INTEGER(I4B) ::     nrho,    ntemp,    nyp
  INTEGER(I4B) :: nse_nrho,nse_ntemp,nse_nyp
  INTEGER(I4B) :: sna_nrho,sna_ntemp,sna_nyp

  REAL(DP) :: precision = 1.0d-9

  character(len=256) eos_filename_stored

! min-max values for tables
  REAL(DP) :: sna_eos_rhomin,  sna_eos_rhomax
  REAL(DP) :: sna_eos_ypmin,   sna_eos_ypmax
  REAL(DP) :: sna_eos_tempmin, sna_eos_tempmax

  REAL(DP) :: nse_eos_rhomin,  nse_eos_rhomax
  REAL(DP) :: nse_eos_ypmin,   nse_eos_ypmax
  REAL(DP) :: nse_eos_tempmin, nse_eos_tempmax

  REAL(DP) :: eos_rhomin,  eos_rhomax
  REAL(DP) :: eos_ypmin,   eos_ypmax
  REAL(DP) :: eos_tempmin, eos_tempmax

  REAL(DP) :: t_max_hack = 320.0d0
  REAL(DP) :: energy_shift = 20.D0
!$threadprivate(energy_shift,t_max_hack)

! number of variables on each table
  INTEGER(I4B) :: nvars     = 28
  INTEGER(I4B) :: sna_nvars = 28
  INTEGER(I4B) :: nse_nvars = 28
  INTEGER(I4B) :: nvars_tab = 28

! allocate table sizes
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: sna_table, nse_table
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: alltables, final_tab

  REAL(DP), DIMENSION(:), ALLOCATABLE ::  sna_yp, sna_logrho, sna_logtemp, &
                                          nse_yp, nse_logrho, nse_logtemp, &
                                          ele_yp, ele_logrho, ele_logtemp, &
                                              yp,     logrho,     logtemp

  ! index variable mapping for sna (nse) tables:
  !  1 -> press
  !  2 -> energy
  !  3 -> entropy
  !  4 -> muhat
  !  5 -> mun
  !  6 -> mup
  !  7 -> dpdn
  !  8 -> dpdt
  !  9 -> dpdy
  ! 10 -> dsdn
  ! 11 -> dsdt
  ! 12 -> dsdy
  ! 13 -> dmuhdn
  ! 14 -> dmuhdt
  ! 15 -> dmuhdy
  ! 16 -> xn
  ! 17 -> xp
  ! 18 -> xa
  ! 19 -> xh
  ! 20 -> ah
  ! 21 -> zh
  ! 22 -> xl
  ! 23 -> al
  ! 24 -> zl

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
  ! 20 -> ahbar
  ! 21 -> zhbar
  ! 22 -> xl
  ! 23 -> albar
  ! 24 -> zlbar

END MODULE Table_Sizes_Mod
