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
MODULE Read_Table_Mod

  USE Table_Sizes_Mod
  USE HDF5
  USE Phase_Space_Input_Mod
  USE Transition_Input_Mod, ONLY : Log10nt_max, Log10nt_min
  USE Tables_Input_Mod, ONLY : only_sna, only_NSE

  IMPLICIT NONE

CONTAINS

  SUBROUTINE READTABLE(eos_filename,flag)
  ! This routine reads the table and initializes
  ! all variables in the module.

    IMPLICIT NONE

    CHARACTER(*) flag
    CHARACTER(*) eos_filename
    CHARACTER(len=100) message

!   HDF5 vars
    INTEGER(HID_T)   :: file_id, dset_id, dspace_id
    INTEGER(HSIZE_T) :: dims1(1), dims2(2), dims3(3)
    INTEGER(I4B) ::  error, rank, accerr
    INTEGER(I4B) ::  i, j, k

    REAL(DP) buffer1,buffer2,buffer3,buffer4

    accerr = 0

    eos_filename_stored = eos_filename

    WRITE (*,"(A8,A6,A12)") "Reading", trim(flag) ," EOS Table"

    CALL h5open_f(error)
    CALL h5fopen_f(trim(adjustl(eos_filename)), H5F_ACC_RDONLY_F, &
                                                    file_id, error)

    WRITE (6,*) 'path to table = ', trim(adjustl(eos_filename))

!   read scalars
    dims1(1)=1
    CALL h5dopen_f(file_id, "pointsrho", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nrho, dims1, error)
    CALL h5dclose_f(dset_id,error)

    IF (error.ne.0) THEN
      STOP "Could not read EOS table file"
    ENDIF

    dims1(1)=1
    CALL h5dopen_f(file_id, "pointstemp", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, ntemp, dims1, error)
    CALL h5dclose_f(dset_id,error)

    IF (error.ne.0) THEN
      STOP "Could not read EOS table file"
    ENDIF

    dims1(1)=1
    CALL h5dopen_f(file_id, "pointsye", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nyp, dims1, error)
    CALL h5dclose_f(dset_id,error)

    IF (error.ne.0) THEN
      STOP "Could not read EOS table file"
    ENDIF

    WRITE (message,"(a25,i5,i5,i5)") "We have nrho ntemp nyp: ", nrho, ntemp, nyp
    WRITE (*,*) message

    IF (ALLOCATED(alltables)) DEALLOCATE(alltables)

    SELECT CASE (flag)
      CASE ("NSE")
        ALLOCATE(alltables(nrho,ntemp,nyp,nse_nvars))
        ALLOCATE(nse_table(nrho,ntemp,nyp,nse_nvars))
        ALLOCATE(nse_logrho(nrho),nse_logtemp(ntemp),nse_yp(nyp))
        nse_nrho  = nrho
        nse_ntemp = ntemp
        nse_nyp   = nyp
      CASE ("SNA")
        ALLOCATE(alltables(nrho,ntemp,nyp,sna_nvars))
        ALLOCATE(sna_table(nrho,ntemp,nyp,sna_nvars))
        ALLOCATE(sna_logrho(nrho),sna_logtemp(ntemp),sna_yp(nyp))
        sna_nrho  = nrho
        sna_ntemp = ntemp
        sna_nyp   = nyp
    END SELECT

  ! index variable mapping:
  !  1 -> press
  !  2 -> energy
  !  3 -> entropy
  !  4 -> muh
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
  ! 20 -> abar
  ! 21 -> zbar
  ! 22 -> xl
  ! 23 -> albar
  ! 24 -> zlbar
  ! 25 -> rad
  ! 26 -> u

    dims3(1)=nrho
    dims3(2)=ntemp
    dims3(3)=nyp

!   get density array
    IF (ALLOCATED(logrho)) DEALLOCATE(logrho)
    ALLOCATE(logrho(nrho))
    dims1(1)=nrho
    CALL h5dopen_f(file_id, "logn", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logrho, dims1, error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
!   check if array bounds are within accepted range
!    TODO: Make this check and similar ones below a subroutine
    IF (only_sna) THEN
      IF (Log10n_min < logrho(1)) THEN
        WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10n_min:"
        WRITE (*,*) "  Minimum density in SNA table higher than "
        WRITE (*,*) "  minimum density set in MERGE code.       "
        WRITE (*,*) "  Log10n_min (MERGE) =  ", Log10n_min
        WRITE (*,*) "  Log10n_min (SNA)   =  ", logrho(1)
        STOP
      ENDIF
      IF (Log10n_max > logrho(dims1(1))) THEN
        WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10n_max:"
        WRITE (*,*) "  Maximum density in SNA table lower than "
        WRITE (*,*) "  maximum density set in MERGE code.       "
        WRITE (*,*) "  Log10n_max (MERGE) =  ", Log10n_max
        WRITE (*,*) "  Log10n_max (SNA)   =  ", logrho(dims1)
        STOP
      ENDIF
    ENDIF

    IF (only_nse) THEN
      IF (Log10n_min < logrho(1)) THEN
        WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10n_min:"
        WRITE (*,*) "  Minimum density in NSE table higher than "
        WRITE (*,*) "  minimum density set in MERGE code.       "
        WRITE (*,*) "  Log10n_min (MERGE) =  ", Log10n_min
        WRITE (*,*) "  Log10n_min (NSE)   =  ", logrho(1)
        STOP
      ENDIF
      IF (Log10n_max > logrho(dims1(1))) THEN
        WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10n_max:"
        WRITE (*,*) "  Maximum density in NSE table lower than "
        WRITE (*,*) "  maximum density set in MERGE code.       "
        WRITE (*,*) "  Log10n_max (MERGE) =  ", Log10n_max
        WRITE (*,*) "  Log10n_max (NSE)   =  ", logrho(dims1(1))
        STOP
      ENDIF
    ENDIF

    IF (.not.only_sna .AND. .not.only_NSE) THEN
      SELECT CASE(flag)
        CASE("NSE")
          IF (Log10n_min < logrho(1)) THEN
            WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10n_min:"
            WRITE (*,*) "  Minimum density in NSE table higher than "
            WRITE (*,*) "  minimum density set in MERGE code.       "
            WRITE (*,*) "  Log10n_min (MERGE) =  ", Log10n_min
            WRITE (*,*) "  Log10n_min (NSE)   =  ", logrho(1)
            STOP
          ENDIF
          IF (Log10nt_max > logrho(dims1(1))) THEN
            WRITE (*,*) "ERROR SETTING MERGE PARAMETERS."
            WRITE (*,*) "  Maximum density in NSE table lower than "
            WRITE (*,*) "  upper threshold density needed for transition."
            WRITE (*,*) "  Log10nt_max (MERGE) =  ", Log10nt_max
            WRITE (*,*) "  Log10n_max  (NSE)   =  ", logrho(dims1(1))
            WRITE (*,*) "  Check User_Guide Section 6.2.3 for discussion"
            STOP
          ENDIF
        CASE("SNA")
          IF (Log10n_max > logrho(dims1(1))) THEN
            WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10n_max:"
            WRITE (*,*) "  Maximum density in SNA table lower than "
            WRITE (*,*) "  maximum density set in MERGE code.       "
            WRITE (*,*) "  Log10n_max (MERGE) =  ", Log10n_max
            WRITE (*,*) "  Log10n_max (SNA)   =  ", logrho(dims1(1))
            STOP
          ENDIF
          IF (Log10nt_min < logrho(1)) THEN
            WRITE (*,*) "ERROR SETTING MERGE PARAMETERS."
            WRITE (*,*) "  Minimum density in SNA table higher than "
            WRITE (*,*) "  lower threshold density needed for transition."
            WRITE (*,*) "  Log10nt_in (MERGE) =  ", Log10nt_min
            WRITE (*,*) "  Log10n_min (SNA)   =  ", logrho(1)
            WRITE (*,*) "  Check User_Guide Section 6.2.3 for discussion"
            STOP
          ENDIF
       END SELECT
    ENDIF

!   get temperature array
    if (ALLOCATED(logtemp)) DEALLOCATE(logtemp)
    ALLOCATE(logtemp(ntemp))
    dims1(1)=ntemp
    CALL h5dopen_f(file_id, "logtemp", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, logtemp, dims1, error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
!   check if array bounds are within accepted range
!    TODO: Make this check and similar ones below a subroutine
    SELECT CASE(flag)
      CASE("SNA")
        IF (Log10T_min < logtemp(1)) THEN
          WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10T_min:"
          WRITE (*,*) "  Minimum temperature in SNA table higher than "
          WRITE (*,*) "  minimum temperature set in MERGE code.       "
          WRITE (*,*) "  Log10T_min (MERGE) =  ", Log10T_min
          WRITE (*,*) "  Log10T_min (SNA)   =  ", logtemp(1)
          STOP
        ENDIF
        IF (Log10T_max > logtemp(dims1(1))) THEN
          WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10T_max:"
          WRITE (*,*) "  Maximum temperature in SNA table lower than "
          WRITE (*,*) "  maximum temperature set in MERGE code.       "
          WRITE (*,*) "  Log10T_max (MERGE) =  ", Log10T_max
          WRITE (*,*) "  Log10T_max (SNA)   =  ", logtemp(dims1(1))
          STOP
        ENDIF
      CASE("NSE")
        IF (Log10T_min < logtemp(1)) THEN
          WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10T_min:"
          WRITE (*,*) "  Minimum temperature in NSE table higher than "
          WRITE (*,*) "  minimum temperature set in MERGE code.       "
          WRITE (*,*) "  Log10T_min (MERGE) =  ", Log10T_min
          WRITE (*,*) "  Log10T_min (NSE)   =  ", logtemp(1)
          STOP
        ENDIF
        IF (Log10T_max > logtemp(dims1(1))) THEN
          WRITE (*,*) "ERROR SETTING MERGE PARAMETER Log10T_max:"
          WRITE (*,*) "  Maximum temperature in NSE table lower than "
          WRITE (*,*) "  maximum temperature set in MERGE code.       "
          WRITE (*,*) "  Log10T_max (MERGE) =  ", Log10T_max
          WRITE (*,*) "  Log10T_max (NSE)   =  ", logtemp(dims1(1))
          STOP
        ENDIF
    END SELECT

!   get proton fraction array
    if (ALLOCATED(yp)) DEALLOCATE(yp)
    ALLOCATE(yp(nyp))
    dims1(1)=nyp
    CALL h5dopen_f(file_id, "ye", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, yp, dims1, error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
!   check if array bounds are within accepted range
!    TODO: Make this check and similar ones below a subroutine

    SELECT CASE(flag)
      CASE("SNA")
        IF (yp_min < yp(1)) THEN
          WRITE (*,*) "ERROR SETTING MERGE PARAMETER yp_min:"
          WRITE (*,*) "  Minimum proton fraction in SNA table higher than "
          WRITE (*,*) "  minimum proton fraction set in MERGE code.       "
          WRITE (*,*) "  yp_min (MERGE) =  ", yp_min
          WRITE (*,*) "  yp_min (SNA)   =  ", yp(1)
          STOP
        ENDIF
        IF (yp_max > yp(dims1(1))) THEN
          WRITE (*,*) "ERROR SETTING MERGE PARAMETER yp_max:"
          WRITE (*,*) "  Maximum proton fraction in SNA table lower than "
          WRITE (*,*) "  maximum proton fraction set in MERGE code.       "
          WRITE (*,*) "  yp_max (MERGE) =  ", yp_max
          WRITE (*,*) "  yp_max (SNA)   =  ", yp(dims1(1))
          STOP
        ENDIF
      CASE("NSE")
        IF (yp_min < yp(1)) THEN
          WRITE (*,*) "ERROR SETTING MERGE PARAMETER yp_min:"
          WRITE (*,*) "  Minimum proton fraction in NSE table higher than "
          WRITE (*,*) "  minimum proton fraction set in MERGE code.       "
          WRITE (*,*) "  yp_min (MERGE) =  ", yp_min
          WRITE (*,*) "  yp_min (NSE)   =  ", yp(1)
          STOP
        ENDIF
        IF (yp_max > yp(dims1(1))) THEN
          WRITE (*,*) "ERROR SETTING MERGE PARAMETER yp_max:"
          WRITE (*,*) "  Maximum proton fraction in NSE table lower than "
          WRITE (*,*) "  maximum proton fraction set in MERGE code.       "
          WRITE (*,*) "  yp_max (MERGE) =  ", yp_max
          WRITE (*,*) "  yp_max (NSE)   =  ", yp(dims1(1))
          STOP
        ENDIF
    END SELECT


! pressure, entropy, energy
    CALL h5dopen_f(file_id, "p", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,1), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "e", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,2), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "s", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,3), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error

! chemical potentials
    CALL h5dopen_f(file_id, "muh", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,4), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "mun", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,5), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "mup", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,6), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error

! pressure derivatives
    CALL h5dopen_f(file_id, "dpdn", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,7), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "dpdt", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,8), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "dpdy", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,9), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error

! entropy derivatives
    CALL h5dopen_f(file_id, "dsdn", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,10), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "dsdt", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,11), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "dsdy", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,12), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error

! chemical potential derivatives
    CALL h5dopen_f(file_id, "dmudn", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,13), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "dmudt", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,14), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "dmudy", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,15), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error

! composition
    CALL h5dopen_f(file_id, "xn", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,16), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "xp", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,17), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "xa", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,18), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error
    CALL h5dopen_f(file_id, "xh", dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,19), dims3,error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error

    select case(flag)
      case("NSE")
        !write (*,*) 'read xl'
        call h5dopen_f(file_id, "xl", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,22),dims3,error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
      ! average light nucleus mass and charge
        !write (*,*) 'read al'
        call h5dopen_f(file_id, "al", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,23),dims3,error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        !write (*,*) 'read zl'
        call h5dopen_f(file_id, "zl", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,24),dims3,error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        ! average heavy nucleus mass and charge
        call h5dopen_f(file_id, "ah", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,20),dims3,error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        call h5dopen_f(file_id, "zh", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,21),dims3,error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        alltables(:,:,:,25:26) = 0.d0
      case("SNA")
        ! average nucleus
        call h5dopen_f(file_id, "abar", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,20),dims3,error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        call h5dopen_f(file_id, "zbar", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,21),dims3,error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        alltables(:,:,:,22:24) = 0.d0
        call h5dopen_f(file_id, "r", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,25),dims3,error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
        call h5dopen_f(file_id, "u", dset_id, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, alltables(:,:,:,26),dims3,error)
        call h5dclose_f(dset_id,error)
        accerr=accerr+error
    end select

    IF (accerr.ne.0) THEN
      STOP "Problem reading EOS table file."
    ENDIF

    CALL h5fclose_f (file_id,error)
    CALL h5close_f (error)

    SELECT CASE (flag)
      CASE ("NSE")
        nse_table   = alltables
        nse_logrho  = logrho
        nse_logtemp = logtemp
        nse_yp      = yp
       ! set min-max values:
        nse_eos_rhomin = 10.0d0**logrho(1)
        nse_eos_rhomax = 10.0d0**logrho(nrho)
        WRITE (*,*)  nse_eos_rhomin, nse_eos_rhomax
        nse_eos_ypmin = yp(1)
        nse_eos_ypmax = yp(nyp)
        WRITE (*,*)  nse_eos_ypmin, nse_eos_ypmax
        nse_eos_tempmin = 10.0d0**logtemp(1)
        nse_eos_tempmax = 10.0d0**logtemp(ntemp)
        WRITE (*,*)  nse_eos_tempmin, nse_eos_tempmax
      CASE("SNA")
        sna_table   = alltables
        sna_logrho  = logrho
        sna_logtemp = logtemp
        sna_yp      = yp
       ! set min-max values:
        sna_eos_rhomin = 10.0d0**logrho(1)
        sna_eos_rhomax = 10.0d0**logrho(nrho)
        WRITE (*,*)  sna_eos_rhomin, sna_eos_rhomax
        sna_eos_ypmin = yp(1)
        sna_eos_ypmax = yp(nyp)
        WRITE (*,*)  sna_eos_ypmin, sna_eos_ypmax
        sna_eos_tempmin = 10.0d0**logtemp(1)
        sna_eos_tempmax = 10.0d0**logtemp(ntemp)
        WRITE (*,*)  sna_eos_tempmin, sna_eos_tempmax
    END SELECT

    WRITE (6,*) "Done reading eos tables", (size(alltables,i), i=1,4)
    WRITE (6,*) "Done reading eos tables", (size(sna_table,i), i=1,4)
    WRITE (6,*) "Done reading eos tables", (size(nse_table,i), i=1,4)

    DEALLOCATE(alltables)

  END SUBROUTINE READTABLE

  SUBROUTINE add_index(string_to_add,index_in_alltables)

    IMPLICIT NONE

    CHARACTER(*) string_to_add
    INTEGER(I4B) :: index_in_alltables

    REAL(DP), allocatable :: temp_table(:,:,:,:)

  ! HDF5 vars
    INTEGER(HID_T) :: file_id,dset_id
    INTEGER(HSIZE_T) :: dims3(3)
    INTEGER(I4B) :: error,accerr

  !set index of newest table member
    index_in_alltables = nvars+1

  !increase table size
    nvars = nvars+1

  !increase size of alltables
    ALLOCATE(temp_table(nrho,ntemp,nyp,nvars))
    temp_table(:,:,:,1:nvars-1) = alltables
    DEALLOCATE(alltables)
    ALLOCATE(alltables(nrho,ntemp,nyp,nvars))
    alltables = temp_table
    DEALLOCATE(temp_table)

  !reopen hdf5 file
    accerr=0

    WRITE (*,*) "Reading EOS Table again to add index",index_in_alltables

    CALL h5open_f(error)

    CALL h5fopen_f (trim(adjustl(eos_filename_stored)), H5F_ACC_RDONLY_F,file_id, error)

    dims3(1)=nrho
    dims3(2)=ntemp
    dims3(3)=nyp
    CALL h5dopen_f(file_id, trim(adjustl(string_to_add)), dset_id, error)
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE,alltables(:,:,:,index_in_alltables), dims3, error)
    CALL h5dclose_f(dset_id,error)
    accerr=accerr+error

    IF (accerr.ne.0) THEN
      STOP "Problem reading EOS table file."
    ENDIF

    CALL h5fclose_f (file_id,error)

    CALL h5close_f (error)

  END SUBROUTINE add_index

END MODULE Read_Table_Mod
