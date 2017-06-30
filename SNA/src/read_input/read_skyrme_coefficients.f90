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

MODULE Read_Skyrme_Coefficients_Mod

  USE Kind_Types_Mod, ONLY : DP
  USE Skyrme_Coefficients_Mod
  USE Make_Tables_Mod
  USE Determine_Nuclear_Properties_Mod, ONLY : SATURATION_DENSITY

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE READ_SKYRME_COEFFICIENTS

    USE Physical_Constants_Mod, Only : ZERO, ONE, TWO, HALF, THREE, PI, R_3_2, &
     Mass_n, Mass_p, Neut_Prot_Mass_Diff, B_alpha, v_alpha, M_alpha, &
     Hbarc_Square, N_Q

    USE Global_Variables_Mod
    USE Print_Parameters_Mod

    IMPLICIT NONE

    REAL(DP) :: Dummy_a, Dummy_b, n_sat, e_bind, stiffening_ratio, &
    stiffening_exponent, stiffening_symmetry
    CHARACTER(LEN=128) :: command

    NAMELIST /TABLE_INPUT/ make_hdf5_table, write_solutions_to_file, & 
                           redefine_print_parameters, &
                           filename_base, output_directory

    NAMELIST /TABLE_OUTPUT/ print_p, print_s, print_e, &
                            print_muh, print_mun, print_mup, &
                            print_dpdn, print_dsdn, print_dmuhdn, &
                            print_dpdt, print_dsdt, print_dmuhdt, &
                            print_dpdy, print_dsdy, print_dmuhdy, &
                            print_abar, print_zbar, print_r, print_u, &
                            print_xn, print_xp, print_xa, print_xh, &
                            print_meff_n, print_meff_p

    NAMELIST /SKYRME_TYPE/ Skyrme_parametrization_type, non_local_terms, &
                           Use_default_constants, High_density_stiffening

    NAMELIST /SKYRME_COEFFICIENTS/ Coeff_a, Coeff_b, Coeff_c, Coeff_d,      &
                                   Coeff_alpha1, Coeff_alpha2, Coeff_delta, &
                                   Coeff_t0, Coeff_t1, Coeff_t2, Coeff_t3,  &
                                   Coeff_x0, Coeff_x1, Coeff_x2, Coeff_x3,  &
                                   Coeff_sigma,                             &
                                   Coeff_qnn, Coeff_qpp, Coeff_qnp, Coeff_qpn

    NAMELIST /SKYRME_STIFFENING/   stiffening_ratio, stiffening_exponent,  &
                                   stiffening_symmetry


!   Set default values for some parameters
    make_hdf5_table  = .TRUE.
    write_solutions_to_file = .FALSE.
    redefine_print_parameters = .FALSE.
    filename_base = 'EOS'
    output_directory = 'EOS'

    Use_default_constants = 0
    High_density_stiffening = .FALSE.
    non_local_terms = 1

    OPEN(10,FILE=trim(adjustl(input_skyrme_filename)), &
         FORM='FORMATTED',STATUS='OLD',ACTION='READ')
    READ(10,NML=TABLE_INPUT)
!   create output directory (may be system/compiler dependent)
    command = ADJUSTL('mkdir ') // ADJUSTL(TRIM(output_directory))
    CALL EXECUTE_COMMAND_LINE (command)
    IF (write_solutions_to_file) THEN
!     create output directories for solutions
      command = ADJUSTL('mkdir ') // ADJUSTL(TRIM(output_directory))// ADJUSTL(TRIM('/SOL'))
      CALL EXECUTE_COMMAND_LINE (command)
      command = ADJUSTL('mkdir ') // ADJUSTL(TRIM(output_directory))// ADJUSTL(TRIM('/NO_SOL'))
      CALL EXECUTE_COMMAND_LINE (command)
    ENDIF
    READ(10,NML=SKYRME_TYPE)

    SELECT CASE (Use_default_constants)
!     Use Default values for neutron and proton masses and alpha particle
!     binding energy
      CASE(0)
        Mass_n  = 939.56540D0
        Mass_p  = 938.27204D0
        Neut_Prot_Mass_Diff = Mass_n - Mass_p
        b_alpha = 28.29552D0 + TWO*Neut_Prot_Mass_Diff
        v_alpha = 24.0D0
      CASE(1)
        Mass_n  = 939.5654D0
        Mass_p  = 939.5654D0
        Neut_Prot_Mass_Diff = 1.29D0
        b_alpha = 28.3D0
        v_alpha = 24.0D0
    END SELECT
    !
    N_Q   = (Mass_n/(TWO*PI*Hbarc_Square))**R_3_2
    M_alpha = 2.D0*(Mass_n + Mass_p) - B_alpha

!   if EoS is stiffened at high densities an extra non-local term is added
    IF (High_density_stiffening) non_local_terms = non_local_terms + 1

!   Allocate standard Skyrme coefficient sizes
    ALLOCATE(Coeff_t3(non_local_terms))
    ALLOCATE(Coeff_x3(non_local_terms))
    ALLOCATE(Coeff_sigma(non_local_terms))

!   Allocate simplified Skyrme coefficient sizes
    ALLOCATE(Coeff_c(non_local_terms))
    ALLOCATE(Coeff_d(non_local_terms))
    ALLOCATE(Coeff_delta(non_local_terms))

!   set all Skyrme coefficients to zero
    Coeff_a = zero
    Coeff_b = zero
    Coeff_c = zero
    Coeff_d = zero
    Coeff_alpha1 = zero
    Coeff_alpha2 = zero
    Coeff_delta  = zero
    Coeff_t0 = zero
    Coeff_t1 = zero
    Coeff_t2 = zero
    Coeff_t3 = zero
    Coeff_x0 = zero
    Coeff_x1 = zero
    Coeff_x2 = zero
    Coeff_x3 = zero
    Coeff_sigma = zero
    Coeff_qnn = zero
    Coeff_qpp = zero
    Coeff_qnp = zero
    Coeff_qpn = zero

!   read Skyrme coefficients
    READ(10,NML=SKYRME_COEFFICIENTS)

!   if Standard Skyrme given then convert to Simplified Skyrme
    SELECT CASE (Skyrme_parametrization_type)
!     Otherwise, convert from standard Skyrme coefficients to simplified ones
      CASE (0)
        Dummy_a = Coeff_t1*(Coeff_x1 + two) + Coeff_t2*(Coeff_x2 + two)
        Dummy_b = half*( Coeff_t2*(two*Coeff_x2 + one) &
                       - Coeff_t1*(two*Coeff_x1 + one) )

        Coeff_alpha1 = (Dummy_a+two*Dummy_b)/8.d0
        Coeff_alpha2 = Dummy_a/8.d0

        Coeff_delta  = Coeff_sigma + one

        Coeff_a     =  Coeff_t0*(one-one*Coeff_x0)/4.d0
        Coeff_b     =  Coeff_t0*(one+two*Coeff_x0)/8.d0
        Coeff_c     =  Coeff_t3*(one - Coeff_x3)/24.d0
        Coeff_d     =  Coeff_t3*(one + two*Coeff_x3)/48.d0

        Coeff_qnn = 3.d0/16.d0*( Coeff_t1*(one-Coeff_x1) &
                               - Coeff_t2*(one+Coeff_x2) )
        Coeff_qnp = 1.d0/ 8.d0*(three*Coeff_t1*(one+Coeff_x1/two) &
                               -      Coeff_t2*(one+Coeff_x2/two) )
        Coeff_qpn = Coeff_qnp
        Coeff_qpp = Coeff_qnn
!     If Skyrme in Simplified form already, set stardard coefficients to zero
!      This is only done for printing purposes.
      CASE(1)
        CONTINUE
      CASE DEFAULT
        STOP 'Skyrme_parametrization_type must be 0 or 1!'
    END SELECT

!   If high density stiffening chosen then compute its terms
!   based on SKYRME_STIFFENING namelist
    SELECT CASE (High_density_stiffening)
      CASE (.FALSE.)
        CONTINUE
      CASE(.TRUE.)
        IF (non_local_terms>2) STOP ' High_density_stiffening only implemented &
                                    for non_local_terms = 1.'
        READ(10,NML=SKYRME_STIFFENING)

        CALL SATURATION_DENSITY(n_sat,e_bind)

        Coeff_delta(2) = stiffening_exponent

        Dummy_a = coeff_delta(1)*(coeff_delta(1)-one)*(coeff_c(1)+coeff_d(1))
        Dummy_b = coeff_delta(2)*(coeff_delta(2)-one)

        Dummy_a = stiffening_ratio*Dummy_a/Dummy_b*&
                  n_sat**(coeff_delta(1)-coeff_delta(2))

        Coeff_c(2) = stiffening_symmetry*Dummy_a
        Coeff_d(2) = (one-stiffening_symmetry)*Dummy_a
    END SELECT

    CALL PRINT_SKYRME_COEFFICIENTS

!   set all print variables to true
    print_p = .TRUE.
    print_s = .TRUE.
    print_e = .TRUE.

    print_muh = .TRUE.
    print_mun = .TRUE.
    print_mup = .TRUE.

    print_dpdn = .TRUE.
    print_dsdn = .TRUE.
    print_dmuhdn = .TRUE.

    print_dpdt = .TRUE.
    print_dsdt = .TRUE.
    print_dmuhdt = .TRUE.

    print_dpdy = .TRUE.
    print_dsdy = .TRUE.
    print_dmuhdy = .TRUE.

    print_abar = .TRUE.
    print_zbar = .TRUE.

    print_r = .TRUE.
    print_u = .TRUE.

    print_xn = .TRUE.
    print_xp = .TRUE.
    print_xa = .TRUE.
    print_xh = .TRUE.

    print_meff_n = .TRUE.
    print_meff_p = .TRUE.

!   If redefining print parameters only print variables set to true
!   in TABLE_OUTPUT namelist
    SELECT CASE (redefine_print_parameters)
      CASE (.FALSE.)
        CONTINUE
      CASE(.TRUE.)
        READ(10,NML=TABLE_OUTPUT)
    END SELECT

    CLOSE(10)

  END SUBROUTINE READ_SKYRME_COEFFICIENTS


  SUBROUTINE PRINT_SKYRME_COEFFICIENTS

    IMPLICIT NONE

    INTEGER(I4B) :: i
    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT   = "(1A10,1ES20.12,1A10)"
    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT_1 = "(1A2,1I1,1A7,1ES20.12,1A18)"
    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT_2 = "(1A6,1I1,1A3,1ES20.12,1A20)"
    CHARACTER(*), PARAMETER :: OUTPUT_FORMAT_3 = "(1A3,1I1,1A6,1ES20.12,1A20)"
    CHARACTER(LEN=128) :: File_Name


    File_Name = ADJUSTL(TRIM(output_directory)) // &
                ADJUSTL(TRIM('/Skyrme_coefficients.dat'))

    OPEN(10,file=File_Name)

!   Write standard Skyrme coefficients
    IF (Skyrme_parametrization_type == 1) then
      WRITE(10,*)
      WRITE(10,*) "Used simplified Skyrme coefficients"
      WRITE(10,*)
    ELSE
      WRITE(10,*)
      WRITE(10,*) "Used standard Skyrme coefficients"
      WRITE(10,*)

      WRITE(10,OUTPUT_FORMAT) 't_0     = ', Coeff_t0, 'MeV fm^3'
      WRITE(10,OUTPUT_FORMAT) 't_1     = ', Coeff_t1, 'MeV fm^5'
      WRITE(10,OUTPUT_FORMAT) 't_2     = ', Coeff_t2, 'MeV fm^5'

      DO i = 1, SIZE(Coeff_t3)
        WRITE(10,OUTPUT_FORMAT_3) 't_3',i,'    = ', Coeff_t3(i), 'MeV fm^(3+3*sigma)'
      ENDDO

      WRITE(10,OUTPUT_FORMAT) 'x_0     = ', Coeff_x0
      WRITE(10,OUTPUT_FORMAT) 'x_1     = ', Coeff_x1
      WRITE(10,OUTPUT_FORMAT) 'x_2     = ', Coeff_x2

      DO i = 1, SIZE(Coeff_x3)
        WRITE(10,OUTPUT_FORMAT_3) 'x_3',i,'    = ', Coeff_x3(i)
      ENDDO

      DO i = 1, SIZE(Coeff_sigma)
        WRITE(10,OUTPUT_FORMAT_2) 'sigma_',i,' = ', Coeff_sigma(i)
      ENDDO
      WRITE(10,*)
      WRITE(10,*) "Convert to simplfied Skyrme coefficients"
      WRITE(10,*)
    ENDIF

!   Write simplified Skyrme coefficients

    WRITE(10,OUTPUT_FORMAT) 'a     = ', Coeff_a, 'MeV fm^3'
    WRITE(10,OUTPUT_FORMAT) 'b     = ', Coeff_b, 'MeV fm^5'

    DO i = 1, SIZE(Coeff_c)
      WRITE(10,OUTPUT_FORMAT_1) 'c_',i,'     = ', Coeff_c(i), 'MeV fm^(3*delta)'
      WRITE(10,OUTPUT_FORMAT_1) 'd_',i,'     = ', Coeff_d(i), 'MeV fm^(3*delta)'
    ENDDO

    DO i = 1, SIZE(Coeff_sigma)
      WRITE(10,OUTPUT_FORMAT_2) 'delta_ ',i,' = ', Coeff_delta(i)
    ENDDO

    WRITE(10,OUTPUT_FORMAT) 'alpha_1 = ', Coeff_alpha1, ' '
    WRITE(10,OUTPUT_FORMAT) 'alpha_2 = ', Coeff_alpha2, ' '

    WRITE(10,OUTPUT_FORMAT) 'q_nn    = ', Coeff_qnn, ' MeV fm^5'
    WRITE(10,OUTPUT_FORMAT) 'q_np    = ', Coeff_qnp, ' MeV fm^5'
    WRITE(10,OUTPUT_FORMAT) 'q_pn    = ', Coeff_qpn, ' MeV fm^5'
    WRITE(10,OUTPUT_FORMAT) 'q_pp    = ', Coeff_qpp, ' MeV fm^5'

  END SUBROUTINE PRINT_SKYRME_COEFFICIENTS


END MODULE Read_Skyrme_Coefficients_Mod
