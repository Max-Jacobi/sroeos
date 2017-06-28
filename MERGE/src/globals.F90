!#############################################################################
!
!   Module that contains parameters and constants

!   author: Andre da Silva Schneider
!   version: 0.0
!   last update: 06/11/2015
!   date:        02/11/2015
!
!#############################################################################

module globals

  use Kind_Types_Mod

  implicit none

  save
! some real numbers from 0. to 10. and pi
  REAL(DP), PARAMETER :: ZERO = 0.0_SP, ONE  = 1.0_SP, TWO = 2.0_SP, THREE = 3.0_SP
  REAL(DP), PARAMETER :: FOUR = 4.0_SP, FIVE = 5.0_SP, SIX = 6.0_SP, TEN = 10.0_SP
  REAL(DP), PARAMETER :: HALF = 0.5_SP, TWOPI2 = TWO*PI*PI, PI2 = PI*PI

! ratios used throughout
  REAL(DP), PARAMETER :: E23 = TWO/THREE, E32 = THREE/TWO, &
                         E35 = THREE/FIVE,E53 = FIVE/THREE, &
                         E52 = FIVE/TWO, E13 = ONE/THREE, &
                         E56 = FIVE/SIX, E92 = THREE*E32, &
                         E715 = 7.d0/15.d0

  CHARACTER(100) :: base, outdir

! constants hbar*c, alpha, masses of proton, neutron and electron (in MeV)
  REAL(DP), PARAMETER :: HC = 197.327_SP,  HC2 = HC*HC, ALPHA = ONE/137.036_SP
  REAL(DP), PARAMETER :: BETA_FAC = (THREE**E53/FIVE**E13)*(PI*ALPHA*HC)**E13

  REAL(DP) :: MEC2 = 0.510999_SP
  REAL(DP) :: MNC2
  REAL(DP) :: MPC2
! volume (fm) and binding energy of alpha particles (MeV)
  REAL(DP) :: v_alpha
  REAL(DP) :: b_alpha

! bulk nuclear matter parameters (may_be changed in bulk_inlist)
  REAL(DP) :: Bind  ! binding energy (MeV)
  REAL(DP) :: n_sat ! saturation density (fm^-3)
  REAL(DP) :: delta ! neutron-proton mass difference (MeV)
  REAL(DP) :: P_S   ! pressure for SNM at saturation density (MeV) (should be zero)
  REAL(DP) :: P_PNM ! pressure for PNM at saturation density (MeV)
  REAL(DP) :: K_0   ! nuclear compressibility (MeV)
  REAL(DP) :: Q_0   ! nuclear skewness coefficient (MeV)
  REAL(DP) :: S_V   ! nuclear symmetry energy parameter (MeV)
  REAL(DP) :: J_S   ! nuclear symmetry energy parameter (MeV) same as S_V
  REAL(DP) :: L_S   ! nuclear symmetry energy density dependence (MeV)
  REAL(DP) :: K_S   ! nuclear symmetry incompressibility (MeV)
  REAL(DP) :: K_TV  ! isospin dependence of nuclear incompressibility (MeV)
  REAL(DP) :: Q_S   ! nuclear symmetry skewness coefficient (MeV)
  REAL(DP) :: S_half! nuclear symmetry energy at half saturation density (MeV)
  REAL(DP) :: A_V   ! Bulk level density parameter (MeV^-1)

! surface parameters
  REAL(DP) :: SIGMA_S ! Surface tension of nuclear matter (MeV fm^-3)
  REAL(DP) :: S_S     ! Surface symmetry energy parameter (MeV^-1)
  REAL(DP) :: A_S     ! Surface level density parameter (MeV^-1)

! from nuclear matter parametrization
  REAL(DP) :: m_star, coeff_a, coeff_b, coeff_c, coeff_d, coeff_e, coeff_f !
  REAL(DP) :: coeff_alpha, coeff_delta, coeff_gamma, K_change, ef_ratio
  REAL(DP) :: coeff_alpha1, coeff_alpha2, qnn, qnp, qpn, qpp, m_alpha

! phase space limits
  REAL(DP) :: ye_min, ye_max, ye_step, ye_incr
  REAL(DP) :: n_min, n_max, n_step, n_incr
  REAL(DP) :: t_min, t_max, t_step, t_incr
  REAL(DP) :: ye_test, n_test, t_test
  REAL(DP) :: x1guess(1), x3guess(3)
  INTEGER(I4B) ::  iye_ini, iye_fin, in_ini, in_fin, it_ini, it_fin

! parameters of EOS
  REAL(DP) :: dens, temp, ye, delta_n, delta_t
  REAL(DP) :: Tcrit, TC0, TC1, TC2, TC3
  REAL(DP) :: exph, expa, facq
  REAL(DP) :: STI, SXI, LOGNNI, LOGNPI, LOGNNO, LOGNPO
  LOGICAL  :: error, test, make_table, guess
  REAL(DP) :: TELE(0:41)
!$omp threadprivate (TELE, dens, temp, ye, delta_n, delta_t)

  REAL(DP) :: ftable(1:200,1:130)
  LOGICAL ::  lgttable(1:200,1:130)

  LOGICAL :: FIX_A0, ASCII
  REAL(DP):: A00

  INTEGER(I4B), DIMENSION(1:3,1000) :: fix_bc_of_NaN
  INTEGER(I4B) :: NaN_count

!  character(80) :: isolist_filename != "src/nse_eos/tables/list_235"
!  character(80) :: isolist_filename = "tables/isolist_bigger"
!  character(80) :: part_filename != "src/nse_eos/tables/partition_table_gr1d"
!  character(80) :: mass_filename != "src/nse_eos/tables/mass_table"
!  character(80) :: ext_name


end module globals
