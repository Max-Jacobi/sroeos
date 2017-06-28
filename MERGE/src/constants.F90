module constants_module

  use Kind_Types_Mod

  implicit none

!---------------------- constants

!from NIST, also defined in lseos.2.7.f for EL_EOS                           
   REAL(DP),parameter :: neutron_mass_MeV  = 939.5654d0 !  neutron mass
   REAL(DP),parameter :: proton_mass_MeV   = 938.2721d0 !   proton mass
   REAL(DP),parameter :: electron_mass_MeV = 0.510999d0 ! electron mass

   REAL(DP),parameter :: mec2 = electron_mass_MeV

   REAL(DP),parameter :: HBC = 197.327d0 
   REAL(DP),parameter :: alpha_coul = 1.d0/137.035999074d0
   REAL(DP),parameter :: clite = 29979245900.0d0

   REAL(DP),parameter :: kb_erg = 1.3806504d-16
   REAL(DP),parameter :: kb_mev = 8.61738568d-11
   REAL(DP),parameter :: erg_to_mev = 6.24150636d5
   REAL(DP),parameter :: mev_to_cgs = 1.60217649d-6
   REAL(DP),parameter :: temp_mev_to_kelvin = 1.1604505d10   
   REAL(DP),parameter :: press_EOS_to_cgs = 1.60217649d33
   REAL(DP),parameter :: energy_EOS_to_cgs = 1.60217649d-6/1.674927211d-24
   REAL(DP),parameter :: entropy_EOS_to_cgs = 1.3806504d-16/1.674927211d-24
   REAL(DP),parameter :: rho_cgs_to_EOS = 5.9704087d-16 !neutron mass 
   REAL(DP),parameter :: press_cgs_to_EOS  = 6.24150964d-34
   REAL(DP),parameter :: energy_cgs_to_EOS = 1.674927211d-24/1.60217649d-6
   REAL(DP),parameter :: entropy_cgs_to_EOS = 1.674927211d-24/1.3806504d-16 
   REAL(DP),parameter :: planck_cgs = 6.626176d-27
   REAL(DP),parameter :: pi1 = 3.14159265358979d0
   REAL(DP),parameter :: avo = 5.9704082443567622d23 ! 1/m_n_cgs
   REAL(DP),parameter :: sac_const = 2.0d0 * pi1 * (hbc)**2
!   REAL(DP),parameter :: avo = 6.0221367d23

   REAL(DP),parameter :: m_n_cgs = neutron_mass_MeV*mev_to_cgs/clite**2
   REAL(DP),parameter :: m_p_cgs = proton_mass_MeV*mev_to_cgs/clite**2

   ! the following constants are used in Luke's mixing script
   !  real*8,parameter :: e_nuc_to_cgs = 1.0d0/1.03654d-18 amu based
   !  m_neutron based
  REAL(DP),parameter :: e_nuc_to_cgs = 1.0d0/1.04541174787614d-18  
  REAL(DP),parameter :: p_nuc_to_cgs = 1.0d0/6.24151d-34


end module constants_module
