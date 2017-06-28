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
MODULE Print_test_Mod

  USE Kind_Types_Mod, ONLY : DP, I4B, LGCL
  USE Physical_Constants_Mod, ONLY : TEN, HALF
  USE Free_Energy_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Print_Test( n, T, Yp, F_o, F_i, F_alpha, F_TR, F_SC, &
                   n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
                   A_heavy, Z_heavy, rad, F, P, S, E, &
                   DF_Dn, DF_Dy, DF_DT, DP_Dn, DP_Dy, DP_DT, &
                   DS_Dn, DS_Dy, DS_DT, DE_Dn, DE_Dy, DE_DT, &
                   mu_no, mu_po, Meff_no, Meff_po, &
                   Dmu_nout_DT, Dmu_nout_Dn, Dmu_nout_Dy, &
                   Dmu_pout_DT, Dmu_pout_Dn, Dmu_pout_Dy)


    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: n, T, Yp, F_o, F_i, F_alpha, F_TR, F_SC, &
                     n_no, n_po, n_ni, n_pi, n_alpha, n_heavy, &
                     A_heavy, Z_heavy, rad, F, P, S, E, &
                     DF_Dn, DF_Dy, DF_DT, DP_Dn, DP_Dy, DP_DT, &
                     DS_Dn, DS_Dy, DS_DT, DE_Dn, DE_Dy, DE_DT, &
                     mu_no, mu_po, Meff_no, Meff_po, &
                     Dmu_nout_DT, Dmu_nout_Dn, Dmu_nout_Dy, &
                     Dmu_pout_DT, Dmu_pout_Dn, Dmu_pout_Dy

    REAL(DP) :: Dmuh_Dn, Dmuh_Dy, Dmuh_DT
!
    Dmuh_Dn = Dmu_nout_Dn - Dmu_pout_Dn
    Dmuh_Dy = Dmu_nout_Dy - Dmu_pout_Dy
    Dmuh_DT = Dmu_nout_DT - Dmu_pout_DT

    WRITE (*,*)
    WRITE (*,603) 'n  =', n, ' fm^-3'
    WRITE (*,603) 'T  =', T, ' MeV  '
    WRITE (*,603) 'Yp = ',  Yp, '     '
    WRITE (*,*)
    WRITE (*,602) 'Helium  = ', n_alpha
    WRITE (*,602) 'Protons = ', n_po
    WRITE (*,602) 'Neutron = ', n_no
    WRITE (*,602) 'Heavy   = ', n_heavy
    WRITE (*,602) 'H_A     = ', A_heavy
    WRITE (*,602) 'H_Z     = ', Z_heavy
    WRITE (*,602) 'H_N     = ', A_heavy - Z_heavy
    WRITE (*,*)
    WRITE (*,602) ' e        =', E/n
    WRITE (*,602) ' f        =', F/n
    WRITE (*,602) ' s        =', S/n
    WRITE (*,602) ' P        =', P
    WRITE (*,602) ' df/dn    =', (DF_DN-F/n)/n
    WRITE (*,602) ' df/dT    =', DF_DT/n
    WRITE (*,602) ' df/dY    =', DF_DY/n
    ! WRITE (*,602) ' d2f/dn2  =', (D2F_DN2-TWO*DF_DN/n+TWO*F_TOTAL/n/n)/n
    ! WRITE (*,602) ' d2f/dT2  =', D2F_DT2/n
    ! WRITE (*,602) ' d2f/dY2  =', D2F_DY2/n
    ! WRITE (*,602) ' d2f/dndT =', (D2F_DTDN-DF_DT/n)/n
    ! WRITE (*,602) ' d2f/dTdY =', D2F_DYDT/n
    ! WRITE (*,602) ' d2f/dYdn =', (D2F_DNDY-DF_DY/n)/n
    ! WRITE (*,602) ' Gamma_s  =', GAMMA_IN
    WRITE (*,602) ' mu_n     =', mu_no
    WRITE (*,602) ' mu_p     =', mu_po
    WRITE (*,602) ' mu_bar   =', mu_no - mu_po
    WRITE (*,*)
    WRITE (*,801)
    WRITE (*,800) 'dP', DP_DN, DP_DY, DP_DT
    WRITE (*,800) 'dS', DS_DN, DS_DY, DS_DT
    WRITE (*,800) 'dF', DF_DN, DF_DY, DF_DT
    WRITE (*,800) 'dmuh',Dmuh_Dn, Dmuh_Dy, Dmuh_DT
    WRITE (*,*)
    WRITE (*,801)
    WRITE (*,800) 'dP', DP_DN, DP_DY, DP_DT
    WRITE (*,800) 'ds', DS_DN/n-S/n**TWO, DS_DY/n, DS_DT/n
    WRITE (*,800) 'de', DE_DN/n-E/n**TWO, DE_DY/n, DE_DT/n
    WRITE (*,800) 'df', DF_DN/n-F/n**TWO, DF_DY/n, DF_DT/n
    WRITE (*,800) 'dmuh',Dmuh_Dn, Dmuh_Dy, Dmuh_DT

    ! WRITE (*,800)
    ! WRITE (*,799) 'n_ni = ', DNNI, 'mu_ni = ', MU_NIN
    ! WRITE (*,799) 'n_pi = ', DNPI, 'mu_pi = ', MU_PIN
    ! WRITE (*,799) 'n_no = ', DNNO, 'mu_no = ', MU_NOUT
    ! WRITE (*,799) 'n_po = ', DNPO, 'mu_po = ', MU_POUT
    ! WRITE (*,800)
    ! WRITE (*,799) 'tau_ni = ', TAU_NIN,  '(m/m*)_ni = ', MNIN_STAR/MNC2
    ! WRITE (*,799) 'tau_pi = ', TAU_PIN,  '(m/m*)_pi = ', MPIN_STAR/MPC2
    ! WRITE (*,799) 'tau_no = ', TAU_NOUT, '(m/m*)_no = ', MNOUT_STAR/MNC2
    ! WRITE (*,799) 'tau_no = ', TAU_POUT, '(m/m*)_po = ', MPOUT_STAR/MPC2
    ! WRITE (*,800)
    ! WRITE (*,799) 'eta_ni = ', ETA_NIN,  'V_ni = ', V_NIN
    ! WRITE (*,799) 'eta_pi = ', ETA_PIN,  'V_pi = ', V_PIN
    ! WRITE (*,799) 'eta_no = ', ETA_NOUT, 'V_no = ', V_NOUT
    ! WRITE (*,799) 'eta_no = ', ETA_POUT, 'V_po = ', V_POUT
    ! WRITE (*,800)
    ! WRITE (*,"(5F14.8)") (DU*F_IN+F_SC)/n, RAD, DNI, DU*F_IN/n, F_SC/n
    ! WRITE (*,*)

602 format (A14,ES22.14)
603 format (A14,ES22.14,A6)
799 format (A12,1es20.12,4x,A12,1es20.12)
800 format (A5,5es20.12)
801 format (13x,'/dn',17x,'/dy',17x,'/dT')

  END SUBROUTINE PRINT_TEST

END MODULE Print_Test_Mod
