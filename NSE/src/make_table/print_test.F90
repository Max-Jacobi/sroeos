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
  USE Physical_Constants_Mod, ONLY : neutron_mass_MeV

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Print_Test( n, t, y, x_n, x_p, x_alpha, x_light, x_heavy, &
                   A_heavy, Z_heavy, A_light, Z_light, P, s, e, &
                   mu_n, mu_p, muh, dmuhdn, dmuhdt, dmuhdy, &
                   dpdn, dpdt, dpdy, dsdn, dsdt, dsdy, dedn, dedt, dedy )


    IMPLICIT NONE

    REAL(DP), INTENT(IN) :: n, t, y, x_n, x_p, x_alpha, x_light, x_heavy, &
                     A_heavy, Z_heavy, A_light, Z_light, P, s, e, &
                     mu_n, mu_p, muh, dmuhdn, dmuhdt, dmuhdy, &
                     dpdn, dpdt, dpdy, dsdn, dsdt, dsdy, dedn, dedt, dedy

    REAL(DP) :: F, DF_DN, DF_DT, DF_DY
    REAL(DP) ::    DS_DN, DS_DT, DS_DY
    REAL(DP) ::    DP_DN, DP_DT, DP_DY

    F = n*(e-T*s)
    DF_DN = F/n + n*(dedn - T*dsdn)
    DF_DT =       n*(dedt - T*dsdt - s)
    DF_DY =       n*(dedy - T*dsdy)

    DS_DN = n*dsdn + s
    DS_DT = n*dsdt
    DS_DY = n*dsdy

    DP_DN = dpdn
    DP_DT = dpdt
    DP_DY = dpdy

    WRITE (*,*)
    WRITE (*,603) 'n  =', n, ' fm^-3'
    WRITE (*,603) 'T  =', T, ' MeV  '
    WRITE (*,603) 'Yp = ',  y, '     '
    WRITE (*,*)
    WRITE (*,602) 'Neutron = ', x_n
    WRITE (*,602) 'Protons = ', x_p
    WRITE (*,602) 'Helium  = ', x_alpha
    WRITE (*,602) 'Light   = ', x_light
    WRITE (*,602) 'Heavy   = ', x_heavy
    WRITE (*,602) 'A_light = ', A_light
    WRITE (*,602) 'Z_light = ', Z_light
    WRITE (*,602) 'N_light = ', A_light - Z_light
    WRITE (*,602) 'A_heavy = ', A_heavy
    WRITE (*,602) 'Z_heavy = ', Z_heavy
    WRITE (*,602) 'N_heavy = ', A_heavy - Z_heavy
    WRITE (*,*)
    WRITE (*,602) ' e        =', e
    WRITE (*,602) ' f        =', e - T*s
    WRITE (*,602) ' s        =', s
    WRITE (*,602) ' P        =', P
    WRITE (*,602) ' df/dn    =', (DF_DN-F/n)/n
    WRITE (*,602) ' df/dT    =', DF_DT/n
    WRITE (*,602) ' df/dY    =', DF_DY/n
    WRITE (*,602) ' mu_n     =', mu_n - neutron_mass_MeV
    WRITE (*,602) ' mu_p     =', mu_p - neutron_mass_MeV
    WRITE (*,602) ' mu_bar   =', muh
    WRITE (*,*)
    WRITE (*,801)
    WRITE (*,800) 'dP', DP_DN, DP_DY, DP_DT
    WRITE (*,800) 'dS', DS_DN, DS_DY, DS_DT
    WRITE (*,800) 'dF', DF_DN, DF_DY, DF_DT
    WRITE (*,800) 'dmu',dmuhdn, dmuhdy, dmuhdt
    WRITE (*,*)
    WRITE (*,801)
    WRITE (*,800) 'dP', dpdn, dpdy, dpdt
    WRITE (*,800) 'ds', dsdn, dsdy, dsdt
    WRITE (*,800) 'de', dedn, dedy, dedt
    WRITE (*,800) 'df', DF_DN/n-F/n**2.d0, DF_DY/n, DF_DT/n
    WRITE (*,800) 'dmu',dmuhdn, dmuhdy, dmuhdt

602 format (A14,ES22.14)
603 format (A14,ES22.14,A6)
799 format (A12,1es20.12,4x,A12,1es20.12)
800 format (A5,5es20.12)
801 format (13x,'/dn',17x,'/dy',17x,'/dT')

  END SUBROUTINE PRINT_TEST

END MODULE Print_Test_Mod
