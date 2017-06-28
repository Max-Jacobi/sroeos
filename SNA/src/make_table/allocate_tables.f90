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
Module Allocate_Mod

  USE Kind_Types_Mod,        ONLY : I4B, DP
  USE Physical_Constants_Mod, ONLY : ZERO
  USE Global_Variables_Mod
  USE Phase_Space_Input_Mod
  USE Allocatable_Tables_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE ALLOCATE_OUTPUT

    IMPLICIT NONE

    INTEGER(I4B) :: i

!   set array dimension along log10(rho), log10(T) and ye
    n_log10n =  n_fin -  n_ini + 1
    n_log10T =  T_fin -  T_ini + 1
    n_Yp     = Yp_fin - Yp_ini + 1

!   allcate log10(rho), log10(T) and ye table dimensions
    allocate(log10n_tab(n_log10n))
    allocate(log10T_tab(n_log10T))
    allocate(yp_tab(n_Yp))

!   determine points where density will be calculated
    do i = n_fin, n_ini, -1 ! log10(n) for rho in 1/fm^3
      log10n_tab(i) = (Log10n_min+dble(i-1)/dble(steps_per_decade_in_n))
    enddo

!   determine points where temperature will be calculated
    do i = t_fin, t_ini, -1 ! log10(T) for T in MeV
      log10t_tab(i) = (Log10T_min+dble(i-1)/dble(steps_per_decade_in_T))
    enddo

!   determine points where electron fraction will be calculated
    do i = yp_ini, yp_fin ! Yp
      yp_tab(i) = yp_min + yp_step*dble(i-1)
    enddo

!   allocate dimensions of saved quantities
!    x = (P,s,E) and y = (n,Yp,T)
!    allocate dx/dy, mu_n, mu_p, mu_hat, A, Z, xn, xp, xa, xh, u, r
    allocate(muh_tab(n_log10n,n_log10T,n_Yp))
    allocate(mun_tab(n_log10n,n_log10T,n_Yp))
    allocate(mup_tab(n_log10n,n_log10T,n_Yp))

    allocate(meff_n_tab(n_log10n,n_log10T,n_Yp))
    allocate(meff_p_tab(n_log10n,n_log10T,n_Yp))

    allocate(p_tab(n_log10n,n_log10T,n_Yp))
    allocate(e_tab(n_log10n,n_log10T,n_Yp))
    allocate(s_tab(n_log10n,n_log10T,n_Yp))

    allocate(dpdn_tab(n_log10n,n_log10T,n_Yp))
    allocate(dsdn_tab(n_log10n,n_log10T,n_Yp))
    allocate(dmudn_tab(n_log10n,n_log10T,n_Yp))

    allocate(dpdt_tab(n_log10n,n_log10T,n_Yp))
    allocate(dsdt_tab(n_log10n,n_log10T,n_Yp))
    allocate(dmudt_tab(n_log10n,n_log10T,n_Yp))

    allocate(dpdy_tab(n_log10n,n_log10T,n_Yp))
    allocate(dsdy_tab(n_log10n,n_log10T,n_Yp))
    allocate(dmudy_tab(n_log10n,n_log10T,n_Yp))

    allocate(xn_tab(n_log10n,n_log10T,n_Yp))
    allocate(xp_tab(n_log10n,n_log10T,n_Yp))
    allocate(xa_tab(n_log10n,n_log10T,n_Yp))
    allocate(xh_tab(n_log10n,n_log10T,n_Yp))

    allocate(abar_tab(n_log10n,n_log10T,n_Yp))
    allocate(zbar_tab(n_log10n,n_log10T,n_Yp))

    allocate(u_tab(n_log10n,n_log10T,n_Yp))
    allocate(r_tab(n_log10n,n_log10T,n_Yp))

    p_tab = zero
    e_tab = zero
    s_tab = zero

    dpdn_tab = zero
    dsdn_tab = zero
    dmudn_tab = zero

    dpdt_tab = zero
    dsdt_tab = zero
    dmudt_tab = zero

    dpdy_tab = zero
    dsdy_tab = zero
    dmudy_tab = zero

    muh_tab  = zero
    mun_tab  = zero
    mup_tab  = zero

    meff_n_tab = zero
    meff_p_tab = zero

    xn_tab   = zero
    xp_tab   = zero
    xa_tab   = zero
    xh_tab   = zero

    abar_tab = zero
    zbar_tab = zero

    u_tab = zero
    r_tab = zero

    pointsrho  = n_log10n
    pointstemp = n_log10T
    pointsye   = n_Yp

  END SUBROUTINE ALLOCATE_OUTPUT

END MODULE Allocate_Mod
