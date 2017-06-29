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
MODULE NSE_TABLE_MOD

  USE Kind_Types_Mod, ONLY : DP,I4B
  USE Physical_Constants_Mod
  USE Table_Sizes_Mod
  USE Phase_Space_Input_Mod
  USE Transition_Input_Mod
  USE READ_EOS
  USE NSE_MOD
  USE ELE_MOD
  USE MERGE_NSE_MOD
  USE WRAP

  IMPLICIT NONE

CONTAINS

  SUBROUTINE NSE_EOS_TABLE

    INTEGER(I4B) :: i,i_n,i_t,i_yp, keytemp, keyerr
    REAL(DP) :: xe, dens, temp, dens_cgs, temp_cgs
    REAL(DP) :: rfeps

    REAL(DP) :: den,tem,ener,pres,entr
    REAL(DP) :: abar,zbar,albar,zlbar,xa,xn,xp,xh,xl,rad,u
    REAL(DP) :: denerdt,dpresdt,dentrdt,denerdd,dpresdd,dentrdd
    REAL(DP) :: denerdy,dpresdy,dentrdy
    REAL(DP) :: mu_p,mu_n,mu_hat,mu_e,mu_nu
    REAL(DP) :: dsdn,dsdt,dsdy,dpdn,dpdt,dpdy,dedn,dedt,dedy
    REAL(DP) :: dpdrhoe,dpderho,gam,fac,cs2

    REAL(DP) :: gam1, etaele, sound

    WRITE (*,*) 'Compute contribution of NSE EOS to final EOS.'

! check if arrays have been ALLOCATEd
    IF (ALLOCATED(yp))        DEALLOCATE(yp)
    IF (ALLOCATED(logtemp))   DEALLOCATE(logtemp)
    IF (ALLOCATED(logrho))    DEALLOCATE(logrho)
    IF (ALLOCATED(alltables)) DEALLOCATE(alltables)

    nyp   = nse_nyp
    nrho  = nse_nrho
    ntemp = nse_ntemp

    ALLOCATE(yp(nyp))
    ALLOCATE(logrho(nrho))
    ALLOCATE(logtemp(ntemp))
    ALLOCATE(alltables(nrho,ntemp,nyp,nvars))

    yp        = nse_yp
    logrho    = nse_logrho
    logtemp   = nse_logtemp
    alltables = nse_table

    eos_rhomin = nse_eos_rhomin
    eos_rhomax = nse_eos_rhomax
    eos_tempmin= nse_eos_tempmin
    eos_tempmax= nse_eos_tempmax
    eos_ypmin  = nse_eos_ypmin
    eos_ypmax  = nse_eos_ypmax

  ! ALLOCATE merge arrays

    ALLOCATE(nse_merge_pres(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_entr(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_ener(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_mu_p(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_mu_n(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_mu_e(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_xa  (n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_xp  (n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_xn  (n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_xh  (n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_abar(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_zbar(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_xl  (n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_r   (n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_u   (n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_albar(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_zlbar(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_dsdn(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_dsdt(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_dsdy(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_dpdn(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_dpdt(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_dpdy(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_dedn(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_dedt(n_fin,t_fin,yp_fin))
    ALLOCATE(nse_merge_dedy(n_fin,t_fin,yp_fin))

    nse_merge_pres = zero
    nse_merge_entr = zero
    nse_merge_ener = zero
    nse_merge_mu_p = zero
    nse_merge_mu_n = zero
    nse_merge_xa   = zero
    nse_merge_xp   = zero
    nse_merge_xn   = zero
    nse_merge_xh   = zero
    nse_merge_abar = zero
    nse_merge_zbar = zero
    nse_merge_xl   = zero
    nse_merge_albar = zero
    nse_merge_zlbar = zero
    nse_merge_r    = zero
    nse_merge_u    = zero
    nse_merge_dsdn = zero
    nse_merge_dsdt = zero
    nse_merge_dsdy = zero
    nse_merge_dpdn = zero
    nse_merge_dpdt = zero
    nse_merge_dpdy = zero
    nse_merge_dedn = zero
    nse_merge_dedt = zero
    nse_merge_dedy = zero

    keytemp = 1

!$OMP PARALLEL DO SCHEDULE(dynamic,1) &
!$OMP PRIVATE(xe,i_t,temp,temp_cgs,i_n,dens,dens_cgs) &
!$OMP PRIVATE(nse_ener,nse_pres,nse_entr,nse_mu_hat,nse_mu_n,nse_mu_p) &
!$OMP PRIVATE(nse_dpresdd,nse_dpresdt,nse_dpresdy) &
!$OMP PRIVATE(nse_dentrdd,nse_dentrdt,nse_dentrdy) &
!$OMP PRIVATE(nse_dmuhdd,nse_dmuhdt,nse_dmuhdy) &
!$OMP PRIVATE(nse_xn,nse_xp,nse_xa,nse_xh,nse_xl) &
!$OMP PRIVATE(nse_abar,nse_zbar,nse_albar,nse_zlbar) &
!$OMP PRIVATE(nse_r,nse_u) &
!$OMP FIRSTPRIVATE(keytemp,keyerr,rfeps) &
!$OMP PRIVATE(abar,zbar,albar,zlbar) &
!$OMP PRIVATE(pres,entr,ener,mu_p,mu_n,mu_hat,mu_e,mu_nu,rad,u) &
!$OMP PRIVATE(xa,xp,xn,xh,xl) &
!$OMP PRIVATE(dsdn,dsdt,dsdy,dpdn,dpdt,dpdy,dedn,dedt,dedy) &
!$OMP PRIVATE(gam,fac,cs2,dpdrhoe,dpderho) &
!$OMP PRIVATE(ele_ener, ele_pres, ele_entr) &
!$OMP PRIVATE(ele_denerdt, ele_dpresdt, ele_dentrdt) &
!$OMP PRIVATE(ele_denerdd, ele_dpresdd, ele_dentrdd) &
!$OMP PRIVATE(ele_denerdy, ele_dpresdy, ele_dentrdy) &
!$OMP PRIVATE(gam1, etaele, sound) &
!$OMP SHARED(final_tab) &
!$OMP SHARED(nse_merge_pres,nse_merge_ener,nse_merge_entr) &
!$OMP SHARED(nse_merge_mu_p,nse_merge_mu_n,nse_merge_mu_e) &
!$OMP SHARED(nse_merge_xa,nse_merge_xn,nse_merge_xp,nse_merge_xh,nse_merge_xl) &
!$OMP SHARED(nse_merge_zbar,nse_merge_abar,nse_merge_zlbar,nse_merge_albar) &
!$OMP SHARED(nse_merge_dsdn,nse_merge_dsdt,nse_merge_dsdy) &
!$OMP SHARED(nse_merge_dpdn,nse_merge_dpdt,nse_merge_dpdy) &
!$OMP SHARED(nse_merge_r,nse_merge_u)
    DO i_yp = yp_ini, yp_fin
      xe = Yp_min + dble(i_Yp-1)*Yp_step
      WRITE (*,*) 'NSE loop. Yp = ', xe
      DO i_t = t_ini, t_fin
        temp = TEN**(Log10T_min+dble(i_t-1)/dble(steps_per_decade_in_T))
        temp_cgs = temp*temp_mev_to_kelvin
        DO i_n = n_ini, n_fin
          dens = TEN**(Log10n_min+dble(i_n-1)/dble(steps_per_decade_in_n))
          dens_cgs = dens/rho_cgs_to_EoS
          IF (log10(dens)>Log10nt_max) CYCLE

          CALL get_eos_full(dens,temp,xe, &
            nse_ener,nse_pres,nse_entr,       &
            nse_mu_hat,nse_mu_n,nse_mu_p,     &
            nse_dpresdd,nse_dpresdt,nse_dpresdy, &
            nse_dentrdd,nse_dentrdt,nse_dentrdy, &
            nse_dmuhdd,nse_dmuhdt,nse_dmuhdy,    &
            nse_xn,nse_xp,nse_xa,nse_xh,nse_xl,  &
            nse_abar,nse_zbar,nse_albar,nse_zlbar, &
            nse_r,nse_u,keytemp,keyerr,rfeps)

            abar = nse_xn+nse_xp+nse_xa/four
            if (nse_abar  > 0.d0) &
              abar = nse_xn+nse_xp+nse_xa/four+nse_xh/nse_abar
            if (nse_albar > 0.d0) &
              abar = nse_xn+nse_xp+nse_xa/four+nse_xl/nse_albar
            if (nse_abar > 0.d0 .and. nse_albar > 0.d0) &
            abar = nse_xn+nse_xp+nse_xa/four+nse_xh/nse_abar+nse_xl/nse_albar

            abar = one/abar
            zbar = xe*abar

!         get lepton+photon part of EoS
          IF (isnan(abar).OR.isnan(zbar).OR.abar<=1.d-10.OR.zbar<=1.d-10) THEN
            WRITE (*,"(A25,15ES14.6)") 'NSE error at (y,T,n):', xe, temp, dens
           WRITE (*,"(A26,15ES14.6)") '     A, Z, Abar, x_h',&
                   nse_xn, nse_xp, nse_xa, nse_xh, nse_xl, & 
                   nse_abar, nse_zbar, nse_albar, nse_zlbar
!           write (*,"(10ES15.7)") '     E, P, S, mu_n, mu_p', & 
!                                  nse_ener,nse_pres,nse_entr,nse_mu_n,nse_mu_p
            CYCLE
          ENDIF

          CALL wrap_timmes(abar,zbar,dens_cgs,temp_cgs, &
                            ele_ener,ele_pres,ele_entr, &
                            ele_denerdt,ele_dpresdt,ele_dentrdt, &
                            ele_denerdd,ele_dpresdd,ele_dentrdd, &
                            ele_denerdy,ele_dpresdy,ele_dentrdy, &
                            gam1,etaele,sound)

! calculate nuclear + electrons + radiation thermo quantities for final table
!         p,s,e
          pres  = nse_pres  + ele_pres*press_cgs_to_EOS
          entr  = nse_entr  + ele_entr*entropy_cgs_to_EOS
          ener  = nse_ener  + ele_ener*energy_cgs_to_EOS

!         mu_p,mu_n,mu_hat,mu_e,mu_nu
          mu_p  = nse_mu_p
          mu_n  = nse_mu_n
          mu_hat= mu_n - mu_p
          mu_e  = temp*etaele + Mass_e
          mu_nu = mu_e - mu_n + mu_p

!         nuclei
          xa   = nse_xa
          xp   = nse_xp
          xn   = nse_xn
          xh   = nse_xh
          abar = nse_abar
          zbar = nse_zbar
          xl   = nse_xl
          albar = nse_albar
          zlbar = nse_zlbar
          rad   = nse_r
          u     = nse_u

!         derivatives
!         ds/dn, ds/dt, ds/dy
          dsdn = nse_dentrdd + ele_dentrdd*entropy_cgs_to_EOS/rho_cgs_to_EOS
          dsdt = nse_dentrdt + ele_dentrdt*entropy_cgs_to_EOS*temp_mev_to_kelvin
          dsdy = nse_dentrdy + ele_dentrdy*entropy_cgs_to_EOS
!         dp/dn, dp/dt, dp/dy
          dpdn = nse_dpresdd + ele_dpresdd*press_cgs_to_EOS/rho_cgs_to_EOS
          dpdt = nse_dpresdt + ele_dpresdt*press_cgs_to_EOS*temp_mev_to_kelvin
          dpdy = nse_dpresdy + ele_dpresdy*press_cgs_to_EOS
!         de/dn, de/dt, de/dy
          dedn = pres/dens/dens + temp*dsdn
          dedt = temp*dsdt
          dedy = - mu_hat + temp*dsdy

!       save fraction of array to be merged later
          IF (log10(dens)<Log10nt_max) THEN
            nse_merge_pres(i_n,i_t,i_yp)  = pres
            nse_merge_entr(i_n,i_t,i_yp)  = entr
            nse_merge_ener(i_n,i_t,i_yp)  = ener

            nse_merge_mu_p(i_n,i_t,i_yp)  = mu_p
            nse_merge_mu_n(i_n,i_t,i_yp)  = mu_n
            nse_merge_mu_e(i_n,i_t,i_yp)  = mu_e

            nse_merge_xa(i_n,i_t,i_yp)    = xa
            nse_merge_xp(i_n,i_t,i_yp)    = xp
            nse_merge_xn(i_n,i_t,i_yp)    = xn
            nse_merge_xh(i_n,i_t,i_yp)    = xh
            nse_merge_abar(i_n,i_t,i_yp)  = abar
            nse_merge_zbar(i_n,i_t,i_yp)  = zbar
            nse_merge_xl(i_n,i_t,i_yp)    = xl
            nse_merge_albar(i_n,i_t,i_yp)  = albar
            nse_merge_zlbar(i_n,i_t,i_yp)  = zlbar

            nse_merge_r(i_n,i_t,i_yp)    = rad
            nse_merge_u(i_n,i_t,i_yp)    = u

            nse_merge_dsdn(i_n,i_t,i_yp)  = dsdn
            nse_merge_dsdt(i_n,i_t,i_yp)  = dsdt
            nse_merge_dsdy(i_n,i_t,i_yp)  = dsdy

            nse_merge_dpdn(i_n,i_t,i_yp)  = dpdn
            nse_merge_dpdt(i_n,i_t,i_yp)  = dpdt
            nse_merge_dpdy(i_n,i_t,i_yp)  = dpdy

!            nse_merge_dedn(i_n,i_t,i_yp)  = dedn
!            nse_merge_dedt(i_n,i_t,i_yp)  = dedt
!            nse_merge_dedy(i_n,i_t,i_yp)  = dedy
        ENDIF

!       Adiabatic index and speed of sound in units of speed of light
          gam = dens/pres*dpdn + temp/dens*(dpdt*dpdt)/pres/dedt
          fac = one + &
          (ener*energy_EOS_to_cgs+clite*clite)*dens_cgs/(pres/press_cgs_to_EOS)
          cs2 = gam/fac
          cs2 = cs2*clite*clite

!         save final table derivatives, adiabatic index and sound speed
          dpdrhoe = (dpdn - dpdt*dedn/dedt)*press_EOS_to_cgs*rho_cgs_to_EoS
          dpderho = (dpdt/dedt)*(press_EOS_to_cgs/energy_EOS_to_cgs)
          dedt    = dedt*energy_EOS_to_cgs

          final_tab(i_n,i_t,i_yp,9:13) = (/dedt,dpdrhoe,dpderho,gam,cs2/)

          final_tab(i_n,i_t,i_yp,14:15) = (/rad,u/)

!         save final table nucleons, light and heavy nuclei properties
          final_tab(i_n,i_t,i_yp,16:24) = &
                                        (/xn,xp,xa,xh,abar,zbar,xl,albar,zlbar/)

!         save final table P,s,e and chemical potentials
          pres = log10(pres) - log10(press_cgs_to_EOS)
          ener = log10(ener+energy_shift)+log10(energy_EOS_to_cgs)

          final_tab(i_n,i_t,i_yp,1:8) = &
                                  (/pres,ener,entr,mu_hat,mu_n,mu_p,mu_e,mu_nu/)

        ENDDO
      ENDDO
    ENDDO
!$END OMP PARALLEL DO

    DEALLOCATE(yp,logrho,logtemp,alltables)

    WRITE (*,*) 'End computation of contribution of NSE EOS to final EOS.'

  END SUBROUTINE NSE_EOS_TABLE

END MODULE NSE_TABLE_MOD
