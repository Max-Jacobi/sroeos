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
MODULE SNA_TABLE_MOD

  USE Kind_Types_Mod, ONLY : DP,I4B
  USE Physical_Constants_Mod
  USE Table_Sizes_Mod
  USE Phase_Space_Input_Mod
  USE Transition_Input_Mod
  USE READ_EOS
  USE SNA_MOD
  USE ELE_MOD
  USE MERGE_SNA_MOD
  USE WRAP

  IMPLICIT NONE

CONTAINS

  SUBROUTINE SNA_EOS_TABLE

    INTEGER(I4B) :: i,i_n,i_t,i_yp, keytemp, keyerr
    REAL(DP) :: xe, dens, temp, dens_cgs, temp_cgs
    REAL(DP) :: rfeps

    REAL(DP) :: den,tem,ener,pres,entr
    REAL(DP) :: abar,zbar,albar,zlbar,xa,xn,xp,xh,xl
    REAL(DP) :: rad,u,meffn,meffp
    REAL(DP) :: denerdt,dpresdt,dentrdt,denerdd,dpresdd,dentrdd
    REAL(DP) :: denerdy,dpresdy,dentrdy
    REAL(DP) :: mu_p,mu_n,mu_hat,mu_e,mu_nu
    REAL(DP) :: dsdn,dsdt,dsdy,dpdn,dpdt,dpdy,dedn,dedt,dedy
    REAL(DP) :: dpdrhoe,dpderho,gam,fac,cs2

    REAL(DP) :: gam1, etaele, sound

    WRITE (*,*) 'Compute contribution of SNA EOS to final EOS.'

! check if arrays have been ALLOCATEd
    IF (ALLOCATED(yp))        DEALLOCATE(yp)
    IF (ALLOCATED(logtemp))   DEALLOCATE(logtemp)
    IF (ALLOCATED(logrho))    DEALLOCATE(logrho)
    IF (ALLOCATED(alltables)) DEALLOCATE(alltables)

    nyp   = sna_nyp
    nrho  = sna_nrho
    ntemp = sna_ntemp

    ALLOCATE(yp(nyp))
    ALLOCATE(logrho(nrho))
    ALLOCATE(logtemp(ntemp))
    ALLOCATE(alltables(nrho,ntemp,nyp,nvars))

    yp        = sna_yp
    logrho    = sna_logrho
    logtemp   = sna_logtemp
    alltables = sna_table

    eos_rhomin = sna_eos_rhomin
    eos_rhomax = sna_eos_rhomax
    eos_tempmin= sna_eos_tempmin
    eos_tempmax= sna_eos_tempmax
    eos_ypmin  = sna_eos_ypmin
    eos_ypmax  = sna_eos_ypmax

  ! ALLOCATE merge arrays

    ALLOCATE(sna_merge_pres(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_entr(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_ener(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_mu_p(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_mu_n(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_mu_e(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_xa  (n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_xp  (n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_xn  (n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_xh  (n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_abar(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_zbar(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_xl  (n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_albar(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_zlbar(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_r   (n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_u   (n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_meffn(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_meffp(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_dsdn(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_dsdt(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_dsdy(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_dpdn(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_dpdt(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_dpdy(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_dedn(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_dedt(n_fin,t_fin,yp_fin))
    ALLOCATE(sna_merge_dedy(n_fin,t_fin,yp_fin))

    sna_merge_pres = zero
    sna_merge_entr = zero
    sna_merge_ener = zero
    sna_merge_mu_p = zero
    sna_merge_mu_n = zero
    sna_merge_xa   = zero
    sna_merge_xp   = zero
    sna_merge_xn   = zero
    sna_merge_xh   = zero
    sna_merge_abar = zero
    sna_merge_zbar = zero
    sna_merge_albar = zero
    sna_merge_zlbar = zero
    sna_merge_r    = zero
    sna_merge_u    = zero
    sna_merge_meffn = zero
    sna_merge_meffp = zero
    sna_merge_dsdn = zero
    sna_merge_dsdt = zero
    sna_merge_dsdy = zero
    sna_merge_dpdn = zero
    sna_merge_dpdt = zero
    sna_merge_dpdy = zero
    sna_merge_dedn = zero
    sna_merge_dedt = zero
    sna_merge_dedy = zero

    keytemp = 1

!$OMP PARALLEL DO SCHEDULE(dynamic,1) &
!$OMP PRIVATE(xe,i_t,temp,temp_cgs,i_n,dens,dens_cgs) &
!$OMP PRIVATE(sna_ener,sna_pres,sna_entr,sna_mu_hat,sna_mu_n,sna_mu_p) &
!$OMP PRIVATE(sna_dpresdd,sna_dpresdt,sna_dpresdy) &
!$OMP PRIVATE(sna_dentrdd,sna_dentrdt,sna_dentrdy) &
!$OMP PRIVATE(sna_dmuhdd,sna_dmuhdt,sna_dmuhdy) &
!$OMP PRIVATE(sna_xn,sna_xp,sna_xa,sna_xh,sna_xl) &
!$OMP PRIVATE(sna_abar,sna_zbar,sna_albar,sna_zlbar) &
!$OMP PRIVATE(sna_r,sna_u,sna_meffn,sna_meffp) &
!$OMP FIRSTPRIVATE(keytemp,keyerr,rfeps) &
!$OMP PRIVATE(abar,zbar,albar,zlbar) &
!$OMP PRIVATE(pres,entr,ener,mu_p,mu_n,mu_hat,mu_e,mu_nu) &
!$OMP PRIVATE(rad,u,meffn,meffp) &
!$OMP PRIVATE(xa,xp,xn,xh,xl) &
!$OMP PRIVATE(dsdn,dsdt,dsdy,dpdn,dpdt,dpdy,dedn,dedt,dedy) &
!$OMP PRIVATE(gam,fac,cs2,dpdrhoe,dpderho) &
!$OMP PRIVATE(ele_ener, ele_pres, ele_entr) &
!$OMP PRIVATE(ele_denerdt, ele_dpresdt, ele_dentrdt) &
!$OMP PRIVATE(ele_denerdd, ele_dpresdd, ele_dentrdd) &
!$OMP PRIVATE(ele_denerdy, ele_dpresdy, ele_dentrdy) &
!$OMP PRIVATE(gam1, etaele, sound) &
!$OMP SHARED(final_tab) &
!$OMP SHARED(sna_merge_pres,sna_merge_ener,sna_merge_entr) &
!$OMP SHARED(sna_merge_mu_p,sna_merge_mu_n,sna_merge_mu_e) &
!$OMP SHARED(sna_merge_xa,sna_merge_xn,sna_merge_xp,sna_merge_xh,sna_merge_xl) &
!$OMP SHARED(sna_merge_zbar,sna_merge_abar,sna_merge_zlbar,sna_merge_albar) &
!$OMP SHARED(sna_merge_dsdn,sna_merge_dsdt,sna_merge_dsdy) &
!$OMP SHARED(sna_merge_dpdn,sna_merge_dpdt,sna_merge_dpdy) &
!$OMP SHARED(sna_merge_r,sna_merge_u,sna_merge_meffn,sna_merge_meffp)
    DO i_yp = yp_ini, yp_fin
      xe = Yp_min + dble(i_Yp-1)*Yp_step
      WRITE (*,*) 'SNA loop. Yp = ', xe
      DO i_t = t_ini, t_fin
        temp = TEN**(Log10T_min+dble(i_t-1)/dble(steps_per_decade_in_T))
        temp_cgs = temp*temp_mev_to_kelvin
        DO i_n = n_ini, n_fin
          dens = TEN**(Log10n_min+dble(i_n-1)/dble(steps_per_decade_in_n))
          dens_cgs = dens/rho_cgs_to_EoS
          IF (log10(dens)<Log10nt_min) CYCLE

          CALL get_eos_full(dens,temp,xe, &
            sna_ener,sna_pres,sna_entr,       &
            sna_mu_hat,sna_mu_n,sna_mu_p,     &
            sna_dpresdd,sna_dpresdt,sna_dpresdy, &
            sna_dentrdd,sna_dentrdt,sna_dentrdy, &
            sna_dmuhdd,sna_dmuhdt,sna_dmuhdy,    &
            sna_xn,sna_xp,sna_xa,sna_xh,sna_xl,  &
            sna_abar,sna_zbar,sna_albar,sna_zlbar, &
            sna_r,sna_u,sna_meffn,sna_meffp,keytemp,keyerr,rfeps)

            abar = sna_xn+sna_xp+sna_xa/four
            if (sna_abar  > 0.d0) &
              abar = sna_xn+sna_xp+sna_xa/four+sna_xh/sna_abar
            if (sna_albar > 0.d0) &
              abar = sna_xn+sna_xp+sna_xa/four+sna_xl/sna_albar
            if (sna_abar > 0.d0 .and. sna_albar > 0.d0) &
            abar = sna_xn+sna_xp+sna_xa/four+sna_xh/sna_abar+sna_xl/sna_albar

            abar = one/abar
            zbar = xe*abar

!         get lepton+photon part of EoS
          IF (isnan(abar).OR.isnan(zbar).OR.abar<=1.d-10.OR.zbar<=1.d-10) THEN
            WRITE (*,"(A26,15ES14.6)") 'Error at (y,T,n) : ', xe, temp, dens
!          WRITE (*,"(A26,15ES14.6)") '     A, Z, Abar, x_h',&
!                                     abar, zbar, sna_abar, sna_xh
!          write (*,"(10ES15.7)") '     E, P, S, mu_n, mu_p', & !                                  sna_ener,sna_pres,sna_entr,sna_mu_n,sna_mu_p
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
          pres  = sna_pres  + ele_pres*press_cgs_to_EOS
          entr  = sna_entr  + ele_entr*entropy_cgs_to_EOS
          ener  = sna_ener  + ele_ener*energy_cgs_to_EOS

!         mu_p,mu_n,mu_hat,mu_e,mu_nu
          mu_p  = sna_mu_p
          mu_n  = sna_mu_n
          mu_hat= mu_n - mu_p
          mu_e  = temp*etaele + Mass_e
          mu_nu = mu_e - mu_n + mu_p

!         nuclei
          xa   = sna_xa
          xp   = sna_xp
          xn   = sna_xn
          xh   = sna_xh
          abar = sna_abar
          zbar = sna_zbar
          xl   = sna_xl
          albar = sna_albar
          zlbar = sna_zlbar
          rad   = sna_r
          u     = sna_u
          meffn = sna_meffn
          meffp = sna_meffp

!         derivatives
!         ds/dn, ds/dt, ds/dy
          dsdn = sna_dentrdd + ele_dentrdd*entropy_cgs_to_EOS/rho_cgs_to_EOS
          dsdt = sna_dentrdt + ele_dentrdt*entropy_cgs_to_EOS*temp_mev_to_kelvin
          dsdy = sna_dentrdy + ele_dentrdy*entropy_cgs_to_EOS
!         dp/dn, dp/dt, dp/dy
          dpdn = sna_dpresdd + ele_dpresdd*press_cgs_to_EOS/rho_cgs_to_EOS
          dpdt = sna_dpresdt + ele_dpresdt*press_cgs_to_EOS*temp_mev_to_kelvin
          dpdy = sna_dpresdy + ele_dpresdy*press_cgs_to_EOS
!         de/dn, de/dt, de/dy
          dedn = pres/dens/dens + temp*dsdn
          dedt = temp*dsdt
          dedy = - mu_hat + temp*dsdy

!       save fraction of array to be merged later
          IF (log10(dens)<Log10nt_max) THEN
            sna_merge_pres(i_n,i_t,i_yp)  = pres
            sna_merge_entr(i_n,i_t,i_yp)  = entr
            sna_merge_ener(i_n,i_t,i_yp)  = ener

            sna_merge_mu_p(i_n,i_t,i_yp)  = mu_p
            sna_merge_mu_n(i_n,i_t,i_yp)  = mu_n
            sna_merge_mu_e(i_n,i_t,i_yp)  = mu_e

            sna_merge_xa(i_n,i_t,i_yp)    = xa
            sna_merge_xp(i_n,i_t,i_yp)    = xp
            sna_merge_xn(i_n,i_t,i_yp)    = xn
            sna_merge_xh(i_n,i_t,i_yp)    = xh
            sna_merge_abar(i_n,i_t,i_yp)  = abar
            sna_merge_zbar(i_n,i_t,i_yp)  = zbar
            sna_merge_xl(i_n,i_t,i_yp)    = zero
            sna_merge_albar(i_n,i_t,i_yp)  = zero
            sna_merge_zlbar(i_n,i_t,i_yp)  = zero

            sna_merge_r(i_n,i_t,i_yp)    = rad
            sna_merge_u(i_n,i_t,i_yp)    = u

            sna_merge_meffn(i_n,i_t,i_yp)    = meffn
            sna_merge_meffp(i_n,i_t,i_yp)    = meffp

            sna_merge_dsdn(i_n,i_t,i_yp)  = dsdn
            sna_merge_dsdt(i_n,i_t,i_yp)  = dsdt
            sna_merge_dsdy(i_n,i_t,i_yp)  = dsdy

            sna_merge_dpdn(i_n,i_t,i_yp)  = dpdn
            sna_merge_dpdt(i_n,i_t,i_yp)  = dpdt
            sna_merge_dpdy(i_n,i_t,i_yp)  = dpdy

!            sna_merge_dedn(i_n,i_t,i_yp)  = dedn
!            sna_merge_dedt(i_n,i_t,i_yp)  = dedt
!            sna_merge_dedy(i_n,i_t,i_yp)  = dedy
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
          dedt = dedt*energy_EOS_to_cgs

          final_tab(i_n,i_t,i_yp,9:13) = (/dedt,dpdrhoe,dpderho,gam,cs2/)

!         save final table nucleons, light and heavy nuclei properties
          final_tab(i_n,i_t,i_yp,16:24) = &
                                        (/xn,xp,xa,xh,abar,zbar,xl,albar,zlbar/)

          final_tab(i_n,i_t,i_yp,14:15) = (/rad,u/)
          final_tab(i_n,i_t,i_yp,25:26) = (/meffn,meffp/)

!         save final table P,s,e and chemical potentials
          pres = log10(pres) - log10(press_cgs_to_EOS)
          ener = log10(ener+energy_shift) + log10(energy_EOS_to_cgs)

          final_tab(i_n,i_t,i_yp,1:8) = &
                                  (/pres,ener,entr,mu_hat,mu_n,mu_p,mu_e,mu_nu/)

        ENDDO
      ENDDO
    ENDDO
!$END OMP PARALLEL DO

    DEALLOCATE(yp,logrho,logtemp,alltables)

    WRITE (*,*) 'End computation of contribution of SNA EOS to final EOS.'

  END SUBROUTINE SNA_EOS_TABLE

END MODULE SNA_TABLE_MOD
