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
MODULE MERGE_TABLE_MOD

  USE Kind_Types_Mod, ONLY : DP,I4B
  USE Physical_Constants_Mod
  USE Table_Sizes_Mod
  USE Phase_Space_Input_Mod
  USE Transition_Input_Mod
  USE SNA_MOD
  USE NSE_MOD
  USE ELE_MOD
  USE MERGE_ELE_MOD
  USE MERGE_NSE_MOD
  USE MERGE_SNA_MOD
  USE WRAP
  USE Tables_Input_Mod 

  IMPLICIT NONE

CONTAINS

  SUBROUTINE MERGE_EOS_TABLE

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

    REAL(DP) :: a_in, a_n, dadn, d2adn2
    REAL(DP) :: sna_free, nse_free, free
    REAL(DP) :: sna_denerdd, sna_denerdy, sna_denerdt
    REAL(DP) :: nse_denerdd, nse_denerdy, nse_denerdt

    WRITE (*,*) 'Merge SNA EOS and NSE EOS.'

    keytemp = 1

!$OMP PARALLEL DO SCHEDULE(dynamic,1) DEFAULT (none) &
!$OMP FIRSTPRIVATE(Yp_ini,Yp_fin,Yp_step,Yp_min) &
!$OMP FIRSTPRIVATE(n_ini,n_fin,steps_per_decade_in_n,Log10n_min) &
!$OMP FIRSTPRIVATE(T_ini,T_fin,steps_per_decade_in_T,Log10T_min) &
!$OMP FIRSTPRIVATE(Log10nt_min,Log10nt_max,energy_shift) &
!$OMP FIRSTPRIVATE(n_transition,n_delta) &
!$OMP SHARED(nse_merge_ener,nse_merge_pres,nse_merge_entr) &
!$OMP SHARED(nse_merge_mu_n,nse_merge_mu_p,nse_merge_mu_e) &
!$OMP SHARED(nse_merge_xn,nse_merge_xp,nse_merge_xa) &
!$OMP SHARED(nse_merge_xl,nse_merge_xh) &
!$OMP SHARED(nse_merge_abar,nse_merge_zbar) &
!$OMP SHARED(nse_merge_albar,nse_merge_zlbar) &
!$OMP SHARED(nse_merge_r,nse_merge_u) &
!$OMP SHARED(nse_merge_meffn,nse_merge_meffp) &
!$OMP SHARED(nse_merge_dsdn,nse_merge_dsdt,nse_merge_dsdy) &
!$OMP SHARED(nse_merge_dpdn,nse_merge_dpdt,nse_merge_dpdy) &
!$OMP SHARED(sna_merge_ener,sna_merge_pres,sna_merge_entr) &
!$OMP SHARED(sna_merge_mu_n,sna_merge_mu_p,sna_merge_mu_e) &
!$OMP SHARED(sna_merge_xn,sna_merge_xp,sna_merge_xa) &
!$OMP SHARED(sna_merge_xl,sna_merge_xh) &
!$OMP SHARED(sna_merge_abar,sna_merge_zbar) &
!$OMP SHARED(sna_merge_albar,sna_merge_zlbar) &
!$OMP SHARED(sna_merge_r,sna_merge_u) &
!$OMP SHARED(sna_merge_meffn,sna_merge_meffp) &
!$OMP SHARED(sna_merge_dsdn,sna_merge_dsdt,sna_merge_dsdy) &
!$OMP SHARED(sna_merge_dpdn,sna_merge_dpdt,sna_merge_dpdy) &
!$OMP PRIVATE(xe,i_t,temp,temp_cgs,i_n,dens,dens_cgs) &
!$OMP PRIVATE(nse_ener,nse_pres,nse_entr,nse_free) &
!$OMP PRIVATE(nse_mu_hat,nse_mu_n,nse_mu_p,nse_mu_e) &
!$OMP PRIVATE(nse_dpresdd,nse_dpresdt,nse_dpresdy) &
!$OMP PRIVATE(nse_dentrdd,nse_dentrdt,nse_dentrdy) &
!$OMP PRIVATE(nse_denerdd,nse_denerdt,nse_denerdy) &
!$OMP PRIVATE(nse_dmuhdd,nse_dmuhdt,nse_dmuhdy) &
!$OMP PRIVATE(nse_xn,nse_xp,nse_xa,nse_xh,nse_xl) &
!$OMP PRIVATE(nse_abar,nse_zbar,nse_albar,nse_zlbar) &
!$OMP PRIVATE(nse_r,nse_u,nse_meffn,nse_meffp) &
!$OMP PRIVATE(sna_ener,sna_pres,sna_entr,sna_free) &
!$OMP PRIVATE(sna_mu_hat,sna_mu_n,sna_mu_p,sna_mu_e) &
!$OMP PRIVATE(sna_dpresdd,sna_dpresdt,sna_dpresdy) &
!$OMP PRIVATE(sna_dentrdd,sna_dentrdt,sna_dentrdy) &
!$OMP PRIVATE(sna_denerdd,sna_denerdt,sna_denerdy) &
!$OMP PRIVATE(sna_dmuhdd,sna_dmuhdt,sna_dmuhdy) &
!$OMP PRIVATE(sna_xn,sna_xp,sna_xa,sna_xh,sna_xl) &
!$OMP PRIVATE(sna_abar,sna_zbar,sna_albar,sna_zlbar) &
!$OMP PRIVATE(sna_r,sna_u,sna_meffn,sna_meffp) &
!$OMP PRIVATE(a_in,a_n,dadn,d2adn2,abar,zbar,albar,zlbar) &
!$OMP PRIVATE(pres,entr,ener,mu_p,mu_n,mu_hat,mu_e,mu_nu) &
!$OMP PRIVATE(rad,u,meffn,meffp,xa,xp,xn,xh,xl) &
!$OMP PRIVATE(dsdn,dsdt,dsdy,dpdn,dpdt,dpdy,dedn,dedt,dedy) &
!$OMP PRIVATE(free,gam,fac,cs2,dpdrhoe,dpderho) &
!$OMP SHARED(final_tab)
    DO i_yp = yp_ini, yp_fin
      xe = Yp_min + dble(i_Yp-1)*Yp_step
      WRITE (6,"(A18,F7.5)") 'MERGE loop. Yp = ', xe
      DO i_t = t_ini, t_fin
        temp = TEN**(Log10T_min+dble(i_t-1)/dble(steps_per_decade_in_T))
        temp_cgs = temp*temp_mev_to_kelvin
        DO i_n = n_ini, n_fin
          dens = TEN**(Log10n_min+dble(i_n-1)/dble(steps_per_decade_in_n))
          dens_cgs = dens/rho_cgs_to_EoS
          IF (log10(dens)>=Log10nt_max .OR. log10(dens)<=Log10nt_min) CYCLE

          ! obtain contribution from nse_table

          nse_pres = nse_merge_pres(i_n,i_t,i_yp)
          nse_entr = nse_merge_entr(i_n,i_t,i_yp)
          nse_ener = nse_merge_ener(i_n,i_t,i_yp)

          nse_mu_p = nse_merge_mu_p(i_n,i_t,i_yp)
          nse_mu_n = nse_merge_mu_n(i_n,i_t,i_yp)
          nse_mu_e = nse_merge_mu_e(i_n,i_t,i_yp)
          nse_mu_hat = nse_mu_n - nse_mu_p

          nse_xa   = nse_merge_xa(i_n,i_t,i_yp)
          nse_xp   = nse_merge_xp(i_n,i_t,i_yp)
          nse_xn   = nse_merge_xn(i_n,i_t,i_yp)
          nse_xh   = nse_merge_xh(i_n,i_t,i_yp)
          nse_abar = nse_merge_abar(i_n,i_t,i_yp)
          nse_zbar = nse_merge_zbar(i_n,i_t,i_yp)
          nse_xl   = nse_merge_xl(i_n,i_t,i_yp)
          nse_albar = nse_merge_albar(i_n,i_t,i_yp)
          nse_zlbar = nse_merge_zlbar(i_n,i_t,i_yp)

          nse_r    = nse_merge_r(i_n,i_t,i_yp)
          nse_u    = nse_merge_u(i_n,i_t,i_yp)

          nse_meffn = nse_merge_meffn(i_n,i_t,i_yp)
          nse_meffp = nse_merge_meffp(i_n,i_t,i_yp)

          nse_dentrdd = nse_merge_dsdn(i_n,i_t,i_yp)
          nse_dentrdt = nse_merge_dsdt(i_n,i_t,i_yp)
          nse_dentrdy = nse_merge_dsdy(i_n,i_t,i_yp)

          nse_dpresdd = nse_merge_dpdn(i_n,i_t,i_yp)
          nse_dpresdt = nse_merge_dpdt(i_n,i_t,i_yp)
          nse_dpresdy = nse_merge_dpdy(i_n,i_t,i_yp)

          ! obtain contribution from sna_table

          sna_pres = sna_merge_pres(i_n,i_t,i_yp)
          sna_entr = sna_merge_entr(i_n,i_t,i_yp)
          sna_ener = sna_merge_ener(i_n,i_t,i_yp)

          sna_mu_p = sna_merge_mu_p(i_n,i_t,i_yp)
          sna_mu_n = sna_merge_mu_n(i_n,i_t,i_yp)
          sna_mu_e = sna_merge_mu_e(i_n,i_t,i_yp)
          sna_mu_hat = sna_mu_n - sna_mu_p

          sna_xa   = sna_merge_xa(i_n,i_t,i_yp)
          sna_xp   = sna_merge_xp(i_n,i_t,i_yp)
          sna_xn   = sna_merge_xn(i_n,i_t,i_yp)
          sna_xh   = sna_merge_xh(i_n,i_t,i_yp)
          sna_abar = sna_merge_abar(i_n,i_t,i_yp)
          sna_zbar = sna_merge_zbar(i_n,i_t,i_yp)
          sna_xl   = sna_merge_xl(i_n,i_t,i_yp)
          sna_albar = sna_merge_albar(i_n,i_t,i_yp)
          sna_zlbar = sna_merge_zlbar(i_n,i_t,i_yp)

          sna_r    = sna_merge_r(i_n,i_t,i_yp)
          sna_u    = sna_merge_u(i_n,i_t,i_yp)

          sna_meffn = sna_merge_meffn(i_n,i_t,i_yp)
          sna_meffp = sna_merge_meffp(i_n,i_t,i_yp)

          sna_dentrdd = sna_merge_dsdn(i_n,i_t,i_yp)
          sna_dentrdt = sna_merge_dsdt(i_n,i_t,i_yp)
          sna_dentrdy = sna_merge_dsdy(i_n,i_t,i_yp)

          sna_dpresdd = sna_merge_dpdn(i_n,i_t,i_yp)
          sna_dpresdt = sna_merge_dpdt(i_n,i_t,i_yp)
          sna_dpresdy = sna_merge_dpdy(i_n,i_t,i_yp)

          ! combine tables
          ! if nt_min < log10(dens) < nt_max combine SNA and NSE EoSs
          ! Assume a(n) = 0.5(1+tanh((log10(n)-n_trans)/n_delta))
          ! Obtain second derivatives da/dn and d²a/dn²
          a_in  = (log10(dens)-n_transition)/n_delta
          a_n   = half*(one+tanh(a_in))
          dadn  = half/cosh(a_in)**two/n_delta/log(ten)/dens
          d2adn2= -dadn/dens*(one+two*tanh(a_in)/log(ten)/n_delta)

          ! combine internal energy
          ener = a_n*sna_ener + (one-a_n)*nse_ener
          ! combine entropy
          entr = a_n*sna_entr + (one-a_n)*nse_entr
          ! combine electron chemical potential
          mu_e = a_n*sna_mu_e + (one-a_n)*nse_mu_e
          ! combine free energy f_T = a(n)*f_LS + (1-a(n))*f_NSE
          sna_free = sna_ener - temp*sna_entr
          nse_free = nse_ener - temp*nse_entr
          ! combine neutron chemical potential
          mu_n = a_n*sna_mu_n + (one-a_n)*nse_mu_n &
                 + dens*dadn*(sna_free-nse_free)
          ! combine proton chemical potential
          mu_p = a_n*sna_mu_p + (one-a_n)*nse_mu_p &
                 + dens*dadn*(sna_free-nse_free)
          ! combine (neutron - proton) chemical potential
          mu_hat = mu_n - mu_p
          free = dens*(a_n*sna_free + (one-a_n)*nse_free)
          ! combine pressure
          pres = a_n*sna_pres+(one-a_n)*nse_pres &
                + dens*dens*dadn*(sna_free-nse_free)
          ! combine pressure derivatives
          dpdn = a_n*sna_dpresdd + (one-a_n)*nse_dpresdd&
                    + (sna_pres-nse_pres)*dadn &
                    + two*dens*(sna_free-nse_free)*dadn &
                    + dens*dens*(sna_free-nse_free)*d2adn2
          dpdy = a_n*sna_dpresdy + (one-a_n)*nse_dpresdy &
                - dens*dens*(sna_mu_hat-nse_mu_hat)*dadn
          dpdt = a_n*sna_dpresdt + (one-a_n)*nse_dpresdt &
                - dens*dens*(sna_entr-nse_entr)*dadn
          ! combine internal energy derivatives
          sna_denerdd = sna_pres/dens/dens + temp*sna_dentrdd
          nse_denerdd = nse_pres/dens/dens + temp*nse_dentrdd
          dedn = a_n*sna_denerdd + (one-a_n)*nse_denerdd &
                + (sna_ener-nse_ener)*dadn

          sna_denerdy = - sna_mu_hat + temp*sna_dentrdy
          nse_denerdy = - nse_mu_hat + temp*nse_dentrdy
          dedy        = a_n*sna_denerdy + (one-a_n)*nse_denerdy

          sna_denerdt = temp*sna_dentrdt
          nse_denerdt = temp*nse_dentrdt
          dedt        = a_n*sna_denerdt + (one-a_n)*nse_denerdt
          ! combine entropy derivatives
          dsdn = a_n*sna_dentrdd + (one-a_n)*nse_dentrdd&
                + (sna_entr-nse_entr)*dadn
          dsdy = a_n*sna_dentrdy + (one-a_n)*nse_dentrdy
          dsdt = a_n*sna_dentrdt + (one-a_n)*nse_dentrdt
          ! combined average nuclear mass number, etc...
          abar  = (a_n*sna_xh + (one-a_n)*nse_xh) &
                / (a_n*sna_xh/sna_abar + (one-a_n)*nse_xh/nse_abar )
          zbar  = (a_n*sna_xh + (one-a_n)*nse_xh) &
                / (a_n*sna_xh/sna_zbar + (one-a_n)*nse_xh/nse_zbar )

          IF (sna_xh == zero .AND. nse_xh == zero) THEN
            abar = 1.d0
            zbar = xe
          ENDIF

          albar = nse_albar 
          zlbar = nse_zlbar

          xa   = a_n*sna_xa + (one-a_n)*nse_xa
          xn   = a_n*sna_xn + (one-a_n)*nse_xn
          xp   = a_n*sna_xp + (one-a_n)*nse_xp
          xh   = a_n*sna_xh + (one-a_n)*nse_xh
          xl   = a_n*sna_xl + (one-a_n)*nse_xl

          rad  = a_n*sna_r
          u    = a_n*sna_u

          meffn  = a_n*sna_meffn + (one-a_n)*nse_meffn
          meffp  = a_n*sna_meffp + (one-a_n)*nse_meffp

!         get lepton+photon part of EoS
          IF (isnan(abar).OR.isnan(zbar).OR.isnan(albar).OR.isnan(zlbar)) THEN
            WRITE (*,"(A25,15ES14.6)") 'MERGE error at (y,T,n):', xe, temp, dens, &
                    abar, sna_xh, sna_abar, nse_xh, nse_abar
            CYCLE
          ENDIF

          IF (abs(1.d0-xn-xp-xa-xl-xh)>1.d-6) THEN
            WRITE (*,"(A25,15ES14.6)") 'Possible abundance error at (y,T,n):', xe, temp, dens
          ENDIF

          IF (xh>1.d-10 .AND. abar < 1.d0) THEN
            WRITE (*,"(A25,15ES14.6)") 'Possible heavy nuclei error at (y,T,n):', xe, temp, dens
          ENDIF

          IF (xl>1.d-10 .AND. albar < 1.d0) THEN
            WRITE (*,"(A25,15ES14.6)") 'Possible light nuclei error at (y,T,n):', xe, temp, dens
          ENDIF

          mu_nu = mu_e - mu_n + mu_p

          dedn = pres/dens/dens + temp*dsdn
          dedt = temp*dsdt
          dedy = - mu_hat + temp*dsdy

!         Adiabatic index and speed of sound in units of speed of light
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
          final_tab(i_n,i_t,i_yp,25:26) = (/meffn,meffp/) 
 
!         save final table nucleons, light and heavy nuclei properties
          final_tab(i_n,i_t,i_yp,16:24) = &
                                        (/xn,xp,xa,xh,abar,zbar,xl,albar,zlbar/)

!         save final table P,s,e and chemical potentials
          pres = log10(pres) - log10(press_cgs_to_EOS)
          ener = log10(ener+energy_shift) + log10(energy_EOS_to_cgs)
          final_tab(i_n,i_t,i_yp,1:8) = &
                                  (/pres,ener,entr,mu_hat,mu_n,mu_p,mu_e,mu_nu/)

        ENDDO
      ENDDO
    ENDDO
!$END OMP PARALLEL DO

    WRITE (*,*) 'End merge of NSE and SNA EOSs.'

    if (allocated(yp))        deallocate(yp)
    if (allocated(logrho))    deallocate(logrho)
    if (allocated(logtemp))   deallocate(logtemp)
    if (allocated(alltables)) deallocate(alltables)

    allocate(yp(yp_fin))
    allocate(logrho(n_fin))
    allocate(logtemp(t_fin))

    write (*,*) 'write yp array'

    do i_yp = yp_ini, yp_fin
      yp(i_yp) = Yp_min + dble(i_Yp-1)*Yp_step
    enddo
    nyp = yp_fin-yp_ini+1

    write (*,*) 'write log10n array'

    do i_n  = n_ini, n_fin
      logrho(i_n)  = Log10n_min+dble(i_n-1)/dble(steps_per_decade_in_n) &
       - log10(rho_cgs_to_EOS)
    enddo
    nrho = n_fin-n_ini+1

    write (*,*) 'write log10temp array'

    do i_t = t_ini, t_fin
      logtemp(i_t) = Log10T_min+dble(i_t-1)/dble(steps_per_decade_in_T)
    enddo

    ntemp = t_fin-t_ini+1

    if (.not. only_NSE .OR. only_SNA) then
      deallocate(sna_merge_pres)
      deallocate(sna_merge_entr)
      deallocate(sna_merge_ener)
      deallocate(sna_merge_mu_p)
      deallocate(sna_merge_mu_n)
      deallocate(sna_merge_xa  )
      deallocate(sna_merge_xp  )
      deallocate(sna_merge_xn  )
      deallocate(sna_merge_xh  )
      deallocate(sna_merge_abar)
      deallocate(sna_merge_zbar)
      deallocate(sna_merge_xl  )
      deallocate(sna_merge_albar)
      deallocate(sna_merge_zlbar)
      deallocate(sna_merge_r   )
      deallocate(sna_merge_u   )
      deallocate(sna_merge_meffn)
      deallocate(sna_merge_meffp)
      deallocate(sna_merge_dsdn)
      deallocate(sna_merge_dsdt)
      deallocate(sna_merge_dsdy)
      deallocate(sna_merge_dpdn)
      deallocate(sna_merge_dpdt)
      deallocate(sna_merge_dpdy)
      deallocate(sna_merge_dedn)
      deallocate(sna_merge_dedt)
      deallocate(sna_merge_dedy)
    endif

    if (.not. only_SNA .OR. only_NSE) then
      deallocate(nse_merge_pres)
      deallocate(nse_merge_entr)
      deallocate(nse_merge_ener)
      deallocate(nse_merge_mu_p)
      deallocate(nse_merge_mu_n)
      deallocate(nse_merge_xa  )
      deallocate(nse_merge_xp  )
      deallocate(nse_merge_xn  )
      deallocate(nse_merge_xh  )
      deallocate(nse_merge_abar)
      deallocate(nse_merge_zbar)
      deallocate(nse_merge_xl  )
      deallocate(nse_merge_albar)
      deallocate(nse_merge_zlbar)
      deallocate(nse_merge_r   )
      deallocate(nse_merge_u   )
      deallocate(nse_merge_meffn)
      deallocate(nse_merge_meffp)
      deallocate(nse_merge_dsdn)
      deallocate(nse_merge_dsdt)
      deallocate(nse_merge_dsdy)
      deallocate(nse_merge_dpdn)
      deallocate(nse_merge_dpdt)
      deallocate(nse_merge_dpdy)
      deallocate(nse_merge_dedn)
      deallocate(nse_merge_dedt)
      deallocate(nse_merge_dedy)
    endif

  END SUBROUTINE MERGE_EOS_TABLE

END MODULE MERGE_TABLE_MOD
