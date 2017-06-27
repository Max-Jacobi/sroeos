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

!
!  This file is originally part of the Timmes EOS available at
!  http://cococubed.asu.edu/code_pages/eos.shtml
!  Reference: F.X. Timmes & D. Arnett, ApJS 125 277 (1999)
!
!  It has been slightly modified by Andre da Silva Schneider
!   to fit SRO_EOS source code
!

module functions

use fdfunctions

implicit none

contains

      subroutine xneroot(mode,den,temp,abar,zbar,ionized, &
                        aa,f,df, &
                        etaele,detadd,detadt,detada,detadz, &
                        etapos,zeff, &
                        xne,dxnedd,dxnedt,dxneda,dxnedz, &
                        xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz, &
                        xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz, &
                        pele,dpeledd,dpeledt,dpeleda,dpeledz, &
                        ppos,dpposdd,dpposdt,dpposda,dpposdz, &
                        pep,dpepdd,dpepdt,dpepda,dpepdz, &
                        eele,deeledd,deeledt,deeleda,deeledz, &
                        epos,deposdd,deposdt,deposda,deposdz, &
                        eep,deepdd,deepdt,deepda,deepdz, &
                        sele,dseledd,dseledt,dseleda,dseledz, &
                        spos,dsposdd,dsposdt,dsposda,dsposdz, &
                        sep,dsepdd,dsepdt,dsepda,dsepdz, &
                        potmult, &
                        eip,deipdd,deipdt,deipda,deipdz, &
                        sip,dsipdd,dsipdt,dsipda,dsipdz)

      include 'implno.dek'
      include 'const.dek'


!c..this routine is called by a root finder to find the degeneracy parameter aa
!c..where the number density from a saha equation equals the number
!c..density as computed by the fermi-dirac integrals.

!c..input is the mode (0 = root find mode, 1 = full calculation),
!c..temperature temp, density dem, average weight abar, average charge zbar,
!c..and the degeneracy parameter (chemical potential/kerg*temp) aa.

!c..everything else is output. see the calling routine for the definitions
!c..of the variables



!c..declare the pass
      integer          mode,ionized,potmult
      double precision den,temp,abar,zbar, &
                       aa,f,df, &
                       etaele,detadd,detadt,detada,detadz, &
                       etapos,zeff, &
                       xne,dxnedd,dxnedt,dxneda,dxnedz, &
                       xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz, &
                       xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz, &
                       pele,dpeledd,dpeledt,dpeleda,dpeledz, &
                       ppos,dpposdd,dpposdt,dpposda,dpposdz, &
                       pep,dpepdd,dpepdt,dpepda,dpepdz, &
                       eele,deeledd,deeledt,deeleda,deeledz, &
                       epos,deposdd,deposdt,deposda,deposdz, &
                       eep,deepdd,deepdt,deepda,deepdz, &
                       sele,dseledd,dseledt,dseleda,dseledz, &
                       spos,dsposdd,dsposdt,dsposda,dsposdz, &
                       sep,dsepdd,dsepdt,dsepda,dsepdz, &
                       eip,deipdd,deipdt,deipda,deipdz, &
                       sip,dsipdd,dsipdt,dsipda,dsipdz




!c..local variables
      double precision kt,beta,beta12,beta32,beta52,  &
                       chi,chiff,chifac,saha,sfac,    &
                       xni,dxnidd,dxnidt,dxnida,dxnidz, &
                       f12,f12eta,f12beta, &
                       f32,f32eta,f32beta, &
                       f52,f52eta,f52beta, &
                       dzeff_deta,dzeffdd,dzeffdt,dzeffda,dzeffdz, &
                       dxne_deta,dxnefer_deta,dxnefer_dbeta, &
                       detap_deta,detap_dbeta, &
                       dxnpfer_detap,dxnpfer_deta,dxnpfer_dbeta, &
                       dxep_deta,dxep_dbeta, &
                       dpele_deta,dpele_dbeta,deele_deta,deele_dbeta, &
                       dppos_detap,dppos_deta,dppos_dbeta, &
                       depos_detap,depos_deta,depos_dbeta, &
                       dsfac_deta,ytot1,zz,y,yy,ww,denion



      double precision xconst,pconst,econst,mecc,dbetadt,safe
      parameter        (xconst  = 2.4883740912221807d30, &
                        pconst  = 1.3581730208282635d24, &
                        econst  = 2.0372595312423953d24, &
                        mecc    = me * clight * clight, &
                        dbetadt = kerg/mecc, &
                        safe    = 0.005d0)


!c..note:
!c..xconst = 8.0d0 * pi * sqrt(2.0d0) * (me/h)**3 * c**3
!c..pconst = xconst * 2.0d0/3.0d0 * me * clight**2
!c..econst = xconst * me * clight**2




!c..some common factors
      ytot1   = 1.0d0/abar
      kt      = kerg * temp
      beta    = kt/mecc
      beta12  = sqrt(beta)
      beta32  = beta * beta12
      xni     = avo * ytot1 * den
      etaele  = aa


!c..ion number density in 1/cm**3,
      xni     = avo * ytot1 * den
      dxnidd  = avo * ytot1
      dxnidt  = 0.0d0
      dxnida  = -xni * ytot1
      dxnidz  = 0.0d0





!c..get the number density of free electrons
!c..saha is the ratio of the ground state to the ionized state
!c..this model is exact for a pure hydrogen composition
!c..denion is a crude pressure ionization model

      if (ionized .eq. 0) then

       denion     = 0.1d0
       chi        = hion * ev2erg * zbar

       chifac     = chi/kt
       yy         = chifac - den/denion

       if (yy .gt. 200.0) then
        saha   = 1.0d90
        f      = 0.0d0
        df     = 1.0d0
        etaele = -100.0d0
        if (mode .eq. 0) return

       else if (yy .lt. -200.0) then
        saha = 0.0d0

       else
        ww   = min(200.0d0,chifac + etaele - den/denion)
        saha = 2.0d0 * exp(ww)
       end if


!c..assume fully ionized
      else
       denion     = 1.0e30
       chi        = 0.0d0
       chifac     = 0.0d0
       saha       = 0.0d0
      end if



!c..the saha factor, effective charge, and
!c..the number density of free electrons

      sfac       = 1.0d0/(1.0d0 + saha)
      dsfac_deta = -sfac*sfac*saha

      zeff       = zbar * sfac
      dzeff_deta = zbar * dsfac_deta

      xne        = xni * zeff
      dxne_deta  = xni * dzeff_deta





!c..get the fermi-dirac integral electron contribution
      call dfermi(0.5d0, etaele, beta, f12, f12eta, f12beta)
      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta)

      zz            = xconst * beta32
      yy            = f12 + beta * f32
      xnefer        = zz * yy
      dxnefer_deta  = zz * (f12eta + beta * f32eta)
      dxnefer_dbeta = xconst * beta12 * (1.5d0 * yy &
                             +  beta * (f12beta + f32 + beta * f32beta))



!c..if the temperature is not too low, get the positron contributions
!c..chemical equilibrium means etaele + etapos = eta_photon = 0.
      etapos        = 0.0d0
      detap_deta    = 0.0d0
      detap_dbeta   = 0.0d0
      xnpfer        = 0.0d0
      dxnpfer_detap = 0.0d0
      dxnpfer_dbeta = 0.0d0

      if (beta .gt. 0.02) then
       etapos      = -aa - 2.0d0/beta
       detap_deta  = -1.0d0
       detap_dbeta = 2.0d0/beta**2
       call dfermi(0.5d0, etapos, beta, f12, f12eta, f12beta)
       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta)
       xnpfer        = zz * (f12 + beta * f32)
       dxnpfer_detap = zz * (f12eta + beta * f32eta)
       dxnpfer_dbeta = xconst * beta12 * (1.5d0 * (f12 + beta * f32) &
                       +  beta * (f12beta + f32 + beta * f32beta))
      end if



!c..charge neutrality means ne_ionizat = ne_elect - ne_posit
      f  = xnefer - xnpfer - xne


!c..derivative of f with eta for newton-like root finders
      df = dxnefer_deta  - dxnpfer_detap * detap_deta  - dxne_deta



!c..if we are in root finder mode, return


      if (mode .eq. 0) return





!c..if we are not in root finder mode, polish off the calculation


!c..all the derivatives are in terms of eta and beta.
!c..we want to convert to temperature, density, abar and zbar derivatives.
!c..so, after the root find above on eta we have: xne = xnefer - xnpfer
!c..taking the derivative of this and solving for the unknown eta derivatives
!c..leads to these expressions:


      dxnpfer_deta  = dxnpfer_detap * detap_deta
      dxnpfer_dbeta = dxnpfer_dbeta + dxnpfer_detap * detap_dbeta
      dxep_deta     = dxnefer_deta  - dxnpfer_deta
      dxep_dbeta    = dxnefer_dbeta - dxnpfer_dbeta


      y      = 1.0d0/(dxep_deta - dxne_deta)

      detadd = (xne/den - dxne_deta/denion)  * y
      detadt = -(dxne_deta*chifac/temp + dxep_dbeta*dbetadt) * y
      detada = -xne/abar * y
      detadz = (xne/zbar + dxne_deta * chifac/zbar)*y



!c..derivatives of the effective charge
      dzeffdd = dzeff_deta * (detadd - 1.0d0/denion)
      dzeffdt = dzeff_deta * (detadt - chifac/temp)
      dzeffda = dzeff_deta * detada
      dzeffdz = sfac + dzeff_deta * (detadz + chifac/zbar)


!c..derivatives of the electron number density
      dxnedd = dxnidd * zeff + xni * dzeffdd
      dxnedt = dxnidt * zeff + xni * dzeffdt
      dxneda = dxnida * zeff + xni * dzeffda
      dxnedz = dxnidz * zeff + xni * dzeffdz



!c..derivatives of the fermi integral electron number densities
      dxneferdd = dxnefer_deta * detadd
      dxneferdt = dxnefer_deta * detadt + dxnefer_dbeta * dbetadt
      dxneferda = dxnefer_deta * detada
      dxneferdz = dxnefer_deta * detadz


!c..derivatives of the fermi integral positron number densities
      dxnpferdd = dxnpfer_deta * detadd
      dxnpferdt = dxnpfer_deta * detadt + dxnpfer_dbeta * dbetadt
      dxnpferda = dxnpfer_deta * detada
      dxnpferdz = dxnpfer_deta * detadz





!c..now get the pressure and energy
!c..for the electrons

      beta52  = beta * beta32
      yy      = pconst * beta52
      zz      = econst * beta52

      call dfermi(1.5d0, etaele, beta, f32, f32eta, f32beta)
      call dfermi(2.5d0, etaele, beta, f52, f52eta, f52beta)

      pele        = yy * (f32 + 0.5d0 * beta * f52)
      dpele_deta  = yy * (f32eta + 0.5d0 * beta * f52eta)
      dpele_dbeta = pconst * beta32 * (2.5d0 * (f32 + 0.5d0*beta* f52) &
                    + beta* (f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta))

      eele        = zz * (f32 + beta * f52)
      deele_deta  = zz * (f32eta + beta * f52eta)
      deele_dbeta = econst * beta32 * (2.5d0 * (f32 + beta * f52) &
                    + beta * (f32beta + f52 + beta * f52beta))


!c..for the positrons
      ppos        = 0.0d0
      dppos_detap = 0.0d0
      dppos_dbeta = 0.0d0
      epos        = 0.0d0
      depos_detap = 0.0d0
      depos_dbeta = 0.0d0

      if (beta .gt. 0.02) then
       call dfermi(1.5d0, etapos, beta, f32, f32eta, f32beta)
       call dfermi(2.5d0, etapos, beta, f52, f52eta, f52beta)

       ppos        = yy * (f32 + 0.5d0 * beta * f52)
       dppos_detap = yy * (f32eta + 0.5d0*beta *f52eta)
       dppos_dbeta = pconst * beta32 * (2.5d0 * (f32 + 0.5d0*beta*f52) &
                     + beta*(f32beta + 0.5d0*f52 + 0.5d0*beta*f52beta))

       epos        = zz * (f32 + beta * f52)
       depos_detap = zz * (f32eta + beta * f52eta)
       depos_dbeta = econst * beta32 * (2.5d0 * (f32 + beta * f52) &
                     + beta * (f32beta + f52 + beta * f52beta))
      end if


!c..derivatives of the electron pressure
      dpeledd = dpele_deta * detadd
      dpeledt = dpele_deta * detadt + dpele_dbeta * dbetadt
      dpeleda = dpele_deta * detada
      dpeledz = dpele_deta * detadz


!c..derivatives of the electron energy
      deeledd = deele_deta * detadd
      deeledt = deele_deta * detadt + deele_dbeta * dbetadt
      deeleda = deele_deta * detada
      deeledz = deele_deta * detadz


!c..derivatives of the positron pressure
      dppos_deta  = dppos_detap * detap_deta
      dppos_dbeta = dppos_dbeta + dppos_detap * detap_dbeta
      dpposdd     = dppos_deta * detadd
      dpposdt     = dppos_deta * detadt + dppos_dbeta * dbetadt
      dpposda     = dppos_deta * detada
      dpposdz     = dppos_deta * detadz


!c..derivatives of the positron energy
      depos_deta  = depos_detap * detap_deta
      depos_dbeta = depos_dbeta + depos_detap * detap_dbeta
      deposdd     = depos_deta * detadd
      deposdt     = depos_deta * detadt + depos_dbeta * dbetadt
      deposda     = depos_deta * detada
      deposdz     = depos_deta * detadz



!c..electron+positron pressure and its derivatives
!c..note: at high temperatures and low densities, dpepdd is very small
!c..and can go negative, so limit it to be positive definite
      pep    = pele    + ppos
      dpepdd = max(dpeledd + dpposdd, 1.0d-30)
      dpepdt = dpeledt + dpposdt
      dpepda = dpeleda + dpposda
      dpepdz = dpeledz + dpposdz


!c..electron+positron thermal energy and its derivatives
      eep    = eele    + epos
      deepdd = deeledd + deposdd
      deepdt = deeledt + deposdt
      deepda = deeleda + deposda
      deepdz = deeledz + deposdz




 114  format(1x,1p5e24.16)


!c..electron entropy in erg/gr/kelvin and its derivatives
      y       = kerg/den

      sele    = ((pele + eele)/kt - etaele*xnefer) * y

      dseledd = ((dpeledd + deeledd)/kt &
                  - detadd*xnefer)*y    &
                  - etaele*dxneferdd*y  &
                  - sele/den

      dseledt = ((dpeledt + deeledt)/kt &
                   - detadt*xnefer     &
                   - etaele*dxneferdt  &
                   - (pele + eele)/(kt*temp))*y

      dseleda = ((dpeleda + deeleda)/kt - detada*xnefer &
                   - etaele*dxneferda)*y

      dseledz = ((dpeledz + deeledz)/kt - detadz*xnefer &
                   - etaele*dxneferdz)*y



!c..positron entropy in erg/gr/kelvin and its derivatives
      spos    = ((ppos + epos)/kt - etapos*xnpfer) * y

      dsposdd = ((dpposdd + deposdd)/kt &
                 - detap_deta*detadd*xnpfer &
                 - etapos*dxnpferdd)*y - spos/den

      dsposdt = ((dpposdt + deposdt)/kt &
                 - (detap_deta*detadt + detap_dbeta*dbetadt)*xnpfer  &
                 - etapos*dxnpferdt  &
                 - (ppos + epos)/(kt*temp))*y

      dsposda = ((dpposda + deposda)/kt &
                 - detap_deta*detada*xnpfer &
                 - etapos*dxnpferda)*y

      dsposdz = ((dpposdz + deposdz)/kt - detap_deta*detadz*xnpfer &
                   - etapos*dxnpferdz)*y


!c..and their sum
      sep     = sele + spos
      dsepdd  = dseledd + dsposdd
      dsepdt  = dseledt + dsposdt
      dsepda  = dseleda + dsposda
      dsepdz  = dseledz + dsposdz



!c..adjust for the rest mass energy of the positrons
      y       = 2.0d0 * mecc
      epos    = epos    + y * xnpfer
      deposdd = deposdd + y * dxnpferdd
      deposdt = deposdt + y * dxnpferdt
      deposda = deposda + y * dxnpferda
      deposdz = deposdz + y * dxnpferdz


!c..and resum
      deepdd = deeledd + deposdd
      deepdt = deeledt + deposdt
      deepda = deeleda + deposda
      deepdz = deeledz + deposdz


!c..convert the electron-positron thermal energy in erg/cm**3 to
!c..a specific thermal energy in erg/gr

      eele    = eele/den
      deeledd = deeledd/den - eele/den
      deeledt = deeledt/den
      deeleda = deeleda/den
      deeledz = deeledz/den

      epos    = epos/den
      deposdd = deposdd/den - epos/den
      deposdt = deposdt/den
      deposda = deposda/den
      deposdz = deposdz/den


!c..and resum
      deepdd = deeledd + deposdd
      deepdt = deeledt + deposdt
      deepda = deeleda + deposda
      deepdz = deeledz + deposdz




!c..and take care of the ionization potential contributions
      if (potmult .eq. 0) then
       eip    = 0.0d0
       deipdd = 0.0d0
       deipdt = 0.0d0
       deipda = 0.0d0
       deipdz = 0.0d0
       sip    = 0.0d0
       dsipdd = 0.0d0
       dsipdt = 0.0d0
       dsipda = 0.0d0
       dsipdz = 0.0d0
      else

       eip    = chi * xne
       deipdd = chi * dxnedd
       deipdt = chi * dxnedt
       deipda = chi * dxneda
       deipdz = chi * dxnedz + hion*ev2erg*xne

!c..the ionization entropy in erg/gr/kelvin and its derivatives
       y      = kerg/den

!c       sip    = (eip/kt - etaele*xne) * y

!c       dsipdd = (deipdd/kt
!c     1            - detadd*xne)*y
!c     2            - etaele*dxnedd*y
!c     3            - sip/den

!c       dsipdt = (deipdt/kt
!c     1             - detadt*xne
!c     2             - etaele*dxnedt
!c     3             - eip/(kt*temp))*y


       sip    = eip/kt * y

       dsipdd = deipdd/kt*y - sip/den

       dsipdt = (deipdt/kt - eip/(kt*temp))*y

       dsipda = deipda/kt*y

       dsipdz = deipdz/kt*y


!c..convert the ionization energy from erg/cm**3 to  erg/gr

       eip    = eip/den
       deipdd = deipdd/den - eip/den
       deipdt = deipdt/den
       deipda = deipda/den
       deipdz = deipdz/den

      end if

      return
      end subroutine xneroot





      subroutine coulomb(den,temp,abar,zbar,&
                         pion,dpiondd,dpiondt,dpionda,dpiondz,&
                         xne,dxnedd,dxnedt,dxneda,dxnedz, plasg,&
                         pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,&
                         ecoul,decouldd,decouldt,decoulda,decouldz,&
                         scoul,dscouldd,dscouldt,dscoulda,dscouldz)
      include 'implno.dek'
      include 'const.dek'


!c..this routine implments coulomb corrections
!c..see yakovlev & shalybkov 1989, uniform background corrections

!c..input


!c..declare the pass
      double precision den,temp,abar,zbar,&
                       pion,dpiondd,dpiondt,dpionda,dpiondz,&
                       xne,dxnedd,dxnedt,dxneda,dxnedz,&
                       plasg,&
                       pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz,&
                       ecoul,decouldd,decouldt,decoulda,decouldz,&
                       scoul,dscouldd,dscouldt,dscoulda,dscouldz


!c..local variables
!c..for the uniform background coulomb correction
      double precision ytot1,kt,ktinv,&
                       s,dsdd,dsdt,dsda,dsdz,sinv,&
                       aele,daeledd,daeledt,daeleda,daeledz,aeleinv,&
                       eplasg,eplasgdd,eplasgdt,eplasgda,eplasgdz,&
                       aion,&
                       lami,inv_lami,lamidd,lamida,lamidz,&
                       plasgdd,plasgdt,plasgda,plasgdz,&
                       x,y,z


      double precision u0,rho6,a1,b1,c1,d1,e1,a2,b2,c2
      parameter        (a1 = -0.898004d0,&
                        b1 =  0.96786d0, &
                        c1 =  0.220703d0,&
                        d1 = -0.86097d0, &
                        e1 =  2.5269d0,  &
                        a2 =  0.29561d0, &
                        b2 =  1.9885d0,  &
                        c2 =  0.288675d0)


!c..various derived constants
      double precision third,forth,fiveth,esqu,forthpi
      parameter        (third   = 1.0d0/3.0d0,&
                        forth   = 4.0d0/3.0d0,&
                        fiveth  = 5.0d0/3.0d0,&
                        esqu    = qe*qe,&
                        forthpi = forth * pi)



!c..common variables
       ytot1   = 1.0d0/abar
       kt      = kerg * temp
       ktinv   = 1.0d0/kt



!c..yakovlev & shalybkov eqs 5, 9 and 10
       s        = forthpi * xne
       dsdd     = forthpi * dxnedd
       dsdt     = forthpi * dxnedt
       dsda     = forthpi * dxneda
       dsdz     = forthpi * dxnedz
       sinv     = 1.0d0/s

!c..electron-sphere radius aele
       aele     = sinv**third
       z        = -third * aele * sinv
       daeledd  = z * dsdd
       daeledt  = z * dsdt
       daeleda  = z * dsda
       daeledz  = z * dsdz
       aeleinv  = 1.0d0/aele

!c..electron coupling parameter eplasg
       eplasg   = esqu * ktinv * aeleinv
       z        = -eplasg * aeleinv
       eplasgdd = z * daeledd
       eplasgdt = z * daeledt - eplasg*ktinv*kerg
       eplasgda = z * daeleda
       eplasgdz = z * daeledz

!c..ion-sphere radius aion
       x        = zbar**third
       aion     = x * aele

!c..ion coupling parameter plasg
       z        = x*x*x*x*x
       plasg    = z * eplasg
       plasgdd  = z * eplasgdd
       plasgdt  = z * eplasgdt
       plasgda  = z * eplasgda
       plasgdz  = z * eplasgdz + fiveth*x*x * eplasg

!c       write(6,*)
!c       write(6,112) aion,aele
!c       write(6,112) plasg,plasgdd,plasgdt,plasgda,plasgdz
!c       write(6,*)



!c..yakovlev & shalybkov 1989 equations 82, 85, 86, 87
       if (plasg .ge. 1.0) then
        x        = plasg**(0.25d0)
        u0       = a1*plasg + b1*x + c1/x + d1
        ecoul    = pion/den * u0
        pcoul    = third * ecoul * den
        scoul    = -avo*ytot1*kerg * &
                    (3.0d0*b1*x - 5.0d0*c1/x &
                   + d1 * (log(plasg) - 1.0d0) - e1)

        y        = avo/abar*kt*(a1 + 0.25d0/plasg*(b1*x - c1/x))
        decouldd = y * plasgdd
        decouldt = y * plasgdt + ecoul/temp
        decoulda = y * plasgda - ecoul/abar
        decouldz = y * plasgdz

        y        = third * den
        dpcouldd = third * ecoul + y*decouldd
        dpcouldt = y * decouldt
        dpcoulda = y * decoulda
        dpcouldz = y * decouldz


        y        = -avo*kerg/(abar*plasg)*(0.75d0*b1*x +1.25d0*c1/x +d1)
        dscouldd = y * plasgdd
        dscouldt = y * plasgdt
        dscoulda = y * plasgda - scoul/abar
        dscouldz = y * plasgdz


!c..yakovlev & shalybkov 1989 equations 102, 103, 104
       else if (plasg .lt. 1.0) then
        x        = plasg*sqrt(plasg)
        y        = plasg**b2
        z        = c2 * x - third * a2 * y
        pcoul    = -pion * z
        ecoul    = 3.0d0 * pcoul/den
        scoul    = -avo/abar*kerg*(c2*x -a2*(b2-1.0d0)/b2*y)

        s        = 1.5d0*c2*x/plasg - third*a2*b2*y/plasg
        dpcouldd = -dpiondd*z - pion*s*plasgdd
        dpcouldt = -dpiondt*z - pion*s*plasgdt
        dpcoulda = -dpionda*z - pion*s*plasgda
        dpcouldz = -dpiondz*z - pion*s*plasgdz

        s        = 3.0d0/den
        decouldd = s * dpcouldd - ecoul/den
        decouldt = s * dpcouldt
        decoulda = s * dpcoulda
        decouldz = s * dpcouldz

        s        = -avo*kerg/(abar*plasg)*(1.5d0*c2*x -a2*(b2-1.0d0)*y)
        dscouldd = s * plasgdd
        dscouldt = s * plasgdt
        dscoulda = s * plasgda - scoul/abar
        dscouldz = s * plasgdz
       end if




      return
      end subroutine coulomb




      subroutine etages(xni,zbar,temp,eta)
      include 'implno.dek'
      include 'const.dek'

!c..this routine makes a damn good guess for the electron degeneracy
!c..parameter eta.
!c..input is the ion number density xni, average charge zbar,
!c..average atomic weigt abar, and temperature temp.
!c..output is a guess at the chemical potential eta


!c..declare the pass
      double precision  xni,zbar,temp,eta

!c..declare
      double precision xne,x,y,z,kt,beta,tmkt,xnefac

      double precision rt2,rt3,rtpi,cpf0,cpf1,cpf2,cpf3,&
                       twoth,fa0,forpi,mecc
      parameter        (rt2     = 1.4142135623730951d0,&
                        rt3     = 1.7320508075688772d0,&
                        rtpi    = 1.7724538509055159d0,&
                        cpf0    = h/(me*clight),&
                        cpf1    = 3.0d0/(8.0d0*pi) * cpf0**3,&
                        cpf2    = 4.0d0/cpf1,&
                        cpf3    = 2.0d0*rt3*rtpi/(rt2*cpf1),&
                        twoth   = 2.0d0/3.0d0,&
                        fa0     = 64.0d0/(9.0d0*pi),&
                        forpi   = 4.0d0 * pi,&
                        mecc    = me * clight * clight)

!c..notes: rt2=sqrt(2)  rt3=sqrt(3)  rtpi=sqrt(pi)


!c..for the purposes of guessing eta, assume full ionization
      xne   = xni * zbar
      kt    = kerg * temp
      beta  = kt/mecc


!c..number density of ionized electrons (c&g 24.354k) and number density at
!c..turning point (c&g 24.354i). if either of these exceed the number density
!c..as given by a saha equation, then pairs are important. set alfa = 1/2.

      if (beta .ge. 1.0) then
       x = cpf2 * beta * beta
      else
       x = cpf3 * beta * (1.0d0 + 0.75d0*beta) * exp(-1.0d0/beta)
      end if
      if (x .ge. xne) then
       eta = -0.5d0


!c..get the dimensionless number density (c&g 24.313), if it is large apply the
!c..formula (c&g 24.309) to get a possible alfa, if not large do a two term
!c..binomial expansion on (c&g 24.309) to estimate alfa.

      else
       z = (xne*cpf1)**twoth
       if (z .ge. 1.0e-6) then
        y = (sqrt(z + 1.0d0) - 1.0d0)/beta
       else
        y = z * (1.0d0 - z * 0.25d0) * 0.5d0/beta
       end if


!c..isolate the constant in front of the number density integral. if it is
!c..small enough run the divine approximation backwards with c&g 24.43. then
!c..join it smoothly with the lower limit.

       x = log10(xne**0.6d0/temp)
       if (x .le. 9.5) then
        z = ((1.0d0 + fa0*beta)*sqrt(1.0d0 + fa0*beta*0.5) - 1.0d0)/fa0
        tmkt    = 2.0d0 * me/h * kt/h
        xnefac  = forpi * tmkt * sqrt(tmkt)
        eta = -log(xnefac*rtpi*(0.5d0+0.75d0*z)/xne)
        if (x .ge. 8.5) eta = eta*(9.5d0-x) + y * (1.0d0 - (9.5d0-x))
       else
        eta = y
       end if
      end if

      return
      end subroutine etages






!c..routine dfermi gets the fermi-dirac functions and their derivaties
!c..routine fdfunc1 forms the integrand of the fermi-dirac functions
!c..routine fdfunc2 same as fdfunc but with the change of variable z**2=x
!c..routine dqleg010 does 10 point gauss-legendre integration  9 fig accuracy
!c..routine dqleg020 does 20 point gauss-legendre integration 14 fig accuracy
!c..routine dqleg040 does 40 point gauss-legendre integration 18 fig accuracy
!c..routine dqleg080 does 80 point gauss-legendre integration 32 fig accuracy
!c..routine dqlag010 does 10 point gauss-laguerre integration  9 fig accuracy
!c..routine dqlag020 does 20 point gauss-laguerre integration 14 fig accuracy
!c..routine dqlag040 does 40 point gauss-laguerre integration 18 fig accuracy
!c..routine dqlag080 does 80 point gauss-laguerre integration 32 fig accuracy



      subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta)
      include 'implno.dek'
!c..
!c..this routine computes the fermi-dirac integrals of
!c..index dk, with degeneracy parameter eta and relativity parameter theta.
!c..input is dk the double precision index of the fermi-dirac function,
!c..eta the degeneracy parameter, and theta the relativity parameter.
!c..the output is fd is computed by applying three 10-point
!c..gauss-legendre and one 10-point gauss-laguerre rules over
!c..four appropriate subintervals. the derivative with respect to eta is
!c..output in fdeta, and the derivative with respct to theta is in fdtheta.
!c..within each subinterval the fd kernel.
!c..
!c..this routine delivers at least 9 figures of accuracy
!c..
!c..reference: j.m. aparicio, apjs 117, 632 1998
!c..
!c..declare
!c..declare
!      external         fdfunc1,fdfunc2
      double precision dk,eta,theta,fd,fdeta,fdtheta,&
                       d,sg,a1,b1,c1,a2,b2,c2,d2,e2,a3,b3,c3,d3,e3,&
                       eta1,xi,xi2,x1,x2,x3,s1,s2,s3,s12,par(3),&
                       res1,dres1,ddres1,res2,dres2,ddres2,&
                       res3,dres3,ddres3,res4,dres4,ddres4


!c   parameters defining the location of the breakpoints for the
!c   subintervals of integration:
      data d   / 3.3609d00 /
      data sg  / 9.1186d-2 /
      data a1  / 6.7774d00 /
      data b1  / 1.1418d00 /
      data c1  / 2.9826d00 /
      data a2  / 3.7601d00 /
      data b2  / 9.3719d-2 /
      data c2  / 2.1063d-2 /
      data d2  / 3.1084d01 /
      data e2  / 1.0056d00 /
      data a3  / 7.5669d00 /
      data b3  / 1.1695d00 /
      data c3  / 7.5416d-1 /
      data d3  / 6.6558d00 /
      data e3  /-1.2819d-1 /


!c   integrand parameters:
      par(1)=dk
      par(2)=eta
      par(3)=theta


!c   definition of xi:
      eta1=sg*(eta-d)
      if (eta1.le.5.d1) then
        xi=log(1.d0+exp(eta1))/sg
      else
        xi=eta-d
      endif
      xi2=xi*xi

!c   definition of the x_i:
      x1=(a1  +b1*xi+c1*   xi2) &
        /(1.d0+c1*xi)
      x2=(a2  +b2*xi+c2*d2*xi2) &
        /(1.d0+e2*xi+c2*   xi2)
      x3=(a3  +b3*xi+c3*d3*xi2) &
        /(1.d0+e3*xi+c3*   xi2)

!c   breakpoints:
      s1=x1-x2
      s2=x1
      s3=x1+x3
      s12=sqrt(s1)

!c   quadrature integrations:

!c 9 significant figure accuracy
!c      call dqleg010(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
!c      call dqleg010(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
!c      call dqleg010(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
!c      call dqlag010(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

!c 14 significant figure accuracy
      call dqleg020(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
      call dqleg020(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
      call dqleg020(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
      call dqlag020(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

!c 18 significant figure accuracy
!c      call dqleg040(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
!c      call dqleg040(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
!c     call dqleg040(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
!c     call dqlag040(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)

!c 32 significant figure accuracy
!c      call dqleg080(fdfunc2, 0.d0,  s12, res1, dres1, ddres1, par,3)
!c      call dqleg080(fdfunc1,   s1,   s2, res2, dres2, ddres2, par,3)
!c      call dqleg080(fdfunc1,   s2,   s3, res3, dres3, ddres3, par,3)
!c      call dqlag080(fdfunc1,   s3, 1.d0, res4, dres4, ddres4, par,3)


!c..sum the contributions
      fd      = res1 + res2 + res3 + res4
      fdeta   = dres1 + dres2 + dres3 + dres4
      fdtheta = ddres1 + ddres2 + ddres3 + ddres4
      return
      end subroutine dfermi





      subroutine dqleg010(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
!c..
!c..10 point gauss-legendre rule for the fermi-dirac function and
!c..its derivatives with respect to eta and theta.
!c..on input f is the name of the subroutine containing the integrand,
!c..a is the lower end point of the interval, b is the higher end point,
!c..par is an array of constant parameters to be passed to subroutine f,
!c..and n is the length of the par array. on output result is the
!c..approximation from applying the 10-point gauss-legendre rule,
!c..dresult is the derivative with respect to eta, and ddresult is the
!c..derivative with respect to theta.
!c..
!c..note: since the number of nodes is even, zero is not an abscissa.
!c..
!c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),&
                       absc1,absc2,center,hlfrun,wg(5),xg(5),&
                       fval1,dfval1,ddfval1,fval2,dfval2,ddfval2

!c the abscissae and weights are given for the interval (-1,1).
!c xg     - abscissae of the 20-point gauss-legendre rule
!c          for half of the usual run (-1,1), i.e.
!c          the positive nodes of the 20-point rule
!c wg     - weights of the 20-point gauss rule.
!c
!c abscissae and weights were evaluated with 100 decimal digit arithmetic.
!
      data xg (  1) /   1.48874338981631210884826001129719984d-1 /
      data xg (  2) /   4.33395394129247190799265943165784162d-1 /
      data xg (  3) /   6.79409568299024406234327365114873575d-1 /
      data xg (  4) /   8.65063366688984510732096688423493048d-1 /
      data xg (  5) /   9.73906528517171720077964012084452053d-1 /

      data wg (  1) /   2.95524224714752870173892994651338329d-1 /
      data wg (  2) /   2.69266719309996355091226921569469352d-1 /
      data wg (  3) /   2.19086362515982043995534934228163192d-1 /
      data wg (  4) /   1.49451349150580593145776339657697332d-1 /
      data wg (  5) /   6.66713443086881375935688098933317928d-2 /


!c           list of major variables
!c           -----------------------
!c
!c           absc   - abscissa
!c           fval*  - function value
!c           result - result of the 10-point gauss formula

      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,5
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end subroutine dqleg010




      subroutine dqleg020(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
!c..
!c..20 point gauss-legendre rule for the fermi-dirac function and
!c..its derivatives with respect to eta and theta.
!c..on input f is the name of the subroutine containing the integrand,
!c..a is the lower end point of the interval, b is the higher end point,
!c..par is an array of constant parameters to be passed to subroutine f,
!c..and n is the length of the par array. on output result is the
!c..approximation from applying the 20-point gauss-legendre rule,
!c..dresult is the derivative with respect to eta, and ddresult is the
!c..derivative with respect to theta.
!c..
!c..note: since the number of nodes is even, zero is not an abscissa.
!c..
!c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),&
                       absc1,absc2,center,hlfrun,wg(10),xg(10),&
                       fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


!c the abscissae and weights are given for the interval (-1,1).
!c xg     - abscissae of the 20-point gauss-legendre rule
!c          for half of the usual run (-1,1), i.e.
!c          the positive nodes of the 20-point rule
!c wg     - weights of the 20-point gauss rule.
!c
!c abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   7.65265211334973337546404093988382110d-2 /
      data xg (  2) /   2.27785851141645078080496195368574624d-1 /
      data xg (  3) /   3.73706088715419560672548177024927237d-1 /
      data xg (  4) /   5.10867001950827098004364050955250998d-1 /
      data xg (  5) /   6.36053680726515025452836696226285936d-1 /
      data xg (  6) /   7.46331906460150792614305070355641590d-1 /
      data xg (  7) /   8.39116971822218823394529061701520685d-1 /
      data xg (  8) /   9.12234428251325905867752441203298113d-1 /
      data xg (  9) /   9.63971927277913791267666131197277221d-1 /
      data xg ( 10) /   9.93128599185094924786122388471320278d-1 /

      data wg (  1) /   1.52753387130725850698084331955097593d-1 /
      data wg (  2) /   1.49172986472603746787828737001969436d-1 /
      data wg (  3) /   1.42096109318382051329298325067164933d-1 /
      data wg (  4) /   1.31688638449176626898494499748163134d-1 /
      data wg (  5) /   1.18194531961518417312377377711382287d-1 /
      data wg (  6) /   1.01930119817240435036750135480349876d-1 /
      data wg (  7) /   8.32767415767047487247581432220462061d-2 /
      data wg (  8) /   6.26720483341090635695065351870416063d-2 /
      data wg (  9) /   4.06014298003869413310399522749321098d-2 /
      data wg ( 10) /   1.76140071391521183118619623518528163d-2 /


!c           list of major variables
!c           -----------------------
!c
!c           absc   - abscissa
!c           fval*  - function value
!c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,10
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end subroutine dqleg020





      subroutine dqleg040(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
!c..
!c..40 point gauss-legendre rule for the fermi-dirac function and
!c..its derivatives with respect to eta and theta.
!c..on input f is the name of the subroutine containing the integrand,
!c..a is the lower end point of the interval, b is the higher end point,
!c..par is an array of constant parameters to be passed to subroutine f,
!c..and n is the length of the par array. on output result is the
!c..approximation from applying the 40-point gauss-legendre rule,
!c..dresult is the derivative with respect to eta, and ddresult is the
!c..derivative with respect to theta.
!c..
!c..note: since the number of nodes is even, zero is not an abscissa.
!c..
!c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),&
                       absc1,absc2,center,hlfrun,wg(20),xg(20),&
                       fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


!c the abscissae and weights are given for the interval (-1,1).
!c xg     - abscissae of the 40-point gauss-legendre rule
!c          for half of the usual run (-1,1), i.e.
!c          the positive nodes of the 40-point rule
!c wg     - weights of the 40-point gauss rule.
!c
!c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.87724175060508219331934440246232946d-2 /
      data xg (  2) /   1.16084070675255208483451284408024113d-1 /
      data xg (  3) /   1.92697580701371099715516852065149894d-1 /
      data xg (  4) /   2.68152185007253681141184344808596183d-1 /
      data xg (  5) /   3.41994090825758473007492481179194310d-1 /
      data xg (  6) /   4.13779204371605001524879745803713682d-1 /
      data xg (  7) /   4.83075801686178712908566574244823004d-1 /
      data xg (  8) /   5.49467125095128202075931305529517970d-1 /
      data xg (  9) /   6.12553889667980237952612450230694877d-1 /
      data xg ( 10) /   6.71956684614179548379354514961494109d-1 /
      data xg ( 11) /   7.27318255189927103280996451754930548d-1 /
      data xg ( 12) /   7.78305651426519387694971545506494848d-1 /
      data xg ( 13) /   8.24612230833311663196320230666098773d-1 /
      data xg ( 14) /   8.65959503212259503820781808354619963d-1 /
      data xg ( 15) /   9.02098806968874296728253330868493103d-1 /
      data xg ( 16) /   9.32812808278676533360852166845205716d-1 /
      data xg ( 17) /   9.57916819213791655804540999452759285d-1 /
      data xg ( 18) /   9.77259949983774262663370283712903806d-1 /
      data xg ( 19) /   9.90726238699457006453054352221372154d-1 /
      data xg ( 20) /   9.98237709710559200349622702420586492d-1 /

      data wg (  1) /   7.75059479784248112637239629583263269d-2 /
      data wg (  2) /   7.70398181642479655883075342838102485d-2 /
      data wg (  3) /   7.61103619006262423715580759224948230d-2 /
      data wg (  4) /   7.47231690579682642001893362613246731d-2 /
      data wg (  5) /   7.28865823958040590605106834425178358d-2 /
      data wg (  6) /   7.06116473912867796954836308552868323d-2 /
      data wg (  7) /   6.79120458152339038256901082319239859d-2 /
      data wg (  8) /   6.48040134566010380745545295667527300d-2 /
      data wg (  9) /   6.13062424929289391665379964083985959d-2 /
      data wg ( 10) /   5.74397690993915513666177309104259856d-2 /
      data wg ( 11) /   5.32278469839368243549964797722605045d-2 /
      data wg ( 12) /   4.86958076350722320614341604481463880d-2 /
      data wg ( 13) /   4.38709081856732719916746860417154958d-2 /
      data wg ( 14) /   3.87821679744720176399720312904461622d-2 /
      data wg ( 15) /   3.34601952825478473926781830864108489d-2 /
      data wg ( 16) /   2.79370069800234010984891575077210773d-2 /
      data wg ( 17) /   2.22458491941669572615043241842085732d-2 /
      data wg ( 18) /   1.64210583819078887128634848823639272d-2 /
      data wg ( 19) /   1.04982845311528136147421710672796523d-2 /
      data wg ( 20) /   4.52127709853319125847173287818533272d-3 /


!c           list of major variables
!c           -----------------------
!c
!c           absc   - abscissa
!c           fval*  - function value
!c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,20
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end subroutine dqleg040





      subroutine dqleg080(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
!c..
!c..80 point gauss-legendre rule for the fermi-dirac function and
!c..its derivatives with respect to eta and theta.
!c..on input f is the name of the subroutine containing the integrand,
!c..a is the lower end point of the interval, b is the higher end point,
!c..par is an array of constant parameters to be passed to subroutine f,
!c..and n is the length of the par array. on output result is the
!c..approximation from applying the 80-point gauss-legendre rule,
!c..dresult is the derivative with respect to eta, and ddresult is the
!c..derivative with respect to theta.
!c..
!c..note: since the number of nodes is even, zero is not an abscissa.
!c..
!c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),&
                       absc1,absc2,center,hlfrun,wg(40),xg(40),&
                       fval1,dfval1,ddfval1,fval2,dfval2,ddfval2


!c the abscissae and weights are given for the interval (-1,1).
!c xg     - abscissae of the 80-point gauss-legendre rule
!c          for half of the usual run (-1,1), i.e.
!c          the positive nodes of the 80-point rule
!c wg     - weights of the 80-point gauss rule.
!c
!c abscissae and weights were evaluated with 100 decimal digit arithmetic.


      data xg (  1) /   1.95113832567939976543512341074545479d-2 /
      data xg (  2) /   5.85044371524206686289933218834177944d-2 /
      data xg (  3) /   9.74083984415845990632784501049369020d-2 /
      data xg (  4) /   1.36164022809143886559241078000717067d-1 /
      data xg (  5) /   1.74712291832646812559339048011286195d-1 /
      data xg (  6) /   2.12994502857666132572388538666321823d-1 /
      data xg (  7) /   2.50952358392272120493158816035004797d-1 /
      data xg (  8) /   2.88528054884511853109139301434713898d-1 /
      data xg (  9) /   3.25664370747701914619112943627358695d-1 /
      data xg ( 10) /   3.62304753499487315619043286358963588d-1 /
      data xg ( 11) /   3.98393405881969227024379642517533757d-1 /
      data xg ( 12) /   4.33875370831756093062386700363181958d-1 /
      data xg ( 13) /   4.68696615170544477036078364935808657d-1 /
      data xg ( 14) /   5.02804111888784987593672750367568003d-1 /
      data xg ( 15) /   5.36145920897131932019857253125400904d-1 /
      data xg ( 16) /   5.68671268122709784725485786624827158d-1 /
      data xg ( 17) /   6.00330622829751743154746299164006848d-1 /
      data xg ( 18) /   6.31075773046871966247928387289336863d-1 /
      data xg ( 19) /   6.60859898986119801735967122844317234d-1 /
      data xg ( 20) /   6.89637644342027600771207612438935266d-1 /
      data xg ( 21) /   7.17365185362099880254068258293815278d-1 /
      data xg ( 22) /   7.44000297583597272316540527930913673d-1 /
      data xg ( 23) /   7.69502420135041373865616068749026083d-1 /
      data xg ( 24) /   7.93832717504605449948639311738454358d-1 /
      data xg ( 25) /   8.16954138681463470371124994012295707d-1 /
      data xg ( 26) /   8.38831473580255275616623043902867064d-1 /
      data xg ( 27) /   8.59431406663111096977192123491656492d-1 /
      data xg ( 28) /   8.78722567678213828703773343639124407d-1 /
      data xg ( 29) /   8.96675579438770683194324071967395986d-1 /
      data xg ( 30) /   9.13263102571757654164733656150947478d-1 /
      data xg ( 31) /   9.28459877172445795953045959075453133d-1 /
      data xg ( 32) /   9.42242761309872674752266004500001735d-1 /
      data xg ( 33) /   9.54590766343634905493481517021029508d-1 /
      data xg ( 34) /   9.65485089043799251452273155671454998d-1 /
      data xg ( 35) /   9.74909140585727793385645230069136276d-1 /
      data xg ( 36) /   9.82848572738629070418288027709116473d-1 /
      data xg ( 37) /   9.89291302499755531026503167136631385d-1 /
      data xg ( 38) /   9.94227540965688277892063503664911698d-1 /
      data xg ( 39) /   9.97649864398237688899494208183122985d-1 /
      data xg ( 40) /   9.99553822651630629880080499094567184d-1 /

      data wg (  1) /   3.90178136563066548112804392527540483d-2 /
      data wg (  2) /   3.89583959627695311986255247722608223d-2 /
      data wg (  3) /   3.88396510590519689317741826687871658d-2 /
      data wg (  4) /   3.86617597740764633270771102671566912d-2 /
      data wg (  5) /   3.84249930069594231852124363294901384d-2 /
      data wg (  6) /   3.81297113144776383442067915657362019d-2 /
      data wg (  7) /   3.77763643620013974897749764263210547d-2 /
      data wg (  8) /   3.73654902387304900267053770578386691d-2 /
      data wg (  9) /   3.68977146382760088391509965734052192d-2 /
      data wg ( 10) /   3.63737499058359780439649910465228136d-2 /
      data wg ( 11) /   3.57943939534160546028615888161544542d-2 /
      data wg ( 12) /   3.51605290447475934955265923886968812d-2 /
      data wg ( 13) /   3.44731204517539287943642267310298320d-2 /
      data wg ( 14) /   3.37332149846115228166751630642387284d-2 /
      data wg ( 15) /   3.29419393976454013828361809019595361d-2 /
      data wg ( 16) /   3.21004986734877731480564902872506960d-2 /
      data wg ( 17) /   3.12101741881147016424428667206035518d-2 /
      data wg ( 18) /   3.02723217595579806612200100909011747d-2 /
      data wg ( 19) /   2.92883695832678476927675860195791396d-2 /
      data wg ( 20) /   2.82598160572768623967531979650145302d-2 /
      data wg ( 21) /   2.71882275004863806744187066805442598d-2 /
      data wg ( 22) /   2.60752357675651179029687436002692871d-2 /
      data wg ( 23) /   2.49225357641154911051178470032198023d-2 /
      data wg ( 24) /   2.37318828659301012931925246135684162d-2 /
      data wg ( 25) /   2.25050902463324619262215896861687390d-2 /
      data wg ( 26) /   2.12440261157820063887107372506131285d-2 /
      data wg ( 27) /   1.99506108781419989288919287151135633d-2 /
      data wg ( 28) /   1.86268142082990314287354141521572090d-2 /
      data wg ( 29) /   1.72746520562693063585842071312909998d-2 /
      data wg ( 30) /   1.58961835837256880449029092291785257d-2 /
      data wg ( 31) /   1.44935080405090761169620745834605500d-2 /
      data wg ( 32) /   1.30687615924013392937868258970563403d-2 /
      data wg ( 33) /   1.16241141207978269164667699954326348d-2 /
      data wg ( 34) /   1.01617660411030645208318503524069436d-2 /
      data wg ( 35) /   8.68394526926085842640945220403428135d-3 /
      data wg ( 36) /   7.19290476811731275267557086795650747d-3 /
      data wg ( 37) /   5.69092245140319864926910711716201847d-3 /
      data wg ( 38) /   4.18031312469489523673930420168135132d-3 /
      data wg ( 39) /   2.66353358951268166929353583166845546d-3 /
      data wg ( 40) /   1.14495000318694153454417194131563611d-3 /


!c           list of major variables
!c           -----------------------
!c
!c           absc   - abscissa
!c           fval*  - function value
!c           result - result of the 20-point gauss formula


      center   = 0.5d0 * (a+b)
      hlfrun   = 0.5d0 * (b-a)
      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,40
        absc1 = center + hlfrun*xg(j)
        absc2 = center - hlfrun*xg(j)
        call f(absc1, par, n, fval1, dfval1, ddfval1)
        call f(absc2, par, n, fval2, dfval2, ddfval2)
        result   = result + (fval1 + fval2)*wg(j)
        dresult  = dresult + (dfval1 + dfval2)*wg(j)
        ddresult = ddresult + (ddfval1 + ddfval2)*wg(j)
      enddo
      result   = result * hlfrun
      dresult  = dresult * hlfrun
      ddresult = ddresult * hlfrun
      return
      end subroutine dqleg080





      subroutine dqlag010(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
!c..
!c..10 point gauss-laguerre rule for the fermi-dirac function.
!c..on input f is the external function defining the integrand
!c..f(x)=g(x)*w(x), where w(x) is the gaussian weight
!c..w(x)=exp(-(x-a)/b) and g(x) a smooth function,
!c..a is the lower end point of the interval, b is the higher end point,
!c..par is an array of constant parameters to be passed to the function f,
!c..and n is the length of the par array. on output result is the
!c..approximation from applying the 10-point gauss-laguerre rule.
!c..since the number of nodes is even, zero is not an abscissa.
!c..
!c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),&
                       absc,wg(10),xg(10),fval,dfval,ddfval


!c the abscissae and weights are given for the interval (0,+inf).
!c xg     - abscissae of the 10-point gauss-laguerre rule
!c wg     - weights of the 10-point gauss rule. since f yet
!c          includes the weight function, the values in wg
!c          are actually exp(xg) times the standard
!c          gauss-laguerre weights
!c
!c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.37793470540492430830772505652711188d-1 /
      data xg (  2) /   7.29454549503170498160373121676078781d-1 /
      data xg (  3) /   1.80834290174031604823292007575060883d00 /
      data xg (  4) /   3.40143369785489951448253222140839067d00 /
      data xg (  5) /   5.55249614006380363241755848686876285d00 /
      data xg (  6) /   8.33015274676449670023876719727452218d00 /
      data xg (  7) /   1.18437858379000655649185389191416139d01 /
      data xg (  8) /   1.62792578313781020995326539358336223d01 /
      data xg (  9) /   2.19965858119807619512770901955944939d01 /
      data xg ( 10) /   2.99206970122738915599087933407991951d01 /

      data wg (  1) /   3.54009738606996308762226891442067608d-1 /
      data wg (  2) /   8.31902301043580738109829658127849577d-1 /
      data wg (  3) /   1.33028856174932817875279219439399369d00 /
      data wg (  4) /   1.86306390311113098976398873548246693d00 /
      data wg (  5) /   2.45025555808301016607269373165752256d00 /
      data wg (  6) /   3.12276415513518249615081826331455472d00 /
      data wg (  7) /   3.93415269556152109865581245924823077d00 /
      data wg (  8) /   4.99241487219302310201148565243315445d00 /
      data wg (  9) /   6.57220248513080297518766871037611234d00 /
      data wg ( 10) /   9.78469584037463069477008663871859813d00 /


!c           list of major variables
!c           -----------------------
!c           absc   - abscissa
!c           fval*  - function value
!c           result - result of the 10-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,10
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end subroutine dqlag010





      subroutine dqlag020(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
!c..
!c..20 point gauss-laguerre rule for the fermi-dirac function.
!c..on input f is the external function defining the integrand
!c..f(x)=g(x)*w(x), where w(x) is the gaussian weight
!c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
!c..a is the lower end point of the interval, b is the higher end point,
!c..par is an array of constant parameters to be passed to the function f,
!c..and n is the length of the par array. on output result is the
!c..approximation from applying the 20-point gauss-laguerre rule.
!c..since the number of nodes is even, zero is not an abscissa.
!c..
!c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),&
                       absc,wg(20),xg(20),fval,dfval,ddfval


!c the abscissae and weights are given for the interval (0,+inf).
!c xg     - abscissae of the 20-point gauss-laguerre rule
!c wg     - weights of the 20-point gauss rule. since f yet
!c          includes the weight function, the values in wg
!c          are actually exp(xg) times the standard
!c          gauss-laguerre weights
!c
!c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   7.05398896919887533666890045842150958d-2 /
      data xg (  2) /   3.72126818001611443794241388761146636d-1 /
      data xg (  3) /   9.16582102483273564667716277074183187d-1 /
      data xg (  4) /   1.70730653102834388068768966741305070d00 /
      data xg (  5) /   2.74919925530943212964503046049481338d00 /
      data xg (  6) /   4.04892531385088692237495336913333219d00 /
      data xg (  7) /   5.61517497086161651410453988565189234d00 /
      data xg (  8) /   7.45901745367106330976886021837181759d00 /
      data xg (  9) /   9.59439286958109677247367273428279837d00 /
      data xg ( 10) /   1.20388025469643163096234092988655158d01 /
      data xg ( 11) /   1.48142934426307399785126797100479756d01 /
      data xg ( 12) /   1.79488955205193760173657909926125096d01 /
      data xg ( 13) /   2.14787882402850109757351703695946692d01 /
      data xg ( 14) /   2.54517027931869055035186774846415418d01 /
      data xg ( 15) /   2.99325546317006120067136561351658232d01 /
      data xg ( 16) /   3.50134342404790000062849359066881395d01 /
      data xg ( 17) /   4.08330570567285710620295677078075526d01 /
      data xg ( 18) /   4.76199940473465021399416271528511211d01 /
      data xg ( 19) /   5.58107957500638988907507734444972356d01 /
      data xg ( 20) /   6.65244165256157538186403187914606659d01 /

      data wg (  1) /   1.81080062418989255451675405913110644d-1 /
      data wg (  2) /   4.22556767878563974520344172566458197d-1 /
      data wg (  3) /   6.66909546701848150373482114992515927d-1 /
      data wg (  4) /   9.15352372783073672670604684771868067d-1 /
      data wg (  5) /   1.16953970719554597380147822239577476d00 /
      data wg (  6) /   1.43135498592820598636844994891514331d00 /
      data wg (  7) /   1.70298113798502272402533261633206720d00 /
      data wg (  8) /   1.98701589079274721410921839275129020d00 /
      data wg (  9) /   2.28663578125343078546222854681495651d00 /
      data wg ( 10) /   2.60583472755383333269498950954033323d00 /
      data wg ( 11) /   2.94978373421395086600235416827285951d00 /
      data wg ( 12) /   3.32539578200931955236951937421751118d00 /
      data wg ( 13) /   3.74225547058981092111707293265377811d00 /
      data wg ( 14) /   4.21423671025188041986808063782478746d00 /
      data wg ( 15) /   4.76251846149020929695292197839096371d00 /
      data wg ( 16) /   5.42172604424557430380308297989981779d00 /
      data wg ( 17) /   6.25401235693242129289518490300707542d00 /
      data wg ( 18) /   7.38731438905443455194030019196464791d00 /
      data wg ( 19) /   9.15132873098747960794348242552950528d00 /
      data wg ( 20) /   1.28933886459399966710262871287485278d01 /


!c           list of major variables
!c           -----------------------
!c           absc   - abscissa
!c           fval*  - function value
!c           result - result of the 20-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,20
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end subroutine dqlag020




      subroutine dqlag040(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
!c..
!c..20 point gauss-laguerre rule for the fermi-dirac function.
!c..on input f is the external function defining the integrand
!c..f(x)=g(x)*w(x), where w(x) is the gaussian weight
!c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
!c..a is the lower end point of the interval, b is the higher end point,
!c..par is an array of constant parameters to be passed to the function f,
!c..and n is the length of the par array. on output result is the
!c..approximation from applying the 20-point gauss-laguerre rule.
!c..since the number of nodes is even, zero is not an abscissa.
!c..
!c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),&
                       absc,wg(40),xg(40),fval,dfval,ddfval


!c the abscissae and weights are given for the interval (0,+inf).
!c xg     - abscissae of the 20-point gauss-laguerre rule
!c wg     - weights of the 20-point gauss rule. since f yet
!c          includes the weight function, the values in wg
!c          are actually exp(xg) times the standard
!c          gauss-laguerre weights
!c
!c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   3.57003943088883851220844712866008554d-2 /
      data xg (  2) /   1.88162283158698516003589346219095913d-1 /
      data xg (  3) /   4.62694281314576453564937524561190364d-1 /
      data xg (  4) /   8.59772963972934922257272224688722412d-1 /
      data xg (  5) /   1.38001082052733718649800032959526559d00 /
      data xg (  6) /   2.02420913592282673344206600280013075d00 /
      data xg (  7) /   2.79336935350681645765351448602664039d00 /
      data xg (  8) /   3.68870267790827020959152635190868698d00 /
      data xg (  9) /   4.71164114655497269361872283627747369d00 /
      data xg ( 10) /   5.86385087834371811427316423799582987d00 /
      data xg ( 11) /   7.14724790810228825068569195197942362d00 /
      data xg ( 12) /   8.56401701758616376271852204208813232d00 /
      data xg ( 13) /   1.01166340484519394068496296563952448d01 /
      data xg ( 14) /   1.18078922940045848428415867043606304d01 /
      data xg ( 15) /   1.36409337125370872283716763606501202d01 /
      data xg ( 16) /   1.56192858933390738372019636521880145d01 /
      data xg ( 17) /   1.77469059500956630425738774954243772d01 /
      data xg ( 18) /   2.00282328345748905296126148101751172d01 /
      data xg ( 19) /   2.24682499834984183513717862289945366d01 /
      data xg ( 20) /   2.50725607724262037943960862094009769d01 /
      data xg ( 21) /   2.78474800091688627207517041404557997d01 /
      data xg ( 22) /   3.08001457394454627007543851961911114d01 /
      data xg ( 23) /   3.39386570849137196090988585862819990d01 /
      data xg ( 24) /   3.72722458804760043283207609906074207d01 /
      data xg ( 25) /   4.08114928238869204661556755816006426d01 /
      data xg ( 26) /   4.45686031753344627071230206344983559d01 /
      data xg ( 27) /   4.85577635330599922809620488067067936d01 /
      data xg ( 28) /   5.27956111872169329693520211373917638d01 /
      data xg ( 29) /   5.73018633233936274950337469958921651d01 /
      data xg ( 30) /   6.21001790727751116121681990578989921d01 /
      data xg ( 31) /   6.72193709271269987990802775518887054d01 /
      data xg ( 32) /   7.26951588476124621175219277242619385d01 /
      data xg ( 33) /   7.85728029115713092805438968334812596d01 /
      data xg ( 34) /   8.49112311357049845427015647096663186d01 /
      data xg ( 35) /   9.17898746712363769923371934806273153d01 /
      data xg ( 36) /   9.93208087174468082501090541654868123d01 /
      data xg ( 37) /   1.07672440639388272520796767611322664d02 /
      data xg ( 38) /   1.17122309512690688807650644123550702d02 /
      data xg ( 39) /   1.28201841988255651192541104389631263d02 /
      data xg ( 40) /   1.42280044469159997888348835359541764d02 /

      data wg (  1) /   9.16254711574598973115116980801374830d-2 /
      data wg (  2) /   2.13420584905012080007193367121512341d-1 /
      data wg (  3) /   3.35718116680284673880510701616292191d-1 /
      data wg (  4) /   4.58540935033497560385432380376452497d-1 /
      data wg (  5) /   5.82068165779105168990996365401543283d-1 /
      data wg (  6) /   7.06495216367219392989830015673016682d-1 /
      data wg (  7) /   8.32026903003485238099112947978349523d-1 /
      data wg (  8) /   9.58878198794443111448122679676028906d-1 /
      data wg (  9) /   1.08727616203054971575386933317202661d00 /
      data wg ( 10) /   1.21746232797778097895427785066560948d00 /
      data wg ( 11) /   1.34969549135676530792393859442394519d00 /
      data wg ( 12) /   1.48425492977684671120561178612978719d00 /
      data wg ( 13) /   1.62144416281182197802316884316454527d00 /
      data wg ( 14) /   1.76159537467676961118424220420981598d00 /
      data wg ( 15) /   1.90507466589479967668299320597279371d00 /
      data wg ( 16) /   2.05228834726171671760199582272947454d00 /
      data wg ( 17) /   2.20369055324509588909828344328140570d00 /
      data wg ( 18) /   2.35979253852320332354037375378901497d00 /
      data wg ( 19) /   2.52117414037643299165313690287422820d00 /
      data wg ( 20) /   2.68849805540884226415950544706374659d00 /
      data wg ( 21) /   2.86252781321044881203476395983104311d00 /
      data wg ( 22) /   3.04415066531151710041043967954333670d00 /
      data wg ( 23) /   3.23440709726353194177490239428867111d00 /
      data wg ( 24) /   3.43452939842774809220398481891602464d00 /
      data wg ( 25) /   3.64599282499408907238965646699490434d00 /
      data wg ( 26) /   3.87058459721651656808475320213444338d00 /
      data wg ( 27) /   4.11049868043282265583582247263951577d00 /
      data wg ( 28) /   4.36846872325406347450808338272945025d00 /
      data wg ( 29) /   4.64795898407446688299303399883883991d00 /
      data wg ( 30) /   4.95344611240989326218696150785562721d00 /
      data wg ( 31) /   5.29084840590073657468737365718858968d00 /
      data wg ( 32) /   5.66820460903297677000730529023263795d00 /
      data wg ( 33) /   6.09679641474342030593376010859198806d00 /
      data wg ( 34) /   6.59310886103999953794429664206294899d00 /
      data wg ( 35) /   7.18249599553689315064429801626699574d00 /
      data wg ( 36) /   7.90666631138422877369310742310586595d00 /
      data wg ( 37) /   8.84089249281034652079125595063026792d00 /
      data wg ( 38) /   1.01408992656211694839094600306940468d01 /
      data wg ( 39) /   1.22100212992046038985226485875881108d01 /
      data wg ( 40) /   1.67055206420242974052468774398573553d01 /


!c           list of major variables
!c           -----------------------
!c           absc   - abscissa
!c           fval*  - function value
!c           result - result of the 20-point gauss formula


      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,40
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end subroutine dqlag040





      subroutine dqlag080(f,a,b,result,dresult,ddresult,par,n)
      include 'implno.dek'
!c..
!c..20 point gauss-laguerre rule for the fermi-dirac function.
!c..on input f is the external function defining the integrand
!c..f(x)=g(x)*w(x), where w(x) is the gaussian weight
!c..w(x)=dexp(-(x-a)/b) and g(x) a smooth function,
!c..a is the lower end point of the interval, b is the higher end point,
!c..par is an array of constant parameters to be passed to the function f,
!c..and n is the length of the par array. on output result is the
!c..approximation from applying the 20-point gauss-laguerre rule.
!c..since the number of nodes is even, zero is not an abscissa.
!c..
!c..declare
      external         f
      integer          j,n
      double precision a,b,result,dresult,ddresult,par(n),&
                      absc,wg(80),xg(80),fval,dfval,ddfval


!c the abscissae and weights are given for the interval (0,+inf).
!c xg     - abscissae of the 20-point gauss-laguerre rule
!c wg     - weights of the 20-point gauss rule. since f yet
!c          includes the weight function, the values in wg
!c          are actually exp(xg) times the standard
!c          gauss-laguerre weights
!c
!c abscissae and weights were evaluated with 100 decimal digit arithmetic.

      data xg (  1) /   1.79604233006983655540103192474016803d-2 /
      data xg (  2) /   9.46399129943539888113902724652172943d-2 /
      data xg (  3) /   2.32622868125867569207706157216349831d-1 /
      data xg (  4) /   4.31992547802387480255786172497770411d-1 /
      data xg (  5) /   6.92828861352021839905702213635446867d-1 /
      data xg (  6) /   1.01523255618947143744625436859935350d00 /
      data xg (  7) /   1.39932768784287277414419051430978382d00 /
      data xg (  8) /   1.84526230383584513811177117769599966d00 /
      data xg (  9) /   2.35320887160926152447244708016140181d00 /
      data xg ( 10) /   2.92336468655542632483691234259732862d00 /
      data xg ( 11) /   3.55595231404613405944967308324638370d00 /
      data xg ( 12) /   4.25122008230987808316485766448577637d00 /
      data xg ( 13) /   5.00944263362016477243367706818206389d00 /
      data xg ( 14) /   5.83092153860871901982127113295605083d00 /
      data xg ( 15) /   6.71598597785131711156550087635199430d00 /
      data xg ( 16) /   7.66499349489177306073418909047823480d00 /
      data xg ( 17) /   8.67833082516770109543442255542661083d00 /
      data xg ( 18) /   9.75641480574293071316550944617366591d00 /
      data xg ( 19) /   1.08996933712878553774361001021489406d01 /
      data xg ( 20) /   1.21086466423656999007054848698315593d01 /
      data xg ( 21) /   1.33837881127786473701629840603833297d01 /
      data xg ( 22) /   1.47256659435085855393358076261838437d01 /
      data xg ( 23) /   1.61348643716624665791658545428990907d01 /
      data xg ( 24) /   1.76120052438144378598635686943586520d01 /
      data xg ( 25) /   1.91577496842412479221729970205674985d01 /
      data xg ( 26) /   2.07727999097920960924419379010489579d01 /
      data xg ( 27) /   2.24579012045404583114095916950877516d01 /
      data xg ( 28) /   2.42138440689586473771922469392447092d01 /
      data xg ( 29) /   2.60414665601655866929390053565435682d01 /
      data xg ( 30) /   2.79416568418594655558233069293692111d01 /
      data xg ( 31) /   2.99153559649009855011270412115737715d01 /
      data xg ( 32) /   3.19635609022089207107748887542636533d01 /
      data xg ( 33) /   3.40873278647261898749834947342860505d01 /
      data xg ( 34) /   3.62877759287814544588031988436216948d01 /
      data xg ( 35) /   3.85660910092922104582563052172908535d01 /
      data xg ( 36) /   4.09235302180312671999095850595544326d01 /
      data xg ( 37) /   4.33614266517312302957826760468219500d01 /
      data xg ( 38) /   4.58811946612788863456266489974878378d01 /
      data xg ( 39) /   4.84843356608331891358737273353563006d01 /
      data xg ( 40) /   5.11724445446070105959889432334907144d01 /
      data xg ( 41) /   5.39472167895544471206210278787572430d01 /
      data xg ( 42) /   5.68104563346362231341248503244102122d01 /
      data xg ( 43) /   5.97640843421099549427295961277471927d01 /
      data xg ( 44) /   6.28101489639264772036272917590288682d01 /
      data xg ( 45) /   6.59508362574560573434640627160792248d01 /
      data xg ( 46) /   6.91884824202362773741980288648237373d01 /
      data xg ( 47) /   7.25255875442633453588389652616568450d01 /
      data xg ( 48) /   7.59648311278641748269449794974796502d01 /
      data xg ( 49) /   7.95090896290888369620572826259980809d01 /
      data xg ( 50) /   8.31614564010536896630429506875848705d01 /
      data xg ( 51) /   8.69252644196156234481165926040448396d01 /
      data xg ( 52) /   9.08041123009407559518411727820318427d01 /
      data xg ( 53) /   9.48018942159474332072071889138735302d01 /
      data xg ( 54) /   9.89228344469405791648019372738036790d01 /
      data xg ( 55) /   1.03171527508039130233047094167345654d02 /
      data xg ( 56) /   1.07552984977539906327607890798975954d02 /
      data xg ( 57) /   1.12072690484128333623930046166211013d02 /
      data xg ( 58) /   1.16736664673503666318157888130801099d02 /
      data xg ( 59) /   1.21551542490952625566863895752110813d02 /
      data xg ( 60) /   1.26524665796515540341570265431653573d02 /
      data xg ( 61) /   1.31664195252120310870089086308006192d02 /
      data xg ( 62) /   1.36979246686936973947570637289463788d02 /
      data xg ( 63) /   1.42480058912161601930826569200455232d02 /
      data xg ( 64) /   1.48178202455004441818652384836007732d02 /
      data xg ( 65) /   1.54086842281798697859417425265596259d02 /
      data xg ( 66) /   1.60221072870095715935268416893010646d02 /
      data xg ( 67) /   1.66598351934053918744521179733712213d02 /
      data xg ( 68) /   1.73239071334249503830906503775056999d02 /
      data xg ( 69) /   1.80167323049032317982430208997701523d02 /
      data xg ( 70) /   1.87411949676963772390490134588021771d02 /
      data xg ( 71) /   1.95008022441532991450390479600599643d02 /
      data xg ( 72) /   2.02998984195074937824807677823714777d02 /
      data xg ( 73) /   2.11439870494836466691484904695542608d02 /
      data xg ( 74) /   2.20402368151735739654044206677763168d02 /
      data xg ( 75) /   2.29983206075680004348410969675844754d02 /
      data xg ( 76) /   2.40319087055841540417597460479219628d02 /
      data xg ( 77) /   2.51615879330499611167444939310973194d02 /
      data xg ( 78) /   2.64213823883199102097696108691435553d02 /
      data xg ( 79) /   2.78766733046004563652014172530611597d02 /
      data xg ( 80) /   2.96966511995651345758852859155703581d02 /

      data wg (  1) /   4.60931031330609664705251321395510083d-2 /
      data wg (  2) /   1.07313007783932752564150320304398860d-1 /
      data wg (  3) /   1.68664429547948111794220457782702406d-1 /
      data wg (  4) /   2.30088089384940054411257181978193282d-1 /
      data wg (  5) /   2.91601302502437964832169318772943752d-1 /
      data wg (  6) /   3.53226753575408236352723125805647046d-1 /
      data wg (  7) /   4.14988177550940466187197686311280092d-1 /
      data wg (  8) /   4.76909792302936241314777025418505661d-1 /
      data wg (  9) /   5.39016218474955374499507656522327912d-1 /
      data wg ( 10) /   6.01332497447190529086765248840739512d-1 /
      data wg ( 11) /   6.63884136396680571849442240727299214d-1 /
      data wg ( 12) /   7.26697163614156688973567296249140514d-1 /
      data wg ( 13) /   7.89798189428428531349793078398788294d-1 /
      data wg ( 14) /   8.53214471438152298354598162431362968d-1 /
      data wg ( 15) /   9.16973983833892698590342900031553302d-1 /
      data wg ( 16) /   9.81105491004005747195060155984218607d-1 /
      data wg ( 17) /   1.04563862580654218147568445663176029d00 /
      data wg ( 18) /   1.11060397300025890771124763259729371d00 /
      data wg ( 19) /   1.17603315841226175056651076519208666d00 /
      data wg ( 20) /   1.24195894449809359279351761817871338d00 /
      data wg ( 21) /   1.30841533303134064261188542845954645d00 /
      data wg ( 22) /   1.37543767574892843813155917093490796d00 /
      data wg ( 23) /   1.44306279387849270398312417207247308d00 /
      data wg ( 24) /   1.51132910758830693847655020559917703d00 /
      data wg ( 25) /   1.58027677653099415830201878723121659d00 /
      data wg ( 26) /   1.64994785280267874116012042819355036d00 /
      data wg ( 27) /   1.72038644781283277182004281452290770d00 /
      data wg ( 28) /   1.79163891476093832891442620527688915d00 /
      data wg ( 29) /   1.86375404864909708435925709028688162d00 /
      data wg ( 30) /   1.93678330603070923513925434327841646d00 /
      data wg ( 31) /   2.01078104701134222912614988175555546d00 /
      data wg ( 32) /   2.08580480238741046429303978512989079d00 /
      data wg ( 33) /   2.16191556924159897378316344048827763d00 /
      data wg ( 34) /   2.23917813882364652373453997447445645d00 /
      data wg ( 35) /   2.31766146114651854068606048043496370d00 /
      data wg ( 36) /   2.39743905144001430514117238638849980d00 /
      data wg ( 37) /   2.47858944444973417756369164455222527d00 /
      data wg ( 38) /   2.56119670357790455335115509222572643d00 /
      data wg ( 39) /   2.64535099306968892850463441000367534d00 /
      data wg ( 40) /   2.73114922289915138861410287131169260d00 /
      data wg ( 41) /   2.81869577775934171703141873747811157d00 /
      data wg ( 42) /   2.90810334368223018934550276777492687d00 /
      data wg ( 43) /   2.99949384839685626832412451829968724d00 /
      data wg ( 44) /   3.09299953469357468116695108353033660d00 /
      data wg ( 45) /   3.18876418994712376429365271501623466d00 /
      data wg ( 46) /   3.28694455975337531998378107012216956d00 /
      data wg ( 47) /   3.38771197960397652334054908762154571d00 /
      data wg ( 48) /   3.49125426598732012281732423782764895d00 /
      data wg ( 49) /   3.59777791769613046096294730174902943d00 /
      data wg ( 50) /   3.70751069001745708341027155659228179d00 /
      data wg ( 51) /   3.82070461965311695152029959430467622d00 /
      data wg ( 52) /   3.93763959771430720676800540657330923d00 /
      data wg ( 53) /   4.05862761338354481597420116187988679d00 /
      data wg ( 54) /   4.18401782381424031850607692334503121d00 /
      data wg ( 55) /   4.31420264929613425820084573217987912d00 /
      data wg ( 56) /   4.44962515053655906604982820155377774d00 /
      data wg ( 57) /   4.59078802263617511042959849148929810d00 /
      data wg ( 58) /   4.73826464598929537394753873505838770d00 /
      data wg ( 59) /   4.89271277966692168696886936743283567d00 /
      data wg ( 60) /   5.05489168534039512820572507135175938d00 /
      data wg ( 61) /   5.22568375594272391089278010166022467d00 /
      data wg ( 62) /   5.40612213379727909853323512340717863d00 /
      data wg ( 63) /   5.59742640184041404016553694158980053d00 /
      data wg ( 64) /   5.80104932137643943530626162455394841d00 /
      data wg ( 65) /   6.01873893878099108768015151514026344d00 /
      data wg ( 66) /   6.25262247491437403092934213480091928d00 /
      data wg ( 67) /   6.50532173517668675787482719663696133d00 /
      data wg ( 68) /   6.78011521200777294201287347980059368d00 /
      data wg ( 69) /   7.08117122025414518776174311916759402d00 /
      data wg ( 70) /   7.41389244615305421421695606226687752d00 /
      data wg ( 71) /   7.78544154841612700386232740339230532d00 /
      data wg ( 72) /   8.20557347814596472333905086100917119d00 /
      data wg ( 73) /   8.68801383996161871469419958058255237d00 /
      data wg ( 74) /   9.25286973415578523923556506201979918d00 /
      data wg ( 75) /   9.93114471840215736008370986534009772d00 /
      data wg ( 76) /   1.07739736414646829405750843522990655d01 /
      data wg ( 77) /   1.18738912465097447081950887710877400d01 /
      data wg ( 78) /   1.34228858497264236139734940154089734d01 /
      data wg ( 79) /   1.59197801616897924449554252200185978d01 /
      data wg ( 80) /   2.14214542964372259537521036186415127d01 /


!c           list of major variables
!c           -----------------------
!c           absc   - abscissa
!c           fval*  - function value
!c           result - result of the 20-point gauss formula

      result   = 0.0d0
      dresult  = 0.0d0
      ddresult = 0.0d0
      do j=1,80
       absc = a+b*xg(j)
       call f(absc, par, n, fval, dfval, ddfval)
       result   = result + fval*wg(j)
       dresult  = dresult + dfval*wg(j)
       ddresult = ddresult + ddfval*wg(j)
      enddo
      result   = result*b
      dresult  = dresult*b
      ddresult = ddresult*b
      return
      end subroutine dqlag080

end module functions
