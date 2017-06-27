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
module eleos

 use functions
 use Physical_Constants_Mod
 use eleos_variables

 implicit none

 contains

      subroutine azbar(xmass,aion,zion,ionmax,ymass,abar,zbar)
      include 'implno.dek'



!c..this routine calculates composition variables for an eos routine

!c..input:
!c..mass fractions     = xmass(1:ionmax)
!c..number of nucleons = aion(1:ionmax)
!c..charge of nucleus  = zion(1:ionmax)
!c..number of isotopes = ionmax

!c..output:
!c..molar abundances        = ymass(1:ionmax),
!c..mean number of nucleons = abar
!c..mean nucleon charge     = zbar


!c..declare
      integer          i,ionmax
      double precision xmass(ionmax),aion(ionmax),zion(ionmax), &
                       ymass(ionmax),abar,zbar,zbarxx,ytot1

      zbarxx  = 0.0d0
      ytot1   = 0.0d0
      do i=1,ionmax
       ymass(i) = xmass(i)/aion(i)
       ytot1    = ytot1 + ymass(i)
       zbarxx   = zbarxx + zion(i) * ymass(i)
      enddo
      abar   = 1.0d0/ytot1
      zbar   = zbarxx * abar
      return
      end subroutine azbar


      subroutine pretty_eos_out(whose)
      include 'implno.dek'
      include 'vector_eos.dek'


!c..writes a pretty output for the eos tester
!c..
!c..declare
      integer     j
      character*7 whose

!c..popular formats
01    format(1x,t2,a,t11,'total',t24,'ion',t34,'electron', &
             t46,'positron',t58,'radiation',t70,'coulomb')
02    format(1x,t2,a,1p6e12.4)
03    format(1x,t2,a6,1pe12.4,t22,a6,1pe12.4, &
               t42,a6,1pe12.4,t62,a6,1pe12.4)



      do j=jlo_eos,jhi_eos


!c..the input
      write(6,03) 'temp =',temp_row(1),'den  =',den_row(1), &
                  'abar =',abar_row(1),'zbar =',zbar_row(1)
      write(6,*) ' '


!c..and the output
!c..first the totals from each of the components
      write(6,01)  whose
      write(6,02) 'pres =', &
                  ptot_row(j),pion_row(j),pele_row(j), &
                  ppos_row(j),prad_row(j),pcou_row(j)
      write(6,02) 'ener =', &
                  etot_row(j),eion_row(j),eele_row(j), &
                  epos_row(j),erad_row(j),ecou_row(j)
      write(6,02) 'entr =', &
                  stot_row(j),sion_row(j),sele_row(j), &
                  spos_row(j),srad_row(j),scou_row(j)


!c..derivatives of the totals with respect to the input variables
      write(6,*)  ' '
      write(6,03) 'dp/dd=',dpd_row(j),'dp/dt=',dpt_row(j), &
                  'dp/da=',dpa_row(j),'dp/dz=',dpz_row(j)
      write(6,03) 'de/dd=',ded_row(j),'de/dt=',det_row(j), &
                  'de/da=',dea_row(j),'de/dz=',dez_row(j)
      write(6,03) 'ds/dd=',dsd_row(j),'ds/dt=',dst_row(j), &
                  'ds/da=',dsa_row(j),'ds/dz=',dsz_row(j)


!c..derivatives of the electron-positron compoenets with
!c..respect to the input variables
      write(6,*) ' '
      write(6,03) 'dpepd=',dpepd_row(j),'dpept=',dpept_row(j), &
                  'dpepa=',dpepa_row(j),'dpepz=',dpepz_row(j)
      write(6,03) 'deepd=',deepd_row(j),'deept=',deept_row(j), &
                  'deepa=',deepa_row(j),'deepz=',deepz_row(j)
      write(6,03) 'dsepd=',dsepd_row(j),'dsept=',dsept_row(j), &
                  'dsepa=',dsepa_row(j),'dsepz=',dsepz_row(j)


!c..the thermodynamic consistency relations, these should all be
!c..at the floating poiint limit of zero
      write(6,*) ' '
      write(6,03) 'maxw1=',dse_row(j),'maxw2=',dpe_row(j), &
                  'maxw3=',dsp_row(j)


!c..number density of electrons, poistrons, matter electrons, and ions
      write(6,03) 'xne  =',xne_row(j),'xnp  =',xnp_row(j), &
                  'xnem =',xnem_row(j),'xni  =',xni_row(j)


!c..derivatibves of the electron number density with
!c..respect to the input variables
      write(6,03) 'dxned=',dxned_row(j),'dxnet=',dxnet_row(j), &
                  'dxnea=',dxnea_row(j),'dxnez=',dxnez_row(j)


!c..electron chemical potential, positron chemical potential
!c..and derivatives of electron chemical potential with respect
!c..to the input variables
      write(6,03) 'eta  =',etaele_row(j),'etap =',etapos_row(j)
      write(6,03) 'detad=',detad_row(j),'detat=',detat_row(j), &
                  'detaa=',detaa_row(j),'detaz=',detaz_row(j)


!c..specific heats, and ratio of electostatic to thermal energy
      write(6,03) 'cp   =',cp_row(j),'cv   =',cv_row(j), &
                  'plasg=',plasg_row(j)

!c..the 3 gammas and the sound speed
      write(6,03) 'gam1 =',gam1_row(j),'gam2 =',gam2_row(j), &
                  'gam3 =',gam3_row(j),'csond=',cs_row(j)
      write(6,*) ' '

      enddo
      return
      end subroutine pretty_eos_out







!c..routine eosfxt computes a stellar eos assuming complete ionozation
!c..routine etages makes a good guess for the chemical potential
!c..routine xneroot gets the thermodynamics of the electrons and positrons


      subroutine eosfxt
      include 'implno.dek'
      include 'const.dek'
      include 'vector_eos.dek'

!c..given a temperature temp [K], density den [g/cm**3], and a composition
!c..characterized by abar (average weight) and zbar (average charge),
!c..this routine returns all the other thermodynamic quantities.

!c..of interest is the pressure [erg/cm**3], specific thermal energy [erg/gr],
!c..the entropy [erg/g/K], with their derivatives with respect to temperature,
!c..density, abar, and zbar.

!c..other quantites such the normalized chemical potential eta (plus its
!c..derivatives), number density of electrons and positron pair (along
!c..with their derivatives), adiabatic indices, specific heats, and
!c..relativistically correct sound speed are also returned.

!c..this routine assumes planckian photons, an ideal gas of ions,
!c..and an electron-positron gas with an arbitrary degree of relativity
!c..and degeneracy. the full fermi-dirac integrals and their derivatives
!c..with respect to eta and beta are computed to machine precision, and
!c..all other derivatives are analytic.

!c..references: cox & giuli (c&g) chapter 24,
!c..            timmes & arnett, apj supp. 125, 277, 1999
!c..            timmes & swesty, apj supp. 126, 501, 2000




!c..a dictionary of terms used:
!c..this routine has now been pipelined.
!c..all the input and output variables are in the file vector_eos.dek.
!c..the vector name is the scaler name appended with an "_row",
!c..for example, temp_row(i), den_row(i), and so on.



!c..input:
!c..temp     = temperature
!c..den      = density
!c..abar     = average number of nucleons per nuclei
!c..zbar     = average number of protons per nuclei


!c..output:

!c..pres     = total pressure
!c..dpresdd  = derivative of total pressure with respect to density
!c..dpresdt  = derivative of total pressure with respect to temperature
!c..dpresda  = derivative of total pressure with respect to abar
!c..dpresdz  = derivative of total pressure with respect to zbar

!c..ener     = total internal energy
!c..denerdd  = derivative of total energy with respect to density
!c..denerdt  = derivative of total energy with respect to temperature
!c..denerda  = derivative of total energy with respect to abar
!c..denerdz  = derivative of total energy with respect to zbar

!c..entr     = total entropy
!c..dentrdd  = derivative of total entropy with respect to density
!c..dentrdt  = derivative of total entropy with respect to temperature
!c..dentrda  = derivative of total entropy with respect to abar
!c..dentrdz  = derivative of total entropy with respect to zbar



!c..prad     = radiation pressure
!c..dpraddd  = derivative of the radiation pressure with density
!c..dpraddt  = derivative of the radiation pressure with temperature
!c..dpradda  = derivative of the radiation pressure with abar
!c..dpraddz  = derivative of the radiation pressure with zbar

!c..erad     = radiation energy
!c..deraddd  = derivative of the radiation energy with density
!c..deraddt  = derivative of the radiation energy with temperature
!c..deradda  = derivative of the radiation energy with abar
!c..deraddz  = derivative of the radiation energy with zbar

!c..srad     = radiation entropy
!c..dsraddd  = derivative of the radiation entropy with density
!c..dsraddt  = derivative of the radiation entropy with temperature
!c..dsradda  = derivative of the radiation entropy with abar
!c..dsraddz  = derivative of the radiation entropy with zbar

!c..radmult  = radiation multiplier (useful for turning radiation off/on)




!c..xni      = number density of ions
!c..dxnidd   = derivative of the ion number density with density
!c..dxnidt   = derivative of the ion number density with temperature
!c..dxnida   = derivative of the ion number density with abar
!c..dxnidz   = derivative of the ion number density with zbar

!c..pion     = ion pressure
!c..dpiondd  = derivative of the ion pressure with density
!c..dpiondt  = derivative of the ion pressure with temperature
!c..dpionda  = derivative of the ion pressure with abar
!c..dpiondz  = derivative of the ion pressure with zbar

!c..eion     = ion energy
!c..deiondd  = derivative of the ion energy with density
!c..deiondt  = derivative of the ion energy with temperature
!c..deionda  = derivative of the ion energy with abar
!c..deiondz  = derivative of the ion energy with zbar

!c..sion     = ion entropy
!c..dsiondd  = derivative of the ion entropy with density
!c..dsiondt  = derivative of the ion entropy with temperature
!c..dsionda  = derivative of the ion entropy with abar
!c..dsiondz  = derivative of the ion entropy with zbar

!c..ionmult  = ion multiplier (useful for turning ions off/on)


!c..etaele   = electron chemical potential
!c..detadd   = derivative of the electron chem potential with density
!c..detadt   = derivative of the electron chem potential with temperature
!c..detada   = derivative of the electron chem potential with abar
!c..detadz   = derivative of the electron chem potential with zbar

!c..etapos   = positron degeneracy parameter




!c..xne       = number density of electrons
!c..dxnedd    = derivative of the electron number density with density
!c..dxnedt    = derivative of the electron number density with temperature
!c..dxneda    = derivative of the electron number density with abar
!c..dxnedz    = derivative of the electron number density with zbar

!c..xnefer    = fermi integral electron number density
!c..dxneferdd = derivative of the fermi electron number density with density
!c..dxneferdt = derivative of the fermi electron number density with temperature
!c..dxneferda = derivative of the fermi electron number density with abar
!c..dxneferdz = derivative of the fermi electron number density with zbar

!c..xnpfer    = fermi integral positron number density
!c..dxnpferdd = derivative of the fermi positron number density with density
!c..dxnpferdt = derivative of the fermi positron number density with temperature
!c..dxnpferda = derivative of the fermi positron number density with abar
!c..dxnpferdz = derivative of the fermi positron number density with zbar

!c..pele      = electron pressure
!c..dpeledd   = derivative of the electron pressure with density
!c..dpeledt   = derivative of the electron pressure with temperature
!c..dpeleda   = derivative of the electron pressure with abar
!c..dpeledz   = derivative of the electron pressure with zbar

!c..eele     = electron energy
!c..deeledd   = derivative of the electron energy with density
!c..deeledt   = derivative of the electron energy with temperature
!c..deeleda   = derivative of the electron energy with abar
!c..deeledz   = derivative of the electron energy with zbar

!c..sele     = electron entropy
!c..dseledd   = derivative of the electron entropy with density
!c..dseledt   = derivative of the electron entropy with temperature
!c..dseleda   = derivative of the electron entropy with abar
!c..dseledz   = derivative of the electron entropy with zbar


!c..ppos     = positron pressure
!c..dpposdd   = derivative of the positron pressure with density
!c..dpposdt   = derivative of the positron pressure with temperature
!c..dpposda   = derivative of the positron pressure with abar
!c..dpposdz   = derivative of the positron pressure with zbar

!c..epos     = electron energy
!c..deposdd   = derivative of the positron energy with density
!c..deposdt   = derivative of the positron energy with temperature
!c..deposda   = derivative of the positron energy with abar
!c..deposdz   = derivative of the positron energy with zbar

!c..spos     = electron entropy
!c..dsposdd   = derivative of the positron entropy with density
!c..dsposdt   = derivative of the positron entropy with temperature
!c..dsposda   = derivative of the positron entropy with abar
!c..dsposdz   = derivative of the positron entropy with zbar

!c..pep      = electron + positron pressure
!c..dpepdd   = derivative of the electron+positron pressure with density
!c..dpepdt   = derivative of the electron+positron pressure with temperature
!c..dpepda   = derivative of the electron+positron pressure with abar
!c..dpepdz   = derivative of the electron+positron pressure with zbar

!c..eep      = electron + ositron energy
!c..deepdd   = derivative of the electron+positron energy with density
!c..deepdt   = derivative of the electron+positron energy with temperature
!c..deepda   = derivative of the electron+positron energy with abar
!c..deepdz   = derivative of the electron+positron energy with zbar

!c..sep      = electron + positron entropy
!c..dsepdd   = derivative of the electron+positron entropy with density
!c..dsepdt   = derivative of the electron+positron entropy with temperature
!c..dsepda   = derivative of the electron+positron entropy with abar
!c..dsepdz   = derivative of the electron+positron entropy with zbar

!c..elemult  = electron multiplier (useful for turning e-e+ off/on)


!c..eip      = ionization potential ennergy
!c..deipdd   = derivative of ionization energy with density
!c..deipdt   = derivative of ionization energy with temperature
!c..deipda   = derivative of ionization energy with abar
!c..deipdz   = derivative of ionization energy with zbar


!c..sip      = ionization potential ennergy
!c..dsipdd   = derivative of ionization energy with density
!c..dsipdt   = derivative of ionization energy with temperature
!c..dsipda   = derivative of ionization energy with abar
!c..dsipdz   = derivative of ionization energy with zbar

!c..potmult  = ionization energy multiplier (useful for turning off ionization additions)



!c..pcoul    = coulomb pressure correction
!c..coulmult = coulomb component multiplier
!c..dpcouldd = derivative of the coulomb pressure with density
!c..dpcouldt = derivative of the coulomb pressure with temperature
!c..dpcoulda = derivative of the coulomb pressure with abar
!c..dpcouldz = derivative of the coulomb pressure with zbar

!c..ecoul    = coulomb energy correction
!c..decouldd = derivative of the coulomb energy with density
!c..decouldt = derivative of the coulomb energy with temperature
!c..decoulda = derivative of the coulomb energy with abar
!c..decouldz = derivative of the coulomb energy with zbar

!c..scoul    = coulomb entropy correction
!c..dscouldd = derivative of the coulomb entropy with density
!c..dscouldt = derivative of the coulomb entropy with temperature
!c..dscoulda = derivative of the coulomb entropy with abar
!c..dscouldz = derivative of the coulomb entropy with zbar


!c..kt       = kerg * temperature
!c..beta     = dimensionless ratio of kerg*temp/me*c^2

!c..chit     = temperature exponent in the pressure equation of state
!c..chid     = density exponent in the pressure equation of state
!c..cv       = specific heat at constant volume
!c..cp       = specific heat at constant pressure
!c..gam1     = first adiabatic exponent
!c..gam2     = second adiabatic exponent
!c..gam3     = third adiabatic exponent
!c..nabad    = adiabatic gradient
!c..sound    = relativistically correct adiabatic sound speed
!c..plasg    = ratio of electrostatic to thermal energy


!c..dse      = thermodynamic consistency check de/dt = t*ds/dt
!c..dpe      = thermodynamic consistency check p = d**2 de/dd + t*dpdt
!c..dsp      = thermodynamic consistency check dp/dt = - d**2 ds/dd

      double precision temp

!c..miscelaneous local variables
      integer          i,j,niter,mode
      double precision kt,ktinv,x,y,z,xx,yy,zz,ages,agesav,agesnew, &
                       ratio,ytot1,f,df,deninv,tempinv


!c..various derived constants
      double precision third,sifac,eostol,fpmin
      parameter        (third  = 1.0d0/3.0d0, &
                        sifac  = 8.6322745944370191d-45, &
                        eostol = 1.0d-13, &
                        fpmin  = 1.0d-14)

!c..note: sifac = h**3/(2.0d0*pi*amu)**1.5d0

      temp = tem

!c..popular format statements for debugging
01    format(1x,5(a,1pe24.16))
02    format(1x,5(a,1pe16.8))
03    format(1x,1p5e16.8)



!c..set the on/off switches
      radmult  = 1
      ionmult  = 0
      ionized  = 0
      elemult  = 1
      coulmult = 0
      potmult  = 0



!c..start pipeline loop
      do j=jlo_eos,jhi_eos

       if (temp_row(j) .le. 0.0) stop 'temp less than 0 in eosfxt'
       if (den_row(j)  .le. 0.0) stop 'den less than 0 in eosfxt'

       temp  = temp_row(j)
       den   = den_row(j)
       abar  = abar_row(j)
       zbar  = zbar_row(j)
       ytot1 = 1.0d0/abar


!c..initialize
       prad     = 0.0d0
       dpraddd  = 0.0d0
       dpraddt  = 0.0d0
       dpradda  = 0.0d0
       dpraddz  = 0.0d0

       erad     = 0.0d0
       deraddd  = 0.0d0
       deraddt  = 0.0d0
       deradda  = 0.0d0
       deraddz  = 0.0d0

       srad     = 0.0d0
       dsraddd  = 0.0d0
       dsraddt  = 0.0d0
       dsradda  = 0.0d0
       dsraddz  = 0.0d0

       xni      = 0.0d0
       dxnidd   = 0.0d0
       dxnidt   = 0.0d0
       dxnida   = 0.0d0
       dxnidz   = 0.0d0

       pion     = 0.0d0
       dpiondd  = 0.0d0
       dpiondt  = 0.0d0
       dpionda  = 0.0d0
       dpiondz  = 0.0d0

       eion     = 0.0d0
       deiondd  = 0.0d0
       deiondt  = 0.0d0
       deionda  = 0.0d0
       deiondz  = 0.0d0

       sion     = 0.0d0
       dsiondd  = 0.0d0
       dsiondt  = 0.0d0
       dsionda  = 0.0d0
       dsiondz  = 0.0d0

       xne      = 0.0d0
       dxnedd   = 0.0d0
       dxnedt   = 0.0d0
       dxneda   = 0.0d0
       dxnedz   = 0.0d0

       etaele   = 0.0d0
       detadd   = 0.0d0
       detadt   = 0.0d0
       detada   = 0.0d0
       detadz   = 0.0d0
       etapos   = 0.0d0

       xnefer    = 0.0d0
       dxneferdd = 0.0d0
       dxneferdt = 0.0d0
       dxneferda = 0.0d0
       dxneferdz = 0.0d0

       xnpfer    = 0.0d0
       dxnpferdd = 0.0d0
       dxnpferdt = 0.0d0
       dxnpferda = 0.0d0
       dxnpferdz = 0.0d0

       pele     = 0.0d0
       dpeledd  = 0.0d0
       dpeledt  = 0.0d0
       dpeleda  = 0.0d0
       dpeledz  = 0.0d0

       eele     = 0.0d0
       deeledd  = 0.0d0
       deeledt  = 0.0d0
       deeleda  = 0.0d0
       deeledz  = 0.0d0

       sele     = 0.0d0
       dseledd  = 0.0d0
       dseledt  = 0.0d0
       dseleda  = 0.0d0
       dseledz  = 0.0d0

       ppos     = 0.0d0
       dpposdd  = 0.0d0
       dpposdt  = 0.0d0
       dpposda  = 0.0d0
       dpeledz  = 0.0d0

       epos     = 0.0d0
       deposdd  = 0.0d0
       deposdt  = 0.0d0
       deposda  = 0.0d0
       deeledz  = 0.0d0

       spos     = 0.0d0
       dsposdd  = 0.0d0
       dsposdt  = 0.0d0
       dsposda  = 0.0d0
       dseledz  = 0.0d0

       pep      = 0.0d0
       dpepdd   = 0.0d0
       dpepdt   = 0.0d0
       dpepda   = 0.0d0
       dpepdz   = 0.0d0

       eep      = 0.0d0
       deepdd   = 0.0d0
       deepdt   = 0.0d0
       deepda   = 0.0d0
       deepdz   = 0.0d0

       sep      = 0.0d0
       dsepdd   = 0.0d0
       dsepdt   = 0.0d0
       dsepda   = 0.0d0
       dsepdz   = 0.0d0

       eip      = 0.0d0
       deipdd   = 0.0d0
       deipdt   = 0.0d0
       deipda   = 0.0d0
       deipdz   = 0.0d0

       sip      = 0.0d0
       dsipdd   = 0.0d0
       dsipdt   = 0.0d0
       dsipda   = 0.0d0
       dsipdz   = 0.0d0

       pcoul    = 0.0d0
       dpcouldd = 0.0d0
       dpcouldt = 0.0d0
       dpcoulda = 0.0d0
       dpcouldz = 0.0d0

       ecoul    = 0.0d0
       decouldd = 0.0d0
       decouldt = 0.0d0
       decoulda = 0.0d0
       decouldz = 0.0d0

       scoul    = 0.0d0
       dscouldd = 0.0d0
       dscouldt = 0.0d0
       dscoulda = 0.0d0
       dscouldz = 0.0d0


       kt      = kerg * temp
       ktinv   = 1.0d0/kt
       deninv  = 1.0d0/den
       tempinv = 1.0d0/temp




!c..radiation section:
       if (radmult .ne. 0) then

!c..pressure in erg/cm**3
        prad    = asol * third * temp * temp * temp * temp
        dpraddd = 0.0d0
        dpraddt = 4.0d0 * prad/temp
        dpradda = 0.0d0
        dpraddz = 0.0d0

!c..energy in erg/gr
        erad    = 3.0d0 * prad * deninv
        deraddd = -erad * deninv
        deraddt = 3.0d0 * dpraddt * deninv
        deradda = 0.0d0
        deraddz = 0.0d0

!c..entropy in erg/g/kelvin
        srad    = (prad*deninv + erad) * tempinv
        dsraddd = (dpraddd*deninv - prad*deninv**2 + deraddd) * tempinv
        dsraddt = (dpraddt*deninv + deraddt - srad)  * tempinv
        dsradda = 0.0d0
        dsraddz = 0.0d0
       end if




!c..ion section:

!c..number density in 1/cm**3,
        xni     = avo * ytot1 * den
        dxnidd  = avo * ytot1
        dxnidt  = 0.0d0
        dxnida  = -xni * ytot1
        dxnidz  = 0.0d0

       if (ionmult .ne. 0) then

!c..pressure in erg/cm**3
        pion    = xni * kt
        dpiondd = dxnidd * kt
        dpiondt = xni * kerg
        dpionda = -pion * ytot1
        dpiondz = 0.0d0

!c.. energy in erg/gr
        eion    = 1.5d0 * pion*deninv
        deiondd = (1.5d0 * dpiondd - eion)*deninv
        deiondt = 1.5d0 * dpiondt*deninv
        deionda = 1.5d0 * dpionda*deninv
        deiondz = 0.0d0


!c..ion degeneracy parameter (c&g 9.60)
        y         = 1.0d0/(abar*kt)
        yy        = y * sqrt(y)
        z         = xni * sifac * yy
        etaion    = log(z)
        xx        = 1.0d0/xni
        detaiondd = dxnidd*xx
        detaiondt = dxnidt*xx - 1.5d0*tempinv
        detaionda = dxnida*xx - 1.5d0*ytot1
        detaiondz = dxnidz*xx


!c..entropy in erg/gr/kelvin
!c..the last term is the usual  etaion * kerg * xni/den
!c..sometimes called the sacker-tetrode equation

        sion    = (eion + pion*deninv)*tempinv - etaion * kerg*avo*ytot1

        dsiondd = (deiondd + dpiondd*deninv - pion*deninv**2)*tempinv  &
                  - detaiondd * kerg * avo*ytot1

        dsiondt = (deiondt + dpiondt*deninv)*tempinv  &
                  - (eion + pion*deninv)*tempinv**2 &
                  - detaiondt * kerg * avo*ytot1

        dsionda = (deionda + dpionda*deninv)*tempinv &
                  - detaionda * kerg * avo*ytot1 &
                  + etaion * kerg * avo * ytot1**2

        dsiondz = 0.0d0
       end if









!c..electron-positron section:
       if (elemult .ne. 0) then

!c..make a good guess at the the electron degeneracy parameter eta
        call etages(xni,zbar,temp,ages)
        agesav = ages


!c..newton-raphson to get the electron/positron quantities
        eosfail = .false.
        mode   = 0
        do i=1,100

         call xneroot(mode,den,temp,abar,zbar,ionized, &
                      ages,f,df, &
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


         if (df .eq. 0.0) goto 11
         ratio   = f/df
         agesnew = ages - ratio
         z       = abs((agesnew - ages)/ages)
         ages    = agesnew
         niter   = i
         if (z .lt. eostol .or. abs(ratio) .le. fpmin) goto 20
        enddo

11      write(6,*)
        write(6,*) 'newton-raphson failed in routine eosfxt'
        write(6,01) 'temp  =',temp/temp_mev_to_kelvin,' den =',den*rho_cgs_to_EOS
        write(6,01) 'z     =',z,' ages=',ages, ' agesav=',agesav
        write(6,01) 'eostol=',eostol
        write(6,01) 'f/df  =',f/df,' f   =',f,    ' df    =',df
        write(6,01) 'fpmin =',fpmin
        write(6,*)
        call flush(6)
        eosfail = .true.
        return
20      continue



!c..with the converged values, get the energy, pressure, and entropy
         mode = 1
         call xneroot(mode,den,temp,abar,zbar,ionized, &
                      ages,f,df, &
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

       end if





!c..coulomb corretions section:
       if (coulmult .ne. 0) then

        call coulomb(den,temp,abar,zbar, &
                     pion,dpiondd,dpiondt,dpionda,dpiondz, &
                     xne,dxnedd,dxnedt,dxneda,dxnedz, &
                     plasg, &
                     pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                     ecoul,decouldd,decouldt,decoulda,decouldz, &
                     scoul,dscouldd,dscouldt,dscoulda,dscouldz)


!c..bomb proof the cooulomb correctins
        x   = prad + pion + pele + ppos + pcoul
        if (x .le. 0.0) then

!c         write(6,*)
!c         write(6,*) 'coulomb corrections are causing a negative pressure'
!c         write(6,*) 'setting all coulomb corrections to zero'
!c         write(6,*)

         pcoul    = 0.0d0
         dpcouldd = 0.0d0
         dpcouldt = 0.0d0
         dpcoulda = 0.0d0
         dpcouldz = 0.0d0
         ecoul    = 0.0d0
         decouldd = 0.0d0
         decouldt = 0.0d0
         decoulda = 0.0d0
         decouldz = 0.0d0
         scoul    = 0.0d0
         dscouldd = 0.0d0
         dscouldt = 0.0d0
         dscoulda = 0.0d0
         dscouldz = 0.0d0
        end if
       end if



!c..sum all the components
       pres    = prad + 0.0d0 + pele + ppos + 0.0d0
       ener    = erad + 0.0d0 + eele + epos + 0.0d0 + 0.0d0
       entr    = srad + 0.0d0 + sele + spos + 0.0d0 + 0.0d0

!       pres    = prad + pion + pele + ppos + pcoul
!       ener    = erad + eion + eele + epos + ecoul + eip
!       entr    = srad + sion + sele + spos + scoul + sip

       dpresdd = dpraddd + 0.0d0 + dpepdd + 0.0d0
       dpresdt = dpraddt + 0.0d0 + dpepdt + 0.0d0
       dpresda = dpradda + 0.0d0 + dpepda + 0.0d0
       dpresdz = dpraddz + 0.0d0 + dpepdz + 0.0d0
       denerdd = deraddd + 0.0d0 + deepdd + 0.0d0 + 0.0d0
       denerdt = deraddt + 0.0d0 + deepdt + 0.0d0 + 0.0d0
       denerda = deradda + 0.0d0 + deepda + 0.0d0 + 0.0d0
       denerdz = deraddz + 0.0d0 + deepdz + 0.0d0 + 0.0d0
       dentrdd = dsraddd + 0.0d0 + dsepdd + 0.0d0 + 0.0d0
       dentrdt = dsraddt + 0.0d0 + dsepdt + 0.0d0 + 0.0d0
       dentrda = dsradda + 0.0d0 + dsepda + 0.0d0 + 0.0d0
       dentrdz = dsraddz + 0.0d0 + dsepdz + 0.0d0 + 0.0d0


!ccc       dpresdd = dpraddd + dpiondd + dpepdd + dpcouldd
!ccc       dpresdt = dpraddt + dpiondt + dpepdt + dpcouldt
!ccc       dpresda = dpradda + dpionda + dpepda + dpcoulda
!ccc       dpresdz = dpraddz + dpiondz + dpepdz + dpcouldz
!ccc
!ccc       denerdd = deraddd + deiondd + deepdd + decouldd + deipdd
!ccc       denerdt = deraddt + deiondt + deepdt + decouldt + deipdt
!ccc       denerda = deradda + deionda + deepda + decoulda + deipda
!ccc       denerdz = deraddz + deiondz + deepdz + decouldz + deipdz
!ccc
!ccc       dentrdd = dsraddd + dsiondd + dsepdd + dscouldd + dsipdd
!ccc       dentrdt = dsraddt + dsiondt + dsepdt + dscouldt + dsipdt
!ccc       dentrda = dsradda + dsionda + dsepda + dscoulda + dsipda
!ccc       dentrdz = dsraddz + dsiondz + dsepdz + dscouldz + dsipdz
!ccc




!c..the temperature and density exponents (c&g 9.81 9.82)
!c..the specific heat at constant volume (c&g 9.92)
!c..the third adiabatic exponent (c&g 9.93)
!c..the first adiabatic exponent (c&g 9.97)
!c..the second adiabatic exponent (c&g 9.105)
!c..the specific heat at constant pressure (c&g 9.98)
!c..and relativistic formula for the sound speed (c&g 14.29)

       zz    = pres/den
       chit  = temp/pres * dpresdt
       chid  = dpresdd/zz
       cv    = denerdt
       x     = zz * chit/(temp * cv)
       gam3  = x + 1.0d0
       gam1  = chit*x + chid
       nabad = x/gam1
       gam2  = 1.0d0/(1.0d0 - nabad)
       cp    = cv * gam1/chid
       z     = 1.0d0 + (ener + clight*clight)/zz
       sound = clight * sqrt(gam1/z)

! gamma should also be:
!       gam1 = dpresdd/zz - dentrdd*den * chit / (dentrdt*temp)


!c..maxwell relations; each is zero if the consistency is perfect
!c..delicate subtraction in very degenerate regions causes roundoff error

       dse = temp*dentrdt/denerdt - 1.0d0

       dpe = (denerdd*den**2 + temp*dpresdt)/pres - 1.0d0

       dsp = -dentrdd*(den**2/dpresdt) - 1.0d0




!c..store this row
        ptot_row(j)   = pres
        dpt_row(j)    = dpresdt
        dpd_row(j)    = dpresdd
        dpa_row(j)    = dpresda
        dpz_row(j)    = dpresdz

        etot_row(j)   = ener
        det_row(j)    = denerdt
        ded_row(j)    = denerdd
        dea_row(j)    = denerda
        dez_row(j)    = denerdz

        stot_row(j)   = entr
        dst_row(j)    = dentrdt
        dsd_row(j)    = dentrdd
        dsa_row(j)    = dentrda
        dsz_row(j)    = dentrdz

        prad_row(j)   = prad
        erad_row(j)   = erad
        srad_row(j)   = srad

        pion_row(j)   = pion
        eion_row(j)   = eion
        sion_row(j)   = sion

        xni_row(j)    = xni

        pele_row(j)   = pele
        ppos_row(j)   = ppos

        dpept_row(j)  = dpepdt
        dpepd_row(j)  = dpepdd
        dpepa_row(j)  = dpepda
        dpepz_row(j)  = dpepdz

        eele_row(j)   = eele
        epos_row(j)   = epos
        deept_row(j)  = deepdt
        deepd_row(j)  = deepdd
        deepa_row(j)  = deepda
        deepz_row(j)  = deepdz

        sele_row(j)   = sele
        spos_row(j)   = spos
        dsept_row(j)  = dsepdt
        dsepd_row(j)  = dsepdd
        dsepa_row(j)  = dsepda
        dsepz_row(j)  = dsepdz

        xnem_row(j)   = xne
        xne_row(j)    = xnefer
        dxnet_row(j)  = dxneferdt + dxnpferdt
        dxned_row(j)  = dxneferdd + dxnpferdd
        dxnea_row(j)  = dxneferda + dxnpferda
        dxnez_row(j)  = dxneferdz + dxnpferdz
        xnp_row(j)    = xnpfer
        zeff_row(j)   = zeff

        etaele_row(j) = etaele
        detat_row(j)  = detadt
        detad_row(j)  = detadd
        detaa_row(j)  = detada
        detaz_row(j)  = detadz
        etapos_row(j) = etapos

        eip_row(j)    = eip
        sip_row(j)    = sip

        pcou_row(j)   = pcoul
        ecou_row(j)   = ecoul
        scou_row(j)   = scoul
        plasg_row(j)  = plasg

        dse_row(j)    = dse
        dpe_row(j)    = dpe
        dsp_row(j)    = dsp

        cv_row(j)     = cv
        cp_row(j)     = cp
        gam1_row(j)   = gam1
        gam2_row(j)   = gam2
        gam3_row(j)   = gam3
        cs_row(j)     = sound


!c..for debugging
!c        crap1_row(j)   = etaele
!c        dcrap1d_row(j) = detadd
!c        dcrap1t_row(j) = detadt
!c        dcrap1a_row(j) = detada
!c        dcrap1z_row(j) = detadz


!c..end of pipeline loop
      enddo
      return
      end subroutine eosfxt

end module eleos
