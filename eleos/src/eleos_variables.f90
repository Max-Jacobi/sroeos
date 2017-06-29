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

module eleos_variables

  use Kind_Types_Mod
  use Physical_Constants_Mod

  implicit none

  save

!c..declare the input
      double precision tem,den,zbar,abar
!$omp threadprivate(tem,den,zbar,abar)

!c..declare everything else
!c..totals
      double precision pres,ener,entr, &
                       dpresdd,dpresdt,dpresda,dpresdz, &
                       denerdd,denerdt,denerda,denerdz, &
                       dentrdd,dentrdt,dentrda,dentrdz
!$omp threadprivate(pres,ener,entr)
!$omp threadprivate(dpresdd,dpresdt,dpresda,dpresdz)
!$omp threadprivate(denerdd,denerdt,denerda,denerdz)
!$omp threadprivate(dentrdd,dentrdt,dentrda,dentrdz)

!c..radiation
      integer          radmult
      double precision prad,erad,srad, &
                       dpraddd,dpraddt,dpradda,dpraddz, &
                       deraddd,deraddt,deradda,deraddz, &
                       dsraddd,dsraddt,dsradda,dsraddz

!$omp threadprivate(radmult,prad,erad,srad)
!$omp threadprivate(dpraddd,dpraddt,dpradda,dpraddz)
!$omp threadprivate(deraddd,deraddt,deradda,deraddz)
!$omp threadprivate(dsraddd,dsraddt,dsradda,dsraddz)

!c..ions
      integer          ionmult
      double precision pion,eion,sion,xni,etaion, &
                       dpiondd,dpiondt,dpionda,dpiondz, &
                       deiondd,deiondt,deionda,deiondz, &
                       dsiondd,dsiondt,dsionda,dsiondz, &
                       dxnidd,dxnidt,dxnida,dxnidz, &
                       detaiondd,detaiondt,detaionda,detaiondz
!$omp threadprivate(ionmult,pion,eion,sion,xni,etaion)
!$omp threadprivate(dpiondd,dpiondt,dpionda,dpiondz)
!$omp threadprivate(deiondd,deiondt,deionda,deiondz)
!$omp threadprivate(dsiondd,dsiondt,dsionda,dsiondz)
!$omp threadprivate(dxnidd,dxnidt,dxnida,dxnidz)
!$omp threadprivate(detaiondd,detaiondt,detaionda,detaiondz)

!c..electron-positrons
      integer          elemult
      double precision etaele,detadd,detadt,detada,detadz, &
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
                       sep,dsepdd,dsepdt,dsepda,dsepdz

!$omp threadprivate(elemult,etaele,detadd,detadt,detada,detadz,etapos,zeff)
!$omp threadprivate(xne,dxnedd,dxnedt,dxneda,dxnedz)
!$omp threadprivate(xnefer,dxneferdd,dxneferdt,dxneferda,dxneferdz)
!$omp threadprivate(xnpfer,dxnpferdd,dxnpferdt,dxnpferda,dxnpferdz)
!$omp threadprivate(pele,dpeledd,dpeledt,dpeleda,dpeledz)
!$omp threadprivate(ppos,dpposdd,dpposdt,dpposda,dpposdz)
!$omp threadprivate(pep,dpepdd,dpepdt,dpepda,dpepdz)
!$omp threadprivate(eele,deeledd,deeledt,deeleda,deeledz)
!$omp threadprivate(epos,deposdd,deposdt,deposda,deposdz)
!$omp threadprivate(eep,deepdd,deepdt,deepda,deepdz)
!$omp threadprivate(sele,dseledd,dseledt,dseleda,dseledz)
!$omp threadprivate(spos,dsposdd,dsposdt,dsposda,dsposdz)
!$omp threadprivate(sep,dsepdd,dsepdt,dsepda,dsepdz)

!c..ionization contributions
      integer          ionized,potmult
      double precision eip,deipdd,deipdt,deipda,deipdz, &
                       sip,dsipdd,dsipdt,dsipda,dsipdz
!$omp threadprivate(ionized,potmult)
!$omp threadprivate(eip,deipdd,deipdt,deipda,deipdz)
!$omp threadprivate(sip,dsipdd,dsipdt,dsipda,dsipdz)

!c..coulomb corrections
      integer          coulmult
      double precision plasg, &
                       pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz, &
                       ecoul,decouldd,decouldt,decoulda,decouldz, &
                       scoul,dscouldd,dscouldt,dscoulda,dscouldz
!$omp threadprivate(coulmult,plasg)
!$omp threadprivate(pcoul,dpcouldd,dpcouldt,dpcoulda,dpcouldz)
!$omp threadprivate(ecoul,decouldd,decouldt,decoulda,decouldz)
!$omp threadprivate(scoul,dscouldd,dscouldt,dscoulda,dscouldz)

!c..various physical quantities based on derivatives
      double precision chit,chid,cv,cp,gam1,gam2,gam3,nabad,sound
!$omp threadprivate(chit,chid,cv,cp,gam1,gam2,gam3,nabad,sound)

!c..for the maxwell relations
      double precision dse,dpe,dsp
!$omp threadprivate(dse,dpe,dsp)

end module eleos_variables
