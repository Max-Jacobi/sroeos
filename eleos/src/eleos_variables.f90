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

!---------------------- constants

! !from NIST, also defined in lseos.2.7.f for EL_EOS
!    real(sp),parameter :: rho_cgs_to_EOS = 5.9704087d-16 !neutron mass
!    real(sp),parameter :: kb_erg = 1.3806504d-16
!    real(sp),parameter :: kb_mev = 8.61738568d-11
!    real(sp),parameter :: mev_to_cgs = 1.60217649d-6
!    real(sp),parameter :: erg_to_mev = 6.24150636d5
!    real(sp),parameter :: temp_mev_to_kelvin = 1.1604505d10
!    real(sp),parameter :: press_EOS_to_cgs = 1.60217649d33
!    real(sp),parameter :: energy_EOS_to_cgs = 1.60217649d-6/1.674927211d-24
!    real(sp),parameter :: entropy_EOS_to_cgs = 1.3806504d-16/1.674927211d-24
!    real(sp),parameter :: press_cgs_to_EOS  = 6.24150964d-34
!    real(sp),parameter :: energy_cgs_to_EOS = 1.674927211d-24/1.60217649d-6
!    real(sp),parameter :: entropy_cgs_to_EOS = 1.674927211d-24/1.3806504d-16
!    real(sp),parameter :: clite = 29979245900.0d0
!    real(sp),parameter :: planck_cgs = 6.626176d-27
!    real(sp),parameter :: pi1 = 3.14159265358979d0
!    real(sp),parameter :: avo = 5.9704082443567622d23 ! 1/m_n_cgs
!    real(sp),parameter :: sac_const = 2.0d0 * pi1 * (197.327053d0)**2
! !   real(sp),parameter :: avo = 6.0221367d23

!    real(sp),parameter :: neutron_mass_MeV = 939.565560d0 !neutron mass
!    real(sp),parameter :: proton_mass_MeV = 938.272d0  ! proton mass
!    real(sp),parameter :: m_n_cgs = neutron_mass_MeV*mev_to_cgs/clite**2
!    real(sp),parameter :: m_p_cgs = proton_mass_MeV*mev_to_cgs/clite**2
!    real(sp),parameter :: HBC = 197.327d0
!    real(sp),parameter :: alpha_coul = 1.d0/137.035999074d0

!    ! the following constants are used in Luke's mixing script
!    !  real*8,parameter :: e_nuc_to_cgs = 1.0d0/1.03654d-18 amu based
!    !  m_neutron based
!   real(sp),parameter :: e_nuc_to_cgs = 1.0d0/1.04541174787614d-18
!   real(sp),parameter :: p_nuc_to_cgs = 1.0d0/6.24151d-34


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
