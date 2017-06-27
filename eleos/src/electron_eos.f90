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
module wrap

  use eleos

  implicit none

  contains

!  I modified this from its original version
!   to only gather the thermodynamical quantities I care about.

  subroutine wrap_timmes(abar,zbar,rho,temp,etot,ptot,stot,&
     dedT,dpdT,dsdT,dedr,dpdr,dsdr,dedy,dpdy,dsdy,gamma,etaele,sound)

! subroutine wrap_timmes(abar,zbar,rho,temp,pele,ppos,&
!     eele,epos,sele,spos,prad,erad,srad,mue,muep,gamma,dedT,dpdT,dpdr,&
!     xne,dxnedd,dxnedt,dxnedz,dxneda,dxpepdt,dxpepdd,dxpepda,dxpepdz, &
!     dxeepdt,dxeepdd,dxeepda,dxeepdz,dxsepdt,dxsepdd,dxsepda,dxsepdz, &
!     dxemudt,dxemudd,dxemuda,dxemudz)

!   abar      is the average ion number
!   zbar      is the average ion charge
!   rho       is the average density in g/cm^3
!   temp      is the temperature in K

!   etot      is the total energy in erg
!   ptot      is the total pressure in
!   stot      is the total entropy

!   dedT      is the derivative of the total internal energy w.r.t. temperature
!   dpdT      is the derivative of the total pressure w.r.t. temperature
!   dsdT      is the derivative of the total entropy w.r.t. temperature

!   dedr      is the derivative of the total internal energy w.r.t. density
!   dpdr      is the derivative of the total pressure w.r.t. density
!   dsdr      is the derivative of the total entropy w.r.t. density

!   dedy      is the derivative of the total internal energy w.r.t. proton fraction
!   dpdy      is the derivative of the total pressure w.r.t. proton fraction
!   dsdy      is the derivative of the total entropy w.r.t. proton fraction

!   gamma     is the first adiabatic coefficient d (log P)/d (log rho)
!   etaele    is the electron-positron chemical potential
!   sound     is the sound speed

  use vector_eos_module
  implicit none

  integer ionmax
  parameter (ionmax=1)
  real*8 xmass(ionmax),ymass(ionmax),&
       aion(ionmax),zion(ionmax),temp,rho
  real*8 abar,zbar,etot,ptot,stot,etaele,sound
  real*8 gamma,dedT,dpdT,dsdT,dedr,dpdr,dsdr,dedy,dpdy,dsdy

  xmass(1) = 1.0d0

  aion(1) = abar
  zion(1) = zbar

  jlo_eos = 1
  jhi_eos = 1

  abar_row(1) = abar
  zbar_row(1) = zbar

  den_row(1) = rho
  temp_row(1) = temp

  call eosfxt

  if(eosfail) then
     write(6,"(1P10E15.6)") rho,temp
     write(6,"(1P10E15.6)") abar,zbar
     stop "eosfxt has failed"
  endif

  etot = etot_row(1)
  ptot = ptot_row(1)
  stot = stot_row(1)

  pele = pele_row(1)
  ppos = ppos_row(1)

!   eele = eele_row(1)
!   epos = epos_row(1)

!   sele = sele_row(1)
!   spos = spos_row(1)

!   prad = prad_row(1)
!   erad = erad_row(1)
!   srad = srad_row(1)

  gamma = gam1_row(1)
  etaele = etaele_row(1)
  sound = cs_row(1)

  dedT = det_row(1)
  dpdT = dpt_row(1)
  dsdT = dst_row(1)

  dedr = ded_row(1)
  dpdr = dpd_row(1)
  dsdr = dsd_row(1)

  dedy = dez_row(1)
  dpdy = dpz_row(1)
  dsdy = dsz_row(1)

  end subroutine wrap_timmes

end module wrap
