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
MODULE READ_EOS

  USE Kind_Types_Mod, ONLY : I4B, DP

  IMPLICIT NONE

CONTAINS
! #########################################################
!
! Copyright C. D. Ott and Evan O'Connor, July 2009
!
!  Modified by A. S. Schneider, May 2016
!  to use with LSEOSF90 code
!
!  Enter with (rho, eps, ye) to retrieve 21
!  NSE quantites below.
!
!  Changed to nuclear units...
!
! UNITS: density              g/cm^3
!        temperature          MeV
!        ye                   number fraction per baryon
!        energy               erg/g
!        pressure             dyn/cm^2
!        chemical potentials  MeV
!        entropy              k_B / baryon
!        cs2                  cm^2/s^2 (not relativistic)
!
! keyerr --> error output; should be 0
! rfeps --> root finding relative accuracy, set around 1.0d-10
! keytemp: 0 -> coming in with rho,eps,ye (solve for temp)
!          1 -> coming in with rho,temperature,ye
!          2 -> coming in with rho,entropy,ye (solve for temp)
!          3 -> coming in with pressure,temp,ye (solve for rho)
!

  ! index variable mapping:
  !  1 -> press
  !  2 -> energy
  !  3 -> entropy
  !  4 -> muh
  !  5 -> mun
  !  6 -> mup
  !  7 -> dpdn
  !  8 -> dpdt
  !  9 -> dpdy
  ! 10 -> dsdn
  ! 11 -> dsdt
  ! 12 -> dsdy
  ! 13 -> dmuhdn
  ! 14 -> dmuhdt
  ! 15 -> dmuhdy
  ! 16 -> xn
  ! 17 -> xp
  ! 18 -> xa
  ! 19 -> xh
  ! 20 -> abar
  ! 21 -> zbar
  ! 22 -> xl
  ! 23 -> albar
  ! 24 -> zlbar


  SUBROUTINE get_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xmuh,xmun,xmup,  &
            xdpdn,xdpdt,xdpdy,xdsdn,xdsdt,xdsdy,xdmuhdn,xdmuhdt,xdmuhdy, &
            xxn,xxp,xxa,xxh,xxl,xabar,xzbar,xalbar,xzlbar,xrad,xu,&
            keytemp,keyerr,rfeps)

    USE Table_Sizes_Mod

    IMPLICIT NONE

    REAL(DP), INTENT(inout)   :: xrho,xtemp,xye
    REAL(DP), INTENT(out)     :: xenr,xprs,xent,xmuh,xmun,xmup
    REAL(DP), INTENT(out)     :: xdpdn,xdpdt,xdpdy
    REAL(DP), INTENT(out)     :: xdsdn,xdsdt,xdsdy
    REAL(DP), INTENT(out)     :: xdmuhdn,xdmuhdt,xdmuhdy
    REAL(DP), INTENT(out)     :: xxa,xxn,xxp,xxh,xxl
    REAL(DP), INTENT(out)     :: xabar,xzbar,xalbar,xzlbar
    REAL(DP), INTENT(out)     :: xrad,xu
    REAL(DP), INTENT(in)      :: rfeps
    INTEGER(I4B), INTENT(in)  :: keytemp
    INTEGER(I4B), INTENT(out) :: keyerr


!   local variables
    REAL(DP)     :: lr,lt,y,xx,xeps,leps,xs,xpressure
    REAL(DP)     :: d1,d2,d3
    REAL(DP)     :: ff(nvars)
    INTEGER(I4B) :: keyerrt = 0
    INTEGER(I4B) :: keyerrr = 0

    if(xrho.gt.eos_rhomax) then
      write (*,*) xrho, eos_rhomax
      stop "nuc_eos: rho > rhomax"
    endif

    if(xrho.lt.eos_rhomin) then
      write (*,*) xrho, eos_rhomin
      stop "nuc_eos: rho < rhomin"
    endif

    if(xye.gt.eos_ypmax) then
      write (*,*) xye, eos_ypmax
      stop "nuc_eos: ye > ypmax"
    endif

    if(xye.lt.eos_ypmin) then
      write (*,*) xye, eos_ypmin
      stop "nuc_eos: ye < ypmin"
    endif

    if(keytemp.eq.1) then
       if(xtemp.gt.eos_tempmax) then
         write (*,*) xtemp, eos_tempmax
         stop "nuc_eos: temp > tempmax"
       endif

       if(xtemp.lt.eos_tempmin) then
         write (*,*) xtemp, eos_tempmin
         stop "nuc_eos: temp < tempmin"
       endif
    endif

    lr = log10(xrho)
    lt = log10(xtemp)
    y = xye
    xeps = xenr + energy_shift
    leps = log10(max(xeps,1.0d0))

    keyerr = 0

    if (keytemp.eq.0) then
       !need to find temperature based on xeps
       call findtemp(lr,lt,y,leps,keyerrt,rfeps)
       if(keyerrt.ne.0) then
          keyerr = keyerrt
          write(*,*) "Did not find temperature", keyerr
       endif
       xtemp = 10.0d0**lt

    elseif (keytemp.eq.2) then
       !need to find temperature based on xent
       xs = xent
       call findtemp_entropy(lr,lt,y,xs,keyerrt,rfeps)
       keyerr = keyerrt
       xtemp = 10.0d0**lt

    elseif (keytemp.eq.3) then
       !need to find rho based on xprs
       xpressure = log10(xprs)
       call findrho_press(lr,lt,y,xpressure,keyerrr,rfeps)
       keyerr = keyerrr
       if(keyerrr.ne.0) then
          write(*,*) "Problem in findrho_press:", keyerr
          keyerr = keyerrr
          return
       endif
       xrho = 10.0d0**lr

    endif

    ! have rho,T,ye; proceed:
    call findall(lr,lt,y,ff)

    !unless we want xprs to be constant (keytemp==3), reset xprs
    if(.not.keytemp.eq.3) then
       xprs = ff(1)
    endif

    !unless we want xenr to be constant (keytemp==0), reset xenr
    if(.not.keytemp.eq.0) then
       xenr = ff(2)
    endif

    !unless we want xent to be constant (keytemp==2), reset xent
    if(.not.keytemp.eq.2) then
       xent = ff(3)
    endif

    xmuh = ff(4)
    xmun = ff(5)
    xmup = ff(6)

    xdpdn = ff(7)
    xdpdt = ff(8)
    xdpdy = ff(9)

    xdsdn = ff(10)
    xdsdt = ff(11)
    xdsdy = ff(12)

    xdmuhdn = ff(13)
    xdmuhdt = ff(14)
    xdmuhdy = ff(15)

    xxn = ff(16)
    xxp = ff(17)
    xxa = ff(18)
    xxh = ff(19)
    xxl = ff(22)

    xabar = ff(20)
    xzbar = ff(21)
    xalbar = ff(23)
    xzlbar = ff(24)

    xrad  = ff(25)
    xu    = ff(26)

  END SUBROUTINE get_eos_full

END MODULE READ_EOS
