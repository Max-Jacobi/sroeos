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

module fdfunctions

	implicit none

	contains

	subroutine fdfunc1(x,par,n,fd,fdeta,fdtheta)
      include 'implno.dek'
!c..
!c..forms the fermi-dirac integrand and its derivatives with eta and theta.
!c..on input x is the integration variable, par(1) is the double precision
!c..index, par(2) is the degeneravy parameter, and par(3) is the relativity
!c..parameter. on output fd is the integrand, fdeta is the derivative
!c..with respect to eta, and fdtheta is the derivative with respect to theta.
!c..
!c..declare
      integer          n
      double precision x,par(n),dk,eta,theta,fd,fdeta,fdtheta,&
                       factor,dxst,denom,denom2,xdk,xdkp1

!c..initialize
      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xdk   = x**dk
      xdkp1 = x * xdk
      dxst  = sqrt(1.0d0 + 0.5d0*x*theta)

!c   avoid overflow in the exponentials at large x
      if ((x-eta) .lt. 1.0d2) then
       factor  = exp(x-eta)
       denom   = factor + 1.0d0
       fd      = xdk * dxst / denom
       fdeta   = fd * factor / denom
       denom2  = 4.0d0 * dxst * denom
       fdtheta = xdkp1 / denom2

      else
       factor   = exp(eta-x)
       fd       = xdk * dxst * factor
       fdeta    = fd
       denom2   = 4.0d0 * dxst
       fdtheta  = xdkp1/denom2 * factor
      endif

      return
      end subroutine fdfunc1




      subroutine fdfunc2(x,par,n,fd,fdeta,fdtheta)
      include 'implno.dek'
!c..
!c..forms the fermi-dirac integrand and its derivatives with eta and theta,
!c..when the z**2=x variable change has been made.
!c..on input x is the integration variable, par(1) is the double precision
!c..index, par(2) is the degeneravy parameter, and par(3) is the relativity
!c..parameter. on output fd is the integrand, fdeta is the derivative
!c..with respect to eta, and fdtheta is the derivative with respect to theta.
!c..
!c..declare
      integer          n
      double precision x,par(n),dk,eta,theta,fd,fdeta,fdtheta,&
                       factor,dxst,denom,denom2,xdk,xdkp1,xsq

      dk    = par(1)
      eta   = par(2)
      theta = par(3)
      xsq   = x * x
      xdk   = x**(2.0d0 * dk + 1.0d0)
      xdkp1 = xsq * xdk
      dxst  = sqrt(1.0d0 + 0.5d0 * xsq * theta)

!c   avoid an overflow in the denominator at large x:
      if ((xsq-eta) .lt. 1.d2) then
       factor  = exp(xsq - eta)
       denom   = factor + 1.0d0
       fd      = 2.0d0 * xdk * dxst/denom
       fdeta   = fd * factor/denom
       denom2  = 4.0d0 * dxst * denom
       fdtheta = 2.0d0 * xdkp1/denom2

      else
       factor  = exp(eta - xsq)
       fd      = 2.0d0 * xdk * dxst * factor
       fdeta   = fd
       denom2  = 4.0d0 * dxst
       fdtheta = 2.0d0 * xdkp1/denom2 * factor
      endif

      return
      end subroutine fdfunc2



end module fdfunctions
