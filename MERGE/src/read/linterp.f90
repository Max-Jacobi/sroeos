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
!    The routine intp3d is a modified version of a routine originally
!    provided by Ewald Mueller of the Max Planck Institute for
!    Astrophysics, Garching, Germany
!
SUBROUTINE intp3d ( x, y, z, f, kt, ft, nx, ny, nz, xt, yt, zt, d1, d2, d3 )

!c
      implicit none
!c
!c---------------------------------------------------------------------
!c
!c     purpose: interpolation of a function of three variables in an
!c              equidistant(!!!) table.
!c
!c     method:  8-point Lagrange linear interpolation formula
!c
!c     x        input vector of first  variable
!c     y        input vector of second variable
!c     z        input vector of third  variable
!c
!c     f        output vector of interpolated function values
!c
!c     kt       vector length of input and output vectors
!c
!c     ft       3d array of tabulated function values
!c     nx       x-dimension of table
!c     ny       y-dimension of table
!c     nz       z-dimension of table
!c     xt       vector of x-coordinates of table
!c     yt       vector of y-coordinates of table
!c     zt       vector of z-coordinates of table
!c
!c     d1       centered derivative of ft with respect to x
!c     d2       centered derivative of ft with respect to y
!c     d3       centered derivative of ft with respect to z
!c     Note that d? only make sense when intp3d is called with kt=1
!c---------------------------------------------------------------------
!c
!c

!c
      integer kt,nx,ny,nz,ktx
      real*8 x(kt),y(kt),z(kt),f(kt)
      real*8 xt(nx),yt(ny),zt(nz)
      real*8 ft(nx,ny,nz)
      real*8 d1,d2,d3
!c
!c
      PARAMETER   (ktx = 400)
      real*8  fh(ktx,8), delx(ktx), dely(ktx), delz(ktx), &
                 a1(ktx), a2(ktx), a3(ktx), a4(ktx), &
                 a5(ktx), a6(ktx), a7(ktx), a8(ktx)

      real*8 dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer n,ix,iy,iz

      IF (kt .GT. ktx)  STOP '***KTX**'
!c
!c
!c------  determine spacing parameters of (equidistant!!!) table
!c
      dx    = (xt(nx) - xt(1)) / DBLE(nx-1)
      dy    = (yt(ny) - yt(1)) / DBLE(ny-1)
      dz    = (zt(nz) - zt(1)) / DBLE(nz-1)
!c
      dxi   = 1.0d0 / dx
      dyi   = 1.0d0 / dy
      dzi   = 1.0d0 / dz
!c
      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi
!c
      dxyzi = dxi * dyi * dzi
!c
!c
!c------- loop over all points to be interpolated
!c
      DO  n = 1, kt
!c
!c------- determine location in (equidistant!!!) table
!c
         ix = 2 + INT( (x(n) - xt(1) - 1.e-10) * dxi )
         iy = 2 + INT( (y(n) - yt(1) - 1.e-10) * dyi )
         iz = 2 + INT( (z(n) - zt(1) - 1.e-10) * dzi )
!c
         ix = MAX( 2, MIN( ix, nx ) )
         iy = MAX( 2, MIN( iy, ny ) )
         iz = MAX( 2, MIN( iz, nz ) )
!c
!c         write(*,*) iy-1,iy,iy+1
!c
!c------- set-up auxiliary arrays for Lagrange interpolation
!c
         delx(n) = xt(ix) - x(n)
         dely(n) = yt(iy) - y(n)
         delz(n) = zt(iz) - z(n)
!c
         fh(n,1) = ft(ix  , iy  , iz  )
         fh(n,2) = ft(ix-1, iy  , iz  )
         fh(n,3) = ft(ix  , iy-1, iz  )
         fh(n,4) = ft(ix  , iy  , iz-1)
         fh(n,5) = ft(ix-1, iy-1, iz  )
         fh(n,6) = ft(ix-1, iy  , iz-1)
         fh(n,7) = ft(ix  , iy-1, iz-1)
         fh(n,8) = ft(ix-1, iy-1, iz-1)
!c
!c------ set up coefficients of the interpolation polynomial and
!c       evaluate function values
!c
         a1(n) = fh(n,1)
         a2(n) = dxi   * ( fh(n,2) - fh(n,1) )
         a3(n) = dyi   * ( fh(n,3) - fh(n,1) )
         a4(n) = dzi   * ( fh(n,4) - fh(n,1) )
         a5(n) = dxyi  * ( fh(n,5) - fh(n,2) - fh(n,3) + fh(n,1) )
         a6(n) = dxzi  * ( fh(n,6) - fh(n,2) - fh(n,4) + fh(n,1) )
         a7(n) = dyzi  * ( fh(n,7) - fh(n,3) - fh(n,4) + fh(n,1) )
         a8(n) = dxyzi * ( fh(n,8) - fh(n,1) + fh(n,2) + fh(n,3) + &
                           fh(n,4) - fh(n,5) - fh(n,6) - fh(n,7) )
!c
         d1 = -a2(n)
         d2 = -a3(n)
         d3 = -a4(n)
         f(n)  = a1(n) +  a2(n) * delx(n) &
                       +  a3(n) * dely(n) &
                       +  a4(n) * delz(n) &
                       +  a5(n) * delx(n) * dely(n) &
                       +  a6(n) * delx(n) * delz(n) &
                       +  a7(n) * dely(n) * delz(n) &
                       +  a8(n) * delx(n) * dely(n) * delz(n)
!c
      ENDDO
!c
      RETURN
END SUBROUTINE
