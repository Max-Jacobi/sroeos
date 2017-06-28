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
!    The routine intp3d_many is a variant of a routine originally 
!    provided by Ewald Mueller of the Max Planck Institute 
!    for Astrophysics, Garching, Germany.
!
SUBROUTINE intp3d_many ( x, y, z, f, kt, ft, nx, ny, nz, nvars, xt, yt, zt)
!
      implicit none
!
!---------------------------------------------------------------------
!
!     purpose: interpolation of a function of three variables in an
!              equidistant(!!!) table.
!
!     method:  8-point Lagrange linear interpolation formula
!
!     x        input vector of first  variable
!     y        input vector of second variable
!     z        input vector of third  variable
!
!     f        output vector of interpolated function values
!
!     kt       vector length of input and output vectors
!
!     ft       3d array of tabulated function values
!     nx       x-dimension of table
!     ny       y-dimension of table
!     nz       z-dimension of table
!     xt       vector of x-coordinates of table
!     yt       vector of y-coordinates of table
!     zt       vector of z-coordinates of table
!
!---------------------------------------------------------------------

      integer kt,nx,ny,nz,iv,nvars
      real*8 :: ft(nx,ny,nz,nvars)

      real*8 x(kt),y(kt),z(kt),f(kt,nvars)
      real*8 xt(nx),yt(ny),zt(nz)
      real*8 d1,d2,d3
!
!
      integer,parameter :: ktx = 1
      real*8  fh(ktx,8,nvars), delx(ktx), dely(ktx), delz(ktx), &
           a1(ktx,nvars), a2(ktx,nvars), a3(ktx,nvars), a4(ktx,nvars), &
           a5(ktx,nvars), a6(ktx,nvars), a7(ktx,nvars), a8(ktx,nvars)

      real*8 dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer n,ix,iy,iz

      IF (kt .GT. ktx)  STOP'***KTX**'
!
!
!------  determine spacing parameters of (equidistant!!!) table
!
      dx    = (xt(nx) - xt(1)) / DBLE(nx-1)
      dy    = (yt(ny) - yt(1)) / DBLE(ny-1)
      dz    = (zt(nz) - zt(1)) / DBLE(nz-1)
!
      dxi   = 1.0d0 / dx
      dyi   = 1.0d0 / dy
      dzi   = 1.0d0 / dz
!
      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi
!
      dxyzi = dxi * dyi * dzi
!
!
!------- loop over all points to be interpolated
!
      do  n = 1, kt
!
!------- determine location in (equidistant!!!) table
!
         ix = 2 + INT( (x(n) - xt(1) - 1.e-10) * dxi )
         iy = 2 + INT( (y(n) - yt(1) - 1.e-10) * dyi )
         iz = 2 + INT( (z(n) - zt(1) - 1.e-10) * dzi )
!
         ix = MAX( 2, MIN( ix, nx ) )
         iy = MAX( 2, MIN( iy, ny ) )
         iz = MAX( 2, MIN( iz, nz ) )
!
!         write(*,*) iy-1,iy,iy+1
!
!------- set-up auxiliary arrays for Lagrange interpolation
!
         delx(n) = xt(ix) - x(n)
         dely(n) = yt(iy) - y(n)
         delz(n) = zt(iz) - z(n)
!
         do iv = 1, nvars
            fh(n,1,iv) = ft(ix  , iy  , iz, iv  )
            fh(n,2,iv) = ft(ix-1, iy  , iz, iv  )
            fh(n,3,iv) = ft(ix  , iy-1, iz, iv  )
            fh(n,4,iv) = ft(ix  , iy  , iz-1, iv)
            fh(n,5,iv) = ft(ix-1, iy-1, iz, iv  )
            fh(n,6,iv) = ft(ix-1, iy  , iz-1, iv)
            fh(n,7,iv) = ft(ix  , iy-1, iz-1, iv)
            fh(n,8,iv) = ft(ix-1, iy-1, iz-1, iv)
!
!------ set up coefficients of the interpolation polynomial and
!       evaluate function values
            !
            a1(n,iv) = fh(n,1,iv)
            a2(n,iv) = dxi   * ( fh(n,2,iv) - fh(n,1,iv) )
            a3(n,iv) = dyi   * ( fh(n,3,iv) - fh(n,1,iv) )
            a4(n,iv) = dzi   * ( fh(n,4,iv) - fh(n,1,iv) )
            a5(n,iv) = dxyi  * ( fh(n,5,iv) - fh(n,2,iv) - fh(n,3,iv) + fh(n,1,iv) )
            a6(n,iv) = dxzi  * ( fh(n,6,iv) - fh(n,2,iv) - fh(n,4,iv) + fh(n,1,iv) )
            a7(n,iv) = dyzi  * ( fh(n,7,iv) - fh(n,3,iv) - fh(n,4,iv) + fh(n,1,iv) )
            a8(n,iv) = dxyzi * ( fh(n,8,iv) - fh(n,1,iv) + fh(n,2,iv) + fh(n,3,iv) + &
                 fh(n,4,iv) - fh(n,5,iv) - fh(n,6,iv) - fh(n,7,iv) )
!
            f(n,iv)  = a1(n,iv) +  a2(n,iv) * delx(n)          &
                 +  a3(n,iv) * dely(n)                         &
                 +  a4(n,iv) * delz(n)                         &
                 +  a5(n,iv) * delx(n) * dely(n)               &
                 +  a6(n,iv) * delx(n) * delz(n)               &
                 +  a7(n,iv) * dely(n) * delz(n)               &
                 +  a8(n,iv) * delx(n) * dely(n) * delz(n)
!
         enddo

      enddo
!

end SUBROUTINE intp3d_many
