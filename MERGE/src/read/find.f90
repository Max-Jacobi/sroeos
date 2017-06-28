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
SUBROUTINE findthis(lr,lt,y,value,array,d1,d2,d3)

  USE Table_Sizes_Mod, ONLY : nrho,ntemp,nyp,logrho,logtemp,yp
  IMPLICIT NONE

  integer rip,rim
  integer tip,tim
  integer yip,yim

  real*8 lr,lt,y,value,d1,d2,d3
  real*8 array(*)

! Ewald's interpolator
  call intp3d(lr,lt,y,value,1,array,nrho,ntemp,nyp,logrho,logtemp,yp,d1,d2,d3)

END SUBROUTINE findthis


SUBROUTINE findall(lr,lt,y,ff)

  USE Table_Sizes_Mod, ONLY : nrho,ntemp,nyp,nvars,logrho,logtemp,yp,alltables
  IMPLICIT NONE

  real*8 ff(nvars)
  real*8 ffx(nvars,1)
  real*8 lr,lt,y
  integer i

! Ewald's interpolator
  call intp3d_many(lr,lt,y,ffx,1,alltables,&
       nrho,ntemp,nyp,nvars,logrho,logtemp,yp)
  ff(:) = ffx(:,1)

END SUBROUTINE findall


SUBROUTINE findall_short(lr,lt,y,ff)

  USE Table_Sizes_Mod, ONLY : nrho,ntemp,nyp,nvars,logrho,logtemp,yp,alltables
  IMPLICIT NONE

  real*8 ffx(8,1)
  real*8 ff(8)
  real*8 lr,lt,y
  integer i
  integer :: nvarsx = 8

! Ewald's interpolator
  call intp3d_many(lr,lt,y,ffx,1,alltables(:,:,:,1:8), &
       nrho,ntemp,nyp,nvarsx,logrho,logtemp,yp)
  ff(:) = ffx(:,1)

END SUBROUTINE findall_short

SUBROUTINE findone(lr,lt,y,ff,index)

  USE Table_Sizes_Mod, ONLY : nrho,ntemp,nyp,nvars,logrho,logtemp,yp,alltables
  IMPLICIT NONE

  real*8 ffx(1,1)
  real*8 ff(1)
  real*8 lr,lt,y
  integer index


! Ewald's interpolator
  call intp3d_many(lr,lt,y,ffx,1,alltables(:,:,:,index), &
       nrho,ntemp,nyp,1,logrho,logtemp,yp)
  ff(:) = ffx(:,1)

END SUBROUTINE findone
