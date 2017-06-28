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
!  This file is originally part of the NLEQSLV of Berend Hasselman available at
!  https://cran.r-project.org/web/packages/nleqslv/index.html
!
!  It has been slightly modified by Andre da Silva Schneider
!   to fit SRO_EOS source code
!
MODULE liqrup_mod

USE Kind_Types_Mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE liqrup(Q,ldq,n,R,ldr,u,v,wk)
      INTEGER(I4B) :: ldq,n,ldr
      REAL(DP) :: Q(ldq,*),R(ldr,*),u(*),v(*),wk(*)

!-----------------------------------------------------------------------------
!
!     Arguments
!
!     Inout  Q       Real(ldq,n)      orthogonal matrix from QR
!     In     ldq     Integer          leading dimension of Q
!     In     n       Integer          order of Q and R
!     Inout  R       Real(ldr,n)      upper triangular matrix R from QR
!     In     ldr     Integer          leading dimension of R
!     In     u       Real(*)          vector u of size n
!     In     v       Real(*)          vector v of size n
!     Out    wk      Real(*)          workspace of size n
!
!     on RETURN
!
!        Q       Q is the matrix with orthonormal columns in a QR
!                decomposition of the matrix B = A + u*v'
!
!        R       R is the upper triangular matrix in a QR
!                decomposition of the matrix B = A + u*v'
!
!     Description
!
!     The matrices Q and R are a QR decomposition of a square matrix
!     A = Q*R.
!     Given Q and R, qrupdt computes a QR decomposition of the rank one
!     modification B = A + u*trans(v) of A. Here u and v are vectors.
!
!     Source : procedure outlined in Dennis & Schnabel (AppENDix A)
!              Algorithm 3.1.4 and 3.4.1a
!              modified (to USE Lapack routines and more)
!
!-----------------------------------------------------------------------------

!     Local variables and functions

      INTEGER(I4B) :: k,i
      REAL(DP) :: ddot

!     calculate wk = trans(Q)*u

      DO i=1,n
         wk(i) = ddot(n,Q(1,i),1,u,1)
      ENDDO

!     zero components wk(n),wk(n-1)...wk(2)
!     and apply rotators to R and Q.

      DO k=n-1,1,-1
         CALL jacrot(wk(k),wk(k+1),k,n,Q,ldq,R,ldr,k)
      ENDDO

!     r(1,1:n) += wk(1)*v(1:n)
      CALL daxpy(n,wk(1),v,1,R(1,1),ldr)

!     R is of upper hessenberg form. Triangularize R.
!      kr argument == k+1 to start applying rotation at column k+1
!      otherwise R(k,k) will be rotated twice and this way it also
!      avoids tiny roundoff errors.

      DO k=1,n-1
         CALL jacrot(R(k,k),R(k+1,k),k,n,Q,ldq,R,ldr,k+1)
      ENDDO

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------------

      SUBROUTINE jacrot(a,b,k,n,Q,ldq,R,ldr,kr)

      REAL(DP) :: a, b
      INTEGER(I4B) :: k,n,ldr,ldq,kr
      REAL(DP) :: Q(ldq,*), R(ldr,*)

!-----------------------------------------------------------------------------
!
!     Arguments
!
!     Inout  a       Real             rotate argument
!     Inout  b       Real             rotate argument to rotate to zero
!     In     k       Integer          row/column number for rotation
!     In     n       Integer          order of Q and R
!     Inout  Q       Real(ldq,n)      orthogonal matrix from QR
!     In     ldq     Integer          leading dimension of Q
!     Inout  R       Real(ldr,n)      upper triangular matrix R from QR
!     In     ldr     Integer          leading dimension of R
!     In     u       Real(*)          vector u of size n
!     In     v       Real(*)          vector v of size n
!     In     kr      Integer          start R rotation in column kr
!                                     (should be k or k+1)
!
!-----------------------------------------------------------------------------

      REAL(DP) :: t
      REAL(DP) :: c,s
      REAL(DP), PARAMETER ::  Rzero=0.0d0

      CALL dlartg(a,b,c,s,t)
      a = t
      b = Rzero
      CALL drot(n-kr+1,R(k,kr),ldr,R(k+1,kr),ldr,c,s)
      CALL drot(n     ,Q(1,k) ,1  ,Q(1,k+1) ,1  ,c,s)

      RETURN
      END SUBROUTINE

END MODULE liqrup_mod
