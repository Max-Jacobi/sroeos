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
MODULE liqrev_mod

USE Kind_Types_Mod
USE lautil_mod

  IMPLICIT NONE

  CONTAINS
!-----------------------------------------------------------------------------

      SUBROUTINE liqrev(n,r,ldr,diag,b,x,sdiag,wrk)
      INTEGER(I4B) :: n,ldr
      REAL(DP) :: r(ldr,*),b(*),x(*),sdiag(*),wrk(*)
      REAL(DP) :: diag

!-----------------------------------------------------------------------------
!
!     Arguments
!
!       In     n      Integer         order of R.
!       Inout  R      Real(ldr,*)     upper triangular matrix R from QR
!                                     unaltered
!                                     strict lower triangle CONTAINS
!                                        transposed strict upper triangle of the upper
!                                        triangular matrix S.
!
!       In     diag   Real            scalar for matrix D
!
!       In     ldr    Integer         leading dimension of the array R.
!       In     b      Real(*)         vector of size n
!
!       Out    x      Real(*)         vector of size n
!                                     on output CONTAINS least squares solution
!                                     of the system R*x = b, D*x = 0.
!
!       Out    sdiag  Real(*)         vector of size n, containing the
!                                     diagonal elements of the upper
!                                     triangular matrix S.
!
!       Out    wrk    Real(*)         workspace of size n.
!
!     Description
!
!     Given an n by n upper triangular matrix R, a diagonal matrix D with positive entries
!     and an n-vector b, determine an x which solves the system
!
!         |R*x| = |b|
!         |D*x| = |0|
!
!     in the least squares sense where D=diag*I.
!     This routine can be USEd for two different purposes.
!     The first is to provide a method of slightly modifying a singular or ill-conditioned matrix.
!     The second is for calculating a least squares solution to the above problem within
!     the context of e.g. a Levenberg-Marquardt algorithm combined with a More-Hebden algorithm
!     to determine a value of D (diagonal mu) such that x has a predetermined 2-norm.
!
!     The routine could also be USEd when the matrix R from the QR decomposition of a Jacobian
!     is ill-conditioned (or singular). Then it is difficult to calculate a Newton step
!     accurately (Dennis and Schnabel). D&S advise perturbing trans(J)*J with a positive
!     diagonal matrix.
!
!     The idea is  to solve (J^T * J + mu*I)x=b where mu is a small positive number.
!     Calculation of mu must be done in the CALLing routine.
!     Using a QR decomposition of J solving this system
!     is equivalent solving (R^T*R + mu*I)x=b, where R comes from the QR decomposition.
!     Solving this system is equivalent to solving the above least squares problem with the
!     elements of the matrix D set to sqrt(mu) which should be done in the CALLing routine.
!
!     On output the routine also provides an upper triangular matrix S such that
!     (see description of arguments above for the details)
!
!         (trans(R)*R + D*D) = trans(S)*S .
!
!     Method USEd here is described in
!     Nocedal and Wright, 2006, Numerical Optimization, Springer, ISBN 978-0-387-30303-1
!     page 258--261 (second edition)
!-----------------------------------------------------------------------------

      INTEGER(I4B) :: j,k
      REAL(DP) :: bj,c,s,sum,temp
      REAL(DP) :: ddot
      REAL(DP), PARAMETER :: Rzero = 0.0d0

!     copy R and b to preserve input and initialise S.
!     Save the diagonal elements of R in wrk.
!     Beware: the algorithm operates on an upper triangular matrix,
!     which is stored in lower triangle of R.
!
      DO j=1,n
         CALL dcopy(n-j+1,r(j,j),ldr,r(j,j),1)
         wrk(j) = r(j,j)
      ENDDO
      CALL dcopy(n,b,1,x,1)

!     eliminate the diagonal matrix D using givens rotations.
!     Nocedal method: start at the bottom right
!     at END of loop R CONTAINS diagonal of S
!     save in sdiag and restore original diagonal of R

      DO j=n,1,-1

!        initialise the row of D to be eliminated

         CALL nuzero(n-j+1,sdiag(j))
         sdiag(j) = diag

!        the transformations to eliminate the row of D

         bj = Rzero
         DO k=j,n

!           determine a givens rotation which eliminates the
!           appropriate element in the current row of D.
!           accumulate the transformation in the row of S.

!           eliminate the diagonal element in row j of D
!           this generates fill-in in columns [j+1 .. n] of row j of D
!           successively eliminate the fill-in with givens rotations
!           for R[j+1,j+1] and D[j,j+1].
!           rows of R have been copied into the columns of R initially (see above)
!           perform all operations on those columns to preserve the original R

            IF (sdiag(k) .ne. Rzero) THEN

               CALL nuvgiv(r(k,k),sdiag(k),c,s)
               IF( k .lt. n ) THEN
                   CALL drot(n-k,r(k+1,k),1,sdiag(k+1),1,c,s)
               ENDIF

!              compute the modified element of (b,0).

               temp =  c*x(k) + s*bj
               bj   = -s*x(k) + c*bj
               x(k) = temp

            ENDIF

         ENDDO

      ENDDO

!     retrieve diagonal of S from diagonal of R
!     restore original diagonal of R

      DO k=1,n
         sdiag(k) = r(k,k)
         r(k,k) = wrk(k)
      ENDDO

!     x now CONTAINS modified b
!     solve trans(S)*x = x
!     still to be done: guard against division by 0 to be absolutely safe
!     CALL dblepr('liqrev sdiag', 12, sdiag, n)
      x(n) = x(n) / sdiag(n)
      DO j=n-1,1,-1
         sum  = ddot(n-j,r(j+1,j),1,x(j+1),1)
         x(j) = (x(j) - sum)/sdiag(j)
      ENDDO

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE dtrstt(S,ldr,n,sdiag,x)
      INTEGER(I4B) ::ldr, n
      REAL(DP) :: S(ldr,*), sdiag(*), x(*)
      INTEGER(I4B) :: j
      REAL(DP) :: sum, ddot

!     solve S*x = x where x is the result from SUBROUTINE liqrev
!     S is a lower triangular matrix with diagonal entries in sdiag()
!     and here it is in the lower triangular part of R as RETURNed by liqrev

      x(1) = x(1) / sdiag(1)
      DO j=2,n
         sum  = ddot(j-1,S(j,1),n,x,1)
         x(j) = (x(j) - sum)/sdiag(j)
      ENDDO

      RETURN
      END SUBROUTINE

! ----------------------------------------------------------------------

      SUBROUTINE nuvgiv(x,y,c,s)
      REAL(DP) :: x,y,c,s

!     Parameters
!
!     Inout   x     Real       x input / c*x+s*y on output
!     Inout   y     Real       y input / 0       on output
!     Out     c     Real       c of tranformation (cosine)
!     Out     s     Real       s of tranformation (  sine)
!
!     Description
!
!     Nuvgiv calculates the givens rotator
!
!             |  c   s |
!         G = |        |
!             | -s   c |
!
!     with  c*c+s*s=1
!
!     for which G * | x | = | z |
!                   | y |   | 0 |
!
!     resulting in
!
!            c * x + s * y = z
!           -s * x + c * y = 0   ==>  s/c = y/x or c/s = x/y
!
!     Use Lapack dlartg routine
!     RETURN c and s and the modified x and y

      REAL(DP) :: t

      REAL(DP), PARAMETER :: Rzero=0.0d0

      CALL dlartg(x,y,c,s,t)
      x = t
      y = Rzero
      RETURN
      END SUBROUTINE

END MODULE liqrev_mod
