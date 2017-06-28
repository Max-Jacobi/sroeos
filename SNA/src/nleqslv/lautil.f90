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
MODULE lautil_mod

USE Kind_Types_Mod, ONLY : DP, I4B

IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------------

      SUBROUTINE liqrfa(a, lda, n, tau, work, wsiz, info)
      INTEGER(I4B) ::  lda, n, wsiz, info
      REAL(DP) ::  a(lda,*), tau(*), work(*)

!-------------------------------------------------------------
!
!     QR decomposition of A(n,n)
!
!     Arguments
!
!      Inout A        Real(Lda,n)    Matrix to transform.
!      In    lda      Integer        Leading dimension of A
!      In    n        Integer        number of rows/cols A
!      Out   tau      Real(n)        Information for recovering
!      Out   work     Real(*)        workspace
!      In    wsiz     Integer        size of work()

!     Lapack blocked QR (much faster for larger n)
!
!-------------------------------------------------------------

      CALL dgeqrf(n,n,A,lda,tau,work,wsiz,info)

      RETURN
      END SUBROUTINE

!=============================================================

      SUBROUTINE liqsiz(n,wrksiz)
      INTEGER(I4B) ::  n, wrksiz

!-------------------------------------------------------------------------
!     Query the size of the double precision work array required
!     for optimal operation of the Lapack QR routines
!-------------------------------------------------------------------------

      REAL(DP) :: A(1), work(1)
      INTEGER(I4B) :: lwork, info

      lwork = -1
      CALL dgeqrf(n,n,A,n,work,work,lwork,info)
      IF( info .ne. 0 ) THEN
          wrksiz = -1
      ELSE
          wrksiz = int(work(1))
      ENDIF

      RETURN
      END SUBROUTINE

!=============================================================

      SUBROUTINE liqrqt(a, lda, n, tau, qty, work, wsiz, info)
      INTEGER(I4B) ::  lda, n, wsiz, info
      REAL(DP) :: a(lda,*), tau(*), qty(*), work(*)

!-------------------------------------------------------------
!      Arguments
!
!      In    A     Real(Lda, n)    QR decomposition
!      In    Lda   Integer         Leading dimension A
!      In    n     Integer         Number of rows/columns in A
!      In    tau   Integer         HoUSEholder constants from QR
!      Inout qty   Real(n)         On input, vector y
!                                  On output, trans(Q)*y
!      Out   work  Real(*)         workspace
!      In    wsiz  Integer         size of work()
!
!     Liqrqt calculates trans(Q)*y from the QR decomposition
!
!     Lapack blocked
!-------------------------------------------------------------

      CALL dormqr('L','T',n,1,n,A,lda,tau,qty,n,work,wsiz,info)

      RETURN
      END SUBROUTINE

!=============================================================

      SUBROUTINE liqrqq(q,ldq,tau,n,work,wsiz,info)
      INTEGER(I4B) ::  n, ldq, wsiz, info
      REAL(DP) :: q(ldq,*),tau(*),work(*)

!     Arguments
!
!     Inout  Q     Real(ldq,*)     On input, QR decomposition
!                                    On output, the full Q
!     In     ldq   Integer         leading dimension of Q
!     In     tau   Real(n)         HoUSEholder constants of
!                                     the QR decomposition
!     In     n     Integer         number of rows/columns in Q
!     Out    work  Real(n)         workspace of length n
!     In     wsiz  Integer         size of work()
!
!     Generate Q from QR decomposition Liqrfa (dgeqr2)
!
!     Lapack blocked
!-------------------------------------------------------------

      CALL dorgqr(n,n,n,q,ldq,tau,work,wsiz,info)

      RETURN
      END SUBROUTINE

!-----------------------------------------------------------------------------

      SUBROUTINE nuzero(n,x)
      INTEGER(I4B) ::  n
      REAL(DP) :: x(*)

!     Parameters:
!
!     In    n        Integer           Number of elements.
!     In    x        Real(*)           Vector of reals.
!
!     Description:
!
!     Nuzero sets all elements of x to 0.
!     Does nothing when n <= 0

      REAL(DP), PARAMETER :: Rzero = 0.d0
      INTEGER(I4B) ::  i

      DO i=1,n
         x(i) = Rzero
      ENDDO

      RETURN
      END SUBROUTINE

END MODULE lautil_mod
