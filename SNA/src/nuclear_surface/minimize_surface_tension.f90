!    This file is part of SRO_EOS.
!
!    SRO_EOS is free software: you can redistribute it AND/or modIFy
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
!    along with SRO_EOS.  IF not, see <http://www.gnu.org/licenses/>.
!
MODULE Minimize_surface_tension_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP
  USE Physical_Constants_Mod, ONLY : ZERO, ONE, TWO, HALF, TEN
  USE Surface_Tension_Mod, ONLY : Surface_Tension
  USE nwnleq_mod
  USE lautil_mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE MINIMIZE_SURFACE_TENSION (n_dim, x_sol, x_dens, &
                                   prot_frac, temperature, residue, lgt)
!   if n = 1 determine symmetric nuclear matter properties
!            => neutron and proton densities assumed to be the same
!   if n = 3 determine neutron rich matter properties
!            => proton radius set at zero
!   lgt = .true. IF ...
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: N_DIM
    REAL(DP), DIMENSION(N_DIM), INTENT(INOUT) :: x_sol, residue
    REAL(DP), INTENT(IN) :: prot_frac, temperature
    REAL(DP), DIMENSION(N_DIM) :: x1, x2, r
    REAL(DP), DIMENSION(N_DIM,N_DIM) :: g
    REAL(DP), DIMENSION(4) :: x_dens
    LOGICAL, INTENT(OUT) :: lgt
    REAL(DP) :: log10_dens_neut_in, log10_dens_prot_in
    REAL(DP) :: log10_dens_neut_out, log10_dens_prot_out
    INTEGER(I4B) :: FLAG, retry
    REAL(DP) :: RES1, RES2

!   parameters for non-linear equation solver "nleqslv"
    INTEGER(I4B) :: N, MAXIT, JACFLG(1:4),  SIZE
    INTEGER(I4B) :: METHOD, GLOBAL, XSCALM, LDR, LRWORK
    INTEGER(I4B) :: NJCNT, NFCNT, ITER, TERMCD, QRWSIZ, OUTOPT(2)
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ICDWRK
    REAL(DP) :: XTOL, FTOL, BTOL, CNDTOL
    REAL(DP) :: STEPMX, DELTA, SIGMA, TRACE, DSUB
    REAL(DP), DIMENSION(:), ALLOCATABLE :: RJAC, RWORK, RCDWRK
    REAL(DP), DIMENSION(:), ALLOCATABLE :: QRWORK, SCALEX

!   if n = 1, then x_sol = (neut/prot diffuseness)
!   if n = 3, then x_sol = (neut half width radius,
!                           neut diffuseness, prot diffuseness)
!   if any x_sol(*) > 4.d0 do not even bother finding solution
!    may take too long or subroutine nwnleq may hang
    lgt = .false.
    IF (N_DIM==1 .and. x_sol(1)>4.d0) then
      LGT =.false.
      RETURN
    ENDIF

    IF (N_DIM==3) then
      if (x_sol(1)>4.d0.or.x_sol(2)>4.D0.or.x_sol(3)>4.d0) then
        LGT =.false.
        RETURN
      ENDIF
    ENDIF

    log10_dens_neut_in  = x_dens(1)
    log10_dens_prot_in  = x_dens(2)
    log10_dens_neut_out = x_dens(3)
    log10_dens_prot_out = x_dens(4)

    retry = 0

!   set parameters for non-linear equation solver
    N = N_DIM
    MAXIT = 250 ; JACFLG(1:4) = (/0,-1,-1, 1/) ; OUTOPT(1:2) = (/0,1/)
    METHOD = 0; GLOBAL = 4; XSCALM = 1; LDR = N; LRWORK = 9*N
    XTOL = 1.D-16; FTOL = 1.D-8; BTOL = 1.D-5; CNDTOL = 1.D-7
    STEPMX = -1.D0; DELTA = -1.D0; SIGMA = 0.5D0; TRACE = 1.D0; DSUB = 0.D0
    SIZE = (METHOD+1)*N*N

10  CONTINUE

    ALLOCATE(RJAC(SIZE),RWORK(9*N),RCDWRK(3*N),ICDWRK(N),SCALEX(N))
    CALL liqsiz(n,qrwsiz)
    ALLOCATE(qrwork(qrwsiz))

    RJAC = ZERO ; RWORK = ZERO; RCDWRK = ZERO ; QRWORK = ZERO ; ICDWRK = 0
    SCALEX = ONE

!   Set initial guesses for
!    symmetric nuclear matter
    IF (N_DIM == 1) THEN
      IF (x2(1)==ZERO) X2(1) = 0.3D0
!       if (xsol(1)/=zero) x2 = xsol
    ENDIF
!    asymmetric nuclear matter
    IF (N_DIM == 3) THEN
!   Assuming a Woods-Saxon density profile
!   Set initial guess for neutron radius thickness
      x2(1) = 2.d0*(HALF-prot_frac)
!   Set initial guess for neutron and proton diffuseness parameters
      x2(2) = MAX(HALF,Temperature/20.D0)*(half/prot_frac)
      x2(3) = MAX(HALF,Temperature/20.D0)
!     if some initial guess already given, use it
      IF (sum(x_sol(1:3))/=zero) x2 = x_sol
    ENDIF

    x1 = ZERO
!   solve equation for equilibrium inside AND outside nuclei
    CALL nwnleq(x2,n,scalex,maxit,jacflg,xtol,ftol,btol,cndtol,method,global, &
                xscalm,stepmx,delta,sigma,rjac,ldr,rwork,lrwork, &
                rcdwrk,icdwrk,qrwork,qrwsiz,&
                jacob_surface_tension_derivatives,surface_tension_derivatives,&
                outopt,x1,r,g,njcnt,nfcnt,iter,termcd)
    deallocate(RJAC, RWORK, RCDWRK, ICDWRK, QRWORK, SCALEX)
!   chech whether solution x1 found is actually a solution to eq being solved.
!   sotemimes output for x1 is not a solution, but a point where nlwleq stalled.
    flag = 0
    CALL surface_tension_derivatives( x1, r, n, flag )
    residue = r

    lgt = .FALSE.
!   IF x1 is a solution THEN set as so
    IF (dot_product(r,r)<1.d-12) THEN
      x_sol = x1
      lgt = .TRUE.
      RETURN
    ENDIF

!   sometimes solution is not found with one method of nleqslv
!   but a different method can find the solution
    IF (retry<3) THEN
      retry = retry+1
!        IF (retry==1) nlglobal = 0
!        IF (retry==2) nlglobal = 2
!        IF (retry==2) nlglobal = 4
      IF (retry==1) global = 5
      IF (retry==2) global = 6
      GOTO 10
    ENDIF

    RETURN

CONTAINS

  SUBROUTINE surface_tension_derivatives ( x, s, m, flag )
!
!
!
    IMPLICIT NONE

    INTEGER (I4B), INTENT(IN) :: m, flag
    REAL(DP), INTENT(IN)  :: x(m)
    REAL(DP), INTENT(OUT) :: s(m)
    REAL(DP) :: eps, deps
    REAL(DP) :: gamma_surface_dn, gamma_surface_up
    REAL(DP) :: d_gamma_d_r_neut, d_gamma_d_a_neut, d_gamma_d_a_prot
    ! neutron/protn radius and diffuseness
    REAL(DP) :: r_neut, r_prot, a_neut, a_prot

    EPS = 1.D-2

!   for symmetric nuclear matter minimize surface tension
!    w.r.t. neut/prot diffuseness
    IF (N==1) THEN
      r_neut = zero ; r_prot = zero ; a_neut = x(1); a_prot = x(1)
      deps = MAX(eps,x(1)*eps)
      a_neut = x(1)-deps
      a_prot = a_neut
      gamma_surface_dn = &
      Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                  log10_dens_neut_in, log10_dens_prot_in, &
                  log10_dens_neut_out, log10_dens_prot_out, Temperature)
      a_neut = x(1)+deps
      a_prot = a_neut
      gamma_surface_up = &
      Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                  log10_dens_neut_in, log10_dens_prot_in, &
                  log10_dens_neut_out, log10_dens_prot_out, Temperature)
      S(1) = (gamma_surface_up - gamma_surface_dn)/(TWO*deps)
      a_neut = x(1)
      a_prot = a_neut
      RETURN
    ELSE
!   for asymmetric nuclear matter minimize surface tension
!    w.r.t. neut radius, neut diffuseness and prot diffuseness
      r_neut = x(1)
      r_prot = zero
      a_neut = x(2)
      a_prot = x(3)
    ENDIF

!   calculate derivative w.r.t. r_neut
   ! can't remember why i made the choice below! :\
    deps = MAX(eps,MAX(x(1),x(3))*eps)
   !deps = MAX(eps,x(1)*eps)
    r_neut = x(1)-deps
    gamma_surface_dn = &
    Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                log10_dens_neut_in, log10_dens_prot_in, &
                log10_dens_neut_out, log10_dens_prot_out, Temperature)
    r_neut = x(1)+deps
    gamma_surface_up = &
    Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                log10_dens_neut_in, log10_dens_prot_in, &
                log10_dens_neut_out, log10_dens_prot_out, Temperature)
    d_gamma_d_r_neut = (gamma_surface_up - gamma_surface_dn)/(TWO*deps)
    S(1) = d_gamma_d_r_neut
!   reset r_neut to input value
    r_neut = x(1)

!   calculate derivative w.r.t. a_neut
    deps = MAX(eps,x(2)*eps)
    a_neut = x(2)-deps
    gamma_surface_dn = &
    Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                log10_dens_neut_in, log10_dens_prot_in, &
                log10_dens_neut_out, log10_dens_prot_out, Temperature)
    a_neut = x(2)+deps
    gamma_surface_up = &
    Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                log10_dens_neut_in, log10_dens_prot_in, &
                log10_dens_neut_out, log10_dens_prot_out, Temperature)
    d_gamma_d_a_neut = (gamma_surface_up - gamma_surface_dn)/(TWO*deps)
    S(2) = d_gamma_d_a_neut

!   reset a_neut to input value
    a_neut = x(2)

!   calculate derivative w.r.t. r_neut
    deps = MAX(eps,x(3)*eps)
    a_prot = x(3)-deps
    gamma_surface_dn = &
    Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                log10_dens_neut_in, log10_dens_prot_in, &
                log10_dens_neut_out, log10_dens_prot_out, Temperature)
    a_prot = x(3)+deps
    gamma_surface_up = &
    Surface_Tension(r_neut, r_prot, a_neut, a_prot, &
                log10_dens_neut_in, log10_dens_prot_in, &
                log10_dens_neut_out, log10_dens_prot_out, Temperature)
    d_gamma_d_a_prot = (gamma_surface_up - gamma_surface_dn)/(TWO*deps)
    S(3) = d_gamma_d_a_prot

!   reset a_prot to input value
    a_prot = x(3)

    RETURN
  END SUBROUTINE surface_tension_derivatives

  SUBROUTINE jacob_surface_tension_derivatives ( jac, ldr, x1, m )
!
!
!
    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN) :: m, ldr
    REAL(DP), INTENT(IN)  :: x1(m)
    REAL(DP), INTENT(OUT) :: jac(ldr,*)

    jac(1:ldr,1:n) = zero

    RETURN
  END SUBROUTINE jacob_surface_tension_derivatives

END SUBROUTINE MINIMIZE_SURFACE_TENSION

END MODULE Minimize_surface_tension_Mod
