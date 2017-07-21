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
MODULE Nuclear_Surface_Fit_Mod

  USE Kind_Types_Mod, ONLY : I4B, DP
  USE Physical_Constants_Mod, &
      ONLY : ZERO, HALF, ONE, TWO, R_5_3, THREE, FOUR, SIX, TEN
  USE Surface_Properties_Mod
  USE Levenberg_Marquardt_Mod

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Nuclear_Surface_Fit(data_points,TAB_SURF_TENS_VAL)
!
! Use surface tension table to fit a function for it of the form
! f_surface(x,T) = sigma_s*h(T/Tc(x))*(2*2^a+q)/(x^-a+q+(1-x)^-a)
! where h(T/Tc(x)) = (1+(T/Tc(x))²)^p
! Parameters are (p,a,q) where 1<p<3, 2<a<4 and q>0 usually q ~ (10-100).
! We perform the fitting by using an initial guesses
! p = 2, a = 5/3 and q = 10 and try to
! minimize chi² = sum_i ((fit(i)-tab(i))/tab(i))²
!
    INTEGER(I4B), INTENT(IN) :: data_points
    REAL(DP), DIMENSION(10000,3), INTENT(IN) :: TAB_SURF_TENS_VAL
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: integer_work_array
    REAL(DP), DIMENSION(:), ALLOCATABLE :: SURFACE_TENSION_ARRAY
    INTEGER(I4B) :: info, dimensions
    REAL(DP) :: tolerance
    REAL(DP), DIMENSION(:), ALLOCATABLE :: fitting_values

    ALLOCATE(SURFACE_TENSION_ARRAY(data_points))
    ALLOCATE(integer_work_array(data_points))

    SURFACE_TENSION_ARRAY(1:data_points) = TAB_SURF_TENS_VAL(1:data_points,3)

    integer_work_array  = 0
    info = 0

    !write (*,*) data_points

    Surface_sigma = SURFACE_TENSION_ARRAY(1)

    !write (*,*) surface_sigma

!   initial guess for lambda, p, q
    dimensions = 3
    ALLOCATE(fitting_values(dimensions))
    fitting_values(1:dimensions) = (/ THREE, R_5_3, TEN /)
    !, ONE, -ONE, ONE, -ONE /)
    tolerance = 1.d-12

    call lmdif1(SURFACE_FIT, data_points, dimensions, &
                fitting_values, SURFACE_TENSION_ARRAY, tolerance, &
                info, integer_work_array)

    Surface_lambda = fitting_values(1)
    Surface_p      = fitting_values(2)
    Surface_q      = fitting_values(3)

!    T_crit_coeffs(1) = fitting_values(4)
!    T_crit_coeffs(2) = fitting_values(5)
!    T_crit_coeffs(3) = fitting_values(6)
!    T_crit_coeffs(4) = fitting_values(7)

    !write (*,*) info
    !write (*,*) fitting_values

CONTAINS

  SUBROUTINE SURFACE_FIT(m, n, x, fvec, iflag)

    USE Make_Tables_Mod, ONLY : output_directory
    USE, INTRINSIC :: IEEE_ARITHMETIC

    IMPLICIT NONE

    INTEGER(I4B), INTENT(IN)  :: m, n
    REAL (DP), INTENT(IN)     :: x(:)
    REAL (DP), INTENT(IN OUT) :: fvec(:)
    INTEGER(I4B), INTENT(IN OUT) :: iflag
    INTEGER(I4B) :: i, j, k
    REAL(DP) :: prot_frac, neutron_excess, temperature
    REAL(DP) :: H, Ta, Tb, Tc, Td, Critical_T
    REAL(DP) :: NUM, DEN, lambda, p, q
    REAL(DP) :: surface_tension, surface_tension_fit
!      ----------
!      calculate the functions at x and return this vector in fvec.
!      ----------
    lambda = X(1)
    P = X(2)
    Q = X(3)

    Ta = T_crit_coeffs(1)
    Tb = T_crit_coeffs(2)
    Tc = T_crit_coeffs(3)
    Td = T_crit_coeffs(4)

    fvec = zero

    OPEN(23,file=trim(adjustl(output_directory))&
    //"/surface_tension_vs_sigma.dat")

    DO k = 1, data_points
      IF (TAB_SURF_TENS_VAL(k,1) < prot_frac) WRITE (23,*)
      prot_frac   = TAB_SURF_TENS_VAL(k,1)
      Temperature = TAB_SURF_TENS_VAL(k,2)
      Surface_tension = TAB_SURF_TENS_VAL(k,3)/Surface_sigma


      neutron_excess = (ONE-TWO*prot_frac)

      Critical_T = T_crit*( Ta + Tb*neutron_excess**TWO  &
                               + Tc*neutron_excess**FOUR &
                               + Td*neutron_excess**SIX )

      IF (Temperature>Critical_T) THEN
        H = ZERO
      ELSE
        H = (ONE - (Temperature/Critical_T)**TWO)**p
      ENDIF

      NUM = TWO*TWO**lambda+q
      DEN = prot_frac**(-lambda)+q+(ONE-prot_frac)**(-lambda)

      surface_tension_fit = (H*NUM/DEN)

      WRITE (23,"(16ES20.12)") prot_frac, Temperature, Surface_tension, &
      Surface_tension_fit, Surface_tension_fit/Surface_tension

      fvec(k) = (surface_tension_fit - surface_tension)/surface_tension

      IF ( Surface_tension == ZERO.or.ieee_is_nan(Surface_tension) ) THEN
        !WRITE (*,*) I, J, H, NUM, DEN
        IF (Surface_tension == ZERO) STOP &
                    'ERROR EVALUATING FVEC(K) AND/OR Surface_tension'
        IF (ieee_is_nan(Surface_tension))  STOP &
                    'ERROR: Surface_tension IS NOT A NUMBER'
      ENDIF

    ENDDO

    CLOSE(23)

    !write (*,*) ' x   = ', x
    !write (*,*) ' f.f = ', dot_product(fvec,fvec)

    RETURN

  END SUBROUTINE SURFACE_FIT

  END SUBROUTINE Nuclear_Surface_Fit

END MODULE Nuclear_Surface_Fit_Mod
