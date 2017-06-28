Schneider-Roberts-Ott Equation of State (SRO EOS)
=================================================

Reference:

Schneider, Roberts, Ott 2017
Submitted to Physical Review C
arXiv: [to be filled in]

License: GPLv3


Introduction:
-------------

This directory (SNA) contains the code that calculates the equation of
state (EOS) of nuclear matter in the single nucleus approximation
(SNA) for a given Skyrme parametrization.

Given a set of densities 'n', temperatures 'T', and proton fractions
'y', the SRO EOS code obtains the composition and thermodynamic
properties of matter in the (n,T,y) grid.

The code is a generalization of the work of Lattimer & Swesty
described in Nucl. Phys. A 535, 331 (1991).


Installation:
-------------

Once compiler, flags, and library dependency folders have been set in
the 'make.inc' file in the SRO EOS root directory, run

  make

This installs the 'test' and 'main' versions of the code.

* 'test' obtains the composition and thermodynamic properties at a
  single (n,T,y) point.

* 'main' obtains the composition and thermodynamic properties for a
  grid of (n,T,y) points.


Basic Usage:
------------

The SRO EOS code depends on two input files:

* The 'test' executable depends on
  'input/test.in'  and 'input/skyrme.in'

* The 'main' executable depends on
  'input/space.in' and 'inputs/skyrme.in'

Note that the 'skyrme.in' file is a symbolic link to one of the
provided Skyrme parameter sets.

Important Note: Some versions of gfortran fail to properly read a
namelist unless an empty line is added after the backslash '/' that
determines the end of each namelist.

**** Running the 'test' code with 'input/test.in' ****

'input/test.in' contains a Fortran namelist 'input_test_list'. It is
to be used with the 'test' executable. It looks like this:

******************************************************************************
&input_test_list
! point in phase space
proton_fraction =  10.0D-02, ! proton fraction for test
temperature     =  1.00D+02, ! temperature in MeV for test
density         =  2.00D-09, ! density in fm^-3 for test
! use initial guesses?
guess                 = .true.,
uniform_sol_guess     = -6.5d0,
non_uniform_sol_guess = -9.0D0, -110.0D0, -2.15D0,
/
*****************************************************************************

You may change any of these input parameters as desired. The code
fails for some values of (n,T,y), which are parameterization
dependent.  The most common failure is that the code fails to find a
non-uniform solution where one should exist. At other times times it
finds a solution that is not the global minimum of the free energy.

[[ERROR HANDLING]]

Based on our experience, safe choices of (n,T,y) are

*  densities n in fm^-3 in the range: -13 < log10(n) < +1;
*  temperatures T in MeV in the range:  -3 < log10(T) < +2.5;
*  proton fractions y in the range: 0.005 < y < 0.7.

The code gets slower the lower the temperature is.

**** Running the 'main' executable with 'input/space.in' ****

'input/space.in' contains a Fortran namelist called 'space_list'.  It
is to be used with the 'main' executable. It looks like this:

*****************************************************************************
&space_list
yp_min      =  0.00499D0,       ! minimum proton fraction
yp_max      =  0.66499D0,       ! maximum proton fraction
steps_in_Yp =  67,              ! steps in Yp every 0.01
! Yp_spacing = 0.01D0,          ! spacing in Yp
! (if used Yp_spacing supersedes option for steps_per_hundredth_in_yp)

log10n_min  = -13.20d0,         ! log10 of minimum density in fm^-3
log10n_max  =   1.0d0,          ! log10 of maximum density in fm^-3
steps_per_decade_in_n  =  30,   ! steps per decade in density
! log10n_spacing = 0.03333333d0,! spacing in log10(n)
! (if used log10n_spacing supersedes option for steps_per_decade_in_n)

log10T_min  =  -3.1d0,          ! log10 of minimum density in fm^-3
log10T_max  =   2.5d0,          ! log10 of maximum density in fm^-3
steps_per_decade_in_T  = 30,    ! steps per decade in density
! log10T_spacing = 0.03333333d0,! spacing in log10(T)
! (supersedes log10T_spacing option for steps_per_decade_in_n)
/
*****************************************************************************

The example 'space_list' above makes an SNA EOS table with

*  67 uniform lin spaces in proton fraction from ye_min to ye_max
* 427 uniform log spaces in density from 10**n_min to 10**n_max
* 169 uniform log spacess in temperature from 10**t_min to 10**t_max

This is a low resolution table. For a high resolution table, we
recommend to double steps_per_hundredth_in_yp, steps_per_decade_in_n,
and steps_per_decade_in_T.


**** Choosing the Skyrme Parametrization in 'input/skyrme.in' ****

'input/skyrme.in' contains a few Fortran namelists that specify the
Skyrme parametrization and other EOS parameters and code control
flags. It is needed for both the 'test' and 'main' executables.

The code will read 'input/skyrme.in', but the file itself is a
symbolic link. By default is set to one of the provided
parametrizations. If a different parametrization is desired, the
symbolic link can be changed to point to one of the other
parameterizations that we provide (or to a user-customized file).

We provide parameter sets for the following parametrizations (as
described in the SRO EOS paper):

KDE0v1, LNS, LS220, LS220*, NRAPR, SKRA, SkT1, SkT1*, Skxs20, SLy4,
and SQMC700.

For details on the output the user can refer to the User Guide.
