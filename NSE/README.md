Schneider-Roberts-Ott Equation of State (SRO EOS)
=================================================

Reference:

Schneider, Roberts, Ott 2017
Submitted to Physical Review C
arXiv: [to be filled in]

License: GPLv3 (see file LICENSE)


This file provides a short description of the SRO EOS NSE code,
installation guidelines, usage, some tips and output.

See the User Guide for more details.


Introduction:
-------------


In this directory you will find the code that calculates the equation
of state (EOS) for an ensamble of nuclei in nuclear statistical
equilibrium (NSE) for a list of nuclides.

Given a set of densities 'n', temperatures 'T', and proton fractions
'y', the NSE code obtains the composition and thermodynamic
properties of matter in the (n,T,y) grid.


Installation:
-------------

Once compiler, flags, and library dependency folders have been set in
the 'make.inc' file in the SRO EOS root directory (one up),
installation should be as simple as typing

    make

This installs the 'test' and 'main' versions of the code.

* 'test' obtains the composition and thermodynamic properties at a
  single (n,T,y) point.

* 'main' obtains the composition and thermodynamic properties for a
  grid of (n,T,y) points.


Usage:
------

The SRO NSE EOS code depends on four (test case) or five (main case) input
files.

Two of these input files, 'input/isotopes.in' and 'input/partition.in', should
not be modified by the user. For details see the User Guide.

The other input file dependencies are

* The 'test' executable depends on
  'input/test.in' and 'input/list.in'

* The 'main' executable depends on
  'input/space.in', 'inputs/list.in', and 'input/output.in'

Obs.: Some versions of gfortran fail to propery read a list unless an empty
line is added after the backslash '/' that determines the end of the list.

**** Running the 'test' code with 'input/test.in' ****

'input/test.in' contains a Fortran namelist 'input_test_list'. It is
to be used with the 'test' executable. It looks like this:

******************************************************************************
&input_test_list
! point in phase space
proton_fraction = 30.0D-02, ! proton fraction for test
temperature     = 1.00D-00, ! temperature in MeV for test
density         = 1.00D-03, ! density in fm^-3 for test
/
******************************************************************************

You may change any of these input parameters as desired. In rare cases
the code fails for some values of (n,T,y), which depends on the list
of nuclei used.  The most common failure is that the code fails to
find a set of isotopes that minimizes the free energy of the
system. To prevent that we start from a high temperature and lower the
temperature slowly until the desired temperature is reached. This
almost always guarantees a solution is found. If not an error message
is output.

[[ERROR HANDLING]]

Based on our experience, safe choices of (n,T,y) are

*  densities n in fm^-3 in the range: -16 < log10(n) < -1;
*  temperatures T in MeV in the range:  -3 < log10(T) < +2.5;
*  proton fractions y in the range: 0.005 < y < 0.7.

The code gets slower as the temperature decreases.

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
log10n_max  =  -1.0d0,          ! log10 of maximum density in fm^-3
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

The example 'space_list' above makes an NSE EOS table with

*  67 uniform lin spaces in proton fraction from ye_min to ye_max
* 367 uniform log spaces in density from 10**n_min to 10**n_max
* 169 uniform log spacess in temperature from 10**t_min to 10**t_max

This is a low resolution table. For a high resolution table, we
recommend to double steps_per_hundredth_in_yp, steps_per_decade_in_n,
and steps_per_decade_in_T.

Details about the output table and debugging info are given in the 'TABLE_INPUT'
namelist in the 'input/output.in' file. It looks like this:

&TABLE_INPUT
  make_hdf5_table = .true.,
  write_solutions_to_file = .false.,
  filename_base   = "NSE_0023",  ! base for HDF5 table name
  output_directory= "NSE_0023",  ! directory for ascii output files
/

The parameter 'write_solutions_to_file' is used if the user desires to perform
some debugging. Its details are found in the User Guide.


**** Choosing the isotope list in 'input/list.in' ****

'input/list.in' contains a list of isotopes to be used when looking for a
in NSE.

The code will read 'input/list.in', but the file itself is a symbolic link.
By default is set to one of the provided lists in the 'tables' directory.
If a different list is desired, the symbolic link can be changed to point to
one of the other lists that we provide (or to a user-customized file).

We provide lists for the isotope sets used in the SRO paper. The list names are
'tables/list_#' where # denotes the number of isotopes in the list.

For details on the output the user can refer to the User Guide.
