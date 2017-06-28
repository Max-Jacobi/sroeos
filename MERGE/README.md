Schneider-Roberts-Ott Equation of State (SRO EOS)
=================================================

Reference:

Schneider, Roberts, Ott 2017, Submitted to Physical Review C, arXiv: [to be filled in]

License: GPLv3


Introduction:
-------------

In this directory you will find the code that merges an SNA and an NSE
table and includes the contributions of electrons, positrons, and
photons to the EOS.


Installation:
-------------

Once compiler, flags, and library dependency folders have been set in
the 'make.inc' file in the SRO EOS root directory, run

  make

This creates the 'merge_eos' executable.


Basic Usage:
------------

The SRO EOS MERGE code depends on three input files: 'input/tables.in',
'input/space.in', and 'input/transition.in'.

**** Running the 'merge_eos' code ****

** The 'input/tables.in' file has the format **

*****************************************************************************
&tables
merge_base = 'SNA_NSE',
sna_eos = '/path/to/SNA/EOS/table',
nse_eos = '/path/to/NSE/EOS/table',
only_sna = .FALSE.,
only_nse = .FALSE.,
/
*****************************************************************************

The 'merge_base' parameter is a base name for the final HDF5
table. 'sna_eos' and 'nse_eos' are, respectively the paths to the HDF5
SNA and NSE EOSs the user desires to merge. 'only_sna' and 'only_nse'
are used if the user desires to generate a pure SNA or NSE table.

** The 'input/space.in' file has the format **

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

Note that the range of this file should be at most equal to the range
of the sum of the SNA and NSE tables being merged. We advise the user
to use an equal or lower resolution in the MERGE code than the
resolutions used to generate the SNA and NSE tables. Please refer to
the User Guide for a more detailed discussion.


** The 'input/transition.in' file has the format **

*****************************************************************************
&transition
  n_transition = -4.00d0,   ! Log10(transition density in fm^-3)
  n_delta      =  0.33d0,   ! transition thickness parameter
  n_tolerance  =  1.d-4,    ! tolerance in transition function
/
*****************************************************************************

This file contains the Fortran namelist 'transition', which details
the parameters 'n_transition' and 'n_delta' that describe where the
transition between the SNA and NSE EOSs takes place and its width, see
Eq. (38) in the SRO paper.  A tolerance 'n_tolerance' is also given.
For details see the User Guide. We recommend values between 1.0d-4 and
1.0d-3 for the tolerance
