Schneider-Roberts-Ott Equation of State (SRO EOS)
=================================================

Reference:

Schneider, Roberts, Ott 2017\\
Submitted to Physical Review C\\
arXiv: [to be filled in]

License: GPLv3 (see file LICENSE)

Documentation: User_Guide/User_Guide.pdf

### Quick Start:

(1) Set up a make.inc file for your machine. Examples are provided,
use a symbolic link to point to your custom file. See the User Guide
for compiler and library dependencies. In short, you'll need Fortran,
C/C++, HDF5, BLAS, and LAPACK.


(2) The SRO EOS comes in three parts:

Directory SNA:

Single nucleus approximation code, implementation of the liquid-drop
model with Skyrme interaction. Creates an SNA only table for the
high-density part of the EOS.

Directory NSE:

Nuclear statistical equilibrium code; free gas, no nuclear
interaction. Creates an NSE only table for the low-density part of the
EOS.

Directory MERGE:

Code that merges SNA and NSE tables and adds electrons, positrons,
and photons.


(3) There are two more directories:

Directory eleos:

The directory 'eleos' contains a copy of the Timmes EOS,
Timmes & Arnett 1999, ApJS, 125:277. The source code in this directory
is accessed by the SNA and MERGE codes, but there is no standalone
executable.

Directory formaline_extract:

The SRO EOS uses the Formaline method by Erik Schnetter to add source
code and input files to the output HDF5 tables. See the User Guide for
a detailed discussion. This directory contains C++ code that will
extract a dataset containing source code or an input file from an
HDF5 table.


(4) Read the User Guide to understand inputs and outputs of the SNA
and NSE codes. Make sure to also understand the input requirements
for the MERGE code. In short, the NSE and SNA tables should have
higher resolution than the desired output resolution for the final,
merged table.


(5) Enjoy!


Feel free to contact us with bug reports and questions at
SROEOS@stellarcollapse.org. Please do note that our resources are
limited and we may not be able to help with basic linux/unix and
compilation questions. We recommend you ask a friendly colleague
at your institution for help with such problems.



