Schneider-Roberts-Ott Equation of State (SRO EOS)
=================================================

Reference:

Schneider, Roberts, Ott 2017
Submitted to Physical Review C
arXiv: [to be filled in]

License: GPLv3 (see file LICENSE)


Formaline
---------

Formaline is used to store source and input files in SRO EOS HDF5
tables. See the User Guide for a detailed discussion of this method
and its rational.

This directory contains C++ code that will extract code and input
files (by dataset name) from a given HDF5 table.

Compile with 'make' (provided ../make.inc is set up correctly).

Run with

  ./extract_from_h5 [HDF5 filename] [Dataset Name]


Note that you can always use h5ls to learn about the contents
of an HDF5 file.

