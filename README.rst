========
LIBDFTD3
========

This is a repackaged version of the `DFTD3 program
<http://www.thch.uni-bonn.de/tc/index.php?section=downloads&subsection=getd3>`_
by S. Grimme and his coworkers.

The original program (V3.1 Rev 1) was downloaded at 2016-04-03. It has been
converted to free format and encapsulated into modules. The source has been
split into two parts:

* A library with the core functionality. This can be directly used by third
  party applications wishing to calculate dispersion with the DFT-D3
  approach.
  
* Additional extensions which are necessary for the command line tool DFTD3 and
  the command line tool itself.


Compilation
===========

Edit the file `make.arch` to reflect your compiler and linker. Then issue
``make``. This will build the following components:

* The library `libdftd3.a` and necessery module files (`*.mod`) in the directory
  `lib/`.

* The command line tool `dftd3` in the directory `prg/`.

* A simple tester for the library (`testapi`( in the directory `test/`. The
  source code of this tester demonstrates how the library can be used.


Credits
=======

When using the library or the DFTD3 tool, please cite:

  S. Grimme, J. Antony, S. Ehrlich and H. Krieg
  J. Chem. Phys, 132 (2010), 154104.
 
If BJ-damping is used 

  S. Grimme, S. Ehrlich and L. Goerigk
  J. Comput. Chem, 32 (2011), 1456-1465.

should be cited as well.


License
=======

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 1, or (at your option) any later version.
