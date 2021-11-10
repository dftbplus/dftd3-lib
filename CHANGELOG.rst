==========
Change Log
==========


0.10
====

Added
-----

* CMake based build system

Fixed
-----

* Corrected entry in the <r^2>/<r^4> data for Ba


0.9.2
=====

Fixed
-----

* Missing initialisation in the gdisp() routine


0.9.1
=====

Changed
-------

* Minor source code clean up


0.9
===

Added
-----

* Implementation of revised damping parameters as desribed in
  D. G. A. Smith, L. A. Burns, K. Patkowski, and C. D. Sherrill
  J. Phys. Chem. Lett., 2016, 7, pp 2197–2203.
  (Functionality should correspond to V3.2 Rev 0 of original dftd3 code.)

Fixed
-----

* Routine dftd3_pbc_dispersion delivers now correct stress tensor.


0.1
===

Added
-----

* Interface (API) for third party codes.

* Tester for the API.

* Default settings for various compilers.


Changed
-------

* Original source file converted to free format using Metcalf's convert
  tool.

* Code packed into modules.

* Core library separated from dftd3 tools.
