module libdftd3_module
  use dftd3_module
  use copyc6_module
  use param_module
  implicit none
  private

  public :: k1, k2, k3
  public :: esym
  public :: copyc6
  public :: pbcrdatomnumber, rdatomnumber, pbcrdcoord, setfuncpar, rdpar
  public :: pbcncoord, getc6
  public :: checkrcov, edisp, pbcadisp, pbcgdisp, pbcwregrad, outg, pbcedisp
  public :: set_criteria
  public :: setr0ab, elem, printoptions, stoprun
  public :: readl
  public :: rdcoord, ncoord, pbccheckrcov, adisp, gdisp, loadoldpar

end module libdftd3_module
