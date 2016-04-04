module dftd3_common
  implicit none

  ! Working precision (double precision)
  integer, parameter :: wp = kind(1.0d0)

  ! global ad hoc parameters
  real(wp), parameter :: k1 = 16.0
  real(wp), parameter :: k2 = 4./3.

  ! reasonable choices are between 3 and 5
  ! this gives smoth curves with maxima around the integer values
  ! k3=3 give for CN=0 a slightly smaller value than computed
  ! for the free atom. This also yields to larger CN for atoms
  ! in larger molecules but with the same chem. environment
  ! which is physically not right
  ! values >5 might lead to bumps in the potential
  real(wp), parameter :: k3 = -4.


  real(wp), parameter :: autoang =0.52917726d0
  real(wp), parameter :: autokcal = 627.509541d0
  real(wp), parameter :: autoev = 27.21138505
  ! J/mol nm^6 - > au
  real(wp), parameter :: c6conv = 1.d-3/2625.4999d0/((0.052917726d0)**6)


contains

    subroutine limit(iat,jat,iadr,jadr)
    integer iat,jat,iadr,jadr,i
    iadr=1
    jadr=1
    i=100
10  if (iat.gt.100) then
      iat=iat-100
      iadr=iadr+1
      goto 10
    end if

    i=100
20  if (jat.gt.100) then
      jat=jat-100
      jadr=jadr+1
      goto 20
    end if

  end subroutine limit

end module dftd3_common
