module param_module
  
  real*8 k1,k2,k3

  ! global ad hoc parameters
  parameter (k1=16.0)
  parameter (k2=4./3.) 

  ! reasonable choices are between 3 and 5
  ! this gives smoth curves with maxima around the integer values
  ! k3=3 give for CN=0 a slightly smaller value than computed
  ! for the free atom. This also yields to larger CN for atoms
  ! in larger molecules but with the same chem. environment
  ! which is physically not right
  ! values >5 might lead to bumps in the potential
  parameter (k3=-4.)

end module param_module
