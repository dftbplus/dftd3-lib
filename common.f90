module common_module
  implicit none

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

end module common_module
