! dftd3 program for computing the dispersion energy and forces from cart
! and atomic numbers as described in
!
! S. Grimme, J. Antony, S. Ehrlich and H. Krieg
! A consistent and accurate ab initio parameterization of density functi
! (DFT-D) for the 94 elements H-Pu
! J. Chem. Phys, 132 (2010), 154104
!
! if BJ-damping is used
! S. Grimme, S. Ehrlich and L. Goerigk, J. Comput. Chem, 32 (2011), 1456
!
! should be cited as well.
!
! Copyright (C) 2009 - 2011 Stefan Grimme, University of Muenster, Germa
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.

! For the GNU General Public License, see <http://www.gnu.org/licenses/>


program dftd3
  use libdftd3_module
  implicit none

  integer maxat,max_elem,maxc
  ! conversion factors
  real*8 autoang,autokcal,c6conv,autoev
  parameter (maxat =20000)
  parameter (max_elem=94)
  ! maximum coordination number references per element
  parameter (maxc =5)
  ! coversion factors
  parameter (autoang =0.52917726d0)
  parameter (autokcal=627.509541d0)
  parameter (autoev = 27.21138505)

  ! DFT-D version
  integer version
  ! number of atoms
  integer n
  ! coordinates in au
  real*8,dimension(:,:), allocatable :: xyz,abc
  ! fixed atoms in geometry opt
  logical fix(maxat)
  ! lattice in au
  real*8 lat(3,3)
  ! gradient
  real*8,dimension(:,:), allocatable :: g
  real*8 g_lat(3,3)
  ! cardinal numbers of elements
  integer,dimension(:), allocatable :: iz
  ! cut-off radii for all element pairs
  real*8 r0ab(max_elem,max_elem)
  ! C6 for all element pairs
  real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
  ! how many different C6 for one element
  integer mxc(max_elem)
  ! C6810
  real*8 c6,c8,c10
  ! coordination numbers of the atoms
  real*8,dimension(:), allocatable :: cn
  ! covalent radii
  real*8 rcov(max_elem)
  ! atomic <r^2>/<r^4> values
  real*8 r2r4(max_elem)
  ! energies
  real*8 e6, e8, e10, e12, disp, e6abc
  ! THE PARAMETERS OF THE METHOD (not all a "free")
  real*8 rs6, rs8, rs10, s6, s18, alp6, alp8, alp10, s42, rs18, alp
  ! printout option
  logical echo
  ! grad ?
  logical grad
  ! analyse results ?
  logical anal
  ! third-order term?
  logical noabc
  ! gradient calctype
  logical numgrad
  ! special parameters
  logical tz
  ! periodic boundary conditions
  logical pbc
  ! repetitions of the unitcell to match the rthr and c_thr
  integer rep_vdw(3),rep_cn(3)
  ! R^2 distance neglect threshold (important for speed in case of large s
  real*8 rthr,rthr2
  ! R^2 distance to cutoff for CN_calculation
  real*8 cn_thr
  ! Integer for assigning only max/min cn C6 (0=normal, 1=min, 2=max)
  ! local and dummy variables
  character*80 atmp,btmp,ctmp,dtmp,etmp,ftmp,func
  integer i,j,k,z,nn,iat,jat,i1,i2
  integer ida(max_elem),ipot
  real*8 x,y,dispr,displ,gdsp,dum,xx(10),dum6(86)
  real*8 dum1,dum2,dum3(3)
  logical ex,pot,fdum
  logical minc6list(max_elem),maxc6list(max_elem),minc6,maxc6

  ! PBE0/def2-QZVP atomic values
  ! data r2r4 /
  ! . 8.0589, 3.4698, 29.0974, 14.8517, 11.8799, 7.8715, 5.5588,
  ! . 4.7566, 3.8025, 3.1036, 26.1552, 17.2304, 17.7210, 12.7442,
  ! . 9.5361, 8.1652, 6.7463, 5.6004, 29.2012, 22.3934, 19.0598,
  ! . 16.8590, 15.4023, 12.5589, 13.4788, 12.2309, 11.2809, 10.5569,
  ! . 10.1428, 9.4907, 13.4606, 10.8544, 8.9386, 8.1350, 7.1251,
  ! . 6.1971, 30.0162, 24.4103, 20.3537, 17.4780, 13.5528, 11.8451,
  ! . 11.0355, 10.1997, 9.5414, 9.0061, 8.6417, 8.9975, 14.0834,
  ! . 11.8333, 10.0179, 9.3844, 8.4110, 7.5152, 32.7622, 27.5708,
  ! . 23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828,
  ! . 15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631,
  ! . 14.5716, 15.8758, 13.8989, 12.4834, 11.4421, 10.2671, 8.3549,
  ! . 7.8496, 7.3278, 7.4820, 13.5124, 11.6554, 10.0959, 9.7340,
  ! . 8.8584, 8.0125, 29.8135, 26.3157, 19.1885, 15.8542, 16.1305,
  ! . 15.6161, 15.1226, 16.1576 /


  ! scale r4/r2 values of the atoms by sqrt(Z)
  ! sqrt is also globally close to optimum
  ! together with the factor 1/2 this yield reasonable
  ! c8 for he, ne and ar. for larger Z, C8 becomes too large
  ! which effectively mimics higher R^n terms neglected due
  ! to stability reasons

  ! r2r4 =sqrt(0.5*r2r4(i)*dfloat(i)**0.5 ) with i=elementnumber
  ! the large number of digits is just to keep the results consistent
  ! with older versions. They should not imply any higher accuracy than
  ! the old values
  data r2r4 / &
      & 2.00734898, 1.56637132, 5.01986934, 3.85379032, 3.64446594, &
      & 3.10492822, 2.71175247, 2.59361680, 2.38825250, 2.21522516, &
      & 6.58585536, 5.46295967, 5.65216669, 4.88284902, 4.29727576, &
      & 4.04108902, 3.72932356, 3.44677275, 7.97762753, 7.07623947, &
      & 6.60844053, 6.28791364, 6.07728703, 5.54643096, 5.80491167, &
      & 5.58415602, 5.41374528, 5.28497229, 5.22592821, 5.09817141, &
      & 6.12149689, 5.54083734, 5.06696878, 4.87005108, 4.59089647, &
      & 4.31176304, 9.55461698, 8.67396077, 7.97210197, 7.43439917, &
      & 6.58711862, 6.19536215, 6.01517290, 5.81623410, 5.65710424, &
      & 5.52640661, 5.44263305, 5.58285373, 7.02081898, 6.46815523, &
      & 5.98089120, 5.81686657, 5.53321815, 5.25477007, 11.02204549, &
      &10.15679528, 9.35167836, 9.06926079, 8.97241155, 8.90092807, &
      & 8.85984840, 8.81736827, 8.79317710, 7.89969626, 8.80588454, &
      & 8.42439218, 8.54289262, 8.47583370, 8.45090888, 8.47339339, &
      & 7.83525634, 8.20702843, 7.70559063, 7.32755997, 7.03887381, &
      & 6.68978720, 6.05450052, 5.88752022, 5.70661499, 5.78450695, &
      & 7.79780729, 7.26443867, 6.78151984, 6.67883169, 6.39024318, &
      & 6.09527958, 11.79156076, 11.10997644, 9.51377795, 8.67197068, &
      & 8.77140725, 8.65402716, 8.53923501, 8.85024712 /

  ! covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009,
  ! values for metals decreased by 10 %
  ! data rcov/
  ! . 0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67
  ! ., 1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54
  ! ., 1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09
  ! ., 1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39
  ! ., 1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26
  ! ., 1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57
  ! ., 1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53
  ! ., 1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32
  ! ., 1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58
  ! ., 1.52, 1.53, 1.54, 1.55 /

  ! these new data are scaled with k2=4./3. and converted a_0 via
  ! autoang=0.52917726d0
  data rcov/ &
      & 0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, &
      & 1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
      & 3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
      & 2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
      & 3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
      & 2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
      & 2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
      & 2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, &
      & 3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
      & 2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
      & 3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
      & 4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
      & 3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
      & 3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
      & 3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
      & 2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
      & 3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
      & 3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
      & 3.82984466, 3.85504098, 3.88023730, 3.90543362 /

  ! k1-k3


  ! write(*,'(94(F12.8,'',''))')r2r4
  ! stop

  ! scale and convert to au
  ! rcov=k2*rcov/autoang
  ! write(*,'(94(F11.8,'',''))')rcov
  ! stop
  ! do i=1,max_elem
  ! dum =0.5*r2r4(i)*dfloat(i)**0.5
  ! store it as sqrt because the geom. av. is taken
  ! r2r4(i)=sqrt(dum)
  ! end do
  ! init
  echo=.true.
  grad=.false.
  pot =.false.
  anal=.false.
  noabc=.true.
  numgrad=.false.
  tz=.false.
  func=' none (read from parameter file)'
  version=0
  pbc=.false.
  fix=.false.
  minc6=.false.
  maxc6=.false.
  minc6list=.false.
  maxc6list=.false.
  fdum=.false.
  ! Cutoff r^2 thresholds for the gradient in bohr^2.
  ! rthr influences N^2 part of the gradient.
  ! rthr2 influences the N^3 part of the gradient. When using
  ! dftd3 in combination with semi-empirical methods or FFs, and large
  ! (>1000 atoms) systems, rthr2 is crucial for speed:
  ! Recommended values are 20^2 to 25^2 bohr.

  ! UR, SE
  rthr=9000.0d0
  rthr2=1600.0d0
  cn_thr=1600.0d0

  ! J/mol nm^6 - > au
  c6conv=1.d-3/2625.4999d0/((0.052917726d0)**6)

  ! set radii
  ! call rdab('~/.r0ab.dat',autoang,max_elem,r0ab)
  call setr0ab(max_elem,autoang,r0ab)

  ! read C6 file by default from $HOME
  ! btmp='~/.c6ab.dat'
  ! inquire(file=btmp,exist=ex)
  ! Muenster directory as second default
  ! if (.not.ex)btmp='/usr/qc/dftd3/c6ab.dat'
  ! call loadc6(btmp,maxc,max_elem,c6ab,mxc)



  ! get coord filename
  call getarg(1,etmp)
  inquire(file=etmp,exist=ex)
  if (.not.ex) call printoptions
  ex=.false.
  ipot=0
  ! options
  do i=1,iargc()
    call getarg(i,ftmp)
    if (index(ftmp,'-h') .ne.0) call printoptions
    if (index(ftmp,'-grad' ).ne.0) grad=.true.
    if (index(ftmp,'-anal' ).ne.0) anal=.true.
    if (index(ftmp,'-noprint').ne.0) echo=.false.
    if (index(ftmp,'-abc' ).ne.0) noabc=.false.
    if (index(ftmp,'-pbc' ).ne.0) pbc=.true.
    if (index(ftmp,'-num' ).ne.0) numgrad=.true.
    if (index(ftmp,'-tz') .ne.0) tz=.true.
    if (index(ftmp,'-old') .ne.0) version=2
    if (index(ftmp,'-zero') .ne.0) version=3
    if (index(ftmp,'-bj') .ne.0) version=4
    if (index(ftmp,'-min') .ne.0) then
      minc6=.true.
      j=0
      do
        call getarg(i+j+1,atmp)
        if (index(atmp,'-').eq.0.and.atmp.ne.'') then
          call elem(atmp,nn)
          if (nn.gt.max_elem.or.nn.lt.1) &
              & call stoprun('Could not recognize min Element')
          minc6list(nn)=.true.
          j=j+1
        else
          exit
        end if
      end do
    end if
    if (index(ftmp,'-max') .ne.0) then
      maxc6=.true.
      k=0
      do
        call getarg(i+k+1,atmp)
        if (index(atmp,'-').eq.0.and.atmp.ne.'') then
          call elem(atmp,nn)
          if (nn.gt.max_elem.or.nn.lt.1) &
              & call stoprun('Could not recognize max Element')
          maxc6list(nn)=.true.
          k=k+1
        else
          exit
        end if
      end do
    end if
    if (index(ftmp,'-pot') .ne.0) then
      pot=.true.
      call getarg(i+1,atmp)
      call readl(atmp,xx,nn)
      ipot=idint(xx(1))
    end if
    if (index(ftmp,'-cnthr') .ne.0) then
      call getarg(i+1,atmp)
      call readl(atmp,xx,nn)
      rthr2=xx(1)
      rthr2=rthr2**2
    end if
    if (index(ftmp,'-func') .ne.0) then
      call getarg(i+1,func)
      ex=.true.
    end if



    if (index(ftmp,'-cutoff') .ne.0) then
      call getarg(i+1,atmp)
      call readl(atmp,xx,nn)
      rthr=xx(1)**2
    end if
    ! if (index(ftmp,'-pot') .ne.0) then
  end do

  ! Check command line input


  if (minc6.and.j.lt.1)then
    call stoprun('No Element given for min/max')
  end if
  if (maxc6.and.k.lt.1)then
    call stoprun('No Element given for min/max')
  end if
  do i=1,max_elem

    if (minc6list(i).and.maxc6list(i)) &
        & call stoprun('Unreasonable min/max input!')
    ! if (minc6list(i)) write(*,*)'min:',i
    ! if (maxc6list(i)) write(*,*)'max:',i
  end do
  ! C6 hard-coded (c6ab.dat not used)
  ! this is alternative to loadc6
  call copyc6(btmp,maxc,max_elem,c6ab,mxc, &
      & minc6,minc6list,maxc6,maxc6list)
  cn_thr=rthr2

  ! write(*,*)'CN(P):',c6ab(15,15,mxc(15),1,2)
  ! write(*,*)'mxc(P):',mxc(15)

  if (pbc) then
    call pbcrdatomnumber(etmp,n)
  else
    call rdatomnumber(etmp,n)
  end if
  ! allocations
  allocate(xyz(3,n))
  allocate(g(3,n))
  allocate(iz(n))
  allocate(cn(n))

  ! reading coordinates and cell in VASP.5.2-format
  ! determing repetitions of unitcell
  if (pbc) then
    call pbcrdcoord(etmp,lat,n,xyz,iz,autoang)
    call set_criteria(rthr,lat,dum3)
    rep_vdw=int(dum3)+1
    call set_criteria(cn_thr,lat,dum3)
    rep_cn=int(dum3)+1
    ! write(*,*)'VDW-cutoff:',sqrt(rthr)*autoang
    ! write(*,*)'CN-cutoff :',sqrt(cn_thr)*autoang
    ! write(*,*)'repvdw:',rep_vdw
    ! write(*,*)'repcn :',rep_cn
    !no pbc
  else
    ! read coordinates, either TM or xmol file
    call rdcoord(etmp,n,xyz,iz,fix,fdum)
    !pbc
  end if
  if (n.lt.1) call stoprun( 'no atoms' )
  if (n.gt.maxat) call stoprun( 'too many atoms' )



  ! the analytical E(3) grad is not available yet
  ! if (grad.and.(.not.noabc))numgrad=.true.

  ! set parameters for functionals
  if (ex) then
    call setfuncpar(func,version,tz,s6,rs6,s18,rs18,alp)
  else
    call rdpar (dtmp,version,s6,s18,rs6,rs18,alp)
  end if

  if (echo)then
    write(*,*)' _________________________________'
    write(*,*)'                                  '
    write(*,*)'|         DFTD3 V3.1 Rev 1        |'
    write(*,*)'| S.Grimme, University Bonn       |'
    write(*,*)'|         October  2015           |'
    write(*,*)'|   see dftd3 -h for options      |'
    write(*,*)' _________________________________'
    write(*,*)
    write(*,*)'Please cite DFT-D3 work done with this code as:'
    write(*,*)'S. Grimme, J. Antony, S. Ehrlich and H. Krieg,'
    write(*,*)'J. Chem. Phys. 132 (2010), 154104'
    write(*,*)'If used with BJ-damping cite also'
    write(*,*)'S. Grimme, S. Ehrlich and L. Goerigk,'
    write(*,*)'J. Comput. Chem. 32 (2011), 1456-1465'
    write(*,*)'For DFT-D2 the reference is'
    write(*,*)'S. Grimme, J. Comput. Chem., 27 (2006), 1787-1799'
    write(*,*)
    write(*,*)' files read :     '
    write(*,*)trim(etmp)
    if (.not.ex)write(*,*)trim(dtmp)
  end if

  if (version.lt.2.or.version.gt.4)stop 'inacceptable version number'

  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  ! all calculations start here
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

  ! CNs for output
  if (pbc) then
    call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,cn_thr)
  else
    call ncoord(n,rcov,iz,xyz,cn,cn_thr)
  end if

  if (version.eq.2)then
    if (echo)write(*,'(''loading DFT-D2 parameters ...'')')
    call loadoldpar(autoang,max_elem,maxc,c6ab,r0ab,dum6)
    ! number of CNs for each element
    mxc=1
    !onvert to au
    c6ab=c6ab*c6conv
  end if


  ! which atoms are present? (for printout)
  if (echo)then
    ida=0
    do i=1,n
      ida(iz(i))=ida(iz(i))+1
    end do
    write(*,'(''C6 coefficients used:'')')
    do i=1,94
      if (ida(i).gt.0)then
        write(*,*) mxc(i),' C6 for element ',i
        do j=1,maxc
          if (c6ab(i,i,j,j,1).gt.0)then
            write(*,'(''Z='',i3,'' CN='',F6.3,5x,''C6(AA)='',F8.2)') &
                & i,c6ab(i,i,j,j,2),c6ab(i,i,j,j,1)
          end if
        end do
      end if
    end do
  end if

  ! output
  if (echo) then
              write(*,'(/''#               XYZ [au]  '',12x,&
     &              '' R0(AA) [Ang.]''2x,&
     &              ''CN'',7x,&
     &              ''C6(AA)     C8(AA)   C10(AA) [au] '')&
     &            ')
    x=0
    btmp=''
    do i=1,n
      z=iz(i)
      call getc6(maxc,max_elem,c6ab,mxc,iz(i),iz(i),cn(i),cn(i),c6)
      do j=1,n
        call getc6(maxc,max_elem,c6ab,mxc,iz(i),iz(j),cn(i),cn(j),dum)
        x=x+dum
      end do
      ! compute C8/C10 for output
      c8 =r2r4(iz(i))**2*3.0d0*c6
      c10=(49.0d0/40.0d0)*c8**2/c6
      dum=0.5*autoang*r0ab(z,z)
      if (version.eq.4)dum=rs6*0.5*autoang*sqrt(c8/c6)
      atmp=' '
      if (fix(i)) then
        atmp='f'
        btmp='f'
      end if
      write(*,'(i4,3F10.5,3x,a2,1x,a1,F7.3,2x,F7.3,3F12.1,L)') &
          & i,xyz(1:3,i),esym(z),atmp, &
          & dum,cn(i), &
          & c6,c8,c10
    end do
    write(*,'(/''molecular C6(AA) [au] = '',F12.2)')x
    if (btmp.eq.'f') then
      write(*,*)' '
      write(*,*)'Caution: Some coordinates fixed &
          &in gradient (marked f, see above).'
      write(*,*)' '
    end if
    if (fdum)then
      write(*,*)'Caution: Dummy atoms found and ignored.'
    end if

  end if


  ! testoutput of radii
  ! do i=1,94
  ! call getc6(maxc,max_elem,c6ab,mxc,i,i,0.d0,0.0d0,c6)
  ! c8 =r2r4(i)**2*3.0d0*c6
  ! write(22,*) i, sqrt(c8/c6)
  ! end do
  ! write(22,*)
  ! do i=1,94
  ! write(22,*) i, r0ab(i,i)
  ! end do
  ! stop

  ! for global ad hoc parameters see
  ! k3 in subroutine getc6, k1 and k2 in subroutine ncoord*
  ! fixed or dependent ones:
  rs8 = rs18
  rs10 = rs18
  alp6 = alp
  alp8 = alp+2.
  alp10= alp8+2.
  ! note: if version=4 (Becke-Johnson), a1=rs6 and a2=rs18
  ! and alp* have no meaning

  !*********************************************************************
  !*********************************************************************
  ! testing code
  ! output of C6=f(CN)
  if (pot.and.ipot.gt.100)then
    x=0
    do i=1,100
      call getc6(maxc,max_elem,c6ab,mxc,ipot-100,ipot-100, &
          & x,x,C6)
      write(2,*) x,c6
      x=x+0.05
    end do
    stop
  end if
  ! Edisp pot curve for testing. Caution: C6 is not constant along R!
  if (pot)then
    write(*,*) 'Computing Edisp potential curve for atom ',ipot
    xyz=0
    iz(1)=ipot
    iz(2)=ipot
    n=2
    xyz(3,2)=1.0/autoang
142 if (pbc) then
      call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab, &
          & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc, &
          & e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)
    else
      call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
          & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr, &
          & e6,e8,e10,e12,e6abc)
    end if
    xyz(3,2)=xyz(3,2)+0.02
    disp=-s6*e6-s18*e8
    write(42,*) xyz(3,2)*autoang,disp*autokcal
    write(43,*) xyz(3,2) ,disp*autokcal
    if (pbc) then
      call pbcncoord(n,rcov,iz,xyz,cn,lat,rep_cn,cn_thr)
    else
      call ncoord(n,rcov,iz,xyz,cn,cn_thr)
    end if
    call getc6(maxc,max_elem,c6ab,mxc,iz(1),iz(2),cn(1),cn(2),c6)
    write(2,*)xyz(3,2)*autoang,-autokcal*c6/xyz(3,2)**6
    if (xyz(3,2).lt.20) goto 142
    write(42,*)
    stop 'pot curve done'
  end if
  ! end testing code
  !*********************************************************************
  !*********************************************************************

  ! check if all parameters have been loaded and are resonable
  do iat=1,n-1
    do jat=iat+1,n
      if (r0ab(iz(jat),iz(iat)).lt.0.1) then
        write(*,*) iat,jat,iz(jat),iz(iat)
        call stoprun( 'radius missing' )
      end if
      if (version.eq.2)then
        c6=c6ab(iz(jat),iz(iat),1,1,1)
      else
        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat), &
            & cn(iat),cn(jat),c6)
      end if
      if (c6.lt.1.d-6) then
        write(*,*) iat,jat,iz(jat),iz(iat),cn(iat),cn(jat)
        call stoprun( 'C6 missing' )
      end if
    end do
  end do

  ! sanity check of read coordniates, based on covalnent radii.
  ! Not omnipotent but better than nothing. S.E. 15.09.2011
  ! call checkcn(n,iz,cn,c6ab,max_elem,maxc)
  if (pbc) then
    call pbccheckrcov(n,iz,rcov,xyz,lat)
  else
    call checkrcov(n,iz,rcov,xyz)
  end if

  !ccccccccccccc
  ! energy call
  !ccccccccccccc
  if (pbc) then
    call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
        & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc, &
        & e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)

  else
    call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
        & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr, &
        & e6,e8,e10,e12,e6abc)
  end if

  e6 = e6 *s6

  ! e6abc= e6abc*s6 ! old and wrong !

  e8 = e8 *s18

  disp =-e6-e8-e6abc

  ! e10 has been tested once again with BJ-damping but has no good effect
  ! e10 = e10 *s18
  ! disp =-e6-e8-e10-e6abc

  ! output
  if (echo) then
    if(version.lt.4)then
      write(*,'(/10x,'' DFT-D V'',i1)') version       
    else
      write(*,'(/10x,'' DFT-D V3(BJ)'')') 
    end if
    write(*,'('' DF '',a50)') func          
    write(*,'('' parameters'')') 
    if(version.eq.2)then
      write(*,'('' s6       :'',f10.4)') s6            
      write(*,'('' alpha6   :'',f10.4)') alp6        
    end if
    if(version.eq.3)then
      write(*,'('' s6       :'',f10.4)') s6            
      write(*,'('' s8       :'',f10.4)') s18           
      write(*,'('' rs6      :'',f10.4)') rs6  
      write(*,'('' rs18     :'',f10.4)') rs18          
      write(*,'('' alpha6   :'',f10.4)') alp6        
      write(*,'('' alpha8   :'',f10.4)') alp8           
      write(*,'('' k1-k3    :'',3f10.4)') k1,k2,k3     
    end if
    if(version.eq.4)then
      write(*,'('' s6       :'',F23.15)') s6            
      write(*,'('' s8       :'',F23.15)') s18           
      write(*,'('' a1       :'',F23.15)') rs6           
      write(*,'('' a2       :'',F23.15)') rs18          
      write(*,'('' k1-k3    :'',3F23.15)') k1,k2,k3     
    end if
    write(*,'('' Cutoff   :'',f10.4,'' a.u.'')') sqrt(rthr) !*autoang
    write(*,'('' CN-Cutoff:'',f10.4,'' a.u.'')')sqrt(cn_thr)!*autoang
    !      if (pbc) then
    !        write(*,'('' Rep_vdw  :'',3I3)') rep_vdw
    !      end if
    write(*,*)
    if (pbc) then
      write(*,'('' Edisp /kcal,au,eV:'',f11.4,X,f12.8,X,f11.7)') &
          & disp*autokcal,disp,disp*autoev
    else
      write(*,'('' Edisp /kcal,au:'',f11.4,f12.8)') disp*autokcal,disp
    end if
    write(*,'(/'' E6    /kcal :'',f11.4)')-e6*autokcal
    if(version.gt.2)then
      write(*,'('' E8    /kcal :'',f11.4)')-e8*autokcal 
      !     write(*,'('' E10   /kcal :'',f11.4)')-e10*autokcal 
      if(.not.noabc) &
          & write(*,'('' E6(ABC) "   :'',2f11.6,F16.12)')-e6abc*autokcal 
      write(*,'('' % E8        :'',f6.2)') -e8/disp/0.01         
      if(.not.noabc) &
          & write(*,'('' % E6(ABC)   :'',f6.2)') -e6abc/disp/0.01        
    end if
  end if
  ! this file for tmer2 read tool
  open(unit=1,file='.EDISP')
  write(1,*) disp
  close(1)
  ! open(unit=87,file='.EABC')
  !
  ! write(87,*) -e6abc
  ! close(87)

  !ccccccccccccccccccccccccc
  ! analyse Edisp pair terms
  !ccccccccccccccccccccccccc
  if (anal) then
    if (pbc) then

      call pbcadisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
          & rs6,rs8,rs10,alp6,alp8,alp10,version,autokcal,autoang, &
          & rthr,rep_vdw,cn_thr,rep_cn,s6,s18,disp*autokcal,lat)
    else
      call adisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
          & rs6,rs8,rs10,alp6,alp8,alp10,version,autokcal, &
          & autoang,rthr,cn_thr,s6,s18,disp*autokcal)
      !pbc
    end if
    !anal
  end if

  !ccccccccccccccccccccccccc
  ! gradient
  !ccccccccccccccccccccccccc
  if (grad)then
    g=0.0d0
    call cpu_time(dum1)
    if (pbc) then
      call pbcgdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab, &
          & rcov,s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,numgrad,&
          & version,g,gdsp,x,g_lat,lat,rep_vdw,rep_cn, &
          & rthr,echo,cn_thr)

    else
      call gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
          & s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,rthr, &
          & numgrad,version,echo,g,gdsp,x,rthr2,fix)
    end if
    call cpu_time(dum2)
    if (echo)write(*,'(''gdisp time  '',f6.1)')dum2-dum1


    ! check if gdisp yields same energy as edisp
    dum=abs((disp-gdsp)/disp)
    !if this check gives compiler errors, replace is with a different NaN ch
    if (ISNAN(dum)) call stoprun('internal NaN-error')
    ! if (dum/=dum) call stoprun('internal NaN-error')
    if (dum.gt.1.d-8) then
      write(*,*) disp,gdsp
      call stoprun('internal error')
      !sanitycheck
    end if
    ! write to energy and gradient files in TM style
    if (pbc) then
      if (echo) then
        write(*,*)'Cartesian gradient written to file dftd3_gradient.'
        write(*,*)'Cartesian cellgradient written &
            & to file dftd3_cellgradient. (a.u.)'
        !echo
      end if
      !*autoev
      g_lat=g_lat
      call pbcwregrad(n,g,g_lat)
      !not pbc
    else
      if (echo) then
        write(*,*) 'Cartesian gradient written to file dftd3_gradient'
        !echo
      end if
      call outg(n,g,'dftd3_gradient')
      ! call wregrad(n,xyz,iz,disp,g)

      !pbc
    end if
    !grad
  end if

  if (echo)write(*,*) 'normal termination of dftd3'

  goto 999
  ! test test tesc test test tesc test test tesc test test tesc test test
  ! test test tesc test test tesc test test tesc test test tesc test test
  ! gradient test 6-7 digits should be the same
  if (pbc) then
    call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab, &
        & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc, &
        & e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)
  else
    call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
        & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr, &
        & e6,e8,e10,e12,e6abc)
  end if

  do i=1,n
    do j=1,3
      xyz(j,i)=xyz(j,i)+0.00001
      if (pbc) then
        call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab, &
            & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc, &
            & e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)
      else
        call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
            & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr, &
            & e6,e8,e10,e12,e6abc)
      end if
      e6 = e6 *s6
      e8 = e8 *s18
      dispr =-e6-e8-e6abc
      xyz(j,i)=xyz(j,i)-0.00002
      if (pbc) then
        call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab, &
            & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc, &
            & e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)
      else
        call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov, &
            & rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,cn_thr, &
            & e6,e8,e10,e12,e6abc)
      end if
      e6 = e6 *s6
      e8 = e8 *s18
      displ =-e6-e8-e6abc
      xyz(j,i)=xyz(j,i)+0.00001
      write(*,'(i4,2E14.6)')i,(dispr-displ)/(0.00002),g(j,i)
    end do
  end do

  if (pbc) then
    if (echo) write(*,*)'Doing numerical stresstensor...'
    allocate(abc(3,n))

    call xyz_to_abc(xyz,abc,lat,n)
    dum1=1.d-5
    if (echo) write(*,*)'step: ',dum1
    do i=1,3
      do j=1,3
        lat(j,i)=lat(j,i)+dum1
        call abc_to_xyz(abc,xyz,lat,n)
        !call edisp...dum1
        call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab, &
            & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc, &
            & e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)

        dispr=-s6*e6-s18*e8-e6abc


        lat(j,i)=lat(j,i)-2*dum1
        call abc_to_xyz(abc,xyz,lat,n)
        !call edisp...dum2
        call pbcedisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab, &
            & rcov,rs6,rs8,rs10,alp6,alp8,alp10,version,noabc, &
            & e6,e8,e10,e12,e6abc,lat,rthr,rep_vdw,cn_thr,rep_cn)

        displ=-s6*e6-s18*e8-e6abc
        dum=(dispr-displ)/(dum1*2.0)

        lat(j,i)=lat(j,i)+dum1
        call abc_to_xyz(abc,xyz,lat,n)

        write(*,'("L"2i1,2E14.6)'),i,j,dum,g_lat(j,i)
        !j
      end do
      !i
    end do



    deallocate(abc)
    !pbc
  end if

  ! test test tesc test test tesc test test tesc test test tesc test test
  ! test test tesc test test tesc test test tesc test test tesc test test

999 deallocate(xyz,g,iz,cn)
end program dftd3
