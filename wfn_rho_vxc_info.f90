!=============================================================================
!
! Utilities:
!
! (1) wfn_rho_vxc_info    Originally By DAS      Last Modified 12/12/2011 (DAS)
!
!     Prints the contents of the header of a WFN, RHO, or VXC file in
!     a human-readable format.
!
!==============================================================================

#include "f_defs.h"

program wfn_rho_vxc_info

  use global_m
  use wfn_rho_vxc_io_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  character*3 :: sheader
  character*32 :: sdate, stime
  character*256 :: infile, usage
  integer :: iflavor, nargs, ii, iat, ik, isym, is

  usage = 'Usage: wfn_rho_vxc_info.x wfn'

! Get file names from command-line arguments

  nargs = iargc()

  if (nargs .ne. 1) then
    call die(usage)
  endif

  call getarg(1, infile)
  call open_file(unit=7,file=TRUNC(infile),form='unformatted',status='old')

  sheader = 'GET'
  iflavor = -1
  call read_binary_header_type(7, sheader, iflavor, kp, gvec, syms, crys, warn = .false., &
    dont_warn_kgrid = .true., sdate = sdate, stime = stime)

  write(6,'(a)') '====== GENERAL =====' 
  write(6,'(a,a)')  'Type: ', sheader
  if(iflavor == 1) then
    write(6,'(a)') 'Flavor: real'
  else
    write(6,'(a)') 'Flavor: complex'
  endif
  write(6,'(a,a)')  'Date created: ', sdate
  write(6,'(a,a)')  'Time created: ', stime
  write(6,'(a,i1)') 'Number of spins: ', kp%nspin

  write(6,'(a)') '====== G-VECTORS =====' 
  write(6,'(a,i8)') 'Number of G-vectors: ', gvec%ng
  write(6,'(a,f12.6)') 'Charge density cutoff: ', gvec%ecutrho
  write(6,'(a,3i8)') 'FFT grid: ', gvec%FFTgrid(1:3)
  if(sheader == 'WFN') then
    write(6,'(a,i8)') 'Max number of wfn G-vectors: ', kp%ngkmax
    write(6,'(a,f12.6)') 'Wavefunction cutoff: ', kp%ecutwfc
  endif

  write(6,'(a)') '====== ATOMS =====' 
  write(6,'(a,i6)') 'Number of atoms: ', crys%nat
  write(6,'(2a10,a25)') 'Index', 'Species', 'Coordinates'
  do iat = 1, crys%nat
    write(6,'(2i10,3f12.6)') iat, crys%atyp(iat), crys%apos(1:3, iat)
  enddo

  write(6,'(a)') '====== LATTICE =====' 
  write(6,'(a,f12.6)') 'Cell volume (real space): ', crys%celvol
  write(6,'(a,f12.6)') 'Lattice constant (real space): ', crys%alat
  write(6,'(a)') 'Lattice vectors (real space):'
  do ii = 1, 3
    write(6,'(10x,3f12.6)') crys%avec(1:3, ii)
  enddo
  write(6,'(a)') 'Metric (real space):'
  do ii = 1, 3
    write(6,'(10x,3f12.6)') crys%adot(1:3, ii)
  enddo

  write(6,'(a,f12.6)') 'Cell volume (reciprocal space): ', crys%recvol
  write(6,'(a,f12.6)') 'Lattice constant (reciprocal space): ', crys%blat
  write(6,'(a)') 'Lattice vectors (reciprocal space):'
  do ii = 1, 3
    write(6,'(10x,3f12.6)') crys%bvec(1:3, ii)
  enddo
  write(6,'(a)') 'Metric (reciprocal space):'
  do ii = 1, 3
    write(6,'(10x,3f12.6)') crys%bdot(1:3, ii)
  enddo

  write(6,'(a)') '====== SYMMETRIES =====' 
  write(6,'(a,i2)') 'Number of symmetries: ', syms%ntran
  if(syms%cell_symmetry == 0) then
    write(6,'(a)') 'Symmetry type: cubic'
  else
    write(6,'(a)') 'Symmetry type: hexagonal'
  endif
  write(6,'(a7,a31,12x,a33)') 'Index', 'Rotation matrix', 'Fractional translations'
  do isym = 1, syms%ntran
    write(6,'(i5,1x,a,2x,3(3i4,2x),3f12.6)') isym, ':', syms%mtrx(1:3, 1:3, isym), syms%tnp(1:3, isym)
  enddo

  if(sheader == 'WFN') then
    write(6,'(a)') '====== K-POINTS =====' 
    write(6,'(a,i8)') 'Number of k-points: ', kp%nrk
    write(6,'(a,i8)') 'Number of bands: ', kp%mnband
    write(6,'(a,3i4)') 'k-grid: ', kp%kgrid(1:3)
    write(6,'(a,3f12.6)') 'k-shifts: ', kp%shift(1:3)
    write(6,'(a)') '[ifmin = lowest occupied band, ifmax = highest occupied band, for each spin]'
    write(6,'(a8,a25,11x,a12,a22)',advance='no') 'Index', 'Coordinates', 'Weight', 'Number of G-vectors'
    if(kp%nspin == 1) then
      write(6,'(2a10)') 'ifmin', 'ifmax'
    else
      write(6,'(4a10)') 'ifmin1', 'ifmax1', 'ifmin2', 'ifmax2'
    endif
    do ik = 1, kp%nrk
      write(6,'(i8,4f12.6,i22,4i10)') ik, kp%rk(1:3, ik), kp%w(ik), kp%ngk(ik), &
        (kp%ifmin(ik, is), kp%ifmax(ik, is), is = 1, kp%nspin)
    enddo

    write(6,'(a)') '====== ENERGIES/OCCUPATIONS ====='
    do is = 1, kp%nspin
      write(6,'(a,i2)') 'Spin ', is
      do ik = 1, kp%nrk
        write(6,'(a,i6)') 'k-point ', ik
        write(6,'(9999999f12.6)') kp%el(:, ik, is)
        write(6,'(9999999f12.6)') kp%occ(:, ik, is)
      enddo
    enddo
  endif

  call dealloc_header_type(sheader, crys, kp)
  call close_file(7)

end program wfn_rho_vxc_info
