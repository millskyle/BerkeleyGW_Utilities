!=============================================================================
!
! Utilities:
!
! (1) convert_old_to_new    Originally By DAS
!
!     Takes binary GWR/GWC (and CD95 and VXC files), either real or complex,
!     in format as of MBPT-1.0, plus input file of other info, and writes corresponding
!     binary WFN (and RHO and VXC) files in new format.
!     See example file "convert_old_to_new.inp" in this directory.
!
!==============================================================================

#include "f_defs.h"

program convert_old_to_new

  use global_m
  use input_utils_m
  use wfn_rho_vxc_io_m
  implicit none

  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  character*3 :: sheader
  character*7 :: sflavor
  character*9 :: symgrp
  character*8 :: latunit, posunit
  character*64 :: str0, str1, str2
  character*256 :: progname, inputfile, inwfn, inrho, invxc, &
   outwfn, outrho, outvxc, str
  integer :: ik, is, ib, ig, iat, iflavor, ii, jj, isym
  integer :: nargs, iter
  real(DP) :: dummycelvol, dummypos(3), dummybdot(3,3), efermi
  real(DP), pointer :: dwfn(:,:)
  complex(DPC), pointer :: zwfn(:,:)
  real(DP) :: qtot, omega_plasma, ksum

! Get file names from command-line arguments

  nargs = iargc()

  call getarg(0, progname)

  if (nargs .ne. 8 .and. nargs .ne. 6 .and. nargs .ne. 4) then
    write(6,'(3a)') 'Usage 1: ', TRUNC(progname), &
      ' Real/Complex convert_old_to_new.inp inwfn outwfn'
    write(6,'(3a)') 'Usage 2: ', TRUNC(progname), &
      ' Real/Complex convert_old_to_new.inp inwfn inrho outwfn outrho'
    write(6,'(3a)') 'Usage 3: ', TRUNC(progname), &
      ' Real/Complex convert_old_to_new.inp inwfn inrho invxc outwfn outrho outvxc'
    stop
  endif

  call getarg(1, sflavor)
  call getarg(2, inputfile)
  call getarg(3, inwfn)

  if(nargs .eq. 8) then
    call getarg(4, inrho)
    call getarg(5, invxc)
    call getarg(6, outwfn)
    call getarg(7, outrho)
    call getarg(8, outvxc)
  else if(nargs .eq. 6) then
    call getarg(4, inrho)
    call getarg(5, outwfn)
    call getarg(6, outrho)
  else
    call getarg(4, outwfn)
  endif

  if (TRUNC(sflavor) .eq. RFLAVOR) then
    iflavor = 1
  else if (TRUNC(sflavor) .eq. CFLAVOR) then
    iflavor = 2
  else
    call die('Unknown flavor ' // sflavor // '.')
  endif

  call open_file(unit=10,file=TRUNC(inputfile),form='formatted',status='old')

  call open_file(unit=21,file=TRUNC(outwfn),form='unformatted',status='replace')

  if(nargs > 4) then
    call open_file(unit=12,file=TRUNC(inrho),form='unformatted',status='old')
    call open_file(unit=22,file=TRUNC(outrho),form='unformatted',status='replace')
    if(nargs > 6) then
      call open_file(unit=13,file=TRUNC(invxc),form='unformatted',status='old')
      call open_file(unit=23,file=TRUNC(outvxc),form='unformatted',status='replace')
    endif
  endif

  write(6,'(a)') 'Using ' // TRUNC(sflavor) // ' flavor and input file ' // TRUNC(inputfile) // ':'
  if(nargs .eq. 8) then
    write(6,'(a)') 'Converting old files ' // TRUNC(inwfn) // ', ' // TRUNC(inrho) // ', ' // TRUNC(invxc)
    write(6,'(a)') 'to new files ' // TRUNC(outwfn) // ', ' // TRUNC(outrho) // ', ' // TRUNC(outvxc)
  else
    write(6,'(a)') 'Converting old file ' // TRUNC(inwfn) // ' to new file ' // TRUNC(outwfn)
  endif

  ! Read input file
  
  read(10,*) gvec%ecutrho, kp%ecutwfc
  read(10,*) symgrp
  read(10,*) latunit
  do jj = 1, 3
    read(10,*) (crys%avec(ii, jj), ii = 1, 3)
  enddo
  read(10,*) posunit
  read(10,*) crys%nat
  SAFE_ALLOCATE(crys%atyp, (crys%nat))
  SAFE_ALLOCATE(crys%apos, (3, crys%nat))
  read(10,*) ((crys%apos(ii, iat), ii = 1, 3), crys%atyp(iat), iat = 1, crys%nat)
    
  ! in round 1, we read and allocate.
  ! in round 2, we read and throw away just to back to the right place.
  do iter = 1, 2

    call open_file(unit=11,file=TRUNC(inwfn),form='unformatted',status='old')

    ! Read inwfn
    
    read(11) ((dummybdot(ii,jj),ii=1,3),jj=1,3)
    read(11) dummycelvol
    
    read(11) syms%ntran
    do isym=1,syms%ntran
      read(11) ((syms%mtrx(ii,jj,isym),ii=1,3),jj=1,3)
    enddo
    do isym=1,syms%ntran
      read(11) (syms%tnp(ii,isym),ii=1,3)
    enddo
    
    read(11) kp%nspin
    read(11) kp%nrk
    read(11) (kp%kgrid(ii),ii=1,3)
    read(11) (kp%shift(ii),ii=1,3)
    if(iter == 1) then
      SAFE_ALLOCATE(kp%w, (kp%nrk))
    endif
    read(11) (kp%w(ik),ik=1,kp%nrk)
    if(iter == 1) then
      SAFE_ALLOCATE(kp%rk, (3,kp%nrk))
    endif
    do ik=1,kp%nrk
      read(11) (kp%rk(ii,ik),ii=1,3)
    enddo
    
    read(11) kp%mnband
    if(iter == 1) then
      SAFE_ALLOCATE(kp%ifmin, (kp%nrk,kp%nspin))
      SAFE_ALLOCATE(kp%ifmax, (kp%nrk,kp%nspin))
    endif
    read(11) ((kp%ifmin(ik,is),ik=1,kp%nrk),is=1,kp%nspin)
    read(11) ((kp%ifmax(ik,is),ik=1,kp%nrk),is=1,kp%nspin)
    
    if(iter == 1) then
      SAFE_ALLOCATE(kp%el, (kp%mnband,kp%nrk,kp%nspin))
    endif
    do is=1,kp%nspin
      do ik=1,kp%nrk
        read(11) (kp%el(ib,ik,is),ib=1,kp%mnband)
      enddo
    enddo

    read(11) (gvec%FFTgrid(ii),ii=1,3)
    read(11) gvec%ng
    if(iter == 1) then
      SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
    endif
    do ig=1,gvec%ng
      read(11) (gvec%components(ii,ig), ii=1,3)
    enddo
    
    ! Read numbers of gvectors on each k-point
    ! throwing away other info

    if(iter == 1) then
      SAFE_ALLOCATE(kp%ngk, (kp%nrk))
      do ik=1,kp%nrk
        read(11) kp%ngk(ik)
        do ig=1,kp%ngk(ik)
          read(11)
        enddo
        do ib=1,kp%mnband
          read(11)
        enddo
      enddo
      call close_file(11)
    endif
  enddo

  write(6,'(3x,"Crystal volume:",f32.14)') dummycelvol
  write(6,'(3x,"Number of G-vectors:",i12)') gvec%ng
  write(6,'(3x,"Number of spins:",i16)') kp%nspin
  write(6,'(3x,"Number of bands:",i16)') kp%mnband
  write(6,'(3x,"Number of k-points:",i13)') kp%nrk
  
  ksum = sum(kp%w(1:kp%nrk))
  if(abs(ksum - 1d0) > TOL_Zero) then
    write(0,'(a,f10.6,a)') 'WARNING: k-point weights sum to ', ksum, '; renormalizing.'
    kp%w(1:kp%nrk) = kp%w(1:kp%nrk) / ksum
  endif
  
  ! paratec shifted grids before r307 have
  ! kgrid = 0, shift = 0 for unknown reasons
  ! This is considered illegal by BGW, so we have to reset it
  if(any(kp%kgrid(1:3) <= 0)) then
    write(0,'(a,3i4)') 'WARNING: illegal value in k-grid. k-grid = ', kp%kgrid(1:3)
    write(0,'(a)') 'Reading kgrid from input file'
    read(10,*) (kp%kgrid(ii),ii=1,3)
    if(any(kp%kgrid(1:3) <= 0)) then
      write(0,*) 'read kgrid = ', kp%kgrid(1:3)
      call die("Read illegal kgrid from input file.")
    endif
  endif
  
  call close_file(10)

  ! Compute some other quantities    
  if (symgrp .eq. 'cubic') then
    syms%cell_symmetry = 0
  elseif (symgrp .eq. 'hexagonal') then
    syms%cell_symmetry = 1
  else
    call die('unknown symmetry group')
  endif
  
  if (latunit .eq. 'bohr') then
    continue
  elseif (latunit .eq. 'angstrom') then
    crys%avec(1:3, 1:3) = crys%avec(1:3, 1:3) / BOHR
  else
    call die('unknown units of the lattice vectors')
  endif
  
  if (posunit .eq. 'bohr') then
    continue
  elseif (posunit .eq. 'angstrom') then
    crys%apos(1:3, 1:crys%nat) = crys%apos(1:3, 1:crys%nat) / BOHR
  elseif (posunit .eq. 'crystal') then
    do iat = 1, crys%nat
      dummypos(1:3) = crys%apos(1:3, iat)
      crys%apos(1:3, iat) = 0
      do ii = 1, 3
        do jj = 1, 3
          crys%apos(jj, iat) = crys%apos(jj, iat) + &
            dummypos(ii) * crys%avec(jj, ii)
        enddo
      enddo
    enddo
  else
    call die('unknown units of the atomic positions')
  endif
  
  crys%alat = sqrt(crys%avec(1, 1)**2 + crys%avec(2, 1)**2 + crys%avec(3, 1)**2)
  crys%avec(1:3, 1:3) = crys%avec(1:3, 1:3) / crys%alat
  crys%apos(1:3, 1:crys%nat) = crys%apos(1:3, 1:crys%nat) / crys%alat
  
  crys%celvol = crys%avec(1, 1) * (crys%avec(2, 2) * crys%avec(3, 3) - crys%avec(2, 3) * crys%avec(3, 2)) - &
    crys%avec(2, 1) * (crys%avec(1, 2) * crys%avec(3, 3) - crys%avec(1, 3) * crys%avec(3, 2)) + &
    crys%avec(3, 1) * (crys%avec(1, 2) * crys%avec(2, 3) - crys%avec(1, 3) * crys%avec(2, 2))
  
  crys%blat = 2.0d0 * PI_D / crys%alat
  crys%bvec(1,1) = (crys%avec(2,2) * crys%avec(3,3) - crys%avec(3,2) * crys%avec(2,3)) / crys%celvol
  crys%bvec(2,1) = (crys%avec(3,2) * crys%avec(1,3) - crys%avec(1,2) * crys%avec(3,3)) / crys%celvol
  crys%bvec(3,1) = (crys%avec(1,2) * crys%avec(2,3) - crys%avec(2,2) * crys%avec(1,3)) / crys%celvol
  crys%bvec(1,2) = (crys%avec(2,3) * crys%avec(3,1) - crys%avec(3,3) * crys%avec(2,1)) / crys%celvol
  crys%bvec(2,2) = (crys%avec(3,3) * crys%avec(1,1) - crys%avec(1,3) * crys%avec(3,1)) / crys%celvol
  crys%bvec(3,2) = (crys%avec(1,3) * crys%avec(2,1) - crys%avec(2,3) * crys%avec(1,1)) / crys%celvol
  crys%bvec(1,3) = (crys%avec(2,1) * crys%avec(3,2) - crys%avec(3,1) * crys%avec(2,2)) / crys%celvol
  crys%bvec(2,3) = (crys%avec(3,1) * crys%avec(1,2) - crys%avec(1,1) * crys%avec(3,2)) / crys%celvol
  crys%bvec(3,3) = (crys%avec(1,1) * crys%avec(2,2) - crys%avec(2,1) * crys%avec(1,2)) / crys%celvol
  
  crys%celvol = abs(crys%celvol) * crys%alat**3
  crys%recvol = (2.0d0 * PI_D)**3 / crys%celvol
  
  do ii=1,3
    do jj=1,3
      crys%adot(jj,ii) = dot_product(crys%avec(1:3,jj), crys%avec(1:3,ii)) * crys%alat**2
    enddo
  enddo
  
  do ii=1,3
    do jj=1,3
      crys%bdot(jj,ii) = dot_product(crys%bvec(1:3,jj), crys%bvec(1:3,ii)) * crys%blat**2
    enddo
  enddo
  
  kp%ngkmax = maxval(kp%ngk(1:kp%nrk))
  
  call find_efermi(.true., efermi, 0d0, kp, kp%mnband, 1, "grid", should_search = .true., should_update = .true., write7 = .false.)
  call calc_qtot(kp, crys%celvol, efermi, qtot, omega_plasma, write7 = .false.)
  
  SAFE_ALLOCATE(kp%occ, (kp%mnband, kp%nrk, kp%nspin))
  do is = 1, kp%nspin
    do ik = 1 , kp%nrk
      do ib = 1, kp%mnband
        if(abs(kp%el(ib, ik, is) - efermi) < TOL_Degeneracy) then
          kp%occ(ib, ik, is) = 0.5
        else
          if(ib < kp%ifmin(ik, is) .or. ib > kp%ifmax(ik, is)) then
            kp%occ(ib, ik, is) = 0
          else
            kp%occ(ib, ik, is) = 1
          endif
        endif
      enddo
    enddo
  enddo
  
  ! Check consistency
  
  if (abs(dummycelvol - crys%celvol) .gt. TOL_Small) then
    write(str1,'(f32.14)') dummycelvol
    write(str2,'(f32.14)') crys%celvol
    str = 'unit cell volume read from file (' // TRUNC(str1) // &
      ') differs from unit cell volume computed from lattice vectors (' // &
      TRUNC(str2) // ')'
    call die(str)
  endif
  
  do ii = 1, 3
    do jj = 1, 3
      if (abs(dummybdot(jj, ii) - crys%bdot(jj, ii)) .gt. TOL_Small) then
        write(str0,'("bdot(",i1,",",i1,")")') jj, ii
        write(str1,'(f32.14)') dummybdot(jj, ii)
        write(str2,'(f32.14)') crys%bdot(jj, ii)
        str = TRUNC(str0) // ' read from file (' // TRUNC(str1) // &
          ') differ from ' // TRUNC(str0) // ' computed from lattice vectors (' // &
          TRUNC(str2) // ')'
        call die(str)
      endif
    enddo
  enddo
  
  write(6,'(a)') "Calculations complete."
      
  ! Finally we can write the header and G-vectors for all the files
  ! The G-vectors will be overwritten later, so we need to do it now
  
  sheader = 'WFN'
  call write_binary_header_type(21, sheader, iflavor, kp, gvec, syms, crys)
  call write_binary_gvectors(21, gvec%ng, gvec%ng, gvec%components)
  
  if(nargs > 4) then
    sheader = 'RHO'
    call write_binary_header_type(22, sheader, iflavor, kp, gvec, syms, crys)
    call write_binary_gvectors(22, gvec%ng, gvec%ng, gvec%components)
    
    if(nargs > 6) then
      sheader = 'VXC'
      call write_binary_header_type(23, sheader, iflavor, kp, gvec, syms, crys)
      call write_binary_gvectors(23, gvec%ng, gvec%ng, gvec%components)
    endif
    
    write(6,'(a)') 'Wrote headers.'
  else
    write(6,'(a)') 'Wrote header.'
  endif

  if(iflavor == 1) then
    SAFE_ALLOCATE(dwfn, (kp%ngkmax, kp%nspin))
  else
    SAFE_ALLOCATE(zwfn, (kp%ngkmax, kp%nspin))
  endif
  
  do ik=1,kp%nrk
    read(11) kp%ngk(ik)
    do ig=1,kp%ngk(ik)
      read(11) (gvec%components(ii,ig), ii=1,3)
    enddo
    call write_binary_gvectors(21, kp%ngk(ik), kp%ngk(ik), gvec%components)

    do ib=1,kp%mnband
      if(iflavor == 1) then
        read(11) ((dwfn(ig,is),ig=1,kp%ngk(ik)),is=1,kp%nspin)
        call write_binary_real_data(21, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
      else
        read(11) ((zwfn(ig,is),ig=1,kp%ngk(ik)),is=1,kp%nspin)
        call write_binary_complex_data(21, kp%ngk(ik), kp%ngkmax, kp%nspin, zwfn)
      endif
    enddo
  enddo

  if(iflavor == 1) then
    SAFE_DEALLOCATE_P(dwfn)
  else
    SAFE_DEALLOCATE_P(zwfn)
  endif

  call close_file(11)
  call close_file(21)

  if(nargs > 4) then

    ! CD95 -> RHO

    call rho_vxc_body(12, 22, iflavor, kp, gvec)
    
    write(6,'(a)') 'Wrote charge density.'
    
    call close_file(12)
    call close_file(22)
    
    ! VXC

    if(nargs > 6) then
      call rho_vxc_body(13, 23, iflavor, kp, gvec)
      
      write(6,'(a)') 'Wrote exchange-correlation potential.'

      call close_file(13)
      call close_file(23)
    endif

  endif

  SAFE_DEALLOCATE_P(gvec%components)

  call dealloc_header_type(sheader, crys, kp)

  write(6,'(3x,"Done",/)')

end program convert_old_to_new


!=============================================================================
subroutine rho_vxc_body(iunit_in, iunit_out, iflavor, kp, gvec)

  use global_m
  use wfn_rho_vxc_io_m
  implicit none

  integer, intent(in) :: iunit_in, iunit_out
  integer, intent(in) :: iflavor
  type(kpoints), intent(in) :: kp
  type(gspace), intent(in) :: gvec

  integer :: nproc_para ! number of divisions of file
  integer :: gs_proc    ! number of G-vectors per division
  integer :: ig, igg, is, ii
  real(DP), pointer :: dwfn(:,:)
  complex(DPC), pointer :: zwfn(:,:)

  if(iflavor == 1) then
    SAFE_ALLOCATE(dwfn, (gvec%ng, kp%nspin))
  else
    SAFE_ALLOCATE(zwfn, (gvec%ng, kp%nspin))
  endif

  read(iunit_in) ! bdate, btime
  read(iunit_in) ! nbands, nspin
  read(iunit_in) nproc_para
  igg = 1
  do ii = 1, nproc_para
    read(iunit_in) gs_proc
    do ig = igg, igg + gs_proc - 1
      read(iunit_in) ! gx, gy, gz
      if(iflavor == 1) then
        read(iunit_in) (dwfn(ig,is),is=1,kp%nspin)
      else
        read(iunit_in) (zwfn(ig,is),is=1,kp%nspin)
      endif
    enddo
    igg = igg + gs_proc
  enddo

  if(iflavor == 1) then
    call write_binary_real_data(iunit_out, gvec%ng, gvec%ng, kp%nspin, dwfn)
  else
    call write_binary_complex_data(iunit_out, gvec%ng, gvec%ng, kp%nspin, zwfn)
  endif

  if(iflavor == 1) then
    SAFE_DEALLOCATE_P(dwfn)
  else
    SAFE_DEALLOCATE_P(zwfn)
  endif

  return
end subroutine rho_vxc_body
