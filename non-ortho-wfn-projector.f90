!=============================================================
! Utility:
! 
! wfn_dotproduct (flavorless). Originally by DAS.
!
! Form the overlap between the bands in two wavefunction files.
! If they are copies of the same file, you can check orthonormality.
! Only bands at corresponding k-points are considered, since
! the overlap is zero by Bloch`s theorem if the k-points differ.
!
! Warning: a known issue is possible errors from gmap when relating
! k-points by symmetry operations.
!
!=============================================================



#include "f_defs.h"

program wfn_dotproduct

  use global_m
  use blas_m
  use find_kpt_match_m
  use gmap_m
  use input_utils_m
  use misc_m
  use sort_m
  use wfn_rho_vxc_io_m
  implicit none

  type(crystal) :: crys1, crys2
  type(symmetry) :: syms1, syms2
  type(kpoints) :: kp1, kp2
  type(gspace) :: gvec1, gvec2, gvec_kpt
  character*3 :: sheader
  character*11 :: inform
  character*256 :: sformat, file_1, file_2, usage
  integer :: ik1, ib1, ik2, ik3, ib2, iflavor, ngkmax, ig, ispin
  integer :: nargs, itqq, kgqq(3), temp
  real(DP), allocatable :: dwfn1(:,:), dwfn2(:,:,:), dph(:)
  complex(DPC), allocatable :: zwfn1(:,:), zwfn2(:,:,:), zph(:)
  integer, allocatable :: gindex1(:), gindex2(:), ind(:), isrt(:), ginv(:)
  logical :: informat
  real(DP) :: doverlap, q_shift(3)
  complex(DPC) :: zoverlap

  usage = 'Usage: wfn_dotproduct.x A|B WFN_1 WFN_2'

  ! Get file names from command-line arguments

  nargs = iargc()

  if (nargs < 3) then
    call die(usage)
  endif

  ! mpiexec from mpich1 may add 4 extra arguments to the list
  if (nargs > 3) then
    write(0,'(a,i3,a)') 'WARNING: ', nargs, ' arguments found, only first 3 being used.'
  endif

  call getarg(1, sformat)
  call getarg(2, file_1)
  call getarg(3, file_2)

  ! Open units

  if (TRUNC(sformat) .eq. 'A') then
    informat = .true.
    inform = 'formatted'
  elseif (TRUNC(sformat) .eq. 'B') then
    informat = .false.
    inform = 'unformatted'
  else  
    call die(usage)
  endif

  ! Read headers

  sheader='WFN'
  iflavor=-1
  write(6,'(3x,a)') 'Reading ' // TRUNC(inform) // ' file ' // TRUNC(file_1)
  call open_file(unit=7,file=TRUNC(file_1),form=inform,status='old')
  call read_header_type(7, informat, sheader, iflavor, kp1, gvec1, syms1, crys1, warn = .false., dont_warn_kgrid = .true.)

  ! must have same flavor as first one
  write(6,'(3x,a)') 'Reading ' // TRUNC(inform) // ' file ' // TRUNC(file_2)
  call open_file(unit=8,file=TRUNC(file_2),form=inform,status='old')
  call read_header_type(8, informat, sheader, iflavor, kp2, gvec2, syms2, crys2, warn = .false., dont_warn_kgrid = .true.)
  call check_header(TRUNC(file_1), kp1, gvec1, syms1, crys1, &
    TRUNC(file_2), kp2, gvec2, syms2, crys2, is_wfn = .true., tolerant = .true.)

  write(6,'(2(a,a,i8,a))') TRUNC(file_1), ' has ', kp1%nrk, ' k-points; ', &
    TRUNC(file_2), ' has ', kp2%nrk, ' k-points'
  if(kp2%nrk > kp1%nrk) then
    write(6,'(a)') "Note: only the second file's k-points will have symmetries applied. Switch the order on the command-line"
    write(6,'(a)') "if you want matches where a symmetry is applied to the first file's k-points."
  endif

  if(all(kp1%kgrid(:) > 0 .and. kp2%kgrid(:) > 0)) then
    q_shift(1:3) = dble(kp2%shift(1:3))/kp2%kgrid(1:3) - dble(kp1%shift(1:3))/kp1%kgrid(1:3)
  else
    q_shift(1:3) = 0d0
  endif
  write(6,'(a,3f12.6)') 'Using q-shift = ', q_shift(1:3)

  ! Read charge density gvectors
  SAFE_ALLOCATE(gvec1%components, (3, gvec1%ng))
  SAFE_ALLOCATE(gvec2%components, (3, gvec2%ng))

  ! use gvec1 as the master list of G-vectors
  call read_gvectors(7, informat, gvec1%ng, gvec1%ng, gvec1%components)
  ! needed to be able to use findvector
  call gvec_index(gvec1)

  ngkmax = max(kp1%ngkmax, kp2%ngkmax)
  if (iflavor .eq. 1) then
    SAFE_ALLOCATE(dwfn1, (gvec1%ng, kp1%nspin))
    SAFE_ALLOCATE(dwfn2, (gvec1%ng, kp2%nspin, kp2%mnband))
    SAFE_ALLOCATE(dph, (gvec1%ng))
  else
    SAFE_ALLOCATE(zwfn1, (gvec2%ng, kp1%nspin))
    SAFE_ALLOCATE(zwfn2, (gvec2%ng, kp2%nspin, kp2%mnband))
    SAFE_ALLOCATE(zph, (gvec1%ng))
  end if
  SAFE_ALLOCATE(gindex1, (ngkmax))
  SAFE_ALLOCATE(gindex2, (ngkmax))
  SAFE_ALLOCATE(ind, (ngkmax))
  SAFE_ALLOCATE(gvec_kpt%components, (3, ngkmax))
  SAFE_ALLOCATE(ginv, (gvec1%ng))
  SAFE_ALLOCATE(gvec1%ekin, (gvec1%ng))
  SAFE_ALLOCATE(isrt, (gvec1%ng))

  if(kp1%nspinor > 1) then
    call die("nspinor > 1 not implemented")
  endif

  do ik1 = 1, kp1%nrk
    write(6,'(/,a,i6,a,3f12.6)') 'k-point ', ik1, ': ', kp1%rk(1:3, ik1)
    call find_kpt_match(kp2, syms2, kp1%rk(1:3, ik1) + q_shift(1:3), ik2, itqq, kgqq)
    if(ik2 == 0) then
      write(6,'(a)') 'No match in ' // TRUNC(file_2)
      call read_gvectors(7, informat, kp1%ngk(ik1), ngkmax, gvec_kpt%components, dont_read = .not. informat)
    else
      write(6,'(a,i6,3a,3f12.6)') 'Matches k-point ', ik2, ' in ', TRUNC(file_2), ':', kp2%rk(1:3, ik2)
      write(6,'(a,i3,a,3i3)') 'with symmetry operation ', itqq, ' and Umklapp ', kgqq
      call read_gvectors(7, informat, kp1%ngk(ik1), ngkmax, gvec_kpt%components)
      write(6,'(2(a10,a16,4x),a10,6x,a20)',advance='no') 'band 1', 'energy 1 ', 'band 2', 'energy 2 ', 'spin ', 'Re overlap '
      if(iflavor == 2) write(6,'(a20)',advance='no') 'Im overlap '
      write(6,'(a20)') '|overlap|^2 '

      do ig = 1, kp1%ngk(ik1)
        call findvector(gindex1(ig), gvec_kpt%components(:, ig), gvec1)
        if (gindex1(ig) ==  0)  call die('could not match G-vector')
      enddo
    endif

    if(ik2 > 0) then
      if(kp1%ngk(ik1) /= kp2%ngk(ik2)) then
        call die("Internal error: ngk mismatch")
      endif

      rewind(8)
      call dealloc_header_type(sheader, crys2, kp2)
      call read_header_type(8, informat, sheader, iflavor, kp2, gvec2, syms2, crys2, warn = .false., dont_warn_kgrid = .true.)
      call check_header(TRUNC(file_1), kp1, gvec1, syms1, crys1, &
        TRUNC(file_2), kp2, gvec2, syms2, crys2, is_wfn = .true., tolerant = .true.)
      call read_gvectors(8, informat, gvec2%ng, gvec2%ng, gvec2%components, dont_read = .not. informat)

      do ik3 = 1, ik2 - 1
        call read_gvectors(8, informat, kp2%ngk(ik3), ngkmax, gvec_kpt%components, dont_read = .not. informat)
        do ib2 = 1, kp2%mnband
          if(iflavor == 1) then
            call read_real_data(8, informat, kp2%ngk(ik3), gvec2%ng, kp2%nspin, &
              dwfn2(:,:,ib2), dont_read = .not. informat)
          else
            call read_complex_data(8, informat, kp2%ngk(ik3), gvec2%ng, kp2%nspin*kp2%nspinor, &
              zwfn2(:,:,ib2), dont_read = .not. informat)
          endif
        enddo
      enddo

      ! now we are at the right point in the file
      call read_gvectors(8, informat, kp2%ngk(ik2), ngkmax, gvec_kpt%components)

      ginv = 0
      do ig = 1, kp2%ngk(ik2)
        call findvector(gindex2(ig), gvec_kpt%components(:, ig), gvec1)
        if (gindex2(ig) ==  0)  call die('could not match G-vector')
        ginv(gindex2(ig)) = ig
      enddo

      call kinetic_energies(gvec1, crys1%bdot, gvec1%ekin, qvec = kp2%rk(:, ik2))
      call sortrx_D(gvec1%ng, gvec1%ekin, isrt, gvec = gvec1%components)

      ! FHJ: Keep the following line. It prevents a compiler bug with sunf90
      ! which changes the value of iflavor.
      iflavor = iflavor
      if(iflavor == 1) then
        call gmap(gvec1, syms1, kp2%ngk(ik2), itqq, kgqq, isrt, ginv, ind, dph, .true.)
      else
        call gmap(gvec1, syms1, kp2%ngk(ik2), itqq, kgqq, isrt, ginv, ind, zph, .true.)
      endif

      if(any(ind(1:kp2%ngk(ik2)) == 0)) call die("ind array from gmap has a zero")

      do ib2 = 1, kp2%mnband
        if(iflavor == 1) then
          call read_real_data(8, informat, kp2%ngk(ik2), gvec2%ng, kp2%nspin, &
            dwfn2(:,:,ib2), gindex = gindex2)
        else
          call read_complex_data(8, informat, kp2%ngk(ik2), gvec2%ng, kp2%nspin*kp2%nspinor, &
            zwfn2(:,:,ib2), gindex = gindex2)
        endif
      enddo

      if(iflavor == 1) then
        do ig = 1, kp2%ngk(ik2)
          dwfn2(ind(ig), :, :) = dwfn2(ind(ig), :, :) * dph(ig)
        enddo
      else
        do ig = 1, kp2%ngk(ik2)
          zwfn2(ind(ig), :, :) = zwfn2(ind(ig), :, :) * zph(ig)
        enddo
      endif
    endif

    do ib1 = 1, kp1%mnband
      if(iflavor == 1) then
        call read_real_data(7, informat, kp1%ngk(ik1), gvec1%ng, kp1%nspin, dwfn1, &
          dont_read = (ik2 == 0) .and. .not. informat, gindex = gindex1)
      else
        call read_complex_data(7, informat, kp1%ngk(ik1), gvec1%ng, kp1%nspin*kp1%nspinor, zwfn1, &
          dont_read = (ik2 == 0) .and. .not. informat, gindex = gindex1)
      endif
      if(ik2 == 0) cycle

      do ib2 = 1, kp2%mnband
        if(iflavor == 1) then
          do ispin = 1, kp1%nspin
            doverlap = ddot(gvec1%ng, dwfn1(:, ispin), 1, dwfn2(:, ispin, ib2), 1)
            write(6,'(2(i8,2x,g20.8),i8,8x,2g20.10)') ib1, kp1%el(ib1, ik1, ispin), ib2, kp2%el(ib2, ik2, ispin), &
              ispin, doverlap, dble(doverlap)**2
          enddo
        else
          do ispin = 1, kp1%nspin
            zoverlap = zdotc(gvec1%ng, zwfn1(:, ispin), 1, zwfn2(:, ispin, ib2), 1)
            write(6,'(2(i8,2x,g20.8),i8,8x,3g20.10)') ib1, kp1%el(ib1, ik1, ispin), ib2, kp2%el(ib2, ik2, ispin), &
              ispin, zoverlap, dble(zoverlap)**2 + aimag(zoverlap)**2
          enddo
        endif
      enddo
    enddo
  enddo

  ! Close files

  call close_file(7)
  call close_file(8)

  ! Deallocate arrays

  if (iflavor .eq. 1) then
    SAFE_DEALLOCATE(dwfn1)
    SAFE_DEALLOCATE(dwfn2)
    SAFE_DEALLOCATE(dph)
  else
    SAFE_DEALLOCATE(zwfn1)
    SAFE_DEALLOCATE(zwfn2)
    SAFE_DEALLOCATE(zph)
  endif

  SAFE_DEALLOCATE_P(gvec1%components)
  SAFE_DEALLOCATE_P(gvec2%components)
  SAFE_DEALLOCATE_P(gvec1%index_vec)
  SAFE_DEALLOCATE_P(gvec_kpt%components)
  SAFE_DEALLOCATE_P(gvec1%ekin)
  SAFE_DEALLOCATE(isrt)
  SAFE_DEALLOCATE(ind)
  SAFE_DEALLOCATE(ginv)
  SAFE_DEALLOCATE(gindex1)
  SAFE_DEALLOCATE(gindex2)
  call dealloc_header_type(sheader, crys1, kp1)
  call dealloc_header_type(sheader, crys2, kp2)

end program wfn_dotproduct
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
