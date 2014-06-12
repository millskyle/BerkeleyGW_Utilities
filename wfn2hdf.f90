!===============================================================================
!
! Utilities:
!
! (1) wfn2hdf         Originally By JIM      Last Modified 4/23/2012 (JIM)
!
!     Converts binary WFN to HDF5
!
!===============================================================================

#include "f_defs.h"

program wfn2hdf
  use global_m
  use hdf5
  use wfn_io_hdf5_m
  use wfn_rho_vxc_io_m
  
  implicit none
  
  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  character*3 :: sheader
  character*256 :: infile = 'WFN'
  character*256 :: outfile = 'WFN.h5'
  integer :: ik, ib, iflavor
  real(DP), pointer :: dwfn(:,:)
  complex(DPC), pointer :: zwfn(:,:)
  real(DP), pointer :: wfn(:,:,:)
  integer :: ioffsetk
  integer, allocatable :: wfn_gvec_all(:,:)
  integer :: error
  
!-------------------------------
! JIM: Open HDF interface
#ifdef HDF5
  call h5open_f(error)
#endif
  
  write(*,*) 'Entering wfn2hdf...'
  
  call open_file(unit = 7,file = TRUNC(infile),form = 'unformatted', status = 'old')
  
  write(*,*) '1. Reading header and gvectors from WFN'
  
  sheader = 'GET'
  iflavor = -1
  call read_binary_header_type(7, sheader, iflavor, kp, gvec, syms, crys, dont_warn_kgrid=.true.)
  
  if (sheader .ne. 'WFN') call die("wfn2hdf only workds for binary WFN files")
  
  SAFE_ALLOCATE(gvec%components, (3, gvec%ng))
  SAFE_ALLOCATE(wfn_gvec_all, (3,sum(kp%ngk)))
  
  call read_binary_gvectors(7, gvec%ng, gvec%ng, gvec%components)
  
  if (iflavor == 1) then
    SAFE_ALLOCATE(dwfn, (kp%ngkmax, kp%nspin*kp%nspinor))
    SAFE_ALLOCATE(wfn, (1, kp%ngkmax, kp%nspin*kp%nspinor))
  else
    SAFE_ALLOCATE(zwfn, (kp%ngkmax, kp%nspin*kp%nspinor))
    SAFE_ALLOCATE(wfn, (2, kp%ngkmax, kp%nspin*kp%nspinor))
  endif
  
  call setup_hdf5_wfn_file(TRUNC(outfile), iflavor, kp)
  
  write(*,*) '2. Writing header and gvectors to WFN.h5'
  
  call write_hdf5_header_type(TRUNC(outfile), sheader, iflavor, kp, gvec, syms, crys)
  call write_hdf5_gvectors(TRUNC(outfile), gvec%ng, gvec%ng, gvec%components)
  
  write(*,*) '3. Reading/Writing wavefunctions from WFN to WFN.h5'
  
  ioffsetk = 0
  do ik = 1, kp%nrk
    ! write(*,'(4X,A,I6)') 'kpt', ik
    write(*, *) 'ngk=', kp%ngk(ik)
    call read_binary_gvectors(7, kp%ngk(ik), kp%ngkmax, gvec%components)
    wfn_gvec_all(:, ioffsetk+1:ioffsetk+kp%ngk(ik))=gvec%components(:,1:kp%ngk(ik))
    do ib = 1, kp%mnband
      ! write(*,'(8X,A,I6)') 'band', ib
      if (iflavor == 1) then
        call read_binary_real_data(7, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, dwfn)
        call write_hdf5_band_real(dwfn(1:kp%ngk(ik),:), kp%ngk(ik), ioffsetk, ib-1, 1)
      else
        call read_binary_complex_data(7, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
        call write_hdf5_band_complex(zwfn(1:kp%ngk(ik),:), kp%ngk(ik), ioffsetk, ib-1, 1)
      endif
    enddo
    ioffsetk = ioffsetk + kp%ngk(ik)
  enddo
  call write_hdf5_wfn_gvectors(wfn_gvec_all, sum(kp%ngk))
  
  call dealloc_header_type(sheader, crys, kp)
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE(wfn_gvec_all)
  
  if(iflavor == 1) then
    SAFE_DEALLOCATE_P(dwfn)
  else
    SAFE_DEALLOCATE_P(zwfn)
  endif

!-------------------------------
! JIM: Close HDF interface
#ifdef HDF5
  call h5close_f(error)
#endif
  
  call close_file(7)
end program wfn2hdf
