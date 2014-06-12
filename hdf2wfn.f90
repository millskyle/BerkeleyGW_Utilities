!===============================================================================
!
! Utilities:
!
! (1) hdf2wfn         Originally By JIM      Last Modified 4/23/2012 (JIM)
!
!     Converts HDF5 WFN to ascii
!
!===============================================================================

#include "f_defs.h"

program hdf2wfn
  use global_m
  use hdf5
  use wfn_io_hdf5_m
  use wfn_rho_vxc_io_m
  
  implicit none
  
  type(crystal) :: crys
  type(symmetry) :: syms
  type(kpoints) :: kp
  type(gspace) :: gvec
  character*3 :: sheader = 'WFN'
  character*256 :: infile = 'WFN.h5'
  character*256 :: outfile = 'WFN.h5_ascii'
  integer :: ik, ib, iflavor
  real(DP), pointer :: dwfn(:,:)
  complex(DPC), pointer :: zwfn(:,:)
  real(DP), pointer :: wfn(:,:,:)
  integer, pointer :: gvec_all(:,:)
  integer :: ioffsetk
  integer :: error
  
!-------------------------------
! JIM: Open HDF interface
#ifdef HDF5
  call h5open_f(error)
#endif
  
  write(*,*) 'Entering hdf2wfn...'
  
  call open_file(unit=8,file=TRUNC(outfile),form='formatted',status='replace')
  
  write(*,*) '1. Reading/Writing header and gvectors from WFN.h5 to WFN.h5_ascii'
  
  call read_hdf5_header_type(infile, sheader, iflavor, kp, gvec, syms, crys)
  call write_format_header_type(8, sheader, iflavor, kp, gvec, syms, crys)
  
  SAFE_ALLOCATE(gvec%components, (3,gvec%ng))
  SAFE_ALLOCATE(gvec_all, (3,sum(kp%ngk)))
  
  if (iflavor == 1) then
    SAFE_ALLOCATE(dwfn, (kp%ngkmax, kp%nspin*kp%nspinor))
    SAFE_ALLOCATE(wfn, (1, sum(kp%ngk), kp%mnband))
    
  else
    SAFE_ALLOCATE(zwfn, (kp%ngkmax, kp%nspin*kp%nspinor))
    SAFE_ALLOCATE(wfn, (1, sum(kp%ngk), 1))
  endif
  
  call read_hdf5_gvectors(infile, gvec%ng, gvec%ng, gvec%components)
  call write_format_gvectors(8, gvec%ng, gvec%ng, gvec%components)
  
  write(*,*) '2. Reading/Writing header wavefunctions from WFN.h5 to WFN.h5_ascii'
  
  call read_hdf5_wfn_gvectors(gvec_all, sum(kp%ngk))
  
  ioffsetk = 0
  do ik = 1, kp%nrk
    ! write(*,'(4X,A,I6)') 'kpt', ik
    gvec%components(:,:kp%ngk(ik)) = gvec_all(:,ioffsetk+1:ioffsetk+kp%ngk(ik))
    call write_format_gvectors(8, kp%ngk(ik), kp%ngkmax, gvec%components)
    do ib = 1, kp%mnband
      ! write(*,'(8X,A,I6)') 'band', ib
      if (iflavor == 1) then
        call read_hdf5_band_real(dwfn(1:kp%ngk(ik),:), kp%ngk(ik), ioffsetk, ib-1, 1)
        call write_format_real_data(8, kp%ngk(ik), kp%ngkmax, kp%nspin, dwfn)
      else
        call read_hdf5_band_complex(zwfn(1:kp%ngk(ik),:), kp%ngk(ik), ioffsetk, ib-1, 1)
        call write_format_complex_data(8, kp%ngk(ik), kp%ngkmax, kp%nspin*kp%nspinor, zwfn)
      endif
    enddo
    ioffsetk = ioffsetk + kp%ngk(ik)
  enddo
  
  SAFE_DEALLOCATE_P(gvec%components)
  SAFE_DEALLOCATE_P(gvec_all)
  call dealloc_header_type(sheader, crys, kp)
  
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
  
  call close_file(8)
end program hdf2wfn
