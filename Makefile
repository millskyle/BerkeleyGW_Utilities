## Makefile (D. Strubbe Oct 2010)
##

PREFIX=../..
include $(PREFIX)/Common/common-rules.mk

INTERNAL=
#BEGIN_INTERNAL_ONLY
INTERNAL= analyzebz.x wfnreduce.x
ifeq ($(findstring -DHDF5,$(MATHFLAG)),-DHDF5) 
  INTERNAL += wfn2hdf.x hdf2wfn.x
endif 

#END_INTERNAL_ONLY

default: mf_convert.x convert_old_to_new.x degeneracy_check.x wfnmerge.x wfn_rho_vxc_info.x scissors2eqp.x wfn_dotproduct.x fix_occ.x $(INTERNAL)
all: default

mf_convert.x: $(GLOBALOBJS) mf_convert.o $(COMMON)/check_inversion.o $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/sort.o
	$(LINK) $(FOPTS) -o $@ $^
	$(INSTALL_CMD)
	ln -sf $(PWD)/mf_convert_wrapper.sh $(PREFIX)/bin

convert_old_to_new.x: $(GLOBALOBJS) convert_old_to_new.o $(COMMON)/check_inversion.o $(COMMON)/blas.o \
                      $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/input_utils.o $(COMMON)/splines.o $(COMMON)/sort.o
	$(LINK) $(FOPTS) -o $@ $^ $(LAPACKLIB)
	$(INSTALL_CMD)

wfnmerge.x: $(GLOBALOBJS) wfnmerge.o $(COMMON)/check_inversion.o $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/sort.o
	$(LINK) $(FOPTS) -o $@ $^
	$(INSTALL_CMD)

degeneracy_check.x: $(GLOBALOBJS) degeneracy_check.o $(COMMON)/check_inversion.o $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/sort.o
	$(LINK) $(FOPTS) -o $@ $^
	$(INSTALL_CMD)

wfn_rho_vxc_info.x: $(GLOBALOBJS) wfn_rho_vxc_info.o $(COMMON)/check_inversion.o $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/sort.o
	$(LINK) $(FOPTS) -o $@ $^
	$(INSTALL_CMD)

scissors2eqp.x: $(GLOBALOBJS) scissors2eqp.o $(COMMON)/check_inversion.o $(COMMON)/wfn_rho_vxc_io.o \
                $(COMMON)/sort.o $(COMMON)/scissors.o $(COMMON)/splines.o
	$(LINK) $(FOPTS) -o $@ $^
	$(INSTALL_CMD)

wfn_dotproduct.x: $(GLOBALOBJS) wfn_dotproduct.o $(COMMON)/check_inversion.o $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/sort.o \
                  $(COMMON)/blas.o $(COMMON)/find_kpt_match.o $(COMMON)/misc.o $(COMMON)/input_utils.o $(COMMON)/gmap.o
	$(LINK) $(FOPTS) -o $@ $^ $(LAPACKLIB)
	$(INSTALL_CMD)

fix_occ.x: $(GLOBALOBJS) fix_occ.o fix_occ_cp.o $(COMMON)/check_inversion.o $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/sort.o
	$(LINK) $(FOPTS) -o $@ $^
	$(INSTALL_CMD)

#BEGIN_INTERNAL_ONLY
analyzebz.x: $(GLOBALOBJS) analyzebz.o $(COMMON)/blas.o $(COMMON)/misc.o $(COMMON)/check_inversion.o $(COMMON)/wfn_rho_vxc_io.o \
             $(COMMON)/sort.o $(COMMON)/fullbz.o $(COMMON)/subgrp.o $(COMMON)/irrbz.o
	$(LINK) $(FOPTS) -o $@ $^ $(LAPACKLIB)
	$(INSTALL_CMD)

wfnreduce.x: $(GLOBALOBJS) wfnreduce.o $(COMMON)/fftw.o $(COMMON)/check_inversion.o \
	$(COMMON)/wfn_rho_vxc_io.o $(COMMON)/blas.o $(COMMON)/sort.o
	$(LINK) $(FOPTS) -o $@ $^ $(FFTWLIB) $(LAPACKLIB)
	$(INSTALL_CMD)

ifeq ($(findstring -DHDF5,$(MATHFLAG)),-DHDF5)
wfn2hdf.x: $(GLOBALOBJS) wfn2hdf.o $(COMMON)/hdf5_io.o $(COMMON)/wfn_io_hdf5.o $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/check_inversion.o $(COMMON)/sort.o
	$(LINK) $(FOPTS) -o $@ $^ $(HDF5LIB)
	$(INSTALL_CMD)

hdf2wfn.x: $(GLOBALOBJS) hdf2wfn.o $(COMMON)/hdf5_io.o $(COMMON)/wfn_io_hdf5.o $(COMMON)/wfn_rho_vxc_io.o $(COMMON)/check_inversion.o $(COMMON)/sort.o
	$(LINK) $(FOPTS) -o $@ $^ $(HDF5LIB)
	$(INSTALL_CMD)
endif
#END_INTERNAL_ONLY

# dependencies
fix_occ.o: $(COMMON)/wfn_rho_vxc_io_m.mod $(COMMON)/sort_m.mod
#BEGIN_INTERNAL_ONLY
analyzebz.o : $(COMMON)/misc_m.mod $(COMMON)/fullbz_m.mod $(COMMON)/irrbz_m.mod
wfnreduce.o: $(COMMON)/wfn_rho_vxc_io_m.mod $(COMMON)/fftw_m.mod $(COMMON)/blas_m.mod
wfn2hdf.o: $(COMMON)/hdf5_io_m.mod $(COMMON)/wfn_io_hdf5_m.mod $(COMMON)/wfn_rho_vxc_io_m.mod
hdf2wfn.o: $(COMMON)/hdf5_io_m.mod $(COMMON)/wfn_io_hdf5_m.mod $(COMMON)/wfn_rho_vxc_io_m.mod
#END_INTERNAL_ONLY
mf_convert.o convert_old_to_new.o degeneracy_check.o analyzebz.o wfn_rho_vxc_info.o scissors2eqp.o : $(COMMON)/wfn_rho_vxc_io_m.mod
scissors2eqp.o : $(COMMON)/scissors_m.mod
convert_old_to_new.o : $(COMMON)/input_utils_m.mod
wfnmerge.o: $(COMMON)/wfn_rho_vxc_io_m.mod
wfn_dotproduct.o: $(COMMON)/find_kpt_match_m.mod $(COMMON)/misc_m.mod $(COMMON)/blas_m.mod $(COMMON)/input_utils_m.mod $(COMMON)/gmap_m.mod
