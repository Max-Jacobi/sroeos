#    This file is part of SRO_EOS.
#
#    SRO_EOS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    SRO_EOS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with SRO_EOS.  If not, see <http://www.gnu.org/licenses/>.
#
MAKEFILE = Makefile

include ../../make.inc

LABLAS=$(LAPACKLIBS) $(BLASLIBS)
HDF5  =$(HDF5INCS) $(HDF5LIBS)

OBJECTS=obj/kind_types.o \
	obj/physical_constants.o \
	obj/table_sizes.o \
	obj/h5_writer.o \
	obj/fdfunctions.o \
	obj/functions.o \
	obj/vector_eos_module.o \
	obj/eleos_variables.o \
	obj/ls_modified_timmes.o \
	obj/electron_eos.o \
	obj/linterp.o \
	obj/linterp_many.o \
	obj/bisection.o \
	obj/find.o \
	obj/findtemp.o \
	obj/findrho.o \
	obj/input_filenames.o \
	obj/phase_space.o \
	obj/tables.o \
	obj/transition.o \
	obj/read_input.o \
	obj/readtable.o \
	obj/get_eos.o \
	obj/write_table.o \
	obj/sna_parameters.o \
	obj/nse_parameters.o \
	obj/ele_parameters.o \
	obj/merge_sna_parameters.o \
	obj/merge_nse_parameters.o \
	obj/merge_ele_parameters.o \
	obj/sna_tab.o \
	obj/nse_tab.o \
	obj/merge_tab.o \
	obj/merge_tables.o \
	obj/include_src_input.o

all: dummy merge

dummy:
	@echo "\n Make code to merge EOS tables\n"
	@echo "\n"

version:
	@ ((echo "character(80) :: git_version  = '";git rev-parse --short HEAD; echo "'") | tr -d '\n'; echo; \
	(echo "character(80) :: version_date = '";date        ; echo "'") | tr -d '\n'; echo;  \
	(echo "integer       :: nchanged = ";git status --short -uno|wc -l) |tr -d '\n'; echo ) > version.inc


merge_eos:  dummy version ../src.tar.gz
	mkdir -p obj mod
	make -f $(MAKEFILE) target "TARGET=merge_eos" "FC=$(F95)" \
		"FEATURES=" \
		"MODOBJ=$(OBJECTS) obj/merge_eos.o" \
		"LIBS=$(LDFLAGS) $(OPTION) $(HDF5) -lstdc++"


target: $(MODOBJ)
	$(FC) -o ../$(TARGET) $(MODOBJ) $(LIBS) $(MODDIR)

../src.tar.gz: dummy
# This target is for creating the source tar ball that we will store
# in the output hdf5 file. If you add subdirectories / source files,
# make sure to update this!
	cd .. && git rev-parse HEAD > gitinfo \
	&& cp ../make.inc . \
	&& tar -czvhf src.tar.gz src/*.F90 src/read/*.f90 \
	src/read/*.F90 src/read_input/*.F90 \
	src/eleos/src/*.f90 src/eleos/src/*.dek src/*.cc src/Makefile \
	make.inc Makefile README.md gitinfo src/version.inc \
	&& rm make.inc

obj/%.o: %.F90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o: read/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o: read/%.F90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o: read_input/%.F90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/kind_types.o: read/kind_types.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/physical_constants.o: read/physical_constants.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/h5_writer.o: ../src/h5_writer.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/write_table.o: write_table.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/phase_space.o : read_input/phase_space.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/readtable.o : read/readtable.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/%.o: eleos/src/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/merge_tables.o: merge_tables.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/merge_eos.o: merge_eos.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/include_src_input.o: include_src_input.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(HDF5)


# uninstall using "make clean"
clean:
	rm -rf ../merge_eos
	rm -rf obj/*.o
	rm -rf mod/*.mod
