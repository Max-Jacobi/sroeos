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

HDF5  =$(HDF5INCS) $(HDF5LIBS)

SOURCES=obj/kind_types.o \
	obj/physical_constants.o \
	obj/fermi_integrals.o \
	obj/h5_writer.o \
	obj/test_input.o \
	obj/phase_space.o \
	obj/input_filenames.o \
	obj/allocatable_tables.o \
	obj/allocate_tables.o \
	obj/make_tables.o \
	obj/share.o \
	obj/isotopes.o \
	obj/nse_read_nuclear_data.o \
	obj/functions.o \
	obj/read_input.o \
	obj/nse_eos.o \
	obj/solve.o \
	obj/write_to_table.o \
	obj/print_test.o \
	obj/include_src_input.o

OBJS=$(SOURCES:.f90=.o)

all: dummy nse_table nse_test main

dummy:
	@echo "\n Make code to obtain NSE EOS\n"
	@echo "\n"

# nse_table: dummy ../src.tar.gz
# 	mkdir -p obj mod
# 	make -f $(MAKEFILE) target "TARGET=nse_eos" "FC=$(F95)" \
# 		"FEATURES=" \
# 		"MODOBJ=$(OBJS) obj/nse_eos.o obj/nse_table.o" \
# 		"LIBS=$(LDFLAGS) $(OPTION) $(HDF5) -lstdc++"

version:
	@ ((echo "character(80) :: git_version  = '";git rev-parse --short HEAD; echo "'") | tr -d '\n'; echo; \
	(echo "character(80) :: version_date = '";date        ; echo "'") | tr -d '\n'; echo;  \
	(echo "integer       :: nchanged = ";git status --short -uno|wc -l) |tr -d '\n'; echo ) > version.inc

main: dummy version ../src.tar.gz 
	mkdir -p obj mod
	make -f $(MAKEFILE) target "TARGET=main" "FC=$(F95)" \
		"FEATURES=" \
		"MODOBJ=$(OBJS) obj/main.o" \
		"LIBS=$(LDFLAGS) $(OPTION) $(HDF5) -lstdc++"

test: dummy
	mkdir -p obj mod
	make -f $(MAKEFILE) target "TARGET=test" "FC=$(F95)" \
		"FEATURES=" \
		"MODOBJ=$(OBJS) obj/test.o" \
		"LIBS=$(LDFLAGS) $(OPTION) $(HDF5) -lstdc++"

target: $(MODOBJ)
	$(FC) -o ../$(TARGET) $(MODOBJ) $(LIBS) $(MODDIR)

../src.tar.gz: dummy
# This target is for creating the source tar ball that we will store
# in the output hdf5 file. If you add subdirectories / source files,
# make sure to update this.
	cd .. && git rev-parse HEAD > gitinfo \
	&& cp ../make.inc . \
	&& tar -czvhf src.tar.gz src/*.F90 src/modules/*.f90 \
	src/modules/*.F90 src/make_table/*.cc \
	src/make_table/*.F90 src/make_table/*.f90 src/Makefile \
        src/read_input/*.F90 src/version.inc \
	make.inc Makefile README.md gitinfo \
	&& rm make.inc


obj/%.o: %.F90
	$(FC) $(OPTION) $(LIBS) -c $< -o $@ $(MODDIR)

obj/%.o: modules/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o: modules/%.F90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/h5_writer.o: make_table/h5_writer.f90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/write_to_table.o: make_table/write_to_table.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/%.o: make_table/%.F90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o: read_input/%.F90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o: tools/%.F90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/solve.o: modules/solve.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/nse_table.o: nse_table.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/test.o: test.F90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/include_src_input.o : make_table/include_src_input.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(HDF5)



# uninstall using "make clean"
clean:
	rm -rf ../nse_table*
	rm -rf obj/*.o
	rm -rf mod/*.mod
