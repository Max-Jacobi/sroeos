# This file is part of SRO_EOS.
#
# SRO_EOS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SRO_EOS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with SRO_EOS.  If not, see <http://www.gnu.org/licenses/>.
#
MAKEFILE = Makefile

include ../../make.inc

LABLAS=$(LAPACKLIBS) $(BLASLIBS)
HDF5  =$(HDF5INCS) $(HDF5LIBS)

F90_OBJS=\
	obj/nwprnt.o \
	obj/lautil.o \
	obj/nwutil.o \
	obj/liqrev.o \
	obj/limhpar.o \
	obj/nwtrup.o \
	obj/nwpdlg.o \
	obj/nwddlg.o \
	obj/nwglsh.o \
	obj/nwqlsh.o \
	obj/nwclsh.o \
	obj/nwpure.o \
	obj/liqrup.o \
	obj/nwmhlm.o \
	obj/nwbjac.o \
	obj/nwnjac.o \
	obj/nwbrdn.o \
	obj/nwnwtn.o \
	obj/nwnleq.o \
	obj/global_variables.o \
	obj/phase_space_point.o \
	obj/physical_constants.o \
        obj/print_parameters.o \
	obj/unused.o \
	obj/polynomial_fit.o \
	obj/fermi_integrals.o \
	obj/Levenberg_Marquardt.o \
	obj/LU_decomposition.o \
	obj/eleos_variables.o \
	obj/vector_eos_module.o \
	obj/fdfunctions.o \
	obj/functions.o \
	obj/ls_modified_timmes.o \
	obj/electron_eos.o \
	obj/nuclear_properties.o \
	obj/make_tables.o \
	obj/skyrme_coefficients.o \
	obj/Skyrme_bulk.o \
	obj/Skyrme_bulk_density_derivatives.o \
	obj/Skyrme_bulk_eta_derivatives.o \
	obj/Skyrme_bulk_temperature_derivatives.o \
	obj/Skyrme_bulk_observables.o \
	obj/symmetry_energy.o \
	obj/determine_nuclear_properties.o \
	obj/read_skyrme_coefficients.o \
	obj/surface_properties.o \
	obj/find_critical_point.o \
	obj/fit_critical_temperature.o \
	obj/find_surface_equilibrium.o \
	obj/surface_tension.o \
	obj/minimize_surface_tension.o \
	obj/nuclear_surface_tension.o \
	obj/nuclear_surface_fit.o \
	obj/function_delu_subroutine.o \
	obj/function_beta_subroutine.o \
	obj/function_sigma_subroutine.o \
	obj/function_critical_temperature.o \
	obj/function_h_subroutine.o \
	obj/heavy_nuclei_radius.o \
	obj/heavy_nuclei_mass_number.o \
	obj/heavy_nuclei_chemical_potential.o \
	obj/free_energy_heavy_surface.o \
	obj/free_energy_heavy_translation.o \
	obj/free_energy_heavy.o \
	obj/alpha_properties.o \
	obj/surface_observables.o \
	obj/determine_surface_properties.o \
	obj/uniform_total_derivatives.o \
	obj/non_uniform_total_derivatives.o \
	obj/outside_properties.o \
	obj/heavy_properties.o \
	obj/free_energy.o \
	obj/find_uniform_solution.o \
	obj/find_non_uniform_solution.o \
	obj/print_test.o \
	obj/main_output.o \
	obj/get_output.o \
	obj/test_input.o \
	obj/phase_space.o \
	obj/read_input.o \
	obj/h5_writer.o \
	obj/allocatable_tables.o \
	obj/allocate_tables.o \
	obj/save_to_table.o \
	obj/write_to_table.o \
	obj/solve.o \
	obj/include_src_input.o

dummy:
	@ echo "\nUsing definitions from make.inc \n"

version:
	@ ((echo -n "character(80) :: git_version  = '";git rev-parse --short HEAD; echo "'") | tr -d '\n'; echo; \
	(echo -n "character(80) :: version_date = '";date        ; echo "'") | tr -d '\n'; echo;  \
	(echo -n "integer       :: nchanged = ";git status --short -uno|wc -l) |tr -d '\n'; echo ) > version.inc


test: dummy version
	mkdir -p obj mod
	make -f $(MAKEFILE) target "TARGET=test" "FC=$(F95)" \
		"FEATURES=" \
		"MODOBJ=obj/kind_types.o $(F90_OBJS) obj/test.o" \
		"LIBS=$(OPTION) $(LABLAS) $(HDF5) -lstdc++"

main: dummy version ../src.tar.gz
	mkdir -p obj mod
	make -f $(MAKEFILE) target "TARGET=main" "FC=$(F95)" \
		"FEATURES=" \
		"MODOBJ=obj/kind_types.o $(F90_OBJS) obj/main.o" \
		"LIBS=$(OPTION) $(LABLAS) $(HDF5) -lstdc++"

target: $(MODOBJ)
	$(FC) $(LDFLAGS) -o ../$(TARGET) $(MODOBJ) $(LIBS) $(MODDIR) $(OPENMP)

# obj/%.o : ../../eleos/src/%.f90
# 	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

../src.tar.gz: dummy
# This target is for creating the source tar ball that we will store
# in the output hdf5 file. If you add subdirectories / source files,
# make sure to update this.
	cd .. && git rev-parse HEAD > gitinfo \
	&& cp ../make.inc . && tar -czvhf src.tar.gz \
	src/bulk_matter/*.f90 \
	src/read_input/*.f90 \
	src/nuclear_surface/*.f90 \
	src/heavy_nuclei/*.f90 \
	src/make_table/*.f90 \
	src/make_table/*.cc \
	src/modules/*.f90 \
	src/nleqslv/*.f90 \
	src/eleos/src/*.f90 \
	src/eleos/src/*.dek \
	src/main.f90 \
	src/test.f90 \
	src/Makefile \
	src/version.inc \
	make.inc \
	Makefile \
	README \
	gitinfo \
	&& rm make.inc

obj/kind_types.o : modules/kind_types.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o : bulk_matter/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o : read_input/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o : nuclear_surface/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o : heavy_nuclei/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o : make_table/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o : make_table/make_output.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o : modules/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/%.o : nleqslv/%.f90
	$(FC) $(OPTION) $(LABLAS) -c $< -o $@ $(MODDIR)

obj/%.o : eleos/src/%.f90
	$(FC) $(OPTION) -c $< -o $@ $(MODDIR)

obj/critical_point.o : nuclear_surface/critical_point.f90
	$(FC) $(OPTION) $(LABLAS) -c $< -o $@ $(MODDIR)

obj/h5_writer.o: make_table/h5_writer.f90
	$(FC) $(OPTION) -c $< -o $@ $(HDF5) $(MODDIR)

obj/write_to_table.o: make_table/write_to_table.f90
	$(FC) $(OPTION) -c $< -o $@  $(HDF5) $(MODDIR)

obj/include_src_input.o : make_table/include_src_input.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(HDF5)

obj/polynomial_fit.o : nuclear_surface/polynomial_fit.f90
	$(FC) $(OPTION) $(LABLAS) -c $< -o $@  $(MODDIR)

obj/modules/solve.o : solve.f90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/test.o : test.f90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

obj/main.o : main.f90
	$(FC) $(LIBS) -c $< -o $@ $(MODDIR)

clean:
	rm  -rf obj/*.o mod/*.mod ../main ../test ../src.tar.gz
#
