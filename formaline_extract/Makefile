include ../make.inc

all: extract_from_h5 

extract_from_h5: extract_from_h5.cc
	$(CXX) -O0 -g $(LDFLAGS) -o extract_from_h5 \
	extract_from_h5.cc $(HDF5INCS) $(HDF5LIBS)

