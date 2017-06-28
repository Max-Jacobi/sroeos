#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <iterator>
#include <vector>
#include <iostream>
#define FALSE 0
#define SRCNAME "SNA-source-tar-gz"


void add_file_to_h5(hid_t file_id, 
		    char const* filename, 
		    char const* dsetname) {

  std::ifstream input(filename, std::ios::binary );
  std::vector<char> buffer((std::istreambuf_iterator<char>(input)),
			   (std::istreambuf_iterator<char>()));

  herr_t status;

  char* cbuf = &buffer[0];
  hsize_t bufsize = buffer.size();

  // create data type
  hid_t datatype = H5Tcreate(H5T_OPAQUE, 1);
  assert(datatype > 0);

  // create data space
  hid_t dataspace = H5Screate_simple(1, &bufsize, NULL);
  assert (dataspace >= 0);

  // create dataset
  hid_t dataset = H5Dcreate(file_id, dsetname, datatype, dataspace,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert (dataset >= 0);

  // write dataset
  status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    cbuf);
  assert (status >= 0);

  // close stuff
  status = H5Sclose (dataspace);
  assert (status >= 0);

  status = H5Tclose(datatype);
  assert (status >= 0);

  status = H5Dclose (dataset);
  assert (status >= 0);

}

extern "C" void include_src_input_(char* h5filename, 
				   char* input_space_filename,
				   char* input_skyrme_filename) {

  herr_t status;

  // open file
  hid_t file_id = H5Fopen(h5filename, H5F_ACC_RDWR, H5P_DEFAULT);

  add_file_to_h5(file_id, "src.tar.gz", "SNA-src.tar.gz");  
  add_file_to_h5(file_id, input_space_filename, "SNA-space.in");
  add_file_to_h5(file_id, input_skyrme_filename, "SNA-skyrme.in");  

  status = H5Fclose(file_id);
  assert (status >= 0);


  return;
}
