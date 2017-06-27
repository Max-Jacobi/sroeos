#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <cstring>
#include <fstream>
#include <iterator>
#include <vector>
#include <iostream>
#define FALSE 0

struct opstruct {
  char* want_name;
  int mydata;
};


herr_t check_if_there(hid_t loc_id, const char *name, void *opdata)
{
  H5G_stat_t statbuf;
  int *mydata = (int*) opdata;
  opstruct *mystruct = (opstruct*) opdata;

  H5Gget_objinfo(loc_id, name, FALSE, &statbuf);
  switch (statbuf.type) {
  case H5G_GROUP:
    break;
  case H5G_DATASET:
    if (strcmp(mystruct->want_name,name) == 0) {
      mystruct->mydata = 1;
    }
    break;
  case H5G_TYPE:
    break;
  default:
    printf(" Unable to identify an object ");
  }
  return 0;
}


int main( int argc, char *argv[] ) {

  if(argc != 3) {
    std::cout << "usage: " << argv[0] << " [filename] [dataset name]" << std::endl;
    return 1;
  }

  // open file
  herr_t status;
  opstruct opdata;
  opdata.mydata = 0;
  opdata.want_name = argv[2];

  hid_t file_id = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);

  if(file_id < 0) {
    std::cout << "\n" << "ERROR: Could not open HDF5 file!" << std::endl;
    return 1;
  }

  
  // iterate over root group to check if source code already 
  // there
  H5Giterate(file_id, "/", NULL, check_if_there, &opdata);

  if(opdata.mydata) {
    std::cout << "\n The HDF5 file contains the requested dataset! Extracting ..."
	      << std::endl;
  } else {
    std::cout << "\n ERROR: The HDF5 file does not contain a dataset named " << 
      argv[2] << " !" << std::endl;
    std::cout << "Aborting!" << std::endl;
    return 1;
  }

  // open dataset
  hid_t dataset = H5Dopen(file_id,argv[2],H5P_DEFAULT);
  assert (dataset >= 0);

  // get the datatype
  hid_t datatype = H5Dget_type (dataset);
  assert (datatype >= 0);

  // get dataspace
  hid_t dataspace = H5Dget_space(dataset);
  assert (dataspace >= 0);

  // get size of our dataset
  const int ndims = H5Sget_simple_extent_ndims(dataspace);
  hsize_t dims[ndims];
  H5Sget_simple_extent_dims(dataspace, dims, NULL);

  std::vector<char> buffer(dims[0]);

  // read the data
  status = H5Dread (dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    &buffer[0]);
  assert (status >= 0);
  
  // close stuff
  status = H5Sclose (dataspace);
  assert (status >= 0);

  status = H5Tclose(datatype);
  assert (status >= 0);

  status = H5Dclose(dataset);
  assert (status >= 0);

  status = H5Fclose(file_id);
  assert (status >= 0);

  // write the source code out
  std::ofstream output(argv[2], std::ios::binary );
  output.write( (char*)&buffer[0], buffer.size() );
  output.close();

  std::cout << "Dataset written to " << argv[2] << " !" << std::endl;
  std::cout << "Done! :-)" << std::endl;
  
  return 0;
}

