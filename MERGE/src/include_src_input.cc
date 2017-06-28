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
				   char* store_filename,
				   char* dataset_name) {

  herr_t status;

  // open file
  hid_t file_id = H5Fopen(h5filename, H5F_ACC_RDWR, H5P_DEFAULT);

  add_file_to_h5(file_id, store_filename, dataset_name);

  status = H5Fclose(file_id);
  assert (status >= 0);

  return;
}

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




void copy_dataset(hid_t outfile_id, hid_t file_id,
		  char* dataset_name) {

  opstruct opdata;
  opdata.mydata = 0;
  opdata.want_name = dataset_name;

  H5Giterate(file_id, "/", NULL, check_if_there, &opdata);

  // do we have the dataset we need?
  if(opdata.mydata) {

    // get dataset
    hid_t dataset = H5Dopen(file_id,dataset_name,H5P_DEFAULT);
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

    // read the data into buffer
    herr_t status = H5Dread (dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			     &buffer[0]);
    assert (status >= 0);

    // close stuff   
    status = H5Dclose(dataset);
    assert (status >= 0);

    // now copy to outfile_id

    // create output dataset, use previous datatype and dataspace
    dataset = H5Dcreate(outfile_id, dataset_name, datatype, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert (dataset >= 0);

    // write dataset
    status = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		      &buffer[0]);
    assert (status >= 0);

    // close stuff
    status = H5Sclose (dataspace);
    assert (status >= 0);

    status = H5Tclose(datatype);
    assert (status >= 0);

    status = H5Dclose (dataset);
    assert (status >= 0);

  } else {
    fprintf(stderr,"PROBLEM: source/input dataset %s not available!\n",
	    dataset_name);
    fprintf(stderr,"PROBLEM: Check input file(s)! Aborting merge!\n");

    abort();
  }

  return;
}




extern "C" void copy_src_input_(char* nse_h5filename,
				char* sna_h5filename,
				char* merged_h5filename,
				int* have_nse,
				int* have_sna) { 

  herr_t status;

  // open outfile
  hid_t outfile_id = H5Fopen(merged_h5filename, H5F_ACC_RDWR, H5P_DEFAULT);

  // see if we have NSE. In this case do NSE first.
  if(*have_nse) {
    // open file
    hid_t file_id = H5Fopen(nse_h5filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert (file_id);

    copy_dataset(outfile_id,file_id,"NSE-src.tar.gz");
    copy_dataset(outfile_id,file_id,"NSE-space.in");
    copy_dataset(outfile_id,file_id,"NSE-partition.in");
    copy_dataset(outfile_id,file_id,"NSE-isotope-properties.in");
    copy_dataset(outfile_id,file_id,"NSE-isotope-list.in");
    copy_dataset(outfile_id,file_id,"NSE-output-list.in");

    status = H5Fclose(file_id);
    assert (status >= 0);
  }

  // now SNA
  if(*have_sna) {
    hid_t file_id = H5Fopen(sna_h5filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert (file_id);
    copy_dataset(outfile_id,file_id,"SNA-src.tar.gz");
    copy_dataset(outfile_id,file_id,"SNA-space.in");
    copy_dataset(outfile_id,file_id,"SNA-skyrme.in");
    status = H5Fclose(file_id);
    assert (status >= 0);
  }
  

  status = H5Fclose(outfile_id);
  assert (status >= 0);

  return;
}
