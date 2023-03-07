#include "read_binary.h"

int read_from_binary(std::string path, std::string prefix, std::string name, std::vector<double>& x, bool verbose) {
	std::fstream file(path + prefix + "_" + name + ".dat", std::ios::binary | std::ios::in | std::ios::ate);
	if (file) {
		std::fstream::pos_type size = file.tellg();
		if (verbose) std::cout << path + prefix + "_" + name + ".dat" + " is open for reading. Size: " << size << " (" << size / sizeof(double) << " numbers)" << std::endl;
		char* memblock = new char[size];
		file.seekg(0, std::ios::beg);
		file.read(memblock, size);
		file.close();

		if (verbose) std::cout << path + prefix + "_" + name + ".dat" + " has been read into memory.";

		double* values = (double*)memblock;

		x = std::vector<double>(values, values + size / sizeof(double));

		delete[] memblock;

		return 0;
	}
	else {
		if (verbose) std::cout << path + prefix + "_" + name + ".dat" + " not found.";
		return -1;
	}
}

int read_from_binary(std::string path, std::string name, std::vector<double>& x, bool verbose) {
	std::fstream file(path + name + ".dat", std::ios::binary | std::ios::in | std::ios::ate);
	if (file) {
		std::fstream::pos_type size = file.tellg();
		if (verbose) std::cout << path + name + ".dat" + " is open for reading. Size: " << size << " (" << size / sizeof(double) << " numbers)" << std::endl;
		char* memblock = new char[size];
		file.seekg(0, std::ios::beg);
		file.read(memblock, size);
		file.close();

		if (verbose) std::cout << path + name + ".dat" + " has been read into memory.";

		double* values = (double*)memblock;

		x = std::vector<double>(values, values + size / sizeof(double));

		delete[] memblock;

		return 0;
	}
	else {
		if (verbose) std::cout << path + name + ".dat" + " not found.";
		return -1;
	}
}

int read_from_binary(std::string path, std::string prefix, std::string name, std::vector<double>& x) {
	return read_from_binary(path, prefix, name, x, false);
}

int read_from_binary(std::string path, std::string name, std::vector<double>& x) {
	return read_from_binary(path, name, x, false);
}

int read_from_hdf5(std::string filename, std::string groupname, std::string fieldname, std::vector<double> &res) {
	if ((fieldname == "density") || (fieldname == "rho")) fieldname = "Density";
	else if (fieldname == "vx") fieldname = "x-velocity";
	else if (fieldname == "vy") fieldname = "y-velocity";
	else if (fieldname == "vz") fieldname = "z-velocity";
	//std::cout << "poop " << filename << std::endl;
	//std::cout << "poop " << groupname << std::endl;
	//std::cout << "poop " << fieldname << std::endl;
	try {
		if ((fieldname == "absv") || (fieldname == "velocity-magnitude") || (fieldname == "v")) {
			std::vector<double> vx, vy, vz;
			read_from_hdf5(filename, groupname, "x-velocity", vx);
			read_from_hdf5(filename, groupname, "y-velocity", vy);
			read_from_hdf5(filename, groupname, "z-velocity", vz);
			res.resize(0); res.reserve(vx.size());
			for (int i = 0; i < vx.size(); i++) res.push_back(std::sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]));
		}
		else {
			H5File file(filename.c_str(), H5F_ACC_RDONLY);
			//std::cout << "foo1" << std::endl;
			Group group = file.openGroup(groupname.c_str());
			//std::cout << "foo2" << std::endl;
			DataSet dataset = group.openDataSet(fieldname);
			//std::cout << "foo3" << std::endl;
			DataSpace dataspace = dataset.getSpace();
			//std::cout << "foo4" << std::endl;
			int rank = dataspace.getSimpleExtentNdims();
			//std::cout << "foo5 " << rank << std::endl;
			hsize_t dims[rank];
			//std::cout << "foo6" << std::endl;
			int ndims = dataspace.getSimpleExtentDims(dims, NULL);
			//std::cout << "foo7 " << ndims << " " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
			hsize_t offset[rank];
			for (auto & o: offset) o = 0;
			//std::cout << "foo8 " << offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
			dataspace.selectHyperslab(H5S_SELECT_SET, dims, offset);
			//std::cout << "foo9" << std::endl;
			DataSpace memspace(rank, dims);
			int res_size = 1;
			for (int i = 0; i < rank; i++) {
				res_size *= dims[i];
			}
			double *data_out = new double[res_size];
			dataset.read(data_out, PredType::NATIVE_DOUBLE, memspace, dataspace);
			res.resize(0); res.reserve(res_size);
			for (int i = 0; i < res_size; i++) res.push_back(data_out[i]);
			//std::cout << groupname << " " << fieldname << std::endl;
		}
	}
	
	// catch failure caused by the H5File operations
	catch( FileIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSet operations
	catch( DataSetIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error ) {
		error.printErrorStack();
		return -1;
	}
	return 0;
}

int getSize(std::string filename, std::string groupname) {
	std::string fieldname = "Density";
	try {
		H5File file(filename.c_str(), H5F_ACC_RDONLY);
		//std::cout << "foo1" << std::endl;
		Group group = file.openGroup(groupname.c_str());
		//std::cout << "foo2" << std::endl;
		DataSet dataset = group.openDataSet(fieldname);
		//std::cout << "foo3" << std::endl;
		DataSpace dataspace = dataset.getSpace();
		//std::cout << "foo4" << std::endl;
		int rank = dataspace.getSimpleExtentNdims();
		//std::cout << "foo5 " << rank << std::endl;
		hsize_t dims[rank];
		//std::cout << "foo6" << std::endl;
		int ndims = dataspace.getSimpleExtentDims(dims, NULL);
		//std::cout << "foo7 " << ndims << " " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
		//hsize_t offset[rank];
		//for (auto & o: offset) o = 0;
		//std::cout << "foo8 " << offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
		//dataspace.selectHyperslab(H5S_SELECT_SET, dims, offset);
		//std::cout << "foo9" << std::endl;
		//DataSpace memspace(rank, dims);
		int res_size = 1;
		for (int i = 0; i < rank; i++) {
			res_size *= dims[i];
		}
		return res_size;
	}
	
	// catch failure caused by the H5File operations
	catch( FileIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSet operations
	catch( DataSetIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error ) {
		error.printErrorStack();
		return -1;
	}
	return 0;
}

int read_from_hdf5(std::string filename, std::string groupname, std::string fieldname, double res[], int size) {
	if ((fieldname == "density") || (fieldname == "rho")) fieldname = "Density";
	else if (fieldname == "vx") fieldname = "x-velocity";
	else if (fieldname == "vy") fieldname = "y-velocity";
	else if (fieldname == "vz") fieldname = "z-velocity";
	//std::cout << "poop " << filename << std::endl;
	//std::cout << "poop " << groupname << std::endl;
	//std::cout << "poop " << fieldname << std::endl;
	try {
		if ((fieldname == "absv") || (fieldname == "velocity-magnitude") || (fieldname == "v")) {
			double vx[size];
			double vy[size];
			double vz[size];
			read_from_hdf5(filename, groupname, "x-velocity", vx, size);
			read_from_hdf5(filename, groupname, "y-velocity", vy, size);
			read_from_hdf5(filename, groupname, "z-velocity", vz, size);
			for (int i = 0; i < size; i++) res[i] = std::sqrt(vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);
		}
		else {
			H5File file(filename.c_str(), H5F_ACC_RDONLY);
			//std::cout << "foo1" << std::endl;
			Group group = file.openGroup(groupname.c_str());
			//std::cout << "foo2" << std::endl;
			DataSet dataset = group.openDataSet(fieldname);
			//std::cout << "foo3" << std::endl;
			DataSpace dataspace = dataset.getSpace();
			//std::cout << "foo4" << std::endl;
			int rank = dataspace.getSimpleExtentNdims();
			//std::cout << "foo5 " << rank << std::endl;
			hsize_t dims[rank];
			//std::cout << "foo6" << std::endl;
			int ndims = dataspace.getSimpleExtentDims(dims, NULL);
			//std::cout << "foo7 " << ndims << " " << dims[0] << " " << dims[1] << " " << dims[2] << std::endl;
			hsize_t offset[rank];
			for (auto & o: offset) o = 0;
			//std::cout << "foo8 " << offset[0] << " " << offset[1] << " " << offset[2] << std::endl;
			dataspace.selectHyperslab(H5S_SELECT_SET, dims, offset);
			//std::cout << "foo9" << std::endl;
			DataSpace memspace(rank, dims);
			dataset.read(res, PredType::NATIVE_DOUBLE, memspace, dataspace);
		}
	}
	
	// catch failure caused by the H5File operations
	catch( FileIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSet operations
	catch( DataSetIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataSpaceIException error ) {
		error.printErrorStack();
		return -1;
	}

	// catch failure caused by the DataSpace operations
	catch( DataTypeIException error ) {
		error.printErrorStack();
		return -1;
	}
	return 0;
}

std::vector<double> read_from_hdf5(std::string filename, std::string group, std::string fieldname) {
	std::vector<double> res;
	read_from_hdf5(filename, group, fieldname, res);
	return res;
}