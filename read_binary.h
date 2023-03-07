#pragma once
#ifndef __READ__BINARY
#define __READ__BINARY

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "string_pad.h"
#include "H5Cpp.h"

using namespace H5;

int read_from_binary(std::string path, std::string prefix, std::string name, std::vector<double>& x, bool verbose);
int read_from_binary(std::string path, std::string name, std::vector<double>& x, bool verbose);

int read_from_binary(std::string path, std::string prefix, std::string name, std::vector<double>& x);
int read_from_binary(std::string path, std::string name, std::vector<double>& x);

int getSize(std::string filename, std::string groupname);
int read_from_hdf5(std::string filename, std::string groupname, std::string fieldname, double res[], int size);
int read_from_hdf5(std::string filename, std::string group, std::string fieldname, std::vector<double> &res);
std::vector<double> read_from_hdf5(std::string filename, std::string group, std::string fieldname);

#endif
