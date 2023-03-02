#ifndef __PARALLEL__
#define __PARALLEL__

#include <vector>
#include <iostream>
#include <thread>
#include "string_pad.h"

std::vector<std::string> distribute_grids(std::string path_to_sim, int frame, std::vector<std::vector<std::string>> cpu_files, int Nthreads, int &size_filename, int &size_groupname, std::vector<int> &number_of_tasks);

#endif