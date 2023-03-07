#ifndef __PARALLEL__
#define __PARALLEL__

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include "string_pad.h"

void get_cpu_list(std::string path_to_sim, int frame, std::vector<std::vector<std::string>> &res);
std::vector<std::vector<std::string>> get_cpu_list(std::string path_to_sim, int frame);

std::vector<std::string> distribute_grids(std::string path_to_sim, int frame, std::vector<std::vector<std::string>> cpu_files, int Nthreads, int &size_filename, int &size_groupname, std::vector<int> &number_of_tasks);

std::vector<std::pair<std::string, int>> distribute_grids(std::string path_to_sim, const std::vector<std::vector<std::string>> &cpu_files, int rank, int Nthreads);
std::vector<std::pair<std::string, int>> distribute_grids(std::string path_to_sim, int frame, int rank, int Nthreads);

#endif