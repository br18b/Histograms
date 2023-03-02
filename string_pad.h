#ifndef STRING_PAD
#define STRING_PAD

#include <iostream>

std::string string_pad_left(std::string str, std::string padding, int total_length);

std::string string_pad_left(std::string str, int padding, int total_length);

std::string string_pad_left(int str, std::string padding, int total_length);

std::string string_pad_left(int str, int padding, int total_length);

std::string cpu_filename(std::string path_to_sim, int frame, int cpu_number);

void delete_trailing(std::string &str, char trail);

void unpack_grids(std::string &filename, std::string &groupname, const std::string &grids, int filename_length, int groupname_length, int index);

#endif
