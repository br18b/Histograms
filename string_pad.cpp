#include "string_pad.h"

std::string string_pad_left(std::string str, std::string padding, int total_length) {
        int n_str = static_cast<int>(str.length());
        int n_padding = static_cast<int>(padding.length());
        std::string res = str;
        int i = n_str;
        while (i + n_padding <= total_length) {
                res = padding + res;
                i += n_padding;
        }
        return res;
}

std::string string_pad_left(std::string str, int padding, int total_length) {
        return string_pad_left(str, std::to_string(padding), total_length);
}

std::string string_pad_left(int str, std::string padding, int total_length) {
        return string_pad_left(std::to_string(str), padding, total_length);
}

std::string string_pad_left(int str, int padding, int total_length) {
        return string_pad_left(std::to_string(str), std::to_string(padding), total_length);
}

std::string cpu_filename(std::string path_to_sim, int frame, int cpu_number) {
	while (path_to_sim.back() == '/') path_to_sim.pop_back();
	return path_to_sim + "/DD" + string_pad_left(frame, 0, 4) + "/data" + string_pad_left(frame, 0, 4) + ".cpu" + string_pad_left(cpu_number, 0, 4);
}

void delete_trailing(std::string &str, char trail) {
	while (str.back() == trail) str.pop_back();
}

void unpack_grids(std::string &filename, std::string &groupname, const std::string &grids, int filename_length, int groupname_length, int index) {
        int length = filename_length + groupname_length;
        filename = grids.substr(length * index, filename_length);
        groupname = grids.substr(length * index + filename_length, groupname_length);
}