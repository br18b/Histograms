#include "parallel.h"

void get_cpu_list(std::string path_to_sim, int frame, std::vector<std::vector<std::string>> &res) {
	while (path_to_sim.back() == '/') path_to_sim.pop_back();
	std::string hierarchy_filename = path_to_sim + "/DD" + string_pad_left(frame, 0, 4) + "/data" + string_pad_left(frame, 0, 4) + ".hierarchy";
	std::ifstream hierarchy_in(hierarchy_filename);
	if (hierarchy_in) {
		//std::cout << "Read the hierarchy file:" << std::endl << hierarchy_filename << std::endl;
		std::string line, token;
		int last_grid = 0;
		while (std::getline(hierarchy_in, line)) {
			std::istringstream iss(line);
			iss >> token;
			if (token == "Grid") {
				iss >> token; // read "="
				iss >> token; // read Grid number
				last_grid = std::stoi(token);
			}
			else if (token == "BaryonFileName") {
				iss >> token; // read "="
				iss >> token; // read "./DDxxxx/dataxxxx.cpuyyyy"
				token = token.substr(token.size() - 4, token.size() - 1);
				while ((token.size() > 0) && (token[0] == '0')) token.erase(0, 1);
				int cpu_n = 0;
				if (token.size() > 0) cpu_n = std::stoi(token);
				if (res.size() <= cpu_n) res.resize(cpu_n + 1);
				res[cpu_n].push_back("Grid" + string_pad_left(last_grid, 0, 8));
			}
		}
		//std::cout << "Obtained group names:" << std::endl;
		//std::cout << path_to_sim + "/DD" + string_pad_left(frame, 0, 4) + "/data" + string_pad_left(frame, 0, 4) + ".cpu0000/" + res[0][0] << std::endl;
		//std::cout << "..." << std::endl;
		//std::cout << path_to_sim + "/DD" + string_pad_left(frame, 0, 4) + "/data" + string_pad_left(frame, 0, 4) + ".cpu" + string_pad_left(res.size() - 1, 0, 4) + "/" + res.back().back() << std::endl;
	}
	else {
		std::cout << "File" << std::endl << hierarchy_filename << std::endl << "does not exist!" << std::endl;
	}
}

std::vector<std::vector<std::string>> get_cpu_list(std::string path_to_sim, int frame) {
	std::vector<std::vector<std::string>> res;
	get_cpu_list(path_to_sim, frame, res);
	return res;
}

std::vector<std::string> distribute_grids(std::string path_to_sim, int frame, std::vector<std::vector<std::string>> cpu_files, int Nthreads, int &size_filename, int &size_groupname, std::vector<int> &number_of_tasks) {
    std::vector<std::vector<std::pair<int, std::string>>> grids_to_thread;
    // big cloodge
    //cpu_files.resize(100);
    grids_to_thread.resize(Nthreads);
    int j = 0;
    int k = 0;
    std::cout << "starting to create distribution of cpu files" << std::endl;
    while (j < cpu_files.size()) {
        int c = 0;
        while (c < cpu_files[j].size()) {
            grids_to_thread[k].push_back({j,cpu_files[j][c]});
            c++; k++;
            if (k == Nthreads) k = 0;
        }
        j++;
    }
    number_of_tasks.resize(Nthreads);
    for (int k = 0; k < grids_to_thread.size(); k++) number_of_tasks[k] = grids_to_thread[k].size();
    std::cout << "done with it" << std::endl;
    std::cout << "Packaging tasks..." << std::endl;
    std::vector<std::string> grids_to_thread_packaged; grids_to_thread_packaged.resize(grids_to_thread.size());
    size_filename = 0; size_groupname = 0;
    while (path_to_sim.back() == '/') path_to_sim.pop_back();
    for (int rank = 0; rank < grids_to_thread.size(); rank++) {
        for (auto & pair: grids_to_thread[rank]) {
            int cpuFileNumber = pair.first;
            std::string filename = path_to_sim + "/DD" + string_pad_left(frame, 0, 4) + "/data" + string_pad_left(frame, 0, 4) + ".cpu" + string_pad_left(cpuFileNumber, 0, 4);
            std::string groupname = pair.second;
            if (filename.size() > size_filename) size_filename = filename.size();
            if (groupname.size() > size_groupname) size_groupname = groupname.size();
        }
    }
    for (int rank = 0; rank < grids_to_thread.size(); rank++) {
        std::string to_package = "";
        for (auto & pair: grids_to_thread[rank]) {
            int cpuFileNumber = pair.first;
            std::string filename = path_to_sim + "/DD" + string_pad_left(frame, 0, 4) + "/data" + string_pad_left(frame, 0, 4) + ".cpu" + string_pad_left(cpuFileNumber, 0, 4);
            std::string groupname = pair.second;
            while (filename.size() < size_filename) filename.push_back(' ');
            while (groupname.size() < size_groupname) groupname.push_back(' ');
            to_package = to_package + filename + groupname;
        }
        grids_to_thread_packaged[rank] = to_package;
    }
    return grids_to_thread_packaged;
}

std::vector<std::pair<std::string, int>> distribute_grids(std::string path_to_sim, const std::vector<std::vector<std::string>> &cpu_files, int rank, int Nthreads) {
    std::vector<std::pair<std::string, int>> grids_to_thread;
    // big cloodge
    //cpu_files.resize(100);
    int j = 0;
    int k = 0;
    //std::cout << "starting to create distribution of cpu files" << std::endl;
    while (j < cpu_files.size()) {
    //while (j < 100) {
        int c = 0;
        while (c < cpu_files[j].size()) {
            if (k == rank) grids_to_thread.push_back({cpu_files[j][c], j});
            c++; k++;
            if (k == Nthreads) k = 0;
        }
        j++;
    }
    return grids_to_thread;
}

std::vector<std::pair<std::string, int>> distribute_grids(std::string path_to_sim, int frame, int rank, int Nthreads) {
    std::vector<std::vector<std::string>> cpu_files = get_cpu_list(path_to_sim, frame);
    return distribute_grids(path_to_sim, cpu_files, rank, Nthreads);
}