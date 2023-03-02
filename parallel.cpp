#include "parallel.h"

std::vector<std::string> distribute_grids(std::string path_to_sim, int frame, std::vector<std::vector<std::string>> cpu_files, int Nthreads, int &size_filename, int &size_groupname, std::vector<int> &number_of_tasks) {
    std::vector<std::vector<std::pair<int, std::string>>> grids_to_thread;
    // big cloodge
    cpu_files.resize(100);
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