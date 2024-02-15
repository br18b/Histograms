#include "statistics.h"

std::function<double(double, double)> weight_default = [](double rho, double v) { return 1; };

std::function<double(double, double)> v2V = [](double rho, double v) { return v*v; };
std::function<double(double, double)> v2M = [](double rho, double v) { return rho*v*v; };
std::function<double(double, double)> v2E = [](double rho, double v) { return 0.5*rho*v*v*v*v; };
std::function<double(double, double)> density = [](double rho, double v) { return rho; };
std::function<double(double, double)> kinetic = [](double rho, double v) { return 0.5*rho*v*v; };
std::function<double(double, double)> thermal = [](double rho, double v) { return rho * std::log(rho); };
std::function<double(double, double)> v_mag = [](double rho, double v) { return v; };
std::function<double(double, double)> logrho = [](double rho, double v) { return std::log(rho); };
std::function<double(double, double)> logrho2 = [](double rho, double v) { return std::log(rho)*std::log(rho); };
std::function<double(double, double)> rho_logrho = [](double rho, double v) { return rho * std::log(rho); };
std::function<double(double, double)> rho_logrho2 = [](double rho, double v) { return rho * std::log(rho) * std::log(rho); };
std::function<double(double, double)> kin_logrho = [](double rho, double v) { return 0.5 * rho * v * v * std::log(rho); };
std::function<double(double, double)> kin_logrho2 = [](double rho, double v) { return 0.5 * rho * v * v * std::log(rho) * std::log(rho); };
std::function<double(double, double)> fun_sv = [](double rho, double v) { return std::log(rho)*v; };
std::function<double(double, double)> fun_sv2 = [](double rho, double v) { return std::log(rho)*v*v; };
std::function<double(double, double)> fun_sv3 = [](double rho, double v) { return std::log(rho)*v*v*v; };
std::function<double(double, double)> fun_sv4 = [](double rho, double v) { return std::log(rho)*v*v*v*v; };
std::function<double(double, double)> fun_s2v = [](double rho, double v) { return std::log(rho)*std::log(rho)*v; };
std::function<double(double, double)> fun_s2v2 = [](double rho, double v) { return std::log(rho)*std::log(rho)*v*v; };
std::function<double(double, double)> fun_s2v3 = [](double rho, double v) { return std::log(rho)*std::log(rho)*v*v*v; };
std::function<double(double, double)> fun_s2v4 = [](double rho, double v) { return std::log(rho)*std::log(rho)*v*v*v*v; };

std::function<double(std::vector<double>, std::vector<int>)> M1DV = [](std::vector<double> stats, std::vector<int> stat_indices) { return std::sqrt(stats[stat_indices[0]] / 3); };
std::function<double(std::vector<double>, std::vector<int>)> M1DM = [](std::vector<double> stats, std::vector<int> stat_indices) { return std::sqrt(stats[stat_indices[0]] / (3 * stat_indices[1])); };
std::function<double(std::vector<double>, std::vector<int>)> M1DE = [](std::vector<double> stats, std::vector<int> stat_indices) { return std::sqrt(stats[stat_indices[0]] / (3 * stat_indices[1])); };
std::function<double(std::vector<double>, std::vector<int>)> muV = [](std::vector<double> stats, std::vector<int> stat_indices) { return stats[stat_indices[0]]; };
std::function<double(std::vector<double>, std::vector<int>)> sigmaV = [](std::vector<double> stats, std::vector<int> stat_indices) {
        double mu = stats[stat_indices[0]];
        double s2 = stats[stat_indices[1]];
        return std::sqrt(s2 - mu*mu);
};

std::function<double(std::vector<double>, std::vector<int>)> muM = [](std::vector<double> stats, std::vector<int> stat_indices) {
        double srho = stats[stat_indices[1]];
        double rho = stats[stat_indices[0]];
        return srho/rho;
};

std::function<double(std::vector<double>, std::vector<int>)> sigmaM = [](std::vector<double> stats, std::vector<int> stat_indices) {
        double srho = stats[stat_indices[1]];
        double s2rho = stats[stat_indices[2]];
        double rho = stats[stat_indices[0]];
        double muM = srho / rho;
        double s2 = s2rho / rho;
        return std::sqrt(s2 - muM);
};

std::function<double(std::vector<double>, std::vector<int>)> muE = [](std::vector<double> stats, std::vector<int> stat_indices) {
        double skin = stats[stat_indices[1]];
        double kin = stats[stat_indices[0]];
        return skin/kin;
};

std::function<double(std::vector<double>, std::vector<int>)> sigmaE = [](std::vector<double> stats, std::vector<int> stat_indices) {
        double skin = stats[stat_indices[1]];
        double s2kin = stats[stat_indices[2]];
        double kin = stats[stat_indices[0]];
        double muE = skin / kin;
        double s2 = s2kin / kin;
        return std::sqrt(s2 - muE);
};

void resolve_set(std::set<double> &temp_set, double &accumulator, bool override) {
	if ((temp_set.size() > 4096) || override) {
		for (auto & s: temp_set) accumulator += s;
		temp_set.clear();
	}
}

void extract_stats(std::string path_to_sim, std::string path_output, std::vector<int> frames, std::vector<std::function<double(double, double)>> stat_funs, std::vector<std::pair<std::function<double(std::vector<double>, std::vector<int>)>, std::vector<int>>> stat_aggregate, std::vector<std::string> stat_filenames, std::vector<std::string> stat_aggregate_filenames) {
	delete_trailing(path_to_sim, '/'); delete_trailing(path_output, '/');
	std::cout << "boo!" << std::endl;
	std::vector<std::set<double>> stats; stats.resize(stat_funs.size());
	std::vector<double> stats_per_frame; stats_per_frame.resize(stat_funs.size());
	std::vector<double> stats_total; stats_total.resize(stat_funs.size());
	for (int i = 0; i < frames.size(); i++) {
		for (auto & stat: stats_per_frame) stat = 0;
		int frame = frames[i];
		std::vector<std::vector<std::string>> cpu_files = get_cpu_list(path_to_sim, frame);
		std::cout << "Analyzing frame " << frame;
		int n = 0;
		for (int j = 0; j < cpu_files.size(); j++) {
			std::string filename = path_to_sim + "/DD" + string_pad_left(frame, 0, 4) + "/data" + string_pad_left(frame, 0, 4) + ".cpu" + string_pad_left(j, 0, 4);
			std::cout << "Reading from: ";
			if (filename.size() > 32) {
				std::cout << filename.substr(0, 16) << "..." << filename.substr(filename.size() - 16, filename.size() - 1) << std::endl;
			}
			else {
				std::cout << filename << std::endl;
			}
			for (int k = 0; k < cpu_files[j].size(); k++) {
				std::string group = cpu_files[j][k];
				std::cout << "Analyzing: " << group << std::endl;
				std::vector<double> densities = read_from_hdf5(filename, group, "Density");
				std::vector<double> absv = read_from_hdf5(filename, group, "absv");

				for (int l = 0; l < densities.size(); l++) {
					for (int m = 0; m < stat_funs.size(); m++) {
						auto stat_function = stat_funs[m];
						stats[m].insert(stat_function(densities[l], absv[l]));
						resolve_set(stats[m], stats_per_frame[m], 0);
						//std::cout << stats[m].size() << std::endl;
					}
					n++;
				}
				for (int m = 0; m < stat_funs.size(); m++) resolve_set(stats[m], stats_per_frame[m], 1);
			}
		}
		for (int m = 0; m < stat_funs.size(); m++) {
			stats_per_frame[m] /= n;
			std::ofstream stat_out(path_output + "/" + stat_filenames[m]);
			stat_out << stats_per_frame[m];
			stat_out.close();
		}
	}
}

void extract_stats_parallel(std::string path_to_sim, std::string path_output, std::vector<int> frames, std::vector<std::function<double(double, double)>> stat_funs, std::vector<std::pair<std::function<double(std::vector<double>, std::vector<int>)>, std::vector<int>>> stat_aggregate, std::vector<std::string> stat_filenames, std::vector<std::string> stat_aggregate_filenames, int Nthreads, int rank) {
	delete_trailing(path_to_sim, '/'); delete_trailing(path_output, '/');
	//std::cout << "boo!" << std::endl;
	std::vector<std::set<double>> stats; stats.resize(stat_funs.size());
	std::vector<double> stats_per_frame; stats_per_frame.resize(stat_funs.size());
	std::vector<double> stats_total; stats_total.resize(stat_funs.size());
	//std::cout << "Nframes = "  << frames.size() << std::endl;
	for (int i = 0; i < frames.size(); i++) {
		for (auto & stat: stats_per_frame) stat = 0;
		int frame = frames[i];
		std::cout << "Current frame: " << frame << std::endl;
		std::string my_grids; int filename_size, groupname_size;
		int my_grid_size;
		if (rank == 0) {
			//std::cout << "I'm rank zero and I got here!" << std::endl;
			std::vector<std::vector<std::string>> cpu_files;
			cpu_files = get_cpu_list(path_to_sim, frame);
			//std::cout << "Analyzing frame " << frame;
			std::vector<int> number_of_tasks;
			std::vector<std::string> grids_to_thread_packaged = distribute_grids(path_to_sim, frame, cpu_files, Nthreads, filename_size, groupname_size, number_of_tasks);
			my_grid_size = number_of_tasks[0];
			my_grids = grids_to_thread_packaged[0];
			for (int r = 1; r < Nthreads; r++) {
				std::cout << "Sending to r = " << r << std::endl;
				int current_size = number_of_tasks[r];
				char* grids_array = new char[(filename_size + groupname_size) * current_size];
				std::cout << "sizes: " << grids_to_thread_packaged[r].size() << " " << (filename_size + groupname_size) * number_of_tasks[r] << std::endl;
				for (int k = 0; k < grids_to_thread_packaged[r].size(); k++) {
					grids_array[k] = grids_to_thread_packaged[r][k];
				}
				MPI_Send(&current_size, 1, MPI_INT, r, 0, MPI_COMM_WORLD);
				MPI_Send(&filename_size, 1, MPI_INT, r, 1, MPI_COMM_WORLD);
				MPI_Send(&groupname_size, 1, MPI_INT, r, 2, MPI_COMM_WORLD);
				MPI_Send(grids_array, (filename_size + groupname_size) * number_of_tasks[r], MPI_CHAR, r, 3, MPI_COMM_WORLD);
				delete[] grids_array;
			}
			std::cout << "everything has been sent!" << std::endl;
			
			std::cout << "my rank: " << rank << ", my size: " << my_grid_size << std::endl;
			std::cout << rank << " poop: " << my_grids[0] << std::endl;
		}
		//std::cout << "I'm rank " << rank << " and I'm receiving" << std::endl;

		if (rank > 0) {
			MPI_Recv(&my_grid_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::cout << "my rank: " << rank << ", my size: " << my_grid_size << std::endl;
			MPI_Recv(&filename_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//std::cout << "filename size: " << filename_size << std::endl;
			MPI_Recv(&groupname_size, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//std::cout << "groupname size: " << groupname_size << std::endl;
			int length = (filename_size + groupname_size) * my_grid_size;
			char* grids_array = new char[length];
			my_grids.resize(length);
			//std::cout << "poop1" <<std::endl;
			MPI_Recv(grids_array, length, MPI_CHAR, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int k = 0; k < length; k++) {
				my_grids[k] = grids_array[k];
			}
			delete[] grids_array;
			std::cout << rank << " poop: " << my_grids[0] << std::endl;
		}

		int n = 0;
		
		std::cout << my_grid_size << std::endl;

		for (int j = 0; j < my_grid_size; j++) {
			std::string task = my_grids.substr((filename_size + groupname_size) * j, filename_size + groupname_size);
			std::string filename = task.substr(0, filename_size);
			std::string groupname = task.substr(filename_size, groupname_size);
			//std::cout << "poop " << (filename_size + groupname_size) * j << " " << (filename_size + groupname_size) * (j + 1) << std::endl;
			//std::cout << "task " << task << std::endl;
			while(filename.back() == ' ') filename.pop_back();
			while(groupname.back() == ' ') groupname.pop_back();
			if (rank == 1) std::cout << filename << " " << groupname << std::endl;
			//std::cout << filename << " " << groupname << std::endl;
			std::vector<double> densities = read_from_hdf5(filename, groupname, "density");
			std::vector<double> absv = read_from_hdf5(filename, groupname, "absv");
			for (int l = 0; l < densities.size(); l++) {
				for (int m = 0; m < stat_funs.size(); m++) {
					auto stat_function = stat_funs[m];
					stats[m].insert(stat_function(densities[l], absv[l]));
					resolve_set(stats[m], stats_per_frame[m], 0);
					//std::cout << stats[m].size() << std::endl;
				}
				n++;
			}
			for (int m = 0; m < stat_funs.size(); m++) resolve_set(stats[m], stats_per_frame[m], 1);
		}
		/*
		for (int m = 0; m < stat_funs.size(); m++) {
			stats_per_frame[m] /= n; // divide by n but it needs to be re-weighted from each drone anyway
		}*/
		if (rank > 0) {
			for (int i = 0; i < stat_funs.size(); i++) {
				int tag = 3 + rank + Nthreads * i + Nthreads * 1 * 0; // index = i + nx * j + nx * ny * k
				std::cout << "rank: " << rank << ", i: " << i << ", tag: " << tag << std::endl; // i = rank, j = i, k = 0,1 (sending int n or double stats)
				MPI_Send(&n, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
				tag = 3 + rank + Nthreads * i + Nthreads * 1 * 1;
				MPI_Send(&stats_per_frame[i], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		//std::cout << "my rank: " << rank << " " << stat_funs.size() << std::endl;
		if (rank == 0) {
			for (int i = 0; i < stat_funs.size(); i++) {
				std::cout << stat_filenames[i] << " " << stats_per_frame[i] << std::endl;
				int N = n; // we will accummulate total N;
				for (int r = 1; r < Nthreads; r++) {
					double stat_from_other = 0;
					int tag = 3 + r + Nthreads * i + Nthreads * 1 * 0;
					std::cout << "rank: " << r << ", i: " << i << ", tag: " << tag << std::endl;
					MPI_Recv(&n, 1, MPI_INT, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					tag = 3 + r + Nthreads * i + Nthreads * 1 * 1;
					MPI_Recv(&stat_from_other, 1, MPI_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					std::cout << "I'm rank " << rank << ", variable: " << stat_filenames[i] << ", value: " << stat_from_other << std::endl;
					std::cout.flush();
					N += n;
					stats_per_frame[i] += stat_from_other;
				}
				stats_per_frame[i] /= N; // now we divide by the collective N = n0 + n1 + ... + nrank + ... + n(Threads-1)
				std::cout << stat_filenames[i] << " " << stats_per_frame[i] << std::endl;
			}
			
			for (int m = 0; m < stat_funs.size(); m++) {
				std::cout << "outputting" << std::endl;
				std::ofstream stat_out(path_output + "/" + stat_filenames[m]);
				stat_out << stats_per_frame[m];
				stat_out.close();
			}
			
		}
	}
}
/*
void extract_1D_histogram_parallel(std::string path_to_sim, std::string path_output, std::vector<int> frames, std::string fieldname_x, std::string fieldname_y, std::function<double(double, double)> weight_fun, std::function<double(double, double)> field_1D, int depth, std::vector<std::string> merge_1D, std::string outputname_x, std::string outputpostfix_weight, int rank, int Nthreads) {
	delete_trailing(path_to_sim, '/'); delete_trailing(path_output, '/');
	
	double* bounds = new double[2];
	bounds[0] = 1e10; // xmin
	bounds[1] = -1e10; // xmax

	int offset_max = 0;
	int global_offset = 0;
	int offset = 0;
	int tag = 0;

	int field_size = getSize(path_to_sim + "/DD" + string_pad_left(frames[0], 0, 4) + "/data" + string_pad_left(frames[0], 0, 4) + ".cpu0000", "Grid00000001");
	//std::cout << field_size << std::endl;
	double field_x[field_size];
	double field_y[field_size];

	for (int i = 0; i < frames.size(); i++) {
		int frame = frames[i];
		std::vector<std::pair<std::string, int>> my_grids = distribute_grids(path_to_sim, frame, rank, Nthreads);
		
		int n = 0;
		// __________________________________
		// finding bounds of the 1D histogram
		// __________________________________
		for (int j = 0; j < my_grids.size(); j++) {
			int cpu_number = my_grids[j].second;
			std::string filename = cpu_filename(path_to_sim, frame, cpu_number);
			std::string groupname = my_grids[j].first;

			read_from_hdf5(filename, groupname, fieldname_x, field_x, field_size);
			read_from_hdf5(filename, groupname, fieldname_y, field_y, field_size);
			
			//std::cout << "Rank " << rank << " adjusting bounds" << std::endl;
			// ________________________________________________
			// each rank will compute bounds for its own subset
			// ________________________________________________
			for (int i = 0; i < field_size; i++) {
				double val = field_1D(field_x[i], field_y[i]);
				if (bounds[0] > val) bounds[0] = val;
				if (bounds[1] < val) bounds[1] = val;
			}
			//MPI_Barrier(MPI_COMM_WORLD);
			if (rank == 0) std::cout << j + 1 << "/" << my_grids.size() << std::endl;
			//if (rank == 52) std::cout << rank << " " << j + 1 << "/" << my_grid_size << std::endl;
		}
	}
	// ____________________________________________________________
	// all ranks except the root will send their bounds to the root
	// ____________________________________________________________
	if (rank > 0) {
		tag = rank;
		//std::cout << "Rank " << rank << " sending bounds..." << std::endl;
		MPI_Send(bounds, 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
		//std::cout << "I'm rank " << rank << " and sent bounds (" << bounds[0] << ", " << bounds[1] << ") to rank 0 with tag " << tag <<std::endl;
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) std::cout << "(waiting on all cores...)" << std::endl;
	if (rank == 0) {
		for (int r = 1; r < Nthreads; r++) {
			double* bounds_other = new double[2];
			tag = r;
			// ____________________________________________
			// the root will receive bounds from all drones
			// ____________________________________________
			MPI_Recv(bounds_other, 2, MPI_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			//std::cout << "I'm rank 0 and received bounds from " << r << " with tag " << tag << ".  Bounds: (" << bounds_other[0] << ", " << bounds_other[1] << ")" << std::endl;
			// _______________________________________________________________
			// the root calculates the global bounds across the whole dataset
			// _______________________________________________________________
			if (bounds[0] > bounds_other[0]) bounds[0] = bounds_other[0];
			if (bounds[1] < bounds_other[1]) bounds[1] = bounds_other[1];
			delete[] bounds_other;
		}
	}
	global_offset += Nthreads;
	// ____________________________________________
	// the root sends the global bounds to everyone
	// ____________________________________________
	if (rank == 0) {
		for (int r = 1; r < Nthreads; r++) {
			tag = global_offset + r;
			//std::cout << "I'm rank 0 and I'm sending bounds to rank " << r << " with tag " << tag <<std::endl;
			MPI_Send(bounds, 2, MPI_DOUBLE, r, tag, MPI_COMM_WORLD);
		}
	}

	if (rank == 0) std::cout << "(computing global bounds...)" << std::endl;
	if (rank > 0) {
		tag = global_offset + rank;
		//std::cout << "I'm rank " << rank << " and I'm receiving bounds from rank 0" << std::endl;
		MPI_Recv(bounds, 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) std::cout << "xmin = " << bounds[0] << ", xmax = " << bounds[1] << std::endl;
	
	global_offset += Nthreads;
	// ______________________________________
	// now we actually compute the histograms
	// ______________________________________
	
	BinaryTree histogram1D;
	histogram1D.initialize(bounds[0], bounds[1], "lin");
	histogram1D.uniform_divide(depth);

	delete[] bounds;
	for (int i = 0; i < frames.size(); i++) {
		int frame = frames[i];
		std::vector<std::pair<std::string, int>> my_grids = distribute_grids(path_to_sim, frame, rank, Nthreads);

		int n = 0;

		for (int j = 0; j < my_grids.size(); j++) {
			int cpu_number = my_grids[j].second;
			std::string filename = cpu_filename(path_to_sim, frame, cpu_number);
			std::string groupname = my_grids[j].first;
			
			//std::cout << "Rank " << rank << " is filling histogram from " << filename << " " << groupname << std::endl;
			
			read_from_hdf5(filename, groupname, fieldname_x, field_x, field_size);
			read_from_hdf5(filename, groupname, fieldname_y, field_y, field_size);
			
			//if (rank == 0) std::cout << "Filling histograms " << field_x.size() << " " << field_x[0] << " " << field_y[0] << std::endl;

			//histogram1D.load_points(field_x, field_y, field_1D, weight_fun, field_size);

			if (rank == 0) std::cout << j + 1 << "/" << my_grids.size() << std::endl;
			// actual points are stored in BinaryTree within Node1D.weight - public interface.
			// all histograms have the same structure (created with uniform_divide)
			// only thing that needs to be accummulated are the weights from different drones
			// size doesn't need to be sent, as each drone is working with the same size, root included
		}
		//MPI_Barrier(MPI_COMM_WORLD);
	}
	
	//MPI_Barrier(MPI_COMM_WORLD);
	if (rank > 0) {
		//std::cout << "I motherfucking got here" << std::endl;
		tag = global_offset + rank;
		int nodesize = histogram1D.nodes.size();

		double* my_weights = new double[nodesize];
		for (int i = 0; i < nodesize; i++) my_weights[i] = histogram1D.nodes[i].weight;
		MPI_Send(my_weights, nodesize, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
		delete[] my_weights;
	}
	if (rank == 0) std::cout << "The head core is collecting histograms... " << std::endl;
	// ________________________________________________
	// collecting all histograms and combining into one
	// ________________________________________________
	if (rank == 0) {
		std::cout << "Collecting histograms ... " << std::endl;
		for (int r = 1; r < Nthreads; r++) {
			tag = global_offset + r;
			int nodesize = histogram1D.nodes.size(); // same nodesize everywhere, no need to communicate this
			double* drone_weights = new double[nodesize];
			MPI_Recv(drone_weights, nodesize, MPI_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int i = 0; i < nodesize; i++) {
				histogram1D.nodes[i].weight += drone_weights[i];
			}
			delete[] drone_weights;
		}
		// __________________________________
		// everything else is a one man's job
		// __________________________________
		for (int m = 0; m < merge_1D.size(); m++) {
			histogram1D.merge(std::stod(merge_1D[m]));
			histogram1D.save_structure(path_output, outputname_x + "_bintree_frames_" + std::to_string(frames.front()) + "-" + std::to_string(frames.back()) + "_" + outputpostfix_weight + "_" + merge_1D[m] + ".txt");
			histogram1D.save_structure(path_output, outputname_x + "_bintree_" + outputpostfix_weight + "_" + merge_1D[m] + ".txt");
			histogram1D.save_leaves(path_output, outputname_x + "_bins_frames_" + std::to_string(frames.front()) + "-" + std::to_string(frames.back()) + "_" + outputpostfix_weight + "_" + merge_1D[m] + ".txt");
			histogram1D.save_leaves(path_output, outputname_x + "_bins_frames_" + outputpostfix_weight + "_" + merge_1D[m] + ".txt");
			histogram1D.prob_to_pdf(path_output, outputname_x + "_CDF_frames_" + std::to_string(frames.front()) + "-" + std::to_string(frames.back()) + "_" + outputpostfix_weight + "_" + merge_1D[m] + ".txt");
			histogram1D.prob_to_pdf(path_output, outputname_x + "_CDF_frames_" + outputpostfix_weight + "_" + merge_1D[m] + ".txt");
		}
	}
}
*/
void initialize_bounds(long double bounds[], const Parameters &params) {
	for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
		bounds[0 + 2 * h1D] = 1e10;
		bounds[1 + 2 * h1D] = -1e10;
	}
	int bounds_offset = 2 * params.transform_1D.size();

	for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
		bounds[bounds_offset + 0 + 4 * h2D] = 1e10;
		bounds[bounds_offset + 1 + 4 * h2D] = -1e10;
		bounds[bounds_offset + 2 + 4 * h2D] = 1e10;
		bounds[bounds_offset + 3 + 4 * h2D] = -1e10;
	}
}

void extract_stats(long double pars[], const Parameters &params, long double field_x[], long double field_y[], int field_size) {
	std::vector<std::multiset<long double>> stat_set;
	stat_set.resize(2 * params.stat_fun.size());

	for (int k = 0; k < params.stat_fun.size(); k++) {
		const auto &stat = params.stat_fun[k];

		for(int i = 0; i < field_size; i++) {
			long double stat_val = stat(field_x[i], field_y[i]);

			if (stat_val > 0) stat_set[2 * k].insert(stat_val);
			else stat_set[2 * k + 1].insert(-stat_val);

			handle_set(pars[2 * k], stat_set[2 * k]);
			handle_set(pars[2 * k + 1], stat_set[2 * k + 1]);
		}
		
		handle_set(pars[2 * k], stat_set[2 * k], true);
		handle_set(pars[2 * k + 1], stat_set[2 * k + 1], true);
	}
}

void adjust_bounds(long double bounds[], const Parameters &params, long double field_x[], long double field_y[], int field_size) {
	for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
		const auto &field_1D = params.transform_1D[h1D];

		long double &xmin = bounds[0 + 2 * h1D];
		long double &xmax = bounds[1 + 2 * h1D];

		for(int i = 0; i < field_size; i++) {
			double val = field_1D(field_x[i], field_y[i]);
			if (xmin > val) xmin = val;
			if (xmax < val) xmax = val;
		}
	}
			
	int bounds_offset = 2 * params.transform_1D.size();

	for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
		const auto &field_2D_x = params.transform_2D[h2D].first;
		const auto &field_2D_y = params.transform_2D[h2D].second;
				
		long double &xmin = bounds[bounds_offset + 0 + 4 * h2D];
		long double &xmax = bounds[bounds_offset + 1 + 4 * h2D];
		long double &ymin = bounds[bounds_offset + 2 + 4 * h2D];
		long double &ymax = bounds[bounds_offset + 3 + 4 * h2D];

		for(int i = 0; i < field_size; i++) {
			double val_x = field_2D_x(field_x[i], field_y[i]);
			double val_y = field_2D_y(field_x[i], field_y[i]);
			if (xmin > val_x) xmin = val_x;
			if (xmax < val_x) xmax = val_x;
			if (ymin > val_y) ymin = val_y;
			if (ymax < val_y) ymax = val_y;
		}
	}
}

void adjust_main_bounds(long double bounds[], long double bounds_other[], const Parameters &params) {
	for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
		long double &xmin = bounds[0 + 2 * h1D];
		long double &xmax = bounds[1 + 2 * h1D];

		const long double &xmin_other = bounds_other[0 + 2 * h1D];
		const long double &xmax_other = bounds_other[1 + 2 * h1D];

		if (xmin > xmin_other) xmin = xmin_other;
		if (xmax < xmax_other) xmax = xmax_other;
	}
			
	int bounds_offset = 2 * params.transform_1D.size();

	for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
		long double &xmin = bounds[bounds_offset + 0 + 4 * h2D];
		long double &xmax = bounds[bounds_offset + 1 + 4 * h2D];
		long double &ymin = bounds[bounds_offset + 2 + 4 * h2D];
		long double &ymax = bounds[bounds_offset + 3 + 4 * h2D];
			
		const long double &xmin_other = bounds_other[bounds_offset + 0 + 4 * h2D];
		const long double &xmax_other = bounds_other[bounds_offset + 1 + 4 * h2D];
		const long double &ymin_other = bounds_other[bounds_offset + 2 + 4 * h2D];
		const long double &ymax_other = bounds_other[bounds_offset + 3 + 4 * h2D];

		if (xmin > xmin_other) xmin = xmin_other;
		if (xmax < xmax_other) xmax = xmax_other;
		if (ymin > ymin_other) ymin = ymin_other;
		if (ymax < ymax_other) ymax = ymax_other;
	}
}

void initialize_histograms(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, long double bounds[], const Parameters &params, int &weights_count) {
	weights_count = 0;
	for (int w = 0; w < params.weight.size(); w++) {
		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			long double &xmin = bounds[0 + 2 * h1D];
			long double &xmax = bounds[1 + 2 * h1D];
			auto &hist1D = histograms1D[h1D + params.transform_1D.size() * w];

			hist1D.initialize(xmin, xmax, "lin");
			hist1D.uniform_divide(params.initial_depth_1D[h1D]);
			weights_count += hist1D.nodes.size();
		}
			
		int bounds_offset = 2 * params.transform_1D.size();

		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			long double &xmin = bounds[bounds_offset + 0 + 4 * h2D];
			long double &xmax = bounds[bounds_offset + 1 + 4 * h2D];
			long double &ymin = bounds[bounds_offset + 2 + 4 * h2D];
			long double &ymax = bounds[bounds_offset + 3 + 4 * h2D];
			auto &hist2D = histograms2D[h2D + params.transform_2D.size() * w];

			hist2D.initialize(xmin, xmax, ymin, ymax, "lin", "lin");
			hist2D.uniform_divide(params.initial_depth_2D[h2D]);
			weights_count += hist2D.nodes.size();
		}
	}
}

void populate_histograms(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, long double field_x[], long double field_y[], int field_size, const Parameters &params) {
	for (int w = 0; w < params.weight.size(); w++) {
		const auto &weight_fun = params.weight[w];
		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			const auto &field_1D = params.transform_1D[h1D];
			auto &hist1D = histograms1D[h1D + params.transform_1D.size() * w];
			
			hist1D.load_points(field_x, field_y, field_1D, weight_fun, field_size);
		}
		
		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			const auto &field_2D_x = params.transform_2D[h2D].first;
			const auto &field_2D_y = params.transform_2D[h2D].second;
			auto &hist2D = histograms2D[h2D + params.transform_2D.size() * w];

			hist2D.load_points(field_x, field_y, field_2D_x, field_2D_y, weight_fun, field_size);
		}
	}
}

void populate_histograms_single_frame(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, long double field_x[], long double field_y[], int field_size, const Parameters &params) {
	for (int w = 0; w < params.weight.size(); w++) {
		const auto &weight_fun = params.weight[w];

		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			const auto &field_1D = params.transform_1D[h1D];

			for (int m = 0; m < params.merge_fraction_1D.size(); m++) {
				int hist_index = w + params.weight.size() * h1D + params.weight.size() * params.transform_1D.size() * m;
				auto &hist1D = histograms1D[hist_index];
				
				hist1D.load_points(field_x, field_y, field_1D, weight_fun, field_size);
			}
		}
		
		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			const auto &field_2D_x = params.transform_2D[h2D].first;
			const auto &field_2D_y = params.transform_2D[h2D].second;

			for (int m = 0; m < params.merge_fraction_2D.size(); m++) {
				int hist_index = w + params.weight.size() * h2D + params.weight.size() * params.transform_2D.size() * m;
				auto &hist2D = histograms2D[hist_index];

				hist2D.load_points(field_x, field_y, field_2D_x, field_2D_y, weight_fun, field_size);
			}
		}
	}
}

void pack_weights(const std::vector<BinaryTree> &histograms1D, const std::vector<Quadtree> &histograms2D, long double histogram_weights[], const Parameters &params) {
	int weights_offset = 0;
	for (int w = 0; w < params.weight.size(); w++) {
		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			const auto &hist1D = histograms1D[h1D + params.transform_1D.size() * w];
			
			for (int i = 0; i < hist1D.nodes.size(); i++) {
				histogram_weights[weights_offset + i] = hist1D.nodes[i].weight;
			}
			weights_offset += hist1D.nodes.size();
		}
		
		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			const auto &hist2D = histograms2D[h2D + params.transform_2D.size() * w];

			for (int i = 0; i < hist2D.nodes.size(); i++) {
				histogram_weights[weights_offset + i] = hist2D.nodes[i].weight;
			}
			weights_offset += hist2D.nodes.size();
		}
	}
}

void unpack_weights(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, long double histogram_weights[], const Parameters &params) {
	int weights_offset = 0;
	for (int w = 0; w < params.weight.size(); w++) {
		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			auto &hist1D = histograms1D[h1D + params.transform_1D.size() * w];
			
			for (int i = 0; i < hist1D.nodes.size(); i++) {
				hist1D.nodes[i].weight += histogram_weights[weights_offset + i];
			}
			weights_offset += hist1D.nodes.size();
		}
		
		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			auto &hist2D = histograms2D[h2D + params.transform_2D.size() * w];

			for (int i = 0; i < hist2D.nodes.size(); i++) {
				hist2D.nodes[i].weight += histogram_weights[weights_offset + i];
			}
			weights_offset += hist2D.nodes.size();
		}
	}
}

void unpack_weights_clean(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, long double histogram_weights[], const Parameters &params) {
	int weights_offset = 0;
	for (int w = 0; w < params.weight.size(); w++) {
		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			auto &hist1D = histograms1D[h1D + params.transform_1D.size() * w];
			
			for (int i = 0; i < hist1D.nodes.size(); i++) {
				hist1D.nodes[i].weight = histogram_weights[weights_offset + i];
			}
			weights_offset += hist1D.nodes.size();
		}
		
		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			auto &hist2D = histograms2D[h2D + params.transform_2D.size() * w];

			for (int i = 0; i < hist2D.nodes.size(); i++) {
				hist2D.nodes[i].weight = histogram_weights[weights_offset + i];
			}
			weights_offset += hist2D.nodes.size();
		}
	}
}

void pack_weights_single_frame(const std::vector<BinaryTree> &histograms1D, const std::vector<Quadtree> &histograms2D, long double histogram_weights[], const Parameters &params) {
	int weights_offset = 0;
	for (int w = 0; w < params.weight.size(); w++) {
		for (int m = 0; m < params.merge_fraction_1D.size(); m++) {
			for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
				int hist_index = w + params.weight.size() * h1D + params.weight.size() * params.transform_1D.size() * m;
				const auto &hist1D = histograms1D[hist_index];
				
				for (int i = 0; i < hist1D.nodes.size(); i++) {
					histogram_weights[weights_offset + i] = hist1D.nodes[i].weight;
				}
				weights_offset += hist1D.nodes.size();
			}
		}
		
		for (int m = 0; m < params.merge_fraction_2D.size(); m++) {
			for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
				int hist_index = w + params.weight.size() * h2D + params.weight.size() * params.transform_2D.size() * m;
				const auto &hist2D = histograms2D[hist_index];

				for (int i = 0; i < hist2D.nodes.size(); i++) {
					histogram_weights[weights_offset + i] = hist2D.nodes[i].weight;
				}
				weights_offset += hist2D.nodes.size();
			}
		}
	}
}

void unpack_weights_single_frame(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, long double histogram_weights[], const Parameters &params) {
	int weights_offset = 0;
	for (int w = 0; w < params.weight.size(); w++) {
		for (int m = 0; m < params.merge_fraction_1D.size(); m++) {
			for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
				int hist_index = w + params.weight.size() * h1D + params.weight.size() * params.transform_1D.size() * m;
				auto &hist1D = histograms1D[hist_index];
				
				for (int i = 0; i < hist1D.nodes.size(); i++) {
					hist1D.nodes[i].weight += histogram_weights[weights_offset + i];
				}
				weights_offset += hist1D.nodes.size();
			}
		}
		
		for (int m = 0; m < params.merge_fraction_2D.size(); m++) {
			for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
				int hist_index = w + params.weight.size() * h2D + params.weight.size() * params.transform_2D.size() * m;
				auto &hist2D = histograms2D[hist_index];

				for (int i = 0; i < hist2D.nodes.size(); i++) {
					hist2D.nodes[i].weight += histogram_weights[weights_offset + i];
				}
				weights_offset += hist2D.nodes.size();
			}
		}
	}
}

void unpack_weights_clean_single_frame(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, long double histogram_weights[], const Parameters &params) {
	int weights_offset = 0;
	for (int w = 0; w < params.weight.size(); w++) {
		for (int m = 0; m < params.merge_fraction_1D.size(); m++) {
			for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
				int hist_index = w + params.weight.size() * h1D + params.weight.size() * params.transform_1D.size() * m;
				auto &hist1D = histograms1D[hist_index];
				
				for (int i = 0; i < hist1D.nodes.size(); i++) {
					hist1D.nodes[i].weight = histogram_weights[weights_offset + i];
				}
				weights_offset += hist1D.nodes.size();
			}
		}
		
		for (int m = 0; m < params.merge_fraction_2D.size(); m++) {
			for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
				int hist_index = w + params.weight.size() * h2D + params.weight.size() * params.transform_2D.size() * m;
				auto &hist2D = histograms2D[hist_index];

				for (int i = 0; i < hist2D.nodes.size(); i++) {
					hist2D.nodes[i].weight = histogram_weights[weights_offset + i];
				}
				weights_offset += hist2D.nodes.size();
			}
		}
	}
}

void collect_weights(std::vector<BinaryTree> &histograms1D, const std::vector<BinaryTree> &histograms1D_other, std::vector<Quadtree> &histograms2D, const std::vector<Quadtree> &histograms2D_other, const Parameters &params) {
	for (int w = 0; w < params.weight.size(); w++) {
		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			auto &hist1D = histograms1D[h1D + params.transform_1D.size() * w];
			const auto &hist1D_other = histograms1D_other[h1D + params.transform_1D.size() * w];
			
			for (int i = 0; i < hist1D.nodes.size(); i++) {
				hist1D.nodes[i].weight += hist1D_other.nodes[i].weight;
			}
		}

		int weight_offset;
		if (params.transform_1D.size() > 0) weight_offset = histograms1D[0].nodes.size() * params.transform_1D.size() * params.weight.size();
		else weight_offset = 0;
		
		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			auto &hist2D = histograms2D[h2D + params.transform_2D.size() * w];
			const auto &hist2D_other = histograms2D_other[h2D + params.transform_2D.size() * w];

			for (int i = 0; i < hist2D.nodes.size(); i++) {
				hist2D.nodes[i].weight += hist2D_other.nodes[i].weight;
			}
		}
	}
}

void save_histograms_single_frame(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, const Parameters &params, int frame, int rank, int Nthreads) {
	int r = 0;
	for (int w = 0; w < params.weight.size(); w++) {
		const std::string &postfix_weight = params.postfix_weight[w];

		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			const std::string &outputname = params.name_transform_1D[h1D];

			for (int m = 0; m < params.merge_fraction_1D.size(); m++) {
				if (r == rank) {
					int hist_index = w + params.weight.size() * h1D + params.weight.size() * params.transform_1D.size() * m;
					auto &hist1D = histograms1D[hist_index];
					const std::string &merge_fraction_postfix = params.merge_fraction_1D[m];

					if (params.saveTree) hist1D.save_structure(params.path_output, outputname + "_bintree_frame_" + std::to_string(frame) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
					if (params.saveBins) hist1D.save_leaves(params.path_output, outputname + "_bins_frame_" + std::to_string(frame) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
					if (params.saveCDF) hist1D.prob_to_pdf(params.path_output, outputname + "_CDF_frame_" + std::to_string(frame) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
					
				}
				r++;
				if (r == Nthreads) r = 0;
			}
		}
		
		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			const std::string &outputname = params.name_transform_2D[h2D];

			for (int m = 0; m < params.merge_fraction_2D.size(); m++) {
				if (r == rank) {
					int hist_index = w + params.weight.size() * h2D + params.weight.size() * params.transform_2D.size() * m;
					auto &hist2D = histograms2D[hist_index];
					const std::string &merge_fraction_postfix = params.merge_fraction_2D[m];

					if (params.saveTree) hist2D.save_structure(params.path_output, outputname + "_bintree_frame_" + std::to_string(frame) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
					if (params.saveBins) hist2D.save_leaves(params.path_output, outputname + "_bins_frame_" + std::to_string(frame) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
					if (params.saveCDF) hist2D.prob_to_pdf(params.path_output, outputname + "_CDF_frame_" + std::to_string(frame) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
				}
				r++;
				if (r == Nthreads) r = 0;
			}
		}
	}
}

void save_histograms(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, const Parameters &params, int rank, int Nthreads) {
	int r = 0;
	int frame_start = params.frames[0];
	int frame_end = params.frames.back();
	for (int w = 0; w < params.weight.size(); w++) {
		const std::string &postfix_weight = params.postfix_weight[w];
		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			if (r == rank) {
				auto &hist1D = histograms1D[h1D + params.transform_1D.size() * w];
				const std::string &outputname = params.name_transform_1D[h1D];

				if (params.saveTree) hist1D.save_structure(params.path_output, outputname + "_bintree_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_default.txt");
				if (params.saveBins) hist1D.save_leaves(params.path_output, outputname + "_bins_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_default.txt");
				if (params.saveCDF) hist1D.prob_to_pdf(params.path_output, outputname + "_CDF_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_default.txt");
				
				for (int m = 0; m < params.merge_fraction_1D.size(); m++){
					const std::string &merge_fraction_postfix = params.merge_fraction_1D[m];
					double merge_fraction = std::stod(params.merge_fraction_1D[m]);

					hist1D.merge(merge_fraction);

					if (params.saveTree) hist1D.save_structure(params.path_output, outputname + "_bintree_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
					if (params.saveBins) hist1D.save_leaves(params.path_output, outputname + "_bins_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
					if (params.saveCDF) hist1D.prob_to_pdf(params.path_output, outputname + "_CDF_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
				}
			}
			r++;
			if (r == Nthreads) r = 0;
		}
		
		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			if (r == rank) {
				auto &hist2D = histograms2D[h2D + params.transform_2D.size() * w];
				const std::string &outputname = params.name_transform_2D[h2D];

				if (params.saveTree) hist2D.save_structure(params.path_output, outputname + "_bintree_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_default.txt");
				if (params.saveBins) hist2D.save_leaves(params.path_output, outputname + "_bins_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_default.txt");
				if (params.saveCDF) hist2D.prob_to_pdf(params.path_output, outputname + "_CDF_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_default.txt");

				for (int m = 0; m < params.merge_fraction_2D.size(); m++){
					const std::string &merge_fraction_postfix = params.merge_fraction_2D[m];
					double merge_fraction = std::stod(params.merge_fraction_2D[m]);

					hist2D.merge(merge_fraction);

					if (params.saveTree) hist2D.save_structure(params.path_output, outputname + "_bintree_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
					if (params.saveBins) hist2D.save_leaves(params.path_output, outputname + "_bins_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
					if (params.saveCDF) hist2D.prob_to_pdf(params.path_output, outputname + "_CDF_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt");
				}
			}
			r++;
			if (r == Nthreads) r = 0;
		}
	}
}

void load_histograms(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, const Parameters &params, int &weights_count) {
	weights_count = 0;
	int frame_start = params.frames[0];
	int frame_end = params.frames.back();
	for (int w = 0; w < params.weight.size(); w++) {
		const std::string &postfix_weight = params.postfix_weight[w];

		for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
			const std::string &outputname = params.name_transform_1D[h1D];

			for (int m = 0; m < params.merge_fraction_1D.size(); m++) {
				int hist_index = w + params.weight.size() * h1D + params.weight.size() * params.transform_1D.size() * m;
				auto &hist1D = histograms1D[hist_index];
				const std::string &merge_fraction_postfix = params.merge_fraction_1D[m];
				std::string filename = outputname + "_bintree_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt";
				//std::cout << "Opening " << filename << std::endl;

				hist1D.load_tree(params.path_output, filename);
				
				weights_count += hist1D.nodes.size();
			}
		}
		
		for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
			const std::string &outputname = params.name_transform_2D[h2D];

			for (int m = 0; m < params.merge_fraction_2D.size(); m++){
				int hist_index = w + params.weight.size() * h2D + params.weight.size() * params.transform_2D.size() * m;
				auto &hist2D = histograms2D[hist_index];
				const std::string &merge_fraction_postfix = params.merge_fraction_2D[m];
				std::string filename = outputname + "_bintree_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + "_" + postfix_weight + "_" + merge_fraction_postfix + ".txt";
				//std::cout << "Opening " << filename << std::endl;

				hist2D.load_tree(params.path_output, filename);
				
				weights_count += hist2D.nodes.size();
			}
		}
	}
}

void write_bounds(long double bounds[], const Parameters &params) {
	for (int h1D = 0; h1D < params.transform_1D.size(); h1D++) {
		for (int w = 0; w < params.weight.size(); w++) {
			long double &xmin = bounds[0 + 2 * h1D];
			long double &xmax = bounds[1 + 2 * h1D];

			std::cout << params.name_transform_1D[h1D] << ", (" << params.postfix_weight[w] << "): (" << xmin << ", " << xmax << ")" << std::endl;
		}
	}
			
	int offset = 2 * params.transform_1D.size();

	for (int h2D = 0; h2D < params.transform_2D.size(); h2D++) {
		for (int w = 0; w < params.weight.size(); w++) {
			long double &xmin = bounds[offset + 0 + 4 * h2D];
			long double &xmax = bounds[offset + 1 + 4 * h2D];
			long double &ymin = bounds[offset + 2 + 4 * h2D];
			long double &ymax = bounds[offset + 3 + 4 * h2D];

			std::cout << params.name_transform_2D[h2D] << ", (" << params.postfix_weight[w] << "): (" << xmin << ", " << xmax << "), (" << ymin << ", " << ymax << ")" << std::endl;
		}
	}
}

void clear_histogram_weights(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D) {
	for (auto &h1D: histograms1D) {
		for (auto &n: h1D.nodes) n.weight = 0;
	}
	
	for (auto &h2D: histograms2D) {
		for (auto &n: h2D.nodes) n.weight = 0;
	}
}

void extract_histograms(Parameters params, int rank, int Nthreads) {
	if (params.extractGlobal) {
		if (rank == 0) std::cout << "Extracting global histograms from " << params.path << ", #frames: " << params.frames.size() << std::endl;
		delete_trailing(params.path, '/'); delete_trailing(params.path_output, '/');
		long long int Nstat = 0;
		int tag = 0; int tag_offset = 0;
		int hist_1D_count = 0;
		int hist_2D_count = 0;
		if (params.weight.size() == 0) {
			params.weight.push_back(weight_default); params.postfix_weight = { "V" };
		}
		if (params.name_transform_1D.size() == params.transform_1D.size()) {
			hist_1D_count = params.weight.size() * params.transform_1D.size();
		}
		if (params.name_transform_2D.size() == params.transform_2D.size()) {
			hist_2D_count = params.weight.size() * params.transform_2D.size();
		}

		int bounds_size = 2 * params.transform_1D.size() + 4 * params.transform_2D.size();
		int stats_size = 2 * params.stat_fun.size();

		long double bounds[bounds_size];

		int field_size = getSize(params.path + "/DD" + string_pad_left(params.frames[0], 0, 4) + "/data" + string_pad_left(params.frames[0], 0, 4) + ".cpu0000", "Grid00000001");
		long double field_x[field_size];
		long double field_y[field_size];

		long double stats[stats_size];
		long double stats_global[stats_size / 2];
		for (auto & s : stats) s = 0;
		for (auto & s : stats_global) s = 0;

		// __________________________________________ //
		// dealing with bounds                        //
		// at the same time, statistics are extracted //
		// __________________________________________ //

		initialize_bounds(bounds, params);
		for (int i = 0; i < params.single_frames.size(); i++) { // each thread adjusts its own bounds
			int frame = params.single_frames[i];
			std::vector<std::pair<std::string, int>> my_grids = distribute_grids(params.path, frame, rank, Nthreads);
			
			if (rank == 0) std::cout << "Adjusting bounds from frame " << frame << std::endl;

			for (auto & s : stats) s = 0; Nstat = 0;

			for (int j = 0; j < my_grids.size(); j++) {
				int cpu_number = my_grids[j].second;
				std::string filename = cpu_filename(params.path, frame, cpu_number);
				std::string groupname = my_grids[j].first;

				read_from_hdf5(filename, groupname, params.fieldname_x, field_x, field_size);
				read_from_hdf5(filename, groupname, params.fieldname_y, field_y, field_size);

				adjust_bounds(bounds, params, field_x, field_y, field_size);
				
				extract_stats(stats, params, field_x, field_y, field_size); Nstat += field_size;

				if (rank == 0) std::cout << j+1 << "/" << my_grids.size() << std::endl;
			}

			if (rank > 0) {
				tag = tag_offset + rank;
				
				MPI_Send(stats, stats_size, MPI_LONG_DOUBLE, 0, tag, MPI_COMM_WORLD);
				tag_offset += Nthreads; tag = tag_offset + rank;
				MPI_Send(&Nstat, 1, MPI_LONG_LONG_INT, 0, tag, MPI_COMM_WORLD);
				tag_offset += Nthreads;
			}

			if (rank == 0) {				
				for (int r = 1; r < Nthreads; r++) {
					tag = tag_offset + r;
					long double stats_other[stats_size];
					MPI_Recv(stats_other, stats_size, MPI_LONG_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					for (int i = 0; i < stats_size; i++) stats[i] += stats_other[i];
				}

				tag_offset += Nthreads;
				
				for (int r = 1; r < Nthreads; r++) { // the head receives the droneses boundses and adjusts the main bounds
					tag = tag_offset + r;
					long long int Nother;
					//std::cout << "Rank " << rank << " receiving bounds from rank " << r << " with tag " << tag << std::endl;
					MPI_Recv(&Nother, 1, MPI_LONG_LONG_INT, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					Nstat += Nother;
				}

				tag_offset += Nthreads;

				for (int i = 0; i < params.stat_fun.size(); i++) {
					long double param = (stats[2 * i] - stats[2 * i + 1]) / Nstat;
					stats_global[i] += param;
					std::string name = params.stat_name[i];
					std::ofstream out;
					out.open(params.path_output + "/" + name + "_frame_" + std::to_string(frame) + ".txt");
					out << param;
					out.close();
				}
			}
		}

		// _______________________________ //
		// exchanging the bounds and stats //
		// _______________________________ //

		if (rank > 0) { // drones send their individual bounds to the head
			tag = tag_offset + rank;
			//std::cout << "Rank " << rank << " sending bounds with tag " << tag << std::endl;
			MPI_Send(bounds, bounds_size, MPI_LONG_DOUBLE, 0, tag, MPI_COMM_WORLD);
			tag_offset += Nthreads;
		}
		if (rank == 0) {
			for (int r = 1; r < Nthreads; r++) { // the head receives the droneses boundses and adjusts the main bounds
				tag = tag_offset + r;
				long double bounds_other[bounds_size];
				//std::cout << "Rank " << rank << " receiving bounds from rank " << r << " with tag " << tag << std::endl;
				MPI_Recv(bounds_other, bounds_size, MPI_LONG_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				adjust_main_bounds(bounds, bounds_other, params);
			}

			tag_offset += Nthreads;

			for (int r = 1; r < Nthreads; r++) {
				tag = tag_offset + r; // the head re-distributes the (now global) bounds to everyone
				MPI_Send(bounds, bounds_size, MPI_LONG_DOUBLE, r, tag, MPI_COMM_WORLD);
			}

			tag_offset += Nthreads;
			
			write_bounds(bounds, params);

			for (int i = 0; i < params.stat_fun.size(); i++) {
				long double param = stats_global[i] / params.frames.size();
				std::string name = params.stat_name[i];
				int frame_start = params.single_frames[0]; int frame_end = params.single_frames[params.single_frames.size() - 1];
				std::ofstream out;
				out.open(params.path_output + "/" + name + "_frames_" + std::to_string(frame_start) + "-" + std::to_string(frame_end) + ".txt");
				out << param;
				out.close();
			}
		}
		if (rank > 0) { // everyone receives the global bounds
			tag = tag_offset + rank;
			MPI_Recv(bounds, bounds_size, MPI_LONG_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			tag_offset += Nthreads;
		}

		// _______________________________________ //
		// the real task - creating the histograms //
		// _______________________________________ //

		//std::cout << rank << ", bounds: " << bounds[0] << " " << bounds[1] << " " << bounds[2] << " " << bounds[3] << std::endl;

		if ((hist_1D_count > 0) || (hist_2D_count > 0)) {

			int weights_count;
			std::vector<BinaryTree> histograms1D; histograms1D.resize(hist_1D_count);
			std::vector<Quadtree> histograms2D; histograms2D.resize(hist_2D_count);
			initialize_histograms(histograms1D, histograms2D, bounds, params, weights_count);

			if (rank == 0) std::cout << "Size: " << histograms1D.size() << ", " << histograms2D.size() << std::endl;

			long double histogram_weights[weights_count];

			for (int i = 0; i < params.frames.size(); i++) { // each thread populates the bounds from its own subset

				int frame = params.frames[i];
				std::vector<std::pair<std::string, int>> my_grids = distribute_grids(params.path, frame, rank, Nthreads);
				
				if (rank == 0) std::cout << "Populating histograms from frame " << frame << std::endl;

				//std::cout << rank << " ";
				//for (int j = 0; j < my_grids.size(); j++) std::cout << my_grids[j].second << " ";
				//std::cout << std::endl;

				for (int j = 0; j < my_grids.size(); j++) {
					int cpu_number = my_grids[j].second;
					std::string filename = cpu_filename(params.path, frame, cpu_number);
					std::string groupname = my_grids[j].first;

					read_from_hdf5(filename, groupname, params.fieldname_x, field_x, field_size);
					read_from_hdf5(filename, groupname, params.fieldname_y, field_y, field_size);

					populate_histograms(histograms1D, histograms2D, field_x, field_y, field_size, params);
					
					//if (histograms2D[0].nodes[93183].weight > 0.5) std::cout << my_grids[j].second << ", " << histograms2D[0].nodes[93183].weight << std::endl;

					if (rank == 0) std::cout << j+1 << "/" << my_grids.size() << std::endl;
				}
			}

			//std::cout << rank << ": Preparing to exchange histograms, tag offset = " << tag_offset << ", #weights = " << weights_count << std::endl;

			// ____________________________ //
			// saving the global histograms //
			// ____________________________ //

			//std::cout << "rank: " << rank << ", " << histograms2D[0].nodes[93183].weight << ", " << histograms2D[0].nodes[93807].weight << std::endl;
			
			if (rank > 0) {
				pack_weights(histograms1D, histograms2D, histogram_weights, params);

				tag = tag_offset + rank;

				MPI_Send(histogram_weights, weights_count, MPI_LONG_DOUBLE, 0, tag, MPI_COMM_WORLD);
				tag_offset += Nthreads; tag = tag_offset + rank;

				//std::cout << "Rank " << rank << " sent weights to 0" << std::endl;

				MPI_Recv(histogram_weights, weights_count, MPI_LONG_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				unpack_weights_clean(histograms1D, histograms2D, histogram_weights, params);

				//std::cout << "Rank " << rank << " received weights from 0" << std::endl;
					
				tag_offset += Nthreads;
			}

			if (rank == 0) {

				for (int r = 1; r < Nthreads; r++) {
					tag = tag_offset + r;
					MPI_Recv(histogram_weights, weights_count, MPI_LONG_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

					//std::cout << "Rank 0 received weights from " << r << std::endl;

					unpack_weights(histograms1D, histograms2D, histogram_weights, params);
				}
					
				tag_offset += Nthreads;
				
				pack_weights(histograms1D, histograms2D, histogram_weights, params);

				for (int r = 1; r < Nthreads; r++) {
					tag = tag_offset + r;
					MPI_Send(histogram_weights, weights_count, MPI_LONG_DOUBLE, r, tag, MPI_COMM_WORLD);

					//std::cout << "Rank 0 sent weights to " << r << std::endl;
				}

				tag_offset += Nthreads;
			}

			//std::cout << "rank: " << rank << ", " << histograms2D[0].nodes[93183].weight << ", " << histograms2D[0].nodes[93807].weight << std::endl;

			save_histograms(histograms1D, histograms2D, params, rank, Nthreads);

			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
}

void extract_frame_histograms(Parameters params, int rank, int Nthreads) {
	if (params.extractSingle) {
		if (rank == 0) std::cout << "Extracting single frames from " << params.path << ", #frames: " << params.single_frames.size() << std::endl;
		delete_trailing(params.path, '/'); delete_trailing(params.path_output, '/');
		int tag = 0; int tag_offset = 0;
		long long int Nstat = 0;
		int hist_1D_count = 0;
		int hist_2D_count = 0;
		if (params.weight.size() == 0) {
			params.weight.push_back(weight_default); params.postfix_weight = { "V" };
		}
		if (params.name_transform_1D.size() == params.transform_1D.size()) {
			hist_1D_count = params.weight.size() * params.merge_fraction_1D.size() * params.transform_1D.size();
		}
		if (params.name_transform_2D.size() == params.transform_2D.size()) {
			hist_2D_count = params.weight.size() * params.merge_fraction_2D.size() * params.transform_2D.size();
		}

		int field_size = getSize(params.path + "/DD" + string_pad_left(params.single_frames[0], 0, 4) + "/data" + string_pad_left(params.single_frames[0], 0, 4) + ".cpu0000", "Grid00000001");
		if (rank == 0) std::cout << "Field count per frame: " << field_size << std::endl;
		long double field_x[field_size];
		long double field_y[field_size];

		// _________________________ //
		// load the global histogram //
		// _________________________ //

		int weights_count;
		std::vector<BinaryTree> histograms1D; histograms1D.resize(hist_1D_count);
		std::vector<Quadtree> histograms2D; histograms2D.resize(hist_2D_count);
		
		load_histograms(histograms1D, histograms2D, params, weights_count);

		long double histogram_weights[weights_count];

		for (int i = 0; i < params.single_frames.size(); i++) { // each thread populates the bounds from its own subset
			int frame = params.single_frames[i];
			std::vector<std::pair<std::string, int>> my_grids = distribute_grids(params.path, frame, rank, Nthreads);

			clear_histogram_weights(histograms1D, histograms2D);
			
			if (rank == 0) std::cout << "Populating histograms from frame " << frame << std::endl;

			for (int j = 0; j < my_grids.size(); j++) {
				int cpu_number = my_grids[j].second;
				std::string filename = cpu_filename(params.path, frame, cpu_number);
				std::string groupname = my_grids[j].first;

				read_from_hdf5(filename, groupname, params.fieldname_x, field_x, field_size);
				read_from_hdf5(filename, groupname, params.fieldname_y, field_y, field_size);
			
				if (rank == 0) std::cout << j+1 << "/" << my_grids.size() << "...";

				populate_histograms_single_frame(histograms1D, histograms2D, field_x, field_y, field_size, params);

				if (rank == 0) std::cout << "done!" << std::endl;
			}

			// ______________________________________________ //
			// exchanging and saving the per-frame histograms //
			// ______________________________________________ //

			if (rank > 0) {
				tag = tag_offset + rank;
				pack_weights_single_frame(histograms1D, histograms2D, histogram_weights, params);
				MPI_Send(histogram_weights, weights_count, MPI_LONG_DOUBLE, 0, tag, MPI_COMM_WORLD);
				tag_offset += Nthreads;
			}
			if (rank == 0) {
				std::cout << "Collecting single frame histograms..." << frame << std::endl;

				for (int r = 1; r < Nthreads; r++) {
					tag = tag_offset + r;

					MPI_Recv(histogram_weights, weights_count, MPI_LONG_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					unpack_weights_single_frame(histograms1D, histograms2D, histogram_weights, params);
				}

				tag_offset += Nthreads;
				
				pack_weights_single_frame(histograms1D, histograms2D, histogram_weights, params);

				std::cout << "Sending weights back to drones" << std::endl;
				for (int r = 1; r < Nthreads; r++) {
					tag = tag_offset + r;
					MPI_Send(histogram_weights, weights_count, MPI_LONG_DOUBLE, r, tag, MPI_COMM_WORLD);
				}
				tag_offset += Nthreads;
			}
			if (rank > 0) {
				tag = tag_offset + rank;

				MPI_Recv(histogram_weights, weights_count, MPI_LONG_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				unpack_weights_clean_single_frame(histograms1D, histograms2D, histogram_weights, params);
					
				tag_offset += Nthreads;

			}
			// ______________________________________________
			// every thread saves only a subset of histograms
			// ______________________________________________

			save_histograms_single_frame(histograms1D, histograms2D, params, frame, rank, Nthreads);

			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
}