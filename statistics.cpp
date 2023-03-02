#include "statistics.h"

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

void extract_1D_histogram_parallel(std::string path_to_sim, std::string path_output, std::vector<int> frames, std::string fieldname_x, std::string fieldname_y, std::function<double(double, double)> weight_fun, std::function<double(double, double)> field_1D, int depth, std::vector<std::string> merge_1D, std::string outputname_x, std::string outputpostfix_weight, int rank, int Nthreads) {
	delete_trailing(path_to_sim, '/'); delete_trailing(path_output, '/');
	
	for (int i = 0; i < frames.size(); i++) {
		int frame = frames[i];
		std::cout << "Current frame: " << frame << std::endl;
		std::string my_grids; int filename_size, groupname_size;
		int my_grid_size;
		if (rank == 0) {
			std::vector<std::vector<std::string>> cpu_files;
			cpu_files = get_cpu_list(path_to_sim, frame);
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
		}
		if (rank > 0) {
			MPI_Recv(&my_grid_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&filename_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&groupname_size, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			int length = (filename_size + groupname_size) * my_grid_size;
			char* grids_array = new char[length];
			my_grids.resize(length);
			MPI_Recv(grids_array, length, MPI_CHAR, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int k = 0; k < length; k++) {
				my_grids[k] = grids_array[k];
			}
			delete[] grids_array;
		}

		int n = 0;
			
		std::cout << my_grid_size << std::endl;

		double* bounds = new double[2];
		bounds[0] = 1e10; // xmin
		bounds[1] = -1e10; // xmax

		// __________________________________
		// finding bounds of the 1D histogram
		// __________________________________

		for (int j = 0; j < my_grid_size; j++) {
			std::string filename, groupname;
			unpack_grids(filename, groupname, my_grids, filename_size, groupname_size, j);
			delete_trailing(filename, ' '); delete_trailing(groupname, ' ');
			
			if (rank == 0) std::cout << "Rank 0 reading from " << groupname << std::endl;
			
			std::vector<double> field_x = read_from_hdf5(filename, groupname, fieldname_x);
			std::vector<double> field_y = read_from_hdf5(filename, groupname, fieldname_y);
			
			if (rank == 0) std::cout << "Adjusting bounds" << std::endl;

			// ________________________________________________
			// each rank will compute bounds for its own subset
			// ________________________________________________

			for (int i = 0; i < field_x.size(); i++) {
				double val = field_1D(field_x[i], field_y[i]);
				if (bounds[0] > val) bounds[0] = val;
				if (bounds[1] < val) bounds[1] = val;
			}
		}

		// ____________________________________________________________
		// all ranks except the root will send their bounds to the root
		// ____________________________________________________________

		if (rank > 0) {
			int tag = 4 + rank;
			MPI_Send(bounds, 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
			//std::cout << "I'm rank " << rank << " and sent bounds (" << bounds[0] << ", " << bounds[1] << ") to rank 0 with tag " << tag <<std::endl;
		}
		//MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0) {
			for (int r = 1; r < Nthreads; r++) {
				double* bounds_other = new double[2];
				int tag = 4 + r;

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

			// ____________________________________________
			// the root sends the global bounds to everyone
			// ____________________________________________

			for (int r = 1; r < Nthreads; r++) {
				int tag = 4 + Nthreads + r;
				//std::cout << "I'm rank 0 and I'm sending bounds to rank " << r << " with tag " << tag <<std::endl;
				MPI_Send(bounds, 2, MPI_DOUBLE, r, tag, MPI_COMM_WORLD);
			}
			//std::cout << "xmin = " << xmin << ", xmax = " << xmax << std::endl;
		}
		//MPI_Barrier(MPI_COMM_WORLD);

		if (rank > 0) {
			int tag = 4 + Nthreads + rank;
			//std::cout << "I'm rank " << rank << " and I'm receiving bounds from rank 0 with tag " << tag <<std::endl;
			MPI_Recv(bounds, 2, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		
		//std::cout << rank << " xmin = " << xmin << ", xmax = " << xmax << std::endl;

		// ______________________________________
		// now we actually compute the histograms
		// ______________________________________
		
		BinaryTree histogram1D;
		histogram1D.initialize(bounds[0], bounds[1], "lin");
		histogram1D.uniform_divide(depth);

		delete[] bounds;

		for (int j = 0; j < my_grid_size; j++) {
			std::string filename, groupname;
			unpack_grids(filename, groupname, my_grids, filename_size, groupname_size, j);
			delete_trailing(filename, ' '); delete_trailing(groupname, ' ');
			
			if (rank == 0) std::cout << "Rank 0 reading from " << groupname << std::endl;
			
			std::vector<double> field_x = read_from_hdf5(filename, groupname, fieldname_x);
			std::vector<double> field_y = read_from_hdf5(filename, groupname, fieldname_y);

			//std::cout << rank << " " << filename << " " << groupname << std::endl;
			
			if (rank == 0) std::cout << "Filling histograms" << std::endl;

			histogram1D.load_points(field_x, field_y, field_1D, weight_fun);

			// actual points are stored in BinaryTree within Node1D.weight - public interface.
			// all histograms have the same structure (created with uniform_divide)
			// only thing that needs to be accummulated are the weights from different drones
			// size doesn't need to be sent, as each drone is working with the same size, root included
		}
		if (rank > 0) {
			int tag = 4 + 2 * Nthreads + rank;
			int nodesize = histogram1D.nodes.size();
			double* my_weights = new double[nodesize];

			for (int i = 0; i < nodesize; i++) my_weights[i] = histogram1D.nodes[i].weight;
			
			std::cout << "Rank " << rank << " is sending weights... " << my_weights[50] << std::endl;

			MPI_Send(my_weights, nodesize, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);

			delete[] my_weights;
		}
		
		if (rank == 0) std::cout << "The head core is waiting for others... " << std::endl;
		//MPI_Barrier(MPI_COMM_WORLD);

		// ________________________________________________
		// collecting all histograms and combining into one
		// ________________________________________________

		if (rank == 0) {
			std::cout << "Collecting histograms ... " << std::endl;
			for (int r = 1; r < Nthreads; r++) {
				int tag = 4 + 2 * Nthreads + r;
				int nodesize = histogram1D.nodes.size(); // same nodesize everywhere, no need to communicate this
				double* drone_weights = new double[nodesize];
				MPI_Recv(drone_weights, nodesize, MPI_DOUBLE, r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				for (int i = 0; i < nodesize; i++) {
					histogram1D.nodes[i].weight += drone_weights[i];
				}
				delete[] drone_weights;
				std::cout << "Rank " << std::to_string(r) << " collected." << std::endl;
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
}