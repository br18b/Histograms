#include <iostream>
#include <functional>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "mpi.h"
#include "node.h"
#include "load.h"
#include "string_pad.h"
#include "statistics.h"

//using namespace H5;

void extract_fields(std::string path_to_sim, std::string path_output, std::vector<int> frames, std::string fieldname_x, std::string fieldname_y, std::vector< std::function<double(double, double)>> weight_fun, std::vector<std::function<double(double, double)>> field_1D, std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> field_2D, std::vector<std::string> scale_1D, std::vector<std::pair<std::string, std::string>> scale_2D, std::vector<int> subdiv_1D, std::vector<int> subdiv_2D, std::vector<std::string> merge_1D, std::vector<std::string> merge_2D, std::vector<std::string> outputname_x, std::vector<std::string> outputname_z, std::vector<std::string> outputpostfix_weight) {
	std::vector<BinaryTree> data_1D;
	std::vector<Quadtree> data_2D;
	data_1D.resize(std::max(1, (int)weight_fun.size()) * field_1D.size());
	data_2D.resize(std::max(1, (int)weight_fun.size()) * field_2D.size());
	initialize(path_to_sim, frames, fieldname_x, fieldname_y, field_1D, field_2D, scale_1D, scale_2D, data_1D, data_2D);
	for (int w = 0; w < weight_fun.size(); w++) {
		for (int k = 0; k < subdiv_1D.size(); k++) {
			data_1D[w * subdiv_1D.size() + k].uniform_divide(subdiv_1D[k]);
		}
		for (int k = 0; k < subdiv_2D.size(); k++) {
			data_2D[w * subdiv_2D.size() + k].uniform_divide(subdiv_2D[k]);
		}
	}
	load_points(path_to_sim, frames, fieldname_x, fieldname_y, weight_fun, field_1D, field_2D, data_1D, data_2D);
	for (int w = 0; w < weight_fun.size(); w++) {
		for (int k = 0; k < field_1D.size(); k++) {
			for (int m = 0; m < merge_1D.size(); m++) {
				data_1D[w * field_1D.size() + k].merge(std::stod(merge_1D[m]));
				if (frames.front() > frames.back()) {
					data_1D[w * field_1D.size() + k].save_structure(path_output, outputname_x[k] + "_bintree_frames_" + std::to_string(frames.front()) + "-" + std::to_string(frames.back()) + "_" + outputpostfix_weight[w] + "_" + merge_1D[m] + ".txt");
					data_1D[w * field_1D.size() + k].save_leaves(path_output, outputname_x[k] + "_bins_frames_" + std::to_string(frames.front()) + "-" + std::to_string(frames.back()) + "_" + outputpostfix_weight[w] + "_" + merge_1D[m] + ".txt");
					data_1D[w * field_1D.size() + k].prob_to_pdf(path_output, outputname_x[k] + "_CDF_frames_" + std::to_string(frames.front()) + "-" + std::to_string(frames.back()) + "_" + outputpostfix_weight[w] + "_" + merge_1D[m] + ".txt");
				}
				else {
				}
			}
		}
		for (int k = 0; k < field_2D.size(); k++) {
			for (int m = 0; m < merge_2D.size(); m++) {
				data_2D[w * field_2D.size() + k].merge(std::stod(merge_2D[m]));
				if (frames.front() > frames.back()) {
					data_2D[w * field_2D.size() + k].save_structure(path_output, outputname_z[k] + "_bintree_frames_" + std::to_string(frames.front()) + "-" + std::to_string(frames.back()) + "_" + outputpostfix_weight[w] + "_" + merge_2D[m] + ".txt");
					data_2D[w * field_2D.size() + k].save_leaves(path_output, outputname_z[k] + "_bins_frames_" + std::to_string(frames.front()) + "-" + std::to_string(frames.back()) + "_" + outputpostfix_weight[w] + "_" + merge_2D[m] + ".txt");
					data_2D[w * field_2D.size() + k].prob_to_pdf(path_output, outputname_z[k] + "_CDF_frames_" + std::to_string(frames.front()) + "-" + std::to_string(frames.back()) + "_" + outputpostfix_weight[w] + "_" + merge_2D[m] + ".txt");
				}
				else {
					data_2D[w * field_2D.size() + k].save_structure(path_output, outputname_z[k] + "_bintree_frame_" + std::to_string(frames.front()) + "_" + outputpostfix_weight[w] + "_" + merge_2D[m] + ".txt");
					data_2D[w * field_2D.size() + k].save_leaves(path_output, outputname_z[k] + "_bins_frame_" + std::to_string(frames.front()) + "_" + outputpostfix_weight[w] + "_" + merge_2D[m] + ".txt");
					data_2D[w * field_2D.size() + k].prob_to_pdf(path_output, outputname_z[k] + "_CDF_frame_" + std::to_string(frames.front()) + "_" + outputpostfix_weight[w] + "_" + merge_2D[m] + ".txt");
				}
			}
		}
	}
}

std::function<double(double, double, double)> magnitude = [](double x, double y, double z) {
	return std::sqrt(x * x + y * y + z * z);
};

std::function<double(double, double, double)> theta = [](double x, double y, double z) {
	if (x * x + y * y + z * z > 0) {
		return std::acos(z / std::sqrt(x * x + y * y + z * z));
	}
	else return 1e10;
};

std::function<double(double, double, double)> phi = [](double x, double y, double z) {
	if (x * x + y * y > 0) {
		return std::atan2(y, x);
	}
	else return 1e10;
};

std::function<bool(double, double, double)> angle_condition = [](double x, double y, double z) {
	return (x * x + y * y + z * z > 0) && (x * x + y * y > 0);
};

std::function<bool(double, double, double, double, double)> angle_magnitude_condition = [](double x, double y, double z, double min, double max) {
	return (x * x + y * y + z * z > min * min) && (x * x + y * y + z * z <= max * max) && (x * x + y * y > 0);
};

std::function<double(double, double)> weight_volume = [](double x, double y) {
	return 1;
};

std::function<double(double, double)> weight_mass_from_rho = [](double x, double y) {
	return x;
};

std::function<double(double, double)> weight_mass_from_s = [](double x, double y) {
	return std::exp(x);
};

std::function<double(double, double)> weight_kinetic_from_rho = [](double x, double y) {
	return 0.5 * x * y * y;
};

std::function<double(double, double)> weight_kinetic_from_s = [](double x, double y) {
	return 0.5 * std::exp(x) * y * y;
};

std::function<double(double, double)> fun_logrho_from_rho = [](double x, double y) {
	return std::log(x);
};

std::function<double(double, double)> fun_logrho_from_s = [](double x, double y) {
	return x;
};

std::function<double(double, double)> fun_absv = [](double x, double y) {
	return y;
};

std::function<double(double, double)> fun_Helmholtz_from_rho = [](double x, double y) {
	return x * std::log(x);
};

std::function<double(double, double)> fun_Helmholtz_from_s = [](double x, double y) {
	return x * std::exp(x);
};

double fx(const double& x, const double& threshold) {
	if (x >= 0) {
		return std::log10((x + threshold) / threshold);
	}
	else {
		return -std::log10((-x + threshold) / threshold);
	}
}

std::function<double(double, double)> fun_symlog_Helmholtz_from_rho = [](double x, double y) {
	return fx(x * std::log(x), 0.3);
};

std::function<double(double, double)> fun_symlog_Helmholtz_from_s = [](double x, double y) {
	return fx(x * std::exp(x), 0.3);
};

std::function<double(double, double)> fun_shiftlog_Helmholtz_from_rho = [](double rho, double absv) {
	double res;
	if ((rho > 0.367512) && (rho < 0.368248)) {
		res = 0.3068528194400547 + 2 * std::log(rho + 0.3678794411714423 + 1E-16);
	}
	else res = 0.3678794411714424 + rho * std::log(rho);
	/*
	if (res <= 0) {
		std::cout << "problem!" << std::endl;
		std::cout << logrho << std::endl;
		std::cout << logrho * std::exp(logrho) << std::endl;
	}*/
	return std::log(res);
};

std::function<double(double, double)> fun_shiftlog_Helmholtz_from_s = [](double logrho, double absv) {
	double res;
	if ((logrho > -1.005) && (logrho < -0.995)) {
		res = -1.693147180559945 + 2 * std::log(logrho + 1 + 1E-16);
	}
	else res = 0.3678794411714424 + logrho * std::exp(logrho);
	/*
	if (res <= 0) {
		std::cout << "problem!" << std::endl;
		std::cout << logrho << std::endl;
		std::cout << logrho * std::exp(logrho) << std::endl;
	}*/
	return std::log(res);
};

std::function<double(double, double)> fun_kinetic_from_rho = [](double x, double y) {
	return 0.5 * x * y * y;
};

std::function<double(double, double)> fun_kinetic_from_s = [](double x, double y) {
	return 0.5 * std::exp(x) * y * y;
};

std::function<double(double, double)> fun_log_kinetic_from_rho = [](double x, double y) {
	return std::log(0.5 * x * y * y);
};

std::function<double(double, double)> fun_log_kinetic_from_s = [](double x, double y) {
	return std::log(0.5 * y * y) + x;
};

std::function<double(double, double)> fun_density_from_rho = [](double x, double y) {
	return x;
};

std::function<double(double, double)> fun_density_from_s = [](double x, double y) {
	return std::exp(x);
};

int main(int argc, char *argv[]) {
	int Nthreads, rank;
	MPI_Init(&argc, &argv);
    //MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &Nthread);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Nthreads);
	if (rank == 0) std::cout << "Nthreads = " << Nthreads << std::endl;
	//std::cout << "rank = " << rank << std::endl;
	//return 0;
	std::string path = "/scratch/06736/rabatinb/projects/1024_Mach_grid/xi_0_mach12.4/";
	std::string path_out = path + "stats";
	std::string filename = "DD0052/data0052.cpu0000";
	std::string ifile = path + filename;
	std::vector<std::function<double(double, double)>> stats = { density, logrho, rho_logrho, kin_logrho, logrho2, rho_logrho2, kin_logrho2, v2V, v2M, v2E, kinetic, thermal, fun_sv2, fun_s2v2, fun_sv4, fun_s2v4 };
	int index_rho = 0;
	int index_s = 1;
	int index_sM = 2;
	int index_sE = 3;
	int index_s2 = 4;
	int index_s2M = 5;
	int index_s2E = 6;
	int index_v2 = 7;
	int index_v2M = 8;
	int index_v2E = 9;
	int index_kin = 10;

	//auto foo = read_from_hdf5("/scratch/06736/rabatinb/projects/1024_Mach_grid/xi_0_mach12.4/DD0052/data0052.cpu0000", "Grid00000001", "density");
	//return 0;
/*
	extract_stats(path, path_out, { 52 }, stats,
		{ { M1DV, { index_v2 } }, { M1DM, { index_v2M, index_rho } }, { M1DE, { index_v2E, index_kin } },
		{ muV, { index_s } }, {muM, { index_rho, index_sM } }, { muE, { index_kin, index_sE } },
		{ sigmaV, { index_s, index_s2 } }, {sigmaM, { index_rho, index_sM, index_s2M } }, 
		{ sigmaE, { index_kin, index_sE, index_s2E } } }, { "rho", "sV", "sM", "sE", "s2V", "s2M", "s2E", "v2V", "v2M", "v2E", "kinetic", "thermal", "sv2", "s2v2", "sv4", "s2v4" },
		{ "M1DV", "M1DM", "M1DE", "muV", "muM", "muE", "sigmaV", "sigmaM", "sigmaE" });
*/
	//extract_stats_parallel(path, path_out, { 52 }, { density }, { }, { "rho" }, { }, Nthreads, rank);
	//extract_fields(path, path_out, { 52 }, "rho", "v", { weight_volume, weight_mass_from_rho, weight_kinetic_from_rho }, { fun_absv, fun_logrho_from_rho, fun_symlog_Helmholtz_from_rho, fun_log_kinetic_from_rho }, { {fun_logrho_from_rho, fun_absv}, {fun_symlog_Helmholtz_from_rho, fun_log_kinetic_from_rho} }, { "lin" ,"lin" ,"lin" ,"lin" }, { {"lin" ,"lin"} ,{"lin" ,"lin"} }, { 10,10,10,10 }, { 10,10 }, { "1e-4", "1e-3", "1e-2" }, { "1e-7", "1e-6", "1e-5", "1e-4", "1e-3", "1e-2" }, { "absv","logrho","Helmholtz","kinetic" }, { "logrho_absv","Helmholtz_kinetic" }, { "V", "M", "E" });
	extract_1D_histogram_parallel(path, path_out, { 52 }, "rho", "absv", weight_volume, fun_logrho_from_rho, 10, { "1e-4", "1e-3", "1e-2"}, "logrho_" + std::to_string(Nthreads), "V", rank, Nthreads);
}