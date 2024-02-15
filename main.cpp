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
#include "read_binary.h"

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

std::function<double(double, double)> logX = [](double x, double y) {
	return std::log(x);
};

std::function<double(double, double)> log2X = [](double x, double y) {
	double foo = std::log(x);
	return foo * foo;
};

std::function<double(double, double)> idX = [](double x, double y) {
	return x;
};

std::function<double(double, double)> x2 = [](double x, double y) {
	return x * x;
};

std::function<double(double, double)> logY = [](double x, double y) {
	return std::log(x);
};

std::function<double(double, double)> log2Y = [](double x, double y) {
	double foo = std::log(y);
	return std::log(x);
};

std::function<double(double, double)> idY = [](double x, double y) {
	return y;
};

std::function<double(double, double)> y2 = [](double x, double y) {
	return y * y;
};

std::function<double(double, double)> Helmholtz = [](double x, double y) {
	return x * std::log(x);
};

std::function<double(double, double)> Helmholtz2 = [](double x, double y) {
	double foo = x * std::log(x);
	return foo * foo;
};

std::function<double(double, double)> Helmholtz_s = [](double x, double y) {
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

std::function<double(double, double)> EK = [](double x, double y) {
	return 0.5 * x * y * y;
};

std::function<double(double, double)> kinetic2 = [](double x, double y) {
	double foo = 0.5 * x * y * y;
	return foo * foo;
};

std::function<double(double, double)> kinetic_s = [](double x, double y) {
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

int main(int argc, char *argv[]) { // mpirun -np 64 main
	int Nthreads, rank;
	std::string param_filename;
	Parameters params;
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &Nthreads);
	if (rank == 0) {
		std::cout << "Program launched with " << Nthreads;
		if (Nthreads == 1) std::cout << " thread" << std::endl;
		else std::cout << " threads" << std::endl;
	}
	if (argc > 1) param_filename = argv[1];
	if (rank == 0) std::cout << "Initializing parameters (input file " << param_filename << ")" << std::endl;
	params.initialize(param_filename);

		extract_histograms(params, rank, Nthreads);

	params.saveTree = false;
	
	extract_frame_histograms(params, rank, Nthreads);

	MPI_Barrier(MPI_COMM_WORLD);
	return 0;
}