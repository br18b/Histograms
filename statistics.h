#ifndef __STAT__
#define __STAT__

#include <iostream>
#include <vector>
#include <set>
#include <functional>
#include <cmath>
#include "mpi.h"
#include "read_binary.h"
#include "parallel.h"
#include "node.h"

extern std::function<double(double, double)> weight_default;

extern std::function<double(double, double)> v2V;
extern std::function<double(double, double)> v2M;
extern std::function<double(double, double)> v2E;
extern std::function<double(double, double)> density;
extern std::function<double(double, double)> kinetic;
extern std::function<double(double, double)> thermal;
extern std::function<double(double, double)> v_mag;
extern std::function<double(double, double)> logrho;
extern std::function<double(double, double)> logrho2;
extern std::function<double(double, double)> rho_logrho;
extern std::function<double(double, double)> rho_logrho2;
extern std::function<double(double, double)> kin_logrho;
extern std::function<double(double, double)> kin_logrho2;
extern std::function<double(double, double)> fun_sv;
extern std::function<double(double, double)> fun_sv2;
extern std::function<double(double, double)> fun_sv3;
extern std::function<double(double, double)> fun_sv4;
extern std::function<double(double, double)> fun_s2v;
extern std::function<double(double, double)> fun_s2v2;
extern std::function<double(double, double)> fun_s2v3;
extern std::function<double(double, double)> fun_s2v4;

extern std::function<double(std::vector<double>, std::vector<int>)> M1DV;
extern std::function<double(std::vector<double>, std::vector<int>)> M1DM;
extern std::function<double(std::vector<double>, std::vector<int>)> M1DE;
extern std::function<double(std::vector<double>, std::vector<int>)> muV;
extern std::function<double(std::vector<double>, std::vector<int>)> sigmaV;
extern std::function<double(std::vector<double>, std::vector<int>)> muM;
extern std::function<double(std::vector<double>, std::vector<int>)> sigmaM;
extern std::function<double(std::vector<double>, std::vector<int>)> muE;
extern std::function<double(std::vector<double>, std::vector<int>)> sigmaE;

void resolve_set(std::set<double> &temp_set, double &accumulator, bool override);
void extract_stats(std::string path_to_sim, std::string path_output, std::vector<int> frames, std::vector<std::function<double(double, double)>> stat_funs, std::vector<std::pair<std::function<double(std::vector<double>, std::vector<int>)>, std::vector<int>>> stat_aggregate, std::vector<std::string> stat_filenames, std::vector<std::string> stat_aggregate_filenames);
void extract_stats_parallel(std::string path_to_sim, std::string path_output, std::vector<int> frames, std::vector<std::function<double(double, double)>> stat_funs, std::vector<std::pair<std::function<double(std::vector<double>, std::vector<int>)>, std::vector<int>>> stat_aggregate, std::vector<std::string> stat_filenames, std::vector<std::string> stat_aggregate_filenames, int Nthread, int rank);

void extract_1D_histogram_parallel(std::string path_to_sim, std::string path_output, std::vector<int> frames, std::string fieldname_x, std::string fieldname_y, std::function<double(double, double)> weight_fun, std::function<double(double, double)> field_1D, int depth, std::vector<std::string> merge_1D, std::string outputname_x, std::string outputpostfix_weight, int rank, int Nthreads);

void initialize_bounds(double bounds[], const Parameters &params);
void adjust_bounds(double bounds[], const Parameters &params, double field_x[], double field_y[], int field_size);
void adjust_main_bounds(double bounds[], double bounds_other[], const Parameters &params);
void write_bounds(double bounds[], const Parameters &params);
void initialize_histograms(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, double bounds[], const Parameters &params, int &weights_count);
void populate_histograms(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, double field_x[], double field_y[], int field_size, const Parameters &params);
void pack_weights(const std::vector<BinaryTree> &histograms1D, const std::vector<Quadtree> &histograms2D, double histogram_weights[], const Parameters &params);
void unpack_weights(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, double histogram_weights[], const Parameters &params);
void save_histograms_single_frame(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, const Parameters &params, std::string path_output, int frame);
void save_histograms(std::vector<BinaryTree> &histograms1D, std::vector<Quadtree> &histograms2D, const Parameters &params, std::string path_output);
void extract_histograms(std::string path_to_sim, std::string path_output, Parameters params, int rank, int Nthreads);

#endif
