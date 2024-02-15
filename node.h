#pragma once
#ifndef NODE
#define NODE

#include <vector>
#include <functional>
#include <tuple>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iterator>
#include <cmath>
#include "read_binary.h"
#include "string_pad.h"

double symlog(double x, double threshold);
double symlogInv(double x, double threshold);
double natural_log(double x);
double dec_log(double x);
double lin(double x);

int parse_token(std::vector<std::string> &temp, int index, const char &delimiter);
std::vector<std::string> parse_line(std::string line, std::vector<char> delimiters);
std::stringstream stream_line(std::string line, std::vector<char> delimiters);

void parse_function(std::function<double(double, double)> &func, const std::string &name);
std::function<double(double, double)> parse_function(const std::string &name);

class Parameters {
public:
	std::string path, path_output;
	std::vector<int> frames;
	std::vector<int> single_frames;
	std::string fieldname_x, fieldname_y;

	std::vector<std::function<double(double, double)>> stat_fun;
	std::vector<std::string> stat_name;

	std::vector<std::function<double(double, double)>> transform_1D;
	std::vector<std::string> name_transform_1D;
	std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> transform_2D;
	std::vector<std::string> name_transform_2D;
	std::vector<std::function<double(double, double)>> weight;
	std::vector<std::string> postfix_weight;
	std::vector<std::string> merge_fraction_1D;
	std::vector<std::string> merge_fraction_2D;
	std::vector<int> initial_depth_1D;
	std::vector<int> initial_depth_2D;
	bool saveBins, saveCDF, saveTree, saveIndividualFrames, extractGlobal, extractSingle;
	void initialize();
	void initialize(std::string path);
	void addFrame(int f);
	void addSingleFrame(int f);
	void setFrames(int start, int end);
	void setSingleFrames(int start, int end);
	void setFieldnames(std::string fname_x, std::string fname_y);
	
	void addStat(std::function<double(double, double)> function, std::string name);

	void add1Dtransform(std::function<double(double, double)> function, std::string name, int initial_depth);
	void add2Dtransform(std::pair<std::function<double(double, double)>, std::function<double(double, double)>> function, std::string name, int initial_depth);
	void add2Dtransform(std::function<double(double, double)> fun_x, std::function<double(double, double)> fun_y, std::string name, int initial_depth);
	void addWeight(std::function<double(double, double)> function, std::string postfix);
	void setMergeFractions(std::vector<std::string> fractions1D, std::vector<std::string> fractions2D);
};

class Node {
public:
	double umin = 0;
	double umax = 0;
	double vmin = 0;
	double vmax = 0;
	double xmin = 0;
	double xmax = 0;
	double ymin = 0;
	double ymax = 0;
	double weight = 0;
	double pdf = 0;
	double prob = 0;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> weights;
	bool operator<(const Node& other) const noexcept;
	bool operator>(const Node& other) const noexcept;
	bool operator==(const Node& other) const noexcept;
	bool inside(const double&, const double&);
	void set_bounds();
	void set_bounds(const double &, const double &, const double &, const double &);
	const Node child1() const;
	const Node child2() const;
	const Node child3() const;
	const Node child4() const;
	const Node child1_log() const;
	const Node child2_log() const;
	const Node child3_log() const;
	const Node child4_log() const;
	const Node child1_loglinear() const;
	const Node child2_loglinear() const;
	const Node child3_loglinear() const;
	const Node child4_loglinear() const;
	const Node child1(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&) const;
	const Node child1(const double&) const;
	const Node child2(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&) const;
	const Node child2(const double&) const;
	const Node child3(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&) const;
	const Node child3(const double&) const;
	const Node child4(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&) const;
	const Node child4(const double&) const;
	const Node childL(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&) const;
	const Node childL(const double&) const;
	const Node childL_linear() const;
	const Node childL_log() const;
	const Node childR(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&) const;
	const Node childR(const double&) const;
	const Node childT(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&) const;
	const Node childT() const;
	const Node childB(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&) const;
	const Node childB() const;
	const Node childB_linear() const;
	void child1_in_place();
	void child2_in_place();
	void child3_in_place();
	void child4_in_place();
	void child4_in_place_log();
	void child1_in_place(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&);
	void child2_in_place(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&);
	void child3_in_place(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&);
	void child4_in_place(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&);
	void child4_in_place(const double&);
	void child4_in_place_loglinear();
	void childL_in_place(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&);
	void childR_in_place(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&);
	void childR_in_place(const double&);
	void childR_in_place_linear();
	void childR_in_place_log();
	void childT_in_place(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&);
	void childT_in_place();
	void childT_in_place_linear();
	void childB_in_place(double (*) (const double&, const double&), double (*) (const double&), double (*) (const double&, const double&), double (*) (const double&), const double&);
	double d2e() const;
	double center_x() const;
	double center_y() const;
	double center_weighted_x(double threshold) const;
	double center_weighted_y() const;
};

class Node1D {
public:
	int depth = 0;
	long double weight = 0;
	//std::multiset<long double> partial_sum;
	int child[2] = { -1,-1 };
	int parent = -1;
};

class Node2D {
public:
	int depth = 0;
	long double weight = 0;
	//std::multiset<long double> partial_sum;
	int child[4] = { -1,-1,-1,-1 };
	int parent = -1;
};

class Node3D {
public:
	int depth = 0;
	long double weight = 0;
	//std::multiset<long double> partial_sum;
	int child[8] = { -1,-1,-1,-1,-1,-1,-1,-1 };
	int parent = -1;
};

class BinaryTree {
private:
	void output(std::ostream& os, int i, int dep, double x, double dx);
	double symlog(double x);
	double symlog_inv(double x);
	double natural_log(double x);
	double exp(double x);
	double dec_log(double x);
	double pow10(double x);
	double lin(double x);
	void subdivide(int parent_index);
	void cleanup();
	void weigh(int node_index, double x_scaled, double current_depth);
	void write_node(std::ofstream& output, int node_index);
	void write_node(std::ofstream& output, int node_index, const double& total_weight);
	void write_leaf(std::ofstream& output, int node_index, double x_scaled, double current_depth, const double& total_weight);
	void traverse_nodes(std::map<double, double>& pdf_weight, int node_index, double x_scaled, double current_depth, const double& total_weight);
public:
	int depth = 0;
	double threshold = 0.3;
	double xmin, xmax;
	std::vector<Node1D> nodes;

	std::function<double(double)> scale;
	std::function<double(double)> scale_inv;

	friend std::ostream& operator<<(std::ostream& os, BinaryTree& b);

	void set_scales(std::string scale);
	void set_bounds(double XMIN, double XMAX);

	void fix_depth();
	void fix_depth(int index, int current_depth);

	void initialize(double XMIN, double XMAX, std::string scale_X);
	void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::string scale_X);
	void uniform_divide(int DEPTH);

	void add_point(double x, double weight);
	void add_point(double x);

	void load_points(const std::vector<double> &pts);
	void load_points(const std::vector<double> &pts, const std::vector<double> &wt);
	void load_points(const std::vector<double> &x, const std::vector<double> &y, std::function<double(double, double)> &fun);
	void load_points(const const std::vector<double> &x, const const std::vector<double> &y, std::function<double(double, double)> &fun, std::function<double(double, double)> &fun_wt);
	void load_points(long double x[], long double y[], const std::function<double(double, double)> &fun, const std::function<double(double, double)> &fun_wt, int size);
	void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x);
	void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_weight);

	double count_total(int i);
	double count_total();

	void merge(const double& merge_weight, int index);

	void merge(double threshold_weight);

	void weigh();

	void save_structure(std::string path, std::string prefix, std::string filename);
	void save_structure(std::string path_out, std::string filename);
	void save_structure(std::string path_out, std::string filename, bool normalize);

	void save_leaves(std::string path, std::string prefix, std::string filename);
	void save_leaves(std::string path_out, std::string filename);

	void prob_to_pdf(std::string path, std::string prefix, std::string filename);
	void prob_to_pdf(std::string path_out, std::string filename);

	void load_tree(std::string path, std::string prefix, std::string filename, std::string scale);
	void load_tree(std::string path, std::string filename);

	void adjust_bounds(double& XMIN, double& XMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x);

	void adjust_bounds(double& XMIN, double& XMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_w);
};

BinaryTree load_tree(std::string path, std::string prefix, std::string filename, std::string scale);

class Quadtree {
private:
	const int children_per_node = 4;
	void output(std::ostream& os, int i, int dep, double x, double y, double dx, double dy);
	void cleanup();
	void weigh(int node_index, double x_scaled, double y_scaled, double current_depth);
	void write_node(std::ofstream& output, int node_index);
	void write_node(std::ofstream& output, int node_index, const double& total_weight);
	void write_leaf(std::ofstream& output, int node_index, double x_scaled, double y_scaled, double current_depth, const double& total_weight);
	void adjust_bounds_x(double& XMIN, double& XMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x);
	void adjust_bounds_x(double& XMIN, double& XMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_w);
	void adjust_bounds_y(double& YMIN, double& YMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_y);
	void adjust_bounds_y(double& YMIN, double& YMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_y, std::vector<std::string> filename_w);
	void traverse_nodes(std::map<double, double>& pdf_weight, int node_index, double x_scaled, double y_scaled, double current_depth, const double& total_weight);
public:
	int depth = 0;
	double threshold = 0.3;
	double xmin, xmax, ymin, ymax;
	void subdivide(int parent_index);

	std::vector<Node2D> nodes;
	//double (Quadtree::* scale_x)(double x);
	std::function<double(double)> scale_x;

	std::function<double(double)> scale_y;

	friend std::ostream& operator<<(std::ostream& os, Quadtree &qt);

	void set_scales(std::string scale_X, std::string scale_Y);
	void set_scales(std::function<double(double)> scale_X, std::function<double(double)> scale_Y);
	void set_bounds(double XMIN, double XMAX, double YMIN, double YMAX);

	void fix_depth();
	void fix_depth(int index, int current_depth);

	void initialize(double XMIN, double XMAX, double YMIN, double YMAX, std::string scale_X, std::string scale_Y);
	void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::string scale_X, std::string scale_Y);
	void uniform_divide(int DEPTH);

	void add_point(double x, double y, double weight);
	void add_point(double x, double y);

	void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y);
	void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_weight);
	void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::function<double(double, double, double)> fun1, std::function<double(double, double, double)> fun2, std::function<bool(double, double, double)> cond);
	void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::function<double(double, double, double)> fun1, std::function<double(double, double, double)> fun2, std::function<bool(double, double, double, double, double)> cond, double cond_param1, double cond_param2);
	void load_points(long double x[], long double y[], const std::function<double(double, double)> &fun_x, const std::function<double(double, double)> &fun_y, const std::function<double(double, double)> &fun_wt, int size);

	double count_total(int i);
	double count_total();

	void merge(const double& merge_weight, int index);
	void merge(const double& merge_weight, int index, double min_depth);

	void merge(double threshold_weight);
	void merge(double threshold_weight, double min_depth);

	void weigh();

	void save_structure(std::string path_out, std::string filename);
	void save_structure(std::string path_out, std::string filename, bool normalize);

	void save_leaves(std::string path, std::string prefix, std::string filename);
	void save_leaves(std::string path_out, std::string filename);

	void prob_to_pdf(std::string path, std::string prefix, std::string filename);
	void prob_to_pdf(std::string path_out, std::string filename);

	void adjust_bounds(double& XMIN, double& XMAX, double& YMIN, double& YMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y);

	void adjust_bounds(double& XMIN, double& XMAX, double& YMIN, double& YMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_w);

	std::pair<double, double> center(int node_index);

	void center_and_dxdy(int node_index, double& center_x, double& center_y, double& dxdy);

	void load_tree(std::string path, std::string prefix, std::string filename, std::string scale_X, std::string scale_Y);
	void load_tree(std::string path, std::string filename);

	double collect_weight_X(double XMIN, double XMAX, std::vector<std::pair<std::pair<double, double>, double>>& pdf_X);
	double collect_weight_Y(double YMIN, double YMAX, std::vector<std::pair<std::pair<double, double>, double>>& pdf_Y);
	void fill_marginalized(int node_index, double current_depth, double x_scaled, double y_scaled, std::vector<std::pair<std::pair<double, double>, double>>& pdf_rho, std::vector<std::pair<std::pair<double, double>, double>>& pdf_v);
	void fill_marginalized(std::vector<std::pair<std::pair<double, double>, double>>& pdf_rho, std::vector<std::pair<std::pair<double, double>, double>>& pdf_v, std::string scale_X, std::string scale_Y);
};

void LoadNodes(std::vector<Node>& nodes, std::string filename);

void Subdivide(std::vector<Node>&, Node&, const int&);
void Subdivide(std::vector<Node>&, Node&, const int&, double(&) (const double&, const double&), double(&) (const double&), double(&) (const double&, const double&), double(&) (const double&), const double&, const double&);
void Subdivide_energy_save(const std::string&, Node&, const int&, double(&) (const double&, const double&), double(&) (const double&), double(&) (const double&, const double&), double(&) (const double&), const double&, const double&);

void Subdivide(std::vector<Node>&, Node&, const int&, const double&, const double&);
void Subdivide_energy_save(std::ofstream&, Node&, const int&, const int &, const double&, const double&, const double&, const double&);
void Subdivide_energy_nodes(std::vector<Node>& nodes, Node& parent, const int& N, const int& threshold_N, const double& threshold, const double& aspect_ratio, const double& dx, const double& dy);
void Subdivide_linear_save(std::ofstream&, Node&, const int&, const double&);
void Subdivide_loglog_save(std::ofstream&, Node&, const int&, const double&, const double&, const double&);
void Subdivide_log_linear(std::vector<Node>& nodes, Node& parent, const int& N, const int& threshold_N, const double& aspect_ratio, const double& dx, const double& dy);

void DataExport(const std::vector<Node>&, const std::string&);
void DataExport(const std::map<int, Node>&, const std::string&);
void DataExport(const std::multiset<Node>&, const std::string&);
void DataExport(const std::vector<Node>&, std::vector<double>&, const std::string&);

void DataExport(const std::vector<std::pair<double, double>>&, const std::string&);

int find_index(double x, std::pair<std::pair<double, double>, double>& pdf);

void write_leaf(std::ofstream& output, Quadtree& q1, Quadtree& q2, int node_index, double x_scaled, double y_scaled, double current_depth, const double& total_weight);

void save_difference(std::string path, std::string prefix, std::string filename, Quadtree& q1, Quadtree& q2);

void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, Quadtree& q1, BinaryTree& X, BinaryTree& Y);
void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_weight, Quadtree& Q, BinaryTree& X, BinaryTree& Y);
void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_weight1, std::vector<std::string> filename_weight2, Quadtree& Q1, Quadtree& Q2, BinaryTree& X1, BinaryTree& X2, BinaryTree& Y1, BinaryTree& Y2);
void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_weight1, std::vector<std::string> filename_weight2, std::vector<std::string> filename_weight3, Quadtree& Q1, Quadtree& Q2, Quadtree& Q3, BinaryTree& X1, BinaryTree& X2, BinaryTree& X3, BinaryTree& Y1, BinaryTree& Y2, BinaryTree& Y3);
void load_points(std::vector<Quadtree>& Q, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::function<double(double, double, double)> fun1, std::function<double(double, double, double)> fun2, std::function<bool(double, double, double, double, double)> cond, std::vector<std::pair<double, double>> x);

//void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector< std::function<double(double, double)>> weight_fun, std::vector<std::function<double(double, double)>> field_1D, std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> field_2D, std::vector<BinaryTree>& data_1D, std::vector<Quadtree>& data_2D);
//void load_points(std::string path, std::vector<int> frames, std::string fieldname_x, std::string fieldname_y, std::vector< std::function<double(double, double)>> weight_fun, std::vector<std::function<double(double, double)>> field_1D, std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> field_2D, std::vector<BinaryTree>& data_1D, std::vector<Quadtree>& data_2D);

void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::string scale_X, std::string scale_Y, Quadtree& Q, BinaryTree& X, BinaryTree& Y);
void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::string scale_X, std::string scale_Y, Quadtree& Q1, Quadtree& Q2, BinaryTree& X1, BinaryTree& X2, BinaryTree& Y1, BinaryTree& Y2);
void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::string scale_X, std::string scale_Y, Quadtree& Q1, Quadtree& Q2, Quadtree& Q3, BinaryTree& X1, BinaryTree& X2, BinaryTree& X3, BinaryTree& Y1, BinaryTree& Y2, BinaryTree& Y3);
void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::function<double(double, double)>> field_1D, std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> field_2D, std::vector<std::string> scale_x, std::vector<std::pair<std::string, std::string>> scale_z, std::vector<BinaryTree> &data_1D, std::vector<Quadtree> &data_2D);
//void initialize(std::string path_to_sim, std::vector<int> frames, std::string fieldname_x, std::string fieldname_y, std::vector<std::function<double(double, double)>> field_1D, std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> field_2D, std::vector<std::string> scale_1D, std::vector<std::pair<std::string, std::string>> scale_2D, std::vector<BinaryTree>& data_1D, std::vector<Quadtree>& data_2D);

class Octree {
private:
	const int children_per_node = 8;
	void output(std::ostream& os, int i, int dep, double x, double y, double z, double dx, double dy, double dz);
	void subdivide(int parent_index);
	void cleanup();
	void weigh(int node_index, double x_scaled, double y_scaled, double z_scaled, double current_depth);
	void write_node(std::ofstream& output, int node_index);
	void write_node(std::ofstream& output, int node_index, const double& total_weight);
	void write_leaf(std::ofstream& output, int node_index, double x_scaled, double y_scaled, double z_scaled, double current_depth, const double& total_weight);
	//void adjust_bounds_x(double& XMIN, double& XMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x);
	//void adjust_bounds_y(double& YMIN, double& YMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_y);
	//void adjust_bounds_z(double& ZMIN, double& ZMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_z);
	void traverse_nodes(std::map<double, double>& pdf_weight, int node_index, double x_scaled, double y_scaled, double z_scaled, double current_depth, const double& total_weight);

	void slice_check_x(int N3D_index, Quadtree& q, std::map<int, int>& N3D_to_n2D, double x0_scaled, double x_scaled, double y_scaled, double z_scaled, double depth);

public:
	int depth = 0;
	double threshold = 0.3;
	double xmin, xmax, ymin, ymax, zmin, zmax;

	std::vector<Node3D> nodes;
	//double (Quadtree::* scale_x)(double x);
	std::function<double(double)> scale_x;

	std::function<double(double)> scale_y;

	std::function<double(double)> scale_z;

	friend std::ostream& operator<<(std::ostream& os, Octree& ot);

	void set_scales(std::string scale_X, std::string scale_Y, std::string scale_Z);
	void set_bounds(double XMIN, double XMAX, double YMIN, double YMAX, double ZMIN, double ZMAX);

	void fix_depth();
	void fix_depth(int index, int current_depth);

	void initialize(double XMIN, double XMAX, double YMIN, double YMAX, double ZMIN, double ZMAX, std::string scale_X, std::string scale_Y, std::string scale_Z);
	void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::string scale_X, std::string scale_Y, std::string scale_Z);
	void uniform_divide(int DEPTH);

	void add_point(double x, double y, double z, double weight);
	void add_point(double x, double y, double z);

	void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z);
	void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::vector<std::string> filename_weight);
	void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::function<double(double, double, double)> fun1, std::function<double(double, double, double)> fun2, std::function<double(double, double, double)> fun3, std::function<bool(double, double, double)> cond);

	double count_total(int i);
	double count_total();

	void merge(const double& merge_weight, int index);
	void merge(const double& merge_weight, int index, double min_depth);

	void merge(double threshold_weight);
	void merge(double threshold_weight, double min_depth);

	void weigh();

	void save_structure(std::string path, std::string prefix, std::string filename);
	void save_structure(std::string path_out, std::string filename);
	void save_structure(std::string path_out, std::string filename, bool normalize);

	void save_leaves(std::string path, std::string prefix, std::string filename);
	void save_leaves(std::string path_out, std::string filename);

	void prob_to_pdf(std::string path, std::string prefix, std::string filename);
	void prob_to_pdf(std::string path_out, std::string filename);

	//void adjust_bounds(double& XMIN, double& XMAX, double& YMIN, double& YMAX, double& ZMIN, double& ZMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z);

	//void adjust_bounds(double& XMIN, double& XMAX, double& YMIN, double& YMAX, double& ZMIN, double& ZMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::vector<std::string> filename_w);

	//std::pair<double, double> center(int node_index);

	//void center_and_dxdy(int node_index, double& center_x, double& center_y, double& dxdy);

	void load_tree(std::string path, std::string prefix, std::string filename, std::string scale_X, std::string scale_Y, std::string scale_Z);

	Quadtree slice(int coordinate, double value);

	//double collect_weight_X(double XMIN, double XMAX, std::vector<std::pair<std::pair<double, double>, double>>& pdf_X);
	//double collect_weight_Y(double YMIN, double YMAX, std::vector<std::pair<std::pair<double, double>, double>>& pdf_Y);
	//double collect_weight_Y(double ZMIN, double ZMAX, std::vector<std::pair<std::pair<double, double>, double>>& pdf_Z);
	//void fill_marginalized(int node_index, double current_depth, double x_scaled, double y_scaled, double z_scaled, std::vector<std::pair<std::pair<double, double>, double>>& pdf_X, std::vector<std::pair<std::pair<double, double>, double>>& pdf_Y, std::vector<std::pair<std::pair<double, double>, double>>& pdf_Z);
	//void fill_marginalized(std::vector<std::pair<std::pair<double, double>, double>>& pdf_X, std::vector<std::pair<std::pair<double, double>, double>>& pdf_Y, std::vector<std::pair<std::pair<double, double>, double>>& pdf_Z, std::string scale_X, std::string scale_Y, std::string scale_Z);
};

#endif
