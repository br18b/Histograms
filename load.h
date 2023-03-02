#pragma once
#ifndef LOADER
#define LOADER

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include "node.h"
#include "read_binary.h"

void load(std::string, std::vector<double>&);
void load(std::string, std::vector<double>&, const int &);
void load(std::string filename, std::vector<std::pair<std::pair<double, double>, double>>& res);
void load(std::string filename, std::vector<std::pair<double, std::pair<double, double>>>& res);
std::vector<std::pair<std::pair<double, double>, double>> load1Dhistogram(std::string filename);
std::vector<std::pair<double, std::pair<double, double>>> loadSigmas(std::string filename);
std::vector<std::pair<std::pair<double, double>, double>> load1Dhistogram(std::string path, std::string prefix, std::string filename);
std::vector<std::pair<double, std::pair<double, double>>> loadSigmas(std::string path, std::string prefix, std::string filename);

std::vector<std::pair<std::pair<double, double>, double>> loadCombineTwo(const std::vector<std::pair<std::pair<double, double>, double>>& data1, const std::vector<std::pair<std::pair<double, double>, double>>& data2, std::function<double(double, double)> function);

void load(std::string filename, std::vector<std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>>& res);
void load(std::string filename, std::vector<int>& res);
std::vector<int> load(std::string filename);
std::vector<std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> load2Dhistogram(std::string filename);
std::vector<std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> load2Dhistogram(std::string path, std::string prefix, std::string filename);

void load_linear(const std::string&, const std::string&, const std::string&, std::vector<Node>&);
void load_log_linear(const std::string&, const std::string&, const std::string&, std::vector<Node>&);
void load_loglog(const std::string&, const std::string&, const std::string&, std::vector<Node>&);
void load(const std::string&, const std::string&, const std::string&, std::vector<Node>&);
void load_binary(const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, std::vector<Node>&);

void load_Pearson(const std::string&, std::vector<double>& res);
void load_Pearson(const std::string&, const std::string&, const std::string&, std::vector<double>& res);
std::vector<double> load_Pearson(const std::string&, const std::string&, const std::string&);

void load_twopars(const std::string&, std::vector<double>& res);
void load_twopars(const std::string&, const std::string&, const std::string&, std::vector<double>& res);
std::vector<double> load_twopars(const std::string&, const std::string&, const std::string&);

void load_single_par(const std::string&, double& res);
void load_single_par(const std::string&, const std::string&, const std::string&, double& res);
double load_single_par(const std::string&, const std::string&, const std::string&);

void load_pars(const std::string&, std::vector<double>& res);
void load_pars(const std::string&, std::vector<double>& res, int start, int end);
void load_pars(const std::string&, const std::string&, const std::string&, std::vector<double>& res);
void load_pars(const std::string&, const std::string&, const std::string&, std::vector<double>& res, int start, int end);
std::vector<double> load_pars(const std::string&, const std::string&, const std::string&);

#endif
