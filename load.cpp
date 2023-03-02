#include "load.h"

void load(std::string filename, std::vector<double>& res) {
	std::ifstream in;
	in.open(filename);
	std::string line;
	res.resize(0);
	while (std::getline(in, line)) {
		res.push_back(std::stold(line));
	}
	std::cout << filename << " loaded! Elements: " << res.size() << std::endl;
}

void load(std::string filename, std::vector<double>& res, const int& N) {
	std::ifstream in;
	in.open(filename);
	std::string line;
	res.resize(0);
	int n = 0;
	while (std::getline(in, line) && n < N) {
		res.push_back(std::stold(line));
		n++;
	}
	std::cout << filename << " loaded! Elements: " << res.size() << std::endl;
}

void load(std::string filename, std::vector<std::pair<std::pair<double, double>, double>>& res) {
    std::ifstream input(filename);
    if (input.is_open()) {
        std::string line;
        res.resize(0);
        double a, b, c, d;
        while (std::getline(input, line) && (line != "")) {
            std::istringstream ss(line);
            if (ss >> a >> b >> c >> d) res.push_back(std::pair<std::pair<double, double>, double>{ {a, b}, d});
        }
        std::cout << filename << " loaded! Elements: " << res.size() << std::endl;
        input.close();
    }
    else std::cout << filename << " not found!" << std::endl;
}

void load(std::string filename, std::vector<std::pair<double, std::pair<double, double>>>& res) {
    std::ifstream input(filename);
    if (input.is_open()) {
        std::string line;
        res.resize(0);
        double a, b, c, d;
        while (std::getline(input, line) && (line != "")) {
            std::istringstream ss(line);
            if (ss >> a >> b >> c) res.push_back(std::pair<double, std::pair<double, double>> { a, { b, c }});
        }
        std::cout << filename << " loaded! Elements: " << res.size() << std::endl;
        input.close();
    }
    else std::cout << filename << " not found!" << std::endl;
}

void load(std::string filename, std::vector<int>& res) {
    std::ifstream input(filename);
    if (input.is_open()) {
        std::string line;
        res.resize(0);
        int a;
        while (std::getline(input, line) && (line != "")) {
            std::istringstream ss(line);
            if (ss >> a) res.push_back(a);
        }
        std::cout << filename << " loaded! Elements: " << res.size() << std::endl;
        input.close();
    }
    else std::cout << filename << " not found!" << std::endl;
}

std::vector<std::pair<std::pair<double, double>, double>> load1Dhistogram(std::string filename) {
    std::vector<std::pair<std::pair<double, double>, double>> res;
    load(filename, res);
    return res;
}

std::vector<std::pair<double, std::pair<double, double>>> loadSigmas(std::string filename) {
    std::vector<std::pair<double, std::pair<double, double>>> res;
    load(filename, res);
    return res;
}

std::vector<std::pair<std::pair<double, double>, double>> load1Dhistogram(std::string path, std::string prefix, std::string filename) {
    return load1Dhistogram(path + prefix + "_" + filename);
}

std::vector<std::pair<double, std::pair<double, double>>> loadSigmas(std::string path, std::string prefix, std::string filename) {
    return loadSigmas(path + prefix + "_" + filename);
}

std::vector<int> load(std::string filename) {
    std::vector<int> res; load(filename, res);
    return res;
}

std::vector<std::pair<std::pair<double, double>, double>> loadCombineTwo(const std::vector<std::pair<std::pair<double, double>, double>>& data1, const std::vector<std::pair<std::pair<double, double>, double>>& data2, std::function<double(double, double)> function) {
    std::vector<std::pair<std::pair<double, double>, double>> res;
    Interpolation fun1, fun2;
    for (const auto& p : data1) fun1.pairs.insert(std::pair<double, double>{0.5 * (p.first.first + p.first.second), p.second});
    for (const auto& p : data2) fun2.pairs.insert(std::pair<double, double>{0.5 * (p.first.first + p.first.second), p.second});
    for (const auto& p : data1) res.push_back(std::pair<std::pair<double, double>, double>{ {p.first.first, p.first.second}, function(fun1.get_value(0.5 * (p.first.first + p.first.second)), fun2.get_value(0.5 * (p.first.first + p.first.second))) });
    return res;
}

void load(std::string filename, std::vector<std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>>& res) {
    std::ifstream input(filename);
    if (input.is_open()) {
        std::string line;
        res.resize(0);
        double a, b, c, d, e, f;
        while (std::getline(input, line) && (line != "")) {
            std::istringstream ss(line);
            if (ss >> a >> b >> c >> d >> e >> f) res.push_back(std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>{ { {a, c}, { b, d }}, f});
        }
        std::cout << filename << " loaded! Elements: " << res.size() << std::endl;
    }
    else std::cout << filename << " not found!" << std::endl;
}

std::vector<std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> load2Dhistogram(std::string filename) {
    std::vector<std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> res;
    load(filename, res);
    return res;
}

std::vector<std::pair<std::pair<std::pair<double, double>, std::pair<double, double>>, double>> load2Dhistogram(std::string path, std::string prefix, std::string filename) {
    return load2Dhistogram(path + prefix + "_" + filename);
}

void load_linear(const std::string& filename_x, const std::string& filename_y, const std::string& filenameW, std::vector<Node>& nodes) {
    std::ifstream inEK, inET, inW;
    inET.open(filename_x);
    inEK.open(filename_y);
    inW.open(filenameW);
    std::string lineEK, lineET, lineW;
    nodes.resize(0);
    Node node;
    node.xmin = 100;
    node.xmax = -100;
    node.ymin = 100;
    node.ymax = -100;
    node.x.resize(0); node.y.resize(0); node.weights.resize(0);
    double x, y, w;
    std::cout << "Loading " << filename_x << ", " << filename_y << ", " << filenameW << "... (this might take a while) ...";
    while (std::getline(inEK, lineEK)) {
        std::getline(inET, lineET);
        std::getline(inW, lineW);
        x = std::stold(lineET);
        y = std::stold(lineEK);
        if (inW.is_open()) w = std::stod(lineW);
        else w = 1;
        node.x.push_back(x);
        node.y.push_back(y);
        node.weights.push_back(w);
        if (node.xmin > x) node.xmin = x;
        if (node.xmax < x) node.xmax = x;
        if (node.ymin > y) node.ymin = y;
        if (node.ymax < y) node.ymax = y;
        //if ((n % frac) == 0) std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << floor(10000. * n / N) / (100.0) << "%" << "   ";
    }
    std::cout << std::endl << "Data loaded! Elements: " << node.x.size() << std::endl;
    nodes.push_back(node);
}

void load_log_linear(const std::string& filename_x, const std::string& filename_y, const std::string& filenameW, std::vector<Node>& nodes) {
    std::ifstream inEK, inET, inW;
    inEK.open(filename_x);
    inET.open(filename_y);
    inW.open(filenameW);
    std::string lineEK, lineET, lineW;
    nodes.resize(0);
    Node node;
    node.xmin = 100;
    node.xmax = 0;
    node.ymin = 100;
    node.ymax = -100;
    node.x.resize(0); node.y.resize(0); node.weights.resize(0);
    double x, y, w;
    std::cout << "Loading " << filename_x << ", " << filename_y << ", " << filenameW << "... (this might take a while) ...";
    while (std::getline(inEK, lineEK)) {
        std::getline(inET, lineET);
        std::getline(inW, lineW);
        x = std::stold(lineET);
        y = std::stold(lineEK);
        if (inW.is_open()) w = std::stod(lineW);
        else w = 1;
        node.x.push_back(x);
        node.y.push_back(y);
        node.weights.push_back(w);
        if (node.xmin > x) node.xmin = x;
        if (node.xmax < x) node.xmax = x;
        if (node.ymin > y) node.ymin = y;
        if (node.ymax < y) node.ymax = y;
        //if ((n % frac) == 0) std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << floor(10000. * n / N) / (100.0) << "%" << "   ";
    }
    std::cout << std::endl << "Data loaded! Elements: " << node.x.size() << std::endl;
    nodes.push_back(node);
}

void load_loglog(const std::string& filename_x, const std::string& filename_y, const std::string& filenameW, std::vector<Node>& nodes) {
    std::ifstream inX, inY, inW;
    inX.open(filename_x);
    inY.open(filename_y);
    inW.open(filenameW);
    std::string lineX, lineY, lineW;
    nodes.resize(0);
    Node node;
    node.xmin = 100;
    node.xmax = 0;
    node.ymin = 100;
    node.ymax = 0;
    node.x.resize(0); node.y.resize(0); node.weights.resize(0);
    double x, y, w;
    std::cout << "Loading " << filename_x << ", " << filename_y << ", " << filenameW << "... (this might take a while) ...";
    while (std::getline(inX, lineX)) {
        std::getline(inY, lineY);
        std::getline(inW, lineW);
        x = std::stold(lineX);
        y = std::stold(lineY);
        if (inW.is_open()) w = std::stod(lineW);
        else w = 1;
        if ((x > 0) && (y > 0)) {
            node.x.push_back(x);
            node.y.push_back(y);
            node.weights.push_back(w);
            if (node.xmin > x) node.xmin = x;
            if (node.xmax < x) node.xmax = x;
            if (node.ymin > y) node.ymin = y;
            if (node.ymax < y) node.ymax = y;
        }
        //if ((n % frac) == 0) std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << floor(10000. * n / N) / (100.0) << "%" << "   ";
    }
    std::cout << std::endl << "Data loaded! Elements: " << node.x.size() << std::endl;
    nodes.push_back(node);
}

void load(const std::string& filenameX, const std::string& filenameY, const std::string& filenameW, std::vector<Node>& nodes) {
    std::ifstream inX, inY, inW;
    inX.open(filenameX);
    inY.open(filenameY);
    inW.open(filenameW);
    std::string lineX, lineY, lineW;
    nodes.resize(0);
    Node below, above;
    below.xmin = 0;
    below.xmax = 0;
    below.ymin = 1;
    below.ymax = 0;
    above.xmin = 0;
    above.xmax = 0;
    above.ymin = 1;
    above.ymax = 0;
    below.x.resize(0); below.y.resize(0); below.weights.resize(0);
    //below.x.reserve(N); below.y.reserve(N); below.weights.reserve(N);
    above.x.resize(0); above.y.resize(0); above.weights.resize(0);
    //above.x.reserve(N); above.y.reserve(N); above.weights.reserve(N);
    double x, y, w;
    std::cout << "Loading ... (this might take a while) ..." << std::endl;
    while (std::getline(inX, lineX)) {
        std::getline(inY, lineY);
        std::getline(inW, lineW);
        x = std::stod(lineX);
        y = std::stod(lineY);
        if (inW.is_open()) w = std::stod(lineW);
        else w = 1;
        if (y > 0) {
            if (x < 0) {
                below.x.push_back(x);
                below.y.push_back(y);
                below.weights.push_back(w);
            }
            else {
                above.x.push_back(x);
                above.y.push_back(y);
                above.weights.push_back(w);
            }
        }
    }
    below.x.shrink_to_fit(); below.y.shrink_to_fit(); below.weights.shrink_to_fit();
    above.x.shrink_to_fit(); above.y.shrink_to_fit(); above.weights.shrink_to_fit();
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bData loaded! Elements: " << below.x.size() + above.x.size() << std::endl;
    nodes.push_back(below);
    nodes.push_back(above);
}

void load_binary(const std::string& path, const std::string& prefix, const std::string& filenameX, const std::string& filenameY, const std::string& filenameW, std::vector<Node>& nodes) {

    std::vector<double> X, Y, W;
    double x, y, w;
    std::cout << "Loading ... (this might take a while) ..." << std::endl;
    ::read_from_binary(path, prefix, filenameX, X);
    ::read_from_binary(path, prefix, filenameY, Y);
    int W_read = ::read_from_binary(path, prefix, filenameW, W);
    if (W_read == 0) std::cout << " Loaded!" << std::endl;
    else std::cout << " Loaded! Using default weights..." << std::endl;
    nodes.resize(0);
    if (filenameX == "Helmholtz") {
        Node below, above;
        below.xmin = 0;
        below.xmax = 0;
        below.ymin = 1;
        below.ymax = 0;
        above.xmin = 0;
        above.xmax = 0;
        above.ymin = 1;
        above.ymax = 0;
        below.x.resize(0); below.y.resize(0); below.weights.resize(0);
        //below.x.reserve(N); below.y.reserve(N); below.weights.reserve(N);
        above.x.resize(0); above.y.resize(0); above.weights.resize(0);
        for (int i = 0; i < X.size(); i++) {
            x = X[i];
            y = Y[i];
            if (W_read == 0) w = W[i];
            else w = 1;
            if (y > 0) {
                if (x < 0) {
                    below.x.push_back(x);
                    below.y.push_back(y);
                    below.weights.push_back(w);
                }
                else {
                    above.x.push_back(x);
                    above.y.push_back(y);
                    above.weights.push_back(w);
                }
            }
        }
        below.x.shrink_to_fit(); below.y.shrink_to_fit(); below.weights.shrink_to_fit();
        above.x.shrink_to_fit(); above.y.shrink_to_fit(); above.weights.shrink_to_fit();
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bData loaded! Elements: " << below.x.size() + above.x.size() << std::endl;
        nodes.push_back(below);
        nodes.push_back(above);
    }
    else {
        Node main_node;
        main_node.x.resize(0); main_node.y.resize(0); main_node.weights.resize(0);
        for (int i = 0; i < X.size(); i++) {
            x = X[i];
            y = Y[i];
            if (W_read == 0) w = W[i];
            else w = 1;
            if (y > 0) {
                main_node.x.push_back(x);
                main_node.y.push_back(y);
                main_node.weights.push_back(w);
            }
        }
        main_node.x.shrink_to_fit(); main_node.y.shrink_to_fit(); main_node.weights.shrink_to_fit();
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bData loaded! Elements: " << main_node.x.size() << std::endl;
        nodes.push_back(main_node);
    }
}

void load_Pearson(const std::string& filename, std::vector<double>& res) {
    std::ifstream input;
    std::string line;
    input.open(filename);
    if (input.is_open()) {
        std::string line;
        res.resize(0);
        std::string a, b;
        while (std::getline(input, line) && (line != "")) {
            std::istringstream ss(line);
            if (ss >> a >> b) res.push_back(std::stod(b));
        }
        std::cout << filename << " loaded!" << std::endl;
        input.close();
    }
    else std::cout << filename << " not found!" << std::endl;
}

void load_Pearson(const std::string& path, const std::string& prefix, const std::string& filename, std::vector<double>& res) {
    load_Pearson(path + prefix + "_" + filename, res);
}

std::vector<double> load_Pearson(const std::string& path, const std::string& prefix, const std::string& filename) {
    std::vector<double> res;
    load_Pearson(path, prefix, filename, res);
    return res;
}

void load_single_par(const std::string& filename, double& res) {
    std::ifstream input;
    std::string line;
    input.open(filename);
    if (input.is_open()) {
        std::string line;
        std::string a;
        while (std::getline(input, line) && (line != "")) {
            std::istringstream ss(line);
            if (ss >> a) res = std::stod(a);
        }
        std::cout << filename << " loaded!" << std::endl;
        input.close();
    }
    else std::cout << filename << " not found!" << std::endl;
}

void load_single_par(const std::string& path, const std::string& prefix, const std::string& filename, double& res) {
    load_single_par(path + prefix + "_" + filename, res);
}

double load_single_par(const std::string& path, const std::string& prefix, const std::string& filename) {
    double res = 0;
    load_single_par(path, prefix, filename, res);
    return res;
}

void load_twopars(const std::string& filename, std::vector<double>& res) {
    std::ifstream input;
    std::string line;
    input.open(filename);
    if (input.is_open()) {
        std::string line;
        res.resize(0);
        std::string a, b;
        while (std::getline(input, line) && (line != "")) {
            std::istringstream ss(line);
            if (ss >> a >> b) {
                res.push_back(std::stod(a));
                res.push_back(std::stod(b));
            }
        }
        std::cout << filename << " loaded!" << std::endl;
        input.close();
    }
    else std::cout << filename << " not found!" << std::endl;
}

void load_twopars(const std::string& path, const std::string& prefix, const std::string& filename, std::vector<double>& res) {
    load_twopars(path + prefix + "_" + filename, res);
}

std::vector<double> load_twopars(const std::string& path, const std::string& prefix, const std::string& filename) {
    std::vector<double> res;
    load_twopars(path, prefix, filename, res);
    return res;
}

void load_pars(const std::string& filename, std::vector<double>& res) {
    std::ifstream input;
    std::string line;
    input.open(filename);
    int i = 0;
    if (input.is_open()) {
        std::string line;
        //res.resize(0);
        std::string a;
        while (std::getline(input, line) && (line != "")) {
            std::istringstream ss(line);
            while (ss >> a) {
                if (res.size() > i) res[i] = std::stod(a);
                else res.push_back(std::stod(a));
                i++;
            }
        }
        std::cout << filename << " loaded!" << std::endl;
        input.close();
    }
    else std::cout << filename << " not found!" << std::endl;
}

void load_pars(const std::string& filename, std::vector<double>& res, int start, int end) {
    std::ifstream input;
    std::string line;
    input.open(filename);
    int i = 0;
    int j = 0;
    if (input.is_open()) {
        std::string line;
        //res.resize(0);
        std::string a;
        while (std::getline(input, line) && (line != "")) {
            std::istringstream ss(line);
            while (ss >> a) {
                if ((j >= start) && (j < end)) {
                    if (res.size() > i) res[i] = std::stod(a);
                    else res.push_back(std::stod(a));
                    i++;
                }
                j++;
            }
        }
        std::cout << filename << " loaded!" << std::endl;
        input.close();
    }
    else std::cout << filename << " not found!" << std::endl;
}

void load_pars(const std::string& path, const std::string& prefix, const std::string& filename, std::vector<double>& res) {
    load_pars(path + prefix + "_" + filename, res);
}

void load_pars(const std::string& path, const std::string& prefix, const std::string& filename, std::vector<double>& res, int start, int end) {
    load_pars(path + prefix + "_" + filename, res, start, end);
}

std::vector<double> load_pars(const std::string& path, const std::string& prefix, const std::string& filename) {
    std::vector<double> res;
    load_pars(path, prefix, filename, res);
    return res;
}