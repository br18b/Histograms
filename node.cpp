#include "node.h"

double symlog(double x, double threshold) {
    if (x > 0) return std::log10(x / threshold + 1);
    else return -std::log10(-x / threshold + 1);
}

double symlogInv(double x, double threshold) {
    if (x > 0) return threshold * (pow(10, x) - 1);
    else return -threshold * (pow(10, -x) - 1);
}

double natural_log(double x) {
    return std::log(x);
}

double dec_log(double x) {
    return std::log10(x);
}

double lin(double x) {
    return x;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void Parameters::initialize() {
    frames.resize(0);
	transform_1D.resize(0); name_transform_1D.resize(0);
	transform_2D.resize(0); name_transform_2D.resize(0);
	weight.resize(0); postfix_weight.resize(0);
	merge_fraction_1D.resize(0);
	merge_fraction_2D.resize(0);
    initial_depth_1D.resize(0);
    initial_depth_2D.resize(0);
}

void Parameters::setFieldnames(std::string fname_x, std::string fname_y) {
    fieldname_x = fname_x;
    fieldname_y = fname_y;
}

void Parameters::addFrame(int f) {
    frames.push_back(f);
}

void Parameters::setFrames(int start, int end) {
    frames.resize(0);
    if (start <= end) {
        for (int f = start; f <= end; f++) {
            frames.push_back(f);
        }
    }
}

void Parameters::add1Dtransform(std::function<double(double, double)> function, std::string name, int initial_depth) {
    transform_1D.push_back(function);
    name_transform_1D.push_back(name);
    initial_depth_1D.push_back(initial_depth);
}

void Parameters::add2Dtransform(std::pair<std::function<double(double, double)>, std::function<double(double, double)>> function, std::string name, int initial_depth) {
    transform_2D.push_back(function);
    name_transform_2D.push_back(name);
    initial_depth_2D.push_back(initial_depth);
}

void Parameters::add2Dtransform(std::function<double(double, double)> fun_x, std::function<double(double, double)> fun_y, std::string name, int initial_depth) {
    transform_2D.push_back(std::make_pair(fun_x, fun_y));
    name_transform_2D.push_back(name);
    initial_depth_2D.push_back(initial_depth);
}

void Parameters::addWeight(std::function<double(double, double)> function, std::string postfix) {
    weight.push_back(function);
    postfix_weight.push_back(postfix);
}

void Parameters::setMergeFractions(std::vector<std::string> fractions1D, std::vector<std::string> fractions2D) {
    merge_fraction_1D.resize(0);
    merge_fraction_2D.resize(0);
    for (auto &f: fractions1D) merge_fraction_1D.push_back(f);
    for (auto &f: fractions2D) merge_fraction_2D.push_back(f);
}

bool Node::inside(const double& x, const double& y) {
    return (x >= xmin) && (x < xmax) && (y >= ymin) && (y < ymax);
}

bool Node::operator<(const Node& other) const noexcept {
    return pdf < other.pdf;
}

bool Node::operator>(const Node& other) const noexcept {
    return pdf > other.pdf;
}

bool Node::operator==(const Node& other) const noexcept {
    return pdf == other.pdf;
}

void LoadNodes(std::vector<Node>& nodes, std::string filename) {
    std::ifstream input(filename);
    double xmin, xmax, ymin, ymax, weight;
    nodes.resize(0);
    while (input >> xmin >> xmax >> ymin >> ymax >> weight) {
        Node n;
        n.xmin = xmin;
        n.xmax = xmax;
        n.ymin = ymin;
        n.ymax = ymax;
        n.weight = weight;
        n.pdf = weight / n.d2e();
        nodes.push_back(n);
    }
}

void Subdivide(std::vector<Node>& nodes, Node& parent, const int &N) {
    double total_weight = 0;
    for (int i = 0; i < parent.weights.size(); i++) {
        total_weight += parent.weights[i];
    }
    //std::cout << "(" << parent.xmin << ", " << parent.xmax << "), (" << parent.ymin << ", " << parent.ymax << "), tot. weight: " << total_weight << std::endl;
    if (total_weight > 4000) {
        Node child1 = parent.child1();
        Node child2 = parent.child2();
        Node child3 = parent.child3();
        Node child4 = parent.child4();

        child1.x.resize(0); child1.y.resize(0); child1.weights.resize(0);
        child2.x.resize(0); child2.y.resize(0); child2.weights.resize(0);
        child3.x.resize(0); child3.y.resize(0); child3.weights.resize(0);
        child4.x.resize(0); child4.y.resize(0); child4.weights.resize(0);
        while (!parent.x.empty()) {
            if (child1.inside(parent.x.back(), parent.y.back())) {
                child1.x.push_back(parent.x.back());
                child1.y.push_back(parent.y.back());
                child1.weights.push_back(parent.weights.back());
            }
            else if (child2.inside(parent.x.back(), parent.y.back())) {
                child2.x.push_back(parent.x.back());
                child2.y.push_back(parent.y.back());
                child2.weights.push_back(parent.weights.back());
            }
            else if (child3.inside(parent.x.back(), parent.y.back())) {
                child3.x.push_back(parent.x.back());
                child3.y.push_back(parent.y.back());
                child3.weights.push_back(parent.weights.back());
            }
            else {
                child4.x.push_back(parent.x.back());
                child4.y.push_back(parent.y.back());
                child4.weights.push_back(parent.weights.back());
            }
            parent.x.pop_back();
            parent.y.pop_back();
            parent.weights.pop_back();
        }

        Subdivide(nodes, child1, N);
        Subdivide(nodes, child2, N);
        Subdivide(nodes, child3, N);
        Subdivide(nodes, child4, N);
    }
    else {
        parent.weight = total_weight / N;
        nodes.push_back(parent);
    }
}

void Subdivide(std::vector<Node>& nodes, Node& parent, const int& N, double(&fx) (const double&, const double&), double(&fy) (const double&), double(&fxInv) (const double&, const double&), double(&fyInv) (const double&), const double& threshold, const double& aspect_ratio) {
    double total_weight = std::accumulate(parent.weights.begin(), parent.weights.end(), 0.0L);
    double dxProjected = fx(parent.xmax, threshold) - fx(parent.xmin, threshold);
    double dyProjected = aspect_ratio * (fy(parent.ymax) - fy(parent.ymin));
    //if ((total_weight > 5000) || (dxProjected > 0.1) || (dyProjected > 0.1) || (dxProjected > 1.5 * dyProjected) || (dyProjected > 1.5 * dxProjected)) {
    if (((total_weight > 5000) && (dxProjected > 0.01) && (dyProjected > 0.01)) || (dxProjected > 0.4) || (dyProjected > 0.4)) {
        if (dxProjected > 1.5 * dyProjected) { // divide horizontally into two
            Node childL = parent.childL(fx, fy, fxInv, fyInv, threshold);
            //Node childR = parent.childR(fx, fy, fxInv, fyInv, threshold);
            childL.x.resize(0); childL.y.resize(0); childL.weights.resize(0);
            //childR.x.resize(0); childR.y.resize(0); childR.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childL.inside(parent.x[i], parent.y[i])) {
                    childL.x.push_back(parent.x[i]);
                    childL.y.push_back(parent.y[i]);
                    childL.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childR_in_place(fx, fy, fxInv, fyInv, threshold);
            Subdivide(nodes, childL, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide(nodes, parent, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
        }
        else if (dyProjected > 1.5 * dxProjected) { // divide vertically into two
            Node childB = parent.childB(fx, fy, fxInv, fyInv, threshold);
            //Node childT = parent.childT(fx, fy, fxInv, fyInv, threshold);
            childB.x.resize(0); childB.y.resize(0); childB.weights.resize(0);
            //childT.x.resize(0); childT.y.resize(0); childT.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childB.inside(parent.x[i], parent.y[i])) {
                    childB.x.push_back(parent.x[i]);
                    childB.y.push_back(parent.y[i]);
                    childB.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childT_in_place(fx, fy, fxInv, fyInv, threshold);
            Subdivide(nodes, childB, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide(nodes, parent, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
        }
        else {
            Node child1 = parent.child1(fx, fy, fxInv, fyInv, threshold);
            Node child2 = parent.child2(fx, fy, fxInv, fyInv, threshold);
            Node child3 = parent.child3(fx, fy, fxInv, fyInv, threshold);
            //Node child4 = parent.child4(fx, fy, fxInv, fyInv, threshold);

            child1.x.resize(0); child1.y.resize(0); child1.weights.resize(0);
            child2.x.resize(0); child2.y.resize(0); child2.weights.resize(0);
            child3.x.resize(0); child3.y.resize(0); child3.weights.resize(0);
            //child4.x.resize(0); child4.y.resize(0); child4.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (child1.inside(parent.x[i], parent.y[i])) {
                    child1.x.push_back(parent.x[i]);
                    child1.y.push_back(parent.y[i]);
                    child1.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child2.inside(parent.x[i], parent.y[i])) {
                    child2.x.push_back(parent.x[i]);
                    child2.y.push_back(parent.y[i]);
                    child2.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child3.inside(parent.x[i], parent.y[i])) {
                    child3.x.push_back(parent.x[i]);
                    child3.y.push_back(parent.y[i]);
                    child3.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.child4_in_place(fx, fy, fxInv, fyInv, threshold);
            Subdivide(nodes, child1, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide(nodes, child2, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide(nodes, child3, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide(nodes, parent, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
        }
    }
    else {
        parent.weight = total_weight / N;
        parent.weights.resize(0);
        parent.x.resize(0);
        parent.y.resize(0);
        nodes.push_back(parent);
        //std::cout << nodes.size() << " " << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << std::endl;
    }
}

void Subdivide(std::vector<Node>& nodes, Node& parent, const int& N, const double& threshold, const double& aspect_ratio) {
    double total_weight = std::accumulate(parent.weights.begin(), parent.weights.end(), 0.0L);
    double dxProjected = sgn(parent.xmax) * std::log10((std::fabs(parent.xmax) + threshold) / threshold) - sgn(parent.xmin) * std::log10((std::fabs(parent.xmin) + threshold) / threshold);
    double dyProjected = aspect_ratio * (std::log10(parent.ymax / parent.ymin));
    //std::cout << "node, w = " << total_weight;
    //if ((total_weight > 5000) || (dxProjected > 0.1) || (dyProjected > 0.1) || (dxProjected > 1.5 * dyProjected) || (dyProjected > 1.5 * dxProjected)) {
    if (((total_weight > 5000) && (dxProjected > 0.01) && (dyProjected > 0.01)) || (dxProjected > 0.4) || (dyProjected > 0.4)) {
        if (dxProjected > 1.5 * dyProjected) { // divide horizontally into two
            //std::cout << ", dividing horizontally" << std::endl;
            Node childL = parent.childL(threshold);
            //Node childR = parent.childR(fx, fy, fxInv, fyInv, threshold);
            childL.x.resize(0); childL.y.resize(0); childL.weights.resize(0);
            //childR.x.resize(0); childR.y.resize(0); childR.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childL.inside(parent.x[i], parent.y[i])) {
                    childL.x.push_back(parent.x[i]);
                    childL.y.push_back(parent.y[i]);
                    childL.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childR_in_place(threshold);
            Subdivide(nodes, childL, N, threshold, aspect_ratio);
            Subdivide(nodes, parent, N, threshold, aspect_ratio);
        }
        else if (dyProjected > 1.5 * dxProjected) { // divide vertically into two
            //std::cout << ", dividing vertically" << std::endl;
            Node childB = parent.childB();
            //Node childT = parent.childT(fx, fy, fxInv, fyInv, threshold);
            childB.x.resize(0); childB.y.resize(0); childB.weights.resize(0);
            //childT.x.resize(0); childT.y.resize(0); childT.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childB.inside(parent.x[i], parent.y[i])) {
                    childB.x.push_back(parent.x[i]);
                    childB.y.push_back(parent.y[i]);
                    childB.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childT_in_place();
            Subdivide(nodes, childB, N, threshold, aspect_ratio);
            Subdivide(nodes, parent, N, threshold, aspect_ratio);
        }
        else {
            //std::cout << ", dividing into four" << std::endl;
            Node child1 = parent.child1(threshold);
            Node child2 = parent.child2(threshold);
            Node child3 = parent.child3(threshold);

            child1.x.resize(0); child1.y.resize(0); child1.weights.resize(0);
            child2.x.resize(0); child2.y.resize(0); child2.weights.resize(0);
            child3.x.resize(0); child3.y.resize(0); child3.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (child1.inside(parent.x[i], parent.y[i])) {
                    child1.x.push_back(parent.x[i]);
                    child1.y.push_back(parent.y[i]);
                    child1.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child2.inside(parent.x[i], parent.y[i])) {
                    child2.x.push_back(parent.x[i]);
                    child2.y.push_back(parent.y[i]);
                    child2.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child3.inside(parent.x[i], parent.y[i])) {
                    child3.x.push_back(parent.x[i]);
                    child3.y.push_back(parent.y[i]);
                    child3.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.child4_in_place(threshold);
            Subdivide(nodes, child1, N, threshold, aspect_ratio);
            Subdivide(nodes, child2, N, threshold, aspect_ratio);
            Subdivide(nodes, child3, N, threshold, aspect_ratio);
            Subdivide(nodes, parent, N, threshold, aspect_ratio);
        }
    }
    else {
        //std::cout << nodes.size() << " " << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << std::endl;
        parent.weight = total_weight / N;
        parent.weights.resize(0);
        parent.x.resize(0);
        parent.y.resize(0);
        nodes.push_back(parent);
    }
}

void Subdivide_energy_save(std::ofstream& nodes_out, Node& parent, const int& N, const double& threshold, const double& aspect_ratio, const double& dx, const double& dy) {
    double total_weight = std::accumulate(parent.weights.begin(), parent.weights.end(), 0.0L);
    double dxProjected = sgn(parent.xmax) * std::log10((std::fabs(parent.xmax) + threshold) / threshold) - sgn(parent.xmin) * std::log10((std::fabs(parent.xmin) + threshold) / threshold);
    double dyProjected = aspect_ratio * (std::log10(parent.ymax / parent.ymin));
    //std::cout << "node, w = " << total_weight;
    //if ((total_weight > 5000) || (dxProjected > 0.1) || (dyProjected > 0.1) || (dxProjected > 1.5 * dyProjected) || (dyProjected > 1.5 * dxProjected)) {
    if (((total_weight > 4000) && (dxProjected > dx / 1000) && (dyProjected > dy / 1000)) || (dxProjected > dx / 50) || (dyProjected > dy / 50)) {
        if (dxProjected > 1.5 * dyProjected) { // divide horizontally into two
            //std::cout << ", dividing horizontally" << std::endl;
            Node childL = parent.childL(threshold);
            //Node childR = parent.childR(fx, fy, fxInv, fyInv, threshold);
            childL.x.resize(0); childL.y.resize(0); childL.weights.resize(0);
            //childR.x.resize(0); childR.y.resize(0); childR.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childL.inside(parent.x[i], parent.y[i])) {
                    childL.x.push_back(parent.x[i]);
                    childL.y.push_back(parent.y[i]);
                    childL.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childR_in_place(threshold);
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_energy_save(nodes_out, childL, N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_save(nodes_out, parent, N, threshold, aspect_ratio, dx, dy);
        }
        else if (dyProjected > 1.5 * dxProjected) { // divide vertically into two
            //std::cout << ", dividing vertically" << std::endl;
            Node childB = parent.childB();
            //Node childT = parent.childT(fx, fy, fxInv, fyInv, threshold);
            childB.x.resize(0); childB.y.resize(0); childB.weights.resize(0);
            //childT.x.resize(0); childT.y.resize(0); childT.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childB.inside(parent.x[i], parent.y[i])) {
                    childB.x.push_back(parent.x[i]);
                    childB.y.push_back(parent.y[i]);
                    childB.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childT_in_place();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_energy_save(nodes_out, childB, N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_save(nodes_out, parent, N, threshold, aspect_ratio, dx, dy);
        }
        else {
            //std::cout << ", dividing into four" << std::endl;
            Node child1 = parent.child1(threshold);
            Node child2 = parent.child2(threshold);
            Node child3 = parent.child3(threshold);

            child1.x.resize(0); child1.y.resize(0); child1.weights.resize(0);
            child2.x.resize(0); child2.y.resize(0); child2.weights.resize(0);
            child3.x.resize(0); child3.y.resize(0); child3.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (child1.inside(parent.x[i], parent.y[i])) {
                    child1.x.push_back(parent.x[i]);
                    child1.y.push_back(parent.y[i]);
                    child1.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child2.inside(parent.x[i], parent.y[i])) {
                    child2.x.push_back(parent.x[i]);
                    child2.y.push_back(parent.y[i]);
                    child2.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child3.inside(parent.x[i], parent.y[i])) {
                    child3.x.push_back(parent.x[i]);
                    child3.y.push_back(parent.y[i]);
                    child3.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.child4_in_place(threshold);
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_energy_save(nodes_out, child1, N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_save(nodes_out, child2, N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_save(nodes_out, child3, N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_save(nodes_out, parent, N, threshold, aspect_ratio, dx, dy);
        }
    }
    else {
        //std::cout << nodes.size() << " " << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << std::endl;
        parent.weight = total_weight / N;
        parent.weights.resize(0);
        parent.x.resize(0);
        parent.y.resize(0);
        parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
        //nodes.push_back(parent);
        if (total_weight > 0) {
            nodes_out << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << " " << parent.weight << std::endl;
        }
    }
}

void Subdivide_energy_nodes(std::vector<Node>& nodes, Node& parent, const int& N, const int& threshold_N, const double& threshold, const double& aspect_ratio, const double& dx, const double& dy) {
    double total_weight = std::accumulate(parent.weights.begin(), parent.weights.end(), 0.0L);
    double dxProjected = symlog(parent.xmax, threshold) - symlog(parent.xmin, threshold);
    double dyProjected = aspect_ratio * (std::log10(parent.ymax / parent.ymin));
    //std::cout << "node, w = " << total_weight;
    //if ((total_weight > 5000) || (dxProjected > 0.1) || (dyProjected > 0.1) || (dxProjected > 1.5 * dyProjected) || (dyProjected > 1.5 * dxProjected)) {
    if (((total_weight > threshold_N) && (dxProjected > dx / 500) && (dyProjected > dy * aspect_ratio / 500)) || (dxProjected > dx / 10) || (dyProjected > dy * aspect_ratio / 10)) {
        if (dxProjected > 1.5 * dyProjected) { // divide horizontally into two
            //std::cout << ", dividing horizontally" << std::endl;
            Node childL = parent.childL(threshold);
            //Node childR = parent.childR(fx, fy, fxInv, fyInv, threshold);
            childL.x.resize(0); childL.y.resize(0); childL.weights.resize(0);
            //childR.x.resize(0); childR.y.resize(0); childR.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childL.inside(parent.x[i], parent.y[i])) {
                    childL.x.push_back(parent.x[i]);
                    childL.y.push_back(parent.y[i]);
                    childL.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childR_in_place(threshold);
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_energy_nodes(nodes, childL, N, threshold_N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_nodes(nodes, parent, N, threshold_N, threshold, aspect_ratio, dx, dy);
        }
        else if (dyProjected > 1.5 * dxProjected) { // divide vertically into two
            //std::cout << ", dividing vertically" << std::endl;
            Node childB = parent.childB();
            //Node childT = parent.childT(fx, fy, fxInv, fyInv, threshold);
            childB.x.resize(0); childB.y.resize(0); childB.weights.resize(0);
            //childT.x.resize(0); childT.y.resize(0); childT.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childB.inside(parent.x[i], parent.y[i])) {
                    childB.x.push_back(parent.x[i]);
                    childB.y.push_back(parent.y[i]);
                    childB.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childT_in_place();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_energy_nodes(nodes, childB, N, threshold_N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_nodes(nodes, parent, N, threshold_N, threshold, aspect_ratio, dx, dy);
        }
        else {
            //std::cout << ", dividing into four" << std::endl;
            Node child1 = parent.child1(threshold);
            Node child2 = parent.child2(threshold);
            Node child3 = parent.child3(threshold);

            child1.x.resize(0); child1.y.resize(0); child1.weights.resize(0);
            child2.x.resize(0); child2.y.resize(0); child2.weights.resize(0);
            child3.x.resize(0); child3.y.resize(0); child3.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (child1.inside(parent.x[i], parent.y[i])) {
                    child1.x.push_back(parent.x[i]);
                    child1.y.push_back(parent.y[i]);
                    child1.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child2.inside(parent.x[i], parent.y[i])) {
                    child2.x.push_back(parent.x[i]);
                    child2.y.push_back(parent.y[i]);
                    child2.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child3.inside(parent.x[i], parent.y[i])) {
                    child3.x.push_back(parent.x[i]);
                    child3.y.push_back(parent.y[i]);
                    child3.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.child4_in_place(threshold);
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_energy_nodes(nodes, child1, N, threshold_N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_nodes(nodes, child2, N, threshold_N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_nodes(nodes, child3, N, threshold_N, threshold, aspect_ratio, dx, dy);
            Subdivide_energy_nodes(nodes, parent, N, threshold_N, threshold, aspect_ratio, dx, dy);
        }
    }
    else {
        //std::cout << nodes.size() << " " << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << std::endl;
        parent.weight = total_weight / N;
        parent.weights.resize(0);
        parent.x.resize(0);
        parent.y.resize(0);
        parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
        parent.pdf = parent.weight / parent.d2e();
        //nodes.push_back(parent);
        if (parent.pdf > 0) nodes.push_back(parent);
    }
}

void Subdivide_log_linear(std::vector<Node>& nodes, Node& parent, const int& N, const int& threshold_N, const double& aspect_ratio, const double& dx, const double& dy) {
    double total_weight = std::accumulate(parent.weights.begin(), parent.weights.end(), 0.0L);
    double dxProjected = std::log10(parent.xmax / parent.xmin);
    double dyProjected = aspect_ratio * (parent.ymax - parent.ymin);
    //std::cout << "node, w = " << total_weight;
    //if ((total_weight > 5000) || (dxProjected > 0.1) || (dyProjected > 0.1) || (dxProjected > 1.5 * dyProjected) || (dyProjected > 1.5 * dxProjected)) {
    if (((total_weight > threshold_N) && (dxProjected > dx / 500) && (dyProjected > dy * aspect_ratio / 500)) || (dxProjected > dx / 10) || (dyProjected > dy * aspect_ratio / 10)) {
        if (dxProjected > 1.5 * dyProjected) { // divide horizontally into two
            //std::cout << ", dividing horizontally" << std::endl;
            Node childL = parent.childL_log();
            //Node childR = parent.childR(fx, fy, fxInv, fyInv, threshold);
            childL.x.resize(0); childL.y.resize(0); childL.weights.resize(0);
            //childR.x.resize(0); childR.y.resize(0); childR.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childL.inside(parent.x[i], parent.y[i])) {
                    childL.x.push_back(parent.x[i]);
                    childL.y.push_back(parent.y[i]);
                    childL.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childR_in_place_log();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_log_linear(nodes, childL, N, threshold_N, aspect_ratio, dx, dy);
            Subdivide_log_linear(nodes, parent, N, threshold_N, aspect_ratio, dx, dy);
        }
        else if (dyProjected > 1.5 * dxProjected) { // divide vertically into two
            //std::cout << ", dividing vertically" << std::endl;
            Node childB = parent.childB_linear();
            //Node childT = parent.childT(fx, fy, fxInv, fyInv, threshold);
            childB.x.resize(0); childB.y.resize(0); childB.weights.resize(0);
            //childT.x.resize(0); childT.y.resize(0); childT.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childB.inside(parent.x[i], parent.y[i])) {
                    childB.x.push_back(parent.x[i]);
                    childB.y.push_back(parent.y[i]);
                    childB.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childT_in_place_linear();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_log_linear(nodes, childB, N, threshold_N, aspect_ratio, dx, dy);
            Subdivide_log_linear(nodes, parent, N, threshold_N, aspect_ratio, dx, dy);
        }
        else {
            //std::cout << ", dividing into four" << std::endl;
            Node child1 = parent.child1_loglinear();
            Node child2 = parent.child2_loglinear();
            Node child3 = parent.child3_loglinear();

            child1.x.resize(0); child1.y.resize(0); child1.weights.resize(0);
            child2.x.resize(0); child2.y.resize(0); child2.weights.resize(0);
            child3.x.resize(0); child3.y.resize(0); child3.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (child1.inside(parent.x[i], parent.y[i])) {
                    child1.x.push_back(parent.x[i]);
                    child1.y.push_back(parent.y[i]);
                    child1.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child2.inside(parent.x[i], parent.y[i])) {
                    child2.x.push_back(parent.x[i]);
                    child2.y.push_back(parent.y[i]);
                    child2.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child3.inside(parent.x[i], parent.y[i])) {
                    child3.x.push_back(parent.x[i]);
                    child3.y.push_back(parent.y[i]);
                    child3.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.child4_in_place_loglinear();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_log_linear(nodes, child1, N, threshold_N, aspect_ratio, dx, dy);
            Subdivide_log_linear(nodes, child2, N, threshold_N, aspect_ratio, dx, dy);
            Subdivide_log_linear(nodes, child3, N, threshold_N, aspect_ratio, dx, dy);
            Subdivide_log_linear(nodes, parent, N, threshold_N, aspect_ratio, dx, dy);
        }
    }
    else {
        //std::cout << nodes.size() << " " << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << std::endl;
        parent.weight = total_weight / N;
        parent.weights.resize(0);
        parent.x.resize(0);
        parent.y.resize(0);
        parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
        parent.pdf = parent.weight / parent.d2e();
        //nodes.push_back(parent);
        if (parent.pdf > 0) nodes.push_back(parent);
    }
}

void Subdivide_linear_save(std::ofstream& nodes_out, Node& parent, const int& N, const double& aspect_ratio) {
    double total_weight = std::accumulate(parent.weights.begin(), parent.weights.end(), 0.0L);
    double dxProjected = parent.xmax - parent.xmin;
    double dyProjected = aspect_ratio * (parent.ymax - parent.ymin);
    //std::cout << "node, w = " << total_weight;
    //if ((total_weight > 5000) || (dxProjected > 0.1) || (dyProjected > 0.1) || (dxProjected > 1.5 * dyProjected) || (dyProjected > 1.5 * dxProjected)) {
    if (((total_weight > 4000) && (dxProjected > 0.01) && (dyProjected > 0.01)) || (dxProjected > 0.4) || (dyProjected > 0.4)) {
        if (dxProjected > 1.5 * dyProjected) { // divide horizontally into two
            //std::cout << ", dividing horizontally" << std::endl;
            Node childL = parent.childL_linear();
            //Node childR = parent.childR(fx, fy, fxInv, fyInv, threshold);
            childL.x.resize(0); childL.y.resize(0); childL.weights.resize(0);
            //childR.x.resize(0); childR.y.resize(0); childR.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childL.inside(parent.x[i], parent.y[i])) {
                    childL.x.push_back(parent.x[i]);
                    childL.y.push_back(parent.y[i]);
                    childL.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childR_in_place_linear();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_linear_save(nodes_out, childL, N, aspect_ratio);
            Subdivide_linear_save(nodes_out, parent, N, aspect_ratio);
        }
        else if (dyProjected > 1.5 * dxProjected) { // divide vertically into two
            //std::cout << ", dividing vertically" << std::endl;
            Node childB = parent.childB();
            //Node childT = parent.childT(fx, fy, fxInv, fyInv, threshold);
            childB.x.resize(0); childB.y.resize(0); childB.weights.resize(0);
            //childT.x.resize(0); childT.y.resize(0); childT.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childB.inside(parent.x[i], parent.y[i])) {
                    childB.x.push_back(parent.x[i]);
                    childB.y.push_back(parent.y[i]);
                    childB.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childT_in_place();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_linear_save(nodes_out, childB, N, aspect_ratio);
            Subdivide_linear_save(nodes_out, parent, N, aspect_ratio);
        }
        else {
            //std::cout << ", dividing into four" << std::endl;
            Node child1 = parent.child1();
            Node child2 = parent.child2();
            Node child3 = parent.child3();

            child1.x.resize(0); child1.y.resize(0); child1.weights.resize(0);
            child2.x.resize(0); child2.y.resize(0); child2.weights.resize(0);
            child3.x.resize(0); child3.y.resize(0); child3.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (child1.inside(parent.x[i], parent.y[i])) {
                    child1.x.push_back(parent.x[i]);
                    child1.y.push_back(parent.y[i]);
                    child1.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child2.inside(parent.x[i], parent.y[i])) {
                    child2.x.push_back(parent.x[i]);
                    child2.y.push_back(parent.y[i]);
                    child2.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child3.inside(parent.x[i], parent.y[i])) {
                    child3.x.push_back(parent.x[i]);
                    child3.y.push_back(parent.y[i]);
                    child3.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.child4_in_place();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_linear_save(nodes_out, child1, N, aspect_ratio);
            Subdivide_linear_save(nodes_out, child2, N, aspect_ratio);
            Subdivide_linear_save(nodes_out, child3, N, aspect_ratio);
            Subdivide_linear_save(nodes_out, parent, N, aspect_ratio);
        }
    }
    else {
        //std::cout << nodes.size() << " " << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << std::endl;
        parent.weight = total_weight / N;
        parent.weights.resize(0);
        parent.x.resize(0);
        parent.y.resize(0);
        parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
        //nodes.push_back(parent);
        if (total_weight > 0) {
            nodes_out << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << " " << parent.weight << std::endl;
        }
    }
}

void Subdivide_loglog_save(std::ofstream& nodes_out, Node& parent, const int& N, const double& aspect_ratio, const double &dx, const double &dy) {
    double total_weight = std::accumulate(parent.weights.begin(), parent.weights.end(), 0.0L);
    double dxProjected = std::log10(parent.xmax / parent.xmin);
    double dyProjected = aspect_ratio * std::log10(parent.ymax / parent.ymin);
    //std::cout << "node, w = " << total_weight;
    //if ((total_weight > 5000) || (dxProjected > 0.1) || (dyProjected > 0.1) || (dxProjected > 1.5 * dyProjected) || (dyProjected > 1.5 * dxProjected)) {
    if (((total_weight > 10000) && (dxProjected > dx / 100) && (dyProjected > dy / 100)) || (dxProjected > dx / 20) || (dyProjected > dy / 20)) {
    //if (total_weight > 4000) {
        if (dxProjected > 1.5 * dyProjected) { // divide horizontally into two
            //std::cout << ", dividing horizontally" << std::endl;
            Node childL = parent.childL_log();
            //Node childR = parent.childR(fx, fy, fxInv, fyInv, threshold);
            childL.x.resize(0); childL.y.resize(0); childL.weights.resize(0);
            //childR.x.resize(0); childR.y.resize(0); childR.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childL.inside(parent.x[i], parent.y[i])) {
                    childL.x.push_back(parent.x[i]);
                    childL.y.push_back(parent.y[i]);
                    childL.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childR_in_place_log();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_loglog_save(nodes_out, childL, N, aspect_ratio, dx, dy);
            Subdivide_loglog_save(nodes_out, parent, N, aspect_ratio, dx, dy);
        }
        else if (dyProjected > 1.5 * dxProjected) { // divide vertically into two
            //std::cout << ", dividing vertically" << std::endl;
            Node childB = parent.childB();
            //Node childT = parent.childT(fx, fy, fxInv, fyInv, threshold);
            childB.x.resize(0); childB.y.resize(0); childB.weights.resize(0);
            //childT.x.resize(0); childT.y.resize(0); childT.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childB.inside(parent.x[i], parent.y[i])) {
                    childB.x.push_back(parent.x[i]);
                    childB.y.push_back(parent.y[i]);
                    childB.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.childT_in_place();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_loglog_save(nodes_out, childB, N, aspect_ratio, dx, dy);
            Subdivide_loglog_save(nodes_out, parent, N, aspect_ratio, dx, dy);
        }
        else {
            //std::cout << ", dividing into four" << std::endl;
            Node child1 = parent.child1_log();
            Node child2 = parent.child2_log();
            Node child3 = parent.child3_log();

            child1.x.resize(0); child1.y.resize(0); child1.weights.resize(0);
            child2.x.resize(0); child2.y.resize(0); child2.weights.resize(0);
            child3.x.resize(0); child3.y.resize(0); child3.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (child1.inside(parent.x[i], parent.y[i])) {
                    child1.x.push_back(parent.x[i]);
                    child1.y.push_back(parent.y[i]);
                    child1.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child2.inside(parent.x[i], parent.y[i])) {
                    child2.x.push_back(parent.x[i]);
                    child2.y.push_back(parent.y[i]);
                    child2.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child3.inside(parent.x[i], parent.y[i])) {
                    child3.x.push_back(parent.x[i]);
                    child3.y.push_back(parent.y[i]);
                    child3.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
            }
            parent.child4_in_place_log();
            parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            Subdivide_loglog_save(nodes_out, child1, N, aspect_ratio, dx, dy);
            Subdivide_loglog_save(nodes_out, child2, N, aspect_ratio, dx, dy);
            Subdivide_loglog_save(nodes_out, child3, N, aspect_ratio, dx, dy);
            Subdivide_loglog_save(nodes_out, parent, N, aspect_ratio, dx, dy);
        }
    }
    else {
        //std::cout << nodes.size() << " " << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << std::endl;
        parent.weight = total_weight / N;
        parent.weights.resize(0);
        parent.x.resize(0);
        parent.y.resize(0);
        parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
        //nodes.push_back(parent);
        if (total_weight > 0) {
            nodes_out << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << " " << parent.weight << std::endl;
        }
    }
}

void Subdivide_energy_save(const std::string &filename, Node& parent, const int& N, double(&fx) (const double&, const double&), double(&fy) (const double&), double(&fxInv) (const double&, const double&), double(&fyInv) (const double&), const double& threshold, const double& aspect_ratio) {
    double total_weight = std::accumulate(parent.weights.begin(), parent.weights.end(), 0.0L);
    double dxProjected = fx(parent.xmax, threshold) - fx(parent.xmin, threshold);
    double dyProjected = aspect_ratio * (fy(parent.ymax) - fy(parent.ymin));
    //if ((total_weight > 5000) || (dxProjected > 0.1) || (dyProjected > 0.1) || (dxProjected > 1.5 * dyProjected) || (dyProjected > 1.5 * dxProjected)) {
    if ((total_weight > 5000) && (dxProjected > 0.01) && (dyProjected > 0.01)) {
        if (dxProjected > 1.5 * dyProjected) { // divide horizontally into two
            Node childL = parent.childL(fx, fy, fxInv, fyInv, threshold);
            //Node childR = parent.childR(fx, fy, fxInv, fyInv, threshold);
            childL.x.resize(0); childL.y.resize(0); childL.weights.resize(0);
            //childR.x.resize(0); childR.y.resize(0); childR.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childL.inside(parent.x[i], parent.y[i])) {
                    childL.x.push_back(parent.x[i]);
                    childL.y.push_back(parent.y[i]);
                    childL.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
                parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            }
            parent.childR_in_place(fx, fy, fxInv, fyInv, threshold);
            Subdivide_energy_save(filename, childL, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide_energy_save(filename, parent, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
        }
        else if (dyProjected > 1.5 * dxProjected) { // divide vertically into two
            Node childB = parent.childB(fx, fy, fxInv, fyInv, threshold);
            //Node childT = parent.childT(fx, fy, fxInv, fyInv, threshold);
            childB.x.resize(0); childB.y.resize(0); childB.weights.resize(0);
            //childT.x.resize(0); childT.y.resize(0); childT.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (childB.inside(parent.x[i], parent.y[i])) {
                    childB.x.push_back(parent.x[i]);
                    childB.y.push_back(parent.y[i]);
                    childB.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
                parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            }
            parent.childT_in_place(fx, fy, fxInv, fyInv, threshold);
            Subdivide_energy_save(filename, childB, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide_energy_save(filename, parent, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
        }
        else {
            Node child1 = parent.child1(fx, fy, fxInv, fyInv, threshold);
            Node child2 = parent.child2(fx, fy, fxInv, fyInv, threshold);
            Node child3 = parent.child3(fx, fy, fxInv, fyInv, threshold);
            //Node child4 = parent.child4(fx, fy, fxInv, fyInv, threshold);

            child1.x.resize(0); child1.y.resize(0); child1.weights.resize(0);
            child2.x.resize(0); child2.y.resize(0); child2.weights.resize(0);
            child3.x.resize(0); child3.y.resize(0); child3.weights.resize(0);
            //child4.x.resize(0); child4.y.resize(0); child4.weights.resize(0);
            int i = 0;
            while (i < parent.x.size()) {
                if (child1.inside(parent.x[i], parent.y[i])) {
                    child1.x.push_back(parent.x[i]);
                    child1.y.push_back(parent.y[i]);
                    child1.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child2.inside(parent.x[i], parent.y[i])) {
                    child2.x.push_back(parent.x[i]);
                    child2.y.push_back(parent.y[i]);
                    child2.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else if (child3.inside(parent.x[i], parent.y[i])) {
                    child3.x.push_back(parent.x[i]);
                    child3.y.push_back(parent.y[i]);
                    child3.weights.push_back(parent.weights[i]);
                    parent.x[i] = parent.x.back();
                    parent.y[i] = parent.y.back();
                    parent.weights[i] = parent.weights.back();
                    parent.x.pop_back();
                    parent.y.pop_back();
                    parent.weights.pop_back();
                }
                else i++;
                parent.x.shrink_to_fit(); parent.y.shrink_to_fit(); parent.weights.shrink_to_fit();
            }
            parent.child4_in_place(fx, fy, fxInv, fyInv, threshold);
            Subdivide_energy_save(filename, child1, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide_energy_save(filename, child2, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide_energy_save(filename, child3, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
            Subdivide_energy_save(filename, parent, N, fx, fy, fxInv, fyInv, threshold, aspect_ratio);
        }
    }
    else {
        parent.weight = total_weight / N;
        parent.weights.resize(0);
        parent.x.resize(0);
        parent.y.resize(0);
        if (total_weight > 0) {
            std::ofstream out;
            out.open(filename, std::ios_base::app);
            out << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << " " << parent.weight << std::endl;
            out.close();
        }
        //nodes.push_back(parent);
        //std::cout << nodes.size() << " " << parent.xmin << " " << parent.xmax << " " << parent.ymin << " " << parent.ymax << std::endl;
    }
}

void Node::set_bounds() {
    xmin = 1e10;
    xmax = -1e10;
    ymin = 1e10;
    ymax = -1e10;
    for (int i = 0; i < x.size(); i++) {
        if (xmin > x[i]) xmin = x[i];
        if (xmax < x[i]) xmax = x[i];
        if (ymin > y[i]) ymin = y[i];
        if (ymax < y[i]) ymax = y[i];
    }
    std::cout << "Master dimensions: (" << xmin << ", " << xmax << ") x (" << ymin << ", " << ymax << ")" << std::endl;
}

void Node::set_bounds(const double &Xmin, const double &Xmax, const double &Ymin, const double &Ymax) {
    xmin = Xmin;
    xmax = Xmax;
    ymin = Ymin;
    ymax = Ymax;
    for (int i = 0; i < x.size(); i++) {
        if (!inside(x[i], y[i])) {
            x[i] = x.back();
            y[i] = y.back();
            weights[i] = weights.back();
            x.pop_back();
            y.pop_back();
            weights.pop_back();
        }
        else i++;
    }
    std::cout << "Master dimensions: (" << xmin << ", " << xmax << ") x (" << ymin << ", " << ymax << ")" << std::endl;
}

Node const Node::child1() const {
    Node res;
    res.xmin = xmin;
    res.xmax = 0.5 * (xmin + xmax);
    res.ymin = ymin;
    res.ymax = 0.5 * (ymin + ymax);
    return res;
}

Node const Node::child2() const {
    Node res;
    res.xmin = xmin;
    res.xmax = 0.5 * (xmin + xmax);
    res.ymin = 0.5 * (ymin + ymax);
    res.ymax = ymax;
    return res;
}

Node const Node::child3() const {
    Node res;
    res.xmin = 0.5 * (xmin + xmax);
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = 0.5 * (ymin + ymax);
    return res;
}

Node const Node::child4() const {
    Node res;
    res.xmin = 0.5 * (xmin + xmax);
    res.xmax = xmax;
    res.ymin = 0.5 * (ymin + ymax);
    res.ymax = ymax;
    return res;
}

Node const Node::child1_log() const {
    Node res;
    res.xmin = xmin;
    res.xmax = std::sqrt(xmin * xmax);
    res.ymin = ymin;
    res.ymax = std::sqrt(ymin * ymax);
    return res;
}

Node const Node::child2_log() const {
    Node res;
    res.xmin = xmin;
    res.xmax = std::sqrt(xmin * xmax);
    res.ymin = std::sqrt(ymin * ymax);
    res.ymax = ymax;
    return res;
}

Node const Node::child3_log() const {
    Node res;
    res.xmin = std::sqrt(xmin * xmax);
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = std::sqrt(ymin * ymax);
    return res;
}

Node const Node::child4_log() const {
    Node res;
    res.xmin = std::sqrt(xmin * xmax);
    res.xmax = xmax;
    res.ymin = std::sqrt(ymin * ymax);
    res.ymax = ymax;
    return res;
}

Node const Node::child1_loglinear() const {
    Node res;
    res.xmin = xmin;
    res.xmax = std::sqrt(xmin * xmax);
    res.ymin = ymin;
    res.ymax = 0.5 * (ymin + ymax);
    return res;
}

Node const Node::child2_loglinear() const {
    Node res;
    res.xmin = xmin;
    res.xmax = std::sqrt(xmin * xmax);
    res.ymin = 0.5 * (ymin + ymax);
    res.ymax = ymax;
    return res;
}

Node const Node::child3_loglinear() const {
    Node res;
    res.xmin = std::sqrt(xmin * xmax);
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = 0.5 * (ymin + ymax);
    return res;
}

Node const Node::child4_loglinear() const {
    Node res;
    res.xmin = std::sqrt(xmin * xmax);
    res.xmax = xmax;
    res.ymin = 0.5 * (ymin + ymax);
    res.ymax = ymax;
    return res;
}

Node const Node::child1(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) const {
    Node res;
    res.xmin = xmin;
    res.xmax = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    res.ymin = ymin;
    res.ymax = fyInv(0.5 * (fy(ymin) + fy(ymax)));
    return res;
}

Node const Node::child1(const double& threshold) const {
    Node res;
    res.xmin = xmin;
    res.xmax = 0.5 * (sgn(xmin) * std::log10((std::fabs(xmin) + threshold) / threshold) + sgn(xmax) * std::log10((std::fabs(xmax) + threshold) / threshold));
    res.xmax = sgn(res.xmax) * threshold * (pow(10, std::fabs(res.xmax)) - 1);
    res.ymin = ymin;
    res.ymax = std::sqrt(ymin * ymax);
    return res;
}

Node const Node::child2(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) const {
    Node res;
    res.xmin = xmin;
    res.xmax = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    res.ymin = fyInv(0.5 * (fy(ymin) + fy(ymax)));
    res.ymax = ymax;
    return res;
}

Node const Node::child2(const double& threshold) const {
    Node res;
    res.xmin = xmin;
    res.xmax = 0.5 * (sgn(xmin) * std::log10((std::fabs(xmin) + threshold) / threshold) + sgn(xmax) * std::log10((std::fabs(xmax) + threshold) / threshold));
    res.xmax = sgn(res.xmax) * threshold * (pow(10, std::fabs(res.xmax)) - 1);
    res.ymin = std::sqrt(ymin * ymax);
    res.ymax = ymax;
    return res;
}

Node const Node::child3(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) const {
    Node res;
    res.xmin = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = fyInv(0.5 * (fy(ymin) + fy(ymax)));
    return res;
}

Node const Node::child3(const double& threshold) const {
    Node res;
    res.xmin = 0.5 * (sgn(xmin) * std::log10((std::fabs(xmin) + threshold) / threshold) + sgn(xmax) * std::log10((std::fabs(xmax) + threshold) / threshold));
    res.xmin = sgn(res.xmin) * threshold * (pow(10, std::fabs(res.xmin)) - 1);
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = std::sqrt(ymin * ymax);
    return res;
}

Node const Node::child4(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) const {
    Node res;
    res.xmin = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    res.xmax = xmax;
    res.ymin = fyInv(0.5 * (fy(ymin) + fy(ymax)));
    res.ymax = ymax;
    return res;
}

Node const Node::child4(const double& threshold) const {
    Node res;
    res.xmin = 0.5 * (sgn(xmin) * std::log10((std::fabs(xmin) + threshold) / threshold) + sgn(xmax) * std::log10((std::fabs(xmax) + threshold) / threshold));
    res.xmin = sgn(res.xmin) * threshold * (pow(10, std::fabs(res.xmin)) - 1);
    res.xmax = xmax;
    res.ymin = std::sqrt(ymin * ymax);
    res.ymax = ymax;
    return res;
}

Node const Node::childL(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) const {
    Node res;
    res.xmin = xmin;
    res.xmax = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    res.ymin = ymin;
    res.ymax = ymax;
    return res;
}

Node const Node::childL(const double& threshold) const {
    Node res;
    res.xmax = 0.5 * (sgn(xmin) * std::log10((std::fabs(xmin) + threshold) / threshold) + sgn(xmax) * std::log10((std::fabs(xmax) + threshold) / threshold));
    res.xmax = sgn(res.xmax) * threshold * (pow(10, std::fabs(res.xmax)) - 1);
    res.xmin = xmin;
    res.ymin = ymin;
    res.ymax = ymax;
    return res;
}

Node const Node::childL_linear() const {
    Node res;
    res.xmax = 0.5 * (xmax + xmin);
    res.xmin = xmin;
    res.ymin = ymin;
    res.ymax = ymax;
    return res;
}

Node const Node::childL_log() const {
    Node res;
    res.xmax = std::sqrt(xmax * xmin);
    res.xmin = xmin;
    res.ymin = ymin;
    res.ymax = ymax;
    return res;
}

Node const Node::childR(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) const {
    Node res;
    res.xmin = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = ymax;
    return res;
}

Node const Node::childR(const double& threshold) const {
    Node res;
    res.xmin = 0.5 * (sgn(xmin) * std::log10((std::fabs(xmin) + threshold) / threshold) + sgn(xmax) * std::log10((std::fabs(xmax) + threshold) / threshold));
    res.xmin = sgn(res.xmin) * threshold * (pow(10, std::fabs(res.xmin)) - 1);
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = ymax;
    return res;
}

Node const Node::childB(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) const {
    Node res;
    res.xmin = xmin;
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = fyInv(0.5 * (fy(ymin) + fy(ymax)));
    return res;
}

Node const Node::childB() const {
    Node res;
    res.xmin = xmin;
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = std::sqrt(ymin * ymax);
    return res;
}

Node const Node::childB_linear() const {
    Node res;
    res.xmin = xmin;
    res.xmax = xmax;
    res.ymin = ymin;
    res.ymax = 0.5*(ymin + ymax);
    return res;
}

Node const Node::childT(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) const {
    Node res;
    res.xmin = xmin;
    res.xmax = xmax;
    res.ymin = fyInv(0.5 * (fy(ymin) + fy(ymax)));
    res.ymax = ymax;
    return res;
}

Node const Node::childT() const {
    Node res;
    res.xmin = xmin;
    res.xmax = xmax;
    res.ymin = std::sqrt(ymin * ymax);
    res.ymax = ymax;
    return res;
}

void Node::child1_in_place() {
    xmax = 0.5 * (xmin + xmax);
    ymax = 0.5 * (ymin + ymax);
}

void Node::child2_in_place() {
    xmax = 0.5 * (xmin + xmax);
    ymin = 0.5 * (ymin + ymax);
}

void Node::child3_in_place() {
    xmin = 0.5 * (xmin + xmax);
    ymax = 0.5 * (ymin + ymax);
}

void Node::child4_in_place() {
    xmin = 0.5 * (xmin + xmax);
    ymin = 0.5 * (ymin + ymax);
}

void Node::child4_in_place_log() {
    xmin = std::sqrt(xmax * xmin);
    ymin = std::sqrt(ymax * ymin);
}

void Node::child1_in_place(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) {
    xmax = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    ymax = fyInv(0.5 * (fy(ymin) + fy(ymax)));
}

void Node::child2_in_place(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) {
    xmax = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    ymin = fyInv(0.5 * (fy(ymin) + fy(ymax)));
}

void Node::child3_in_place(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) {
    xmin = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    ymax = fyInv(0.5 * (fy(ymin) + fy(ymax)));
}

void Node::child4_in_place(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) {
    xmin = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
    ymin = fyInv(0.5 * (fy(ymin) + fy(ymax)));
}

void Node::child4_in_place(const double& threshold) {
    xmin = 0.5 * (sgn(xmin) * std::log10((std::fabs(xmin) + threshold) / threshold) + sgn(xmax) * std::log10((std::fabs(xmax) + threshold) / threshold));
    xmin = sgn(xmin) * threshold * (pow(10, std::fabs(xmin)) - 1);
    ymin = std::sqrt(ymin * ymax);
}

void Node::child4_in_place_loglinear() {
    xmin = std::sqrt(xmin * xmax);
    ymin = 0.5 * (ymin + ymax);
}

void Node::childL_in_place(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) {
    xmax = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
}

void Node::childR_in_place(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) {
    xmin = fxInv(0.5 * (fx(xmin, threshold) + fx(xmax, threshold)), threshold);
}

void Node::childR_in_place(const double& threshold) {
    xmin = 0.5 * (sgn(xmin) * std::log10((std::fabs(xmin) + threshold) / threshold) + sgn(xmax) * std::log10((std::fabs(xmax) + threshold) / threshold));
    xmin = sgn(xmin) * threshold * (pow(10, std::fabs(xmin)) - 1);
}

void Node::childR_in_place_linear() {
    xmin = 0.5 * (xmin + xmax);
}

void Node::childR_in_place_log() {
    xmin = std::sqrt(xmin * xmax);
}

void Node::childB_in_place(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) {
    ymax = fyInv(0.5 * (fy(ymin) + fy(ymax)));
}

void Node::childT_in_place(double (*fx) (const double&, const double&), double (*fy) (const double&), double (*fxInv) (const double&, const double&), double (*fyInv) (const double&), const double& threshold) {
    ymin = fyInv(0.5 * (fy(ymin) + fy(ymax)));
}

void Node::childT_in_place() {
    ymin = std::sqrt(ymin * ymax);
}

void Node::childT_in_place_linear() {
    ymin = 0.5*(ymin + ymax);
}

double Node::d2e() const {
    return (xmax - xmin) * (ymax - ymin);
}

double Node::center_x() const {
    return 0.5 * (xmin + xmax);
}

double Node::center_y() const {
    return 0.5 * (ymin + ymax);
}

double Node::center_weighted_x(double threshold) const {
    return symlogInv(0.5 * (umin + umax), threshold);
}

double Node::center_weighted_y() const {
    return pow(10, 0.5 * (vmin + vmax));
}

void DataExport(const std::vector<Node>& nodes, const std::string& filename) {
    std::ofstream out;
    out.open(filename);
    std::string line;
    for (auto& node : nodes) {
        out << node.xmin << " " << node.xmax << " " << node.ymin << " " << node.ymax << " " << node.weight << std::endl;
    }
    out.close();
}

void DataExport(const std::multiset<Node>& nodes, const std::string& filename) {
    std::ofstream out;
    out.open(filename);
    std::string line;
    for (auto& node : nodes) {
        out << node.xmin << " " << node.xmax << " " << node.ymin << " " << node.ymax << " " << node.weight << std::endl;
    }
    out.close();
}

void DataExport(const std::map<int, Node>& nodes, const std::string& filename) {
    std::ofstream out;
    out.open(filename);
    std::string line;
    for (auto& node : nodes) {
        out << node.second.xmin << " " << node.second.xmax << " " << node.second.ymin << " " << node.second.ymax << " " << node.second.weight << std::endl;
    }
    out.close();
}

void DataExport(const std::vector<std::pair<double, double>>& contours, const std::string& filename) {
    std::ofstream out;
    out.open(filename);
    std::string line;
    for (auto& contour : contours) {
        out << contour.first << " " << contour.second << std::endl;
    }
    out.close();
}

void DataExport(const std::vector<Node>& nodes, std::vector<double>& weights, const std::string& filename) {
    std::ofstream out;
    out.open(filename);
    std::string line;
    for (int i = 0; i < nodes.size(); i++) {
        out << nodes[i].xmin << " " << nodes[i].xmax << " " << nodes[i].ymin << " " << nodes[i].ymax << " " << weights[i] << std::endl;
    }
    out.close();
}

double BinaryTree::symlog(double x) {
    if (x > 0) return std::log10(x / threshold + 1);
    else return -std::log10(-x / threshold + 1);
}

double BinaryTree::symlog_inv(double x) {
    if (x > 0) return threshold * (pow(10, x) - 1);
    else return -threshold * (pow(10, -x) - 1);
}

double BinaryTree::natural_log(double x) {
    return std::log(x);
}

double BinaryTree::exp(double x) {
    return std::exp(x);
}

double BinaryTree::dec_log(double x) {
    return std::log10(x);
}

double BinaryTree::pow10(double x) {
    return std::pow(10, x);
}

double BinaryTree::lin(double x) {
    return x;
}

void BinaryTree::output(std::ostream& os, int i, int dep, double x, double dx) {
    if (nodes[i].child[0] >= 0) {
        dx /= 2;
        output(os, nodes[i].child[0], dep + 1, x, dx);
        output(os, nodes[i].child[1], dep + 1, x + dx, dx);
    }
    else {
        os << x + dx / 2 << " " << nodes[i].weight << std::endl;
    }
}

void Quadtree::output(std::ostream& os, int i, int dep, double x, double y, double dx, double dy) {
    if (nodes[i].child[0] >= 0) {
        dx /= 2; dy /= 2;
        output(os, nodes[i].child[0], dep + 1, x, y, dx, dy);
        output(os, nodes[i].child[1], dep + 1, x + dx, y, dx, dy);
        output(os, nodes[i].child[2], dep + 1, x + dx, y + dy, dx, dy);
        output(os, nodes[i].child[3], dep + 1, x, y + dy, dx, dy);
    }
    else {
        os << x + dx / 2 << " " << y + dy / 2 << " " << nodes[i].weight << std::endl;
    }
}

void Octree::output(std::ostream& os, int i, int dep, double x, double y, double z, double dx, double dy, double dz) {
    if (nodes[i].child[0] >= 0) {
        dx /= 2; dy /= 2; dz /= 2;
        output(os, nodes[i].child[0], dep + 1, x, y, z, dx, dy, dz);
        output(os, nodes[i].child[1], dep + 1, x + dx, y, z, dx, dy, dz);
        output(os, nodes[i].child[2], dep + 1, x + dx, y + dy, z, dx, dy, dz);
        output(os, nodes[i].child[3], dep + 1, x, y + dy, z, dx, dy, dz);
        output(os, nodes[i].child[4], dep + 1, x, y, z + dz, dx, dy, dz);
        output(os, nodes[i].child[5], dep + 1, x + dx, y, z + dz, dx, dy, dz);
        output(os, nodes[i].child[6], dep + 1, x + dx, y + dy, z + dz, dx, dy, dz);
        output(os, nodes[i].child[7], dep + 1, x, y + dy, z + dz, dx, dy, dz);
    }
    else {
        os << x + dx / 2 << " " << y + dy / 2 << " " << z + dz / 2 << " " << nodes[i].weight << std::endl;
    }
}

std::ostream& operator<<(std::ostream& os, BinaryTree& qt) {
    qt.output(os, 0, 0, qt.xmin, qt.xmax - qt.xmin);
    return os;
}

std::ostream& operator<<(std::ostream& os, Quadtree& qt) {
    qt.output(os, 0, 0, qt.xmin, qt.ymin, qt.xmax - qt.xmin, qt.ymax - qt.ymin);
    return os;
}

void BinaryTree::subdivide(int parent_index) {
    if (nodes.size() > parent_index) {
        Node1D n;
        n.parent = parent_index;
        n.depth = nodes[parent_index].depth + 1;
        for (int i = 0; i < 2; i++) {
            nodes[parent_index].child[i] = nodes.size();
            nodes.push_back(n);
        }
    }
}

void Quadtree::subdivide(int parent_index) {
    if (nodes.size() > parent_index) {
        Node2D n;
        n.parent = parent_index;
        n.depth = nodes[parent_index].depth + 1;
        for (int i = 0; i < children_per_node; i++) {
            nodes[parent_index].child[i] = nodes.size();
            nodes.push_back(n);
        }
    }
}

void Octree::subdivide(int parent_index) {
    if (nodes.size() > parent_index) {
        Node3D n;
        n.parent = parent_index;
        n.depth = nodes[parent_index].depth + 1;
        for (int i = 0; i < children_per_node; i++) {
            nodes[parent_index].child[i] = nodes.size();
            nodes.push_back(n);
        }
    }
}

void BinaryTree::set_scales(std::string scale_X) {
    if ((scale_X == "lin") || (scale_X == "linear")) {
        scale = [](double x) {return 1; };
    }
    else if ((scale_X == "log") || (scale_X == "ln")) {
        scale = [](double x) {return std::exp(-x); };
    }
    else if ((scale_X == "log10") || (scale_X == "lg")) {
        scale = [](double x) {return std::pow(10, -x); };
    }
    else if ((scale_X == "symlog")) {
        scale = [this](double x) {return 1 / (threshold * (std::pow(10, std::abs(x)) - 1)); };
    }
}

void Quadtree::set_scales(std::string scale_X, std::string scale_Y) {
    if ((scale_X == "lin") || (scale_X == "linear")) {
        scale_x = [](double x) {return 1; };
    }
    else if ((scale_X == "log") || (scale_X == "ln")) {
        scale_x = [](double x) {return std::exp(-x); };
    }
    else if ((scale_X == "log10") || (scale_X == "lg")) {
        scale_x = [](double x) {return std::pow(10, -x); };
    }
    else if ((scale_X == "symlog")) {
        scale_x = [this](double x) {return 1 / (threshold * (std::pow(10, std::abs(x)) - 1)); };
    }

    if ((scale_Y == "lin") || (scale_Y == "linear")) {
        scale_y = [](double x) {return 1; };
    }
    else if ((scale_Y == "log") || (scale_Y == "ln")) {
        scale_y = [](double x) {return std::exp(-x); };
    }
    else if ((scale_Y == "log10") || (scale_Y == "lg")) {
        scale_y = [](double x) {return std::pow(10, -x); };
    }
    else if ((scale_Y == "symlog")) {
        scale_y = [this](double x) {return 1 / (threshold * (std::pow(10, std::abs(x)) - 1)); };
    }
}

void Quadtree::set_scales(std::function<double(double)> scale_X, std::function<double(double)> scale_Y) {
    scale_x = scale_X;
    scale_y = scale_Y;
}

void Octree::set_scales(std::string scale_X, std::string scale_Y, std::string scale_Z) {
    if ((scale_X == "lin") || (scale_X == "linear")) {
        scale_x = [](double x) {return lin(x); };
    }
    else if ((scale_X == "log") || (scale_X == "ln")) {
        scale_x = [](double x) {return std::exp(-x); };
    }
    else if ((scale_X == "log10") || (scale_X == "lg")) {
        scale_x = [](double x) {return std::pow(10, -x); };
    }
    else if ((scale_X == "symlog")) {
        scale_x = [this](double x) {return 1 / (threshold * (std::pow(10, std::abs(x)) - 1)); };
    }

    if ((scale_Y == "lin") || (scale_Y == "linear")) {
        scale_y = [](double x) {return 1; };
    }
    else if ((scale_Y == "log") || (scale_Y == "ln")) {
        scale_y = [](double x) {return std::exp(-x); };
    }
    else if ((scale_Y == "log10") || (scale_Y == "lg")) {
        scale_y = [](double x) {return std::pow(10, -x); };
    }
    else if ((scale_Y == "symlog")) {
        scale_y = [this](double x) {return 1 / (threshold * (std::pow(10, std::abs(x)) - 1)); };
    }

    if ((scale_Z == "lin") || (scale_Z == "linear")) {
        scale_z = [](double x) {return 1; };
    }
    else if ((scale_Y == "log") || (scale_Y == "ln")) {
        scale_z = [](double x) {return std::exp(-x); };
    }
    else if ((scale_Y == "log10") || (scale_Y == "lg")) {
        scale_z = [](double x) {return std::pow(10, -x); };
    }
    else if ((scale_Y == "symlog")) {
        scale_z = [this](double x) {return 1 / (threshold * (std::pow(10, std::abs(x)) - 1)); };
    }
}

void BinaryTree::set_bounds(double XMIN, double XMAX) {
    xmin = XMIN;
    xmax = XMAX;
}

void Quadtree::set_bounds(double XMIN, double XMAX, double YMIN, double YMAX) {
    xmin = XMIN;
    xmax = XMAX;
    ymin = YMIN;
    ymax = YMAX;
}

void Octree::set_bounds(double XMIN, double XMAX, double YMIN, double YMAX, double ZMIN, double ZMAX) {
    xmin = XMIN;
    xmax = XMAX;
    ymin = YMIN;
    ymax = YMAX;
    zmin = ZMIN;
    zmax = ZMAX;
}

void BinaryTree::initialize(double XMIN, double XMAX, std::string scale_X) {
    set_scales(scale_X);
    set_bounds(XMIN, XMAX);
    nodes.clear();
    Node1D n;
    n.parent = -1;
    n.depth = 1;
    nodes.push_back(n);
}

void Quadtree::initialize(double XMIN, double XMAX, double YMIN, double YMAX, std::string scale_X, std::string scale_Y) {
    set_scales(scale_X, scale_Y);
    set_bounds(XMIN, XMAX, YMIN, YMAX);
    nodes.clear();
    Node2D n;
    n.parent = -1;
    n.depth = 1;
    nodes.push_back(n);
}

void Octree::initialize(double XMIN, double XMAX, double YMIN, double YMAX, double ZMIN, double ZMAX, std::string scale_X, std::string scale_Y, std::string scale_Z) {
    set_scales(scale_X, scale_Y, scale_Z);
    set_bounds(XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX);
    nodes.clear();
    Node3D n;
    n.parent = -1;
    n.depth = 1;
    nodes.push_back(n);
}

void BinaryTree::initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::string scale_X) {
    std::string backspace = "";
    for (int i = 0; i < 100; i++) backspace = backspace + "\b";
    set_scales(scale_X);
    std::vector<double> x;
    double XMIN = 1e10;
    double YMIN = 1e10;
    double XMAX = -1e10;
    double YMAX = -1e10;
    int N = filename_x.size();
    for (int i = 0; i < N; i++) {
        //std::cout << backspace << "Reading " << filename_x[i] << "...             ";
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << backspace << "Filename " << filename_x[i] << " not found!                      " << std::endl;
            continue;
        }
        bool adjusted = false;
        for (int j = 0; j < x.size(); j++) {
            if (XMIN > x[j]) {
                XMIN = x[j];
                adjusted = true;
            }
            if (XMAX < x[j]) {
                XMAX = x[j];
                adjusted = true;
            }
        }
        if (adjusted) std::cout << backspace << "Bounds adjusted: (" << XMIN << ", " << XMAX << ")              ";
    }
    std::cout << std::endl;
    set_bounds(XMIN, XMAX);
    nodes.clear();
    Node1D n;
    n.parent = -1;
    n.depth = 1;
    nodes.push_back(n);
}

void Quadtree::initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::string scale_X, std::string scale_Y) {
    std::string backspace = "";
    for (int i = 0; i < 100; i++) backspace = backspace + "\b";
    set_scales(scale_X, scale_Y);
    std::vector<double> x, y;
    double XMIN = 1e10;
    double YMIN = 1e10;
    double XMAX = -1e10;
    double YMAX = -1e10;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    for (int i = 0; i < N; i++) {
        std::cout << backspace << "Reading " << filename_x[i] << "...             ";
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        bool adjusted = false;
        for (int j = 0; j < x.size(); j++) {
            if (XMIN > x[j]) {
                XMIN = x[j];
                adjusted = true;
            }
            if (YMIN > y[j]) {
                YMIN = y[j];
                adjusted = true;
            }
            if (XMAX < x[j]) {
                XMAX = x[j];
                adjusted = true;
            }
            if (YMAX < y[j]) {
                YMAX = y[j];
                adjusted = true;
            }
        }
        if (adjusted) std::cout << backspace << "Bounds adjusted: (" << XMIN << ", " << XMAX << "), (" << YMIN << ", " << YMAX << ")                 ";
    }
    std::cout << std::endl;
    set_bounds(XMIN, XMAX, YMIN, YMAX);
    nodes.clear();
    Node2D n;
    n.parent = -1;
    n.depth = 1;
    nodes.push_back(n);
}

void Octree::initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::string scale_X, std::string scale_Y, std::string scale_Z) {
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    set_scales(scale_X, scale_Y, scale_Z);
    std::vector<double> x, y, z;
    double XMIN = 1e10;
    double YMIN = 1e10;
    double ZMIN = 1e10;
    double XMAX = -1e10;
    double YMAX = -1e10;
    double ZMAX = -1e10;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_z.size()) N = filename_z.size();
    for (int i = 0; i < N; i++) {
        std::cout << backspace << "Reading " << filename_x[i] << "...             ";
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_z[i], z) == -1) {
            std::cout << "Filename " << filename_z[i] << " not found!" << std::endl;
            continue;
        }
        bool adjusted = false;
        for (int j = 0; j < x.size(); j++) {
            if (XMIN > x[j]) {
                XMIN = x[j];
                adjusted = true;
            }
            if (YMIN > y[j]) {
                YMIN = y[j];
                adjusted = true;
            }
            if (ZMIN > z[j]) {
                ZMIN = z[j];
                adjusted = true;
            }
            if (XMAX < x[j]) {
                XMAX = x[j];
                adjusted = true;
            }
            if (YMAX < y[j]) {
                YMAX = y[j];
                adjusted = true;
            }
            if (ZMAX < z[j]) {
                ZMAX = z[j];
                adjusted = true;
            }
        }
        if (adjusted) std::cout << backspace << "Bounds adjusted: (" << XMIN << ", " << XMAX << "), (" << YMIN << ", " << YMAX << "), (" << ZMIN << ", " << ZMAX << ")                 ";
    }
    std::cout << std::endl;
    set_bounds(XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX);
    nodes.clear();
    Node3D n;
    n.parent = -1;
    n.depth = 1;
    nodes.push_back(n);
}

void BinaryTree::add_point(double x, double weight) {
    double x_scaled;
    x_scaled = (x - xmin) / (xmax - xmin);

    int candidate_index = 0;
    if ((x_scaled >= 0) && (x_scaled <= 1)) {
        while (nodes[candidate_index].child[0] >= 0) { // checking in CCW order
            if (x_scaled < 0.5) { // down left
                candidate_index = nodes[candidate_index].child[0];
                x_scaled = 2 * x_scaled;
            }
            else if (x_scaled >= 0.5) { // down right
                candidate_index = nodes[candidate_index].child[1];
                x_scaled = 2 * (x_scaled - 0.5);
            }
        }
        nodes[candidate_index].weight += weight;
    }
}

void Quadtree::add_point(double x, double y, double weight) {
    double x_scaled, y_scaled;
    x_scaled = (x - xmin) / (xmax - xmin);
    y_scaled = (y - ymin) / (ymax - ymin);

    int candidate_index = 0;
    if ((x_scaled >= 0) && (x_scaled <= 1) && (y_scaled >= 0) && (y_scaled <= 1)) {
        while (nodes[candidate_index].child[0] >= 0) { // checking in CCW order
            if ((x_scaled < 0.5) && (y_scaled < 0.5)) { // down left
                candidate_index = nodes[candidate_index].child[0];
                x_scaled = 2 * x_scaled;
                y_scaled = 2 * y_scaled;
            }
            else if ((x_scaled > 0.5) && (y_scaled < 0.5)) { // down right
                candidate_index = nodes[candidate_index].child[1];
                x_scaled = 2 * (x_scaled - 0.5);
                y_scaled = 2 * y_scaled;
            }
            else if ((x_scaled > 0.5) && (y_scaled > 0.5)) { // up right
                candidate_index = nodes[candidate_index].child[2];
                x_scaled = 2 * (x_scaled - 0.5);
                y_scaled = 2 * (y_scaled - 0.5);
            }
            else { // up left
                candidate_index = nodes[candidate_index].child[3];
                x_scaled = 2 * x_scaled;
                y_scaled = 2 * (y_scaled - 0.5);
            }
        }
        nodes[candidate_index].weight += weight;
    }
}

void Octree::add_point(double x, double y, double z, double weight) {
    double x_scaled, y_scaled, z_scaled;
    x_scaled = (x - xmin) / (xmax - xmin);
    y_scaled = (y - ymin) / (ymax - ymin);
    z_scaled = (z - zmin) / (zmax - zmin);

    int candidate_index = 0;
    if ((x_scaled >= 0) && (x_scaled <= 1) && (y_scaled >= 0) && (y_scaled <= 1) && (z_scaled >= 0) && (z_scaled <= 1)) {
        while (nodes[candidate_index].child[0] >= 0) { // checking in CCW order
            if ((x_scaled < 0.5) && (y_scaled < 0.5) && (z_scaled < 0.5)) { // down left down
                candidate_index = nodes[candidate_index].child[0];
                x_scaled = 2 * x_scaled;
                y_scaled = 2 * y_scaled;
                z_scaled = 2 * z_scaled;
            }
            else if ((x_scaled >= 0.5) && (y_scaled < 0.5) && (z_scaled < 0.5)) { // down right down
                candidate_index = nodes[candidate_index].child[1];
                x_scaled = 2 * (x_scaled - 0.5);
                y_scaled = 2 * y_scaled;
                z_scaled = 2 * z_scaled;
            }
            else if ((x_scaled >= 0.5) && (y_scaled >= 0.5) && (z_scaled < 0.5)) { // up right down
                candidate_index = nodes[candidate_index].child[2];
                x_scaled = 2 * (x_scaled - 0.5);
                y_scaled = 2 * (y_scaled - 0.5);
                z_scaled = 2 * z_scaled;
            }
            else if ((x_scaled < 0.5) && (y_scaled >= 0.5) && (z_scaled < 0.5)) { // up left down
                candidate_index = nodes[candidate_index].child[3];
                x_scaled = 2 * x_scaled;
                y_scaled = 2 * (y_scaled - 0.5);
                z_scaled = 2 * z_scaled;
            }
            else if ((x_scaled < 0.5) && (y_scaled < 0.5) && (z_scaled >= 0.5)) { // down left up
                candidate_index = nodes[candidate_index].child[4];
                x_scaled = 2 * x_scaled;
                y_scaled = 2 * y_scaled;
                z_scaled = 2 * (z_scaled - 0.5);
            }
            else if ((x_scaled >= 0.5) && (y_scaled < 0.5) && (z_scaled >= 0.5)) { // down right up
                candidate_index = nodes[candidate_index].child[5];
                x_scaled = 2 * (x_scaled - 0.5);
                y_scaled = 2 * y_scaled;
                z_scaled = 2 * (z_scaled - 0.5);
            }
            else if ((x_scaled >= 0.5) && (y_scaled >= 0.5) && (z_scaled >= 0.5)) { // up right up
                candidate_index = nodes[candidate_index].child[6];
                x_scaled = 2 * (x_scaled - 0.5);
                y_scaled = 2 * (y_scaled - 0.5);
                z_scaled = 2 * (z_scaled - 0.5);
            }
            else if ((x_scaled < 0.5) && (y_scaled >= 0.5) && (z_scaled >= 0.5)) { // up left up
                candidate_index = nodes[candidate_index].child[7];
                x_scaled = 2 * x_scaled;
                y_scaled = 2 * (y_scaled - 0.5);
                z_scaled = 2 * (z_scaled - 0.5);
            }
        }
        nodes[candidate_index].weight += weight;
    }
}

void BinaryTree::add_point(double x) {
    add_point(x, 1.0);
}

void Quadtree::add_point(double x, double y) {
    add_point(x, y, 1.0);
}

void Octree::add_point(double x, double y, double z) {
    add_point(x, y, z, 1.0);
}

void BinaryTree::uniform_divide(int DEPTH) {
    nodes.clear();
    Node1D n;
    n.parent = -1;
    n.depth = 1;
    nodes.push_back(n);
    std::queue<int> to_divide;
    to_divide.push(0);
    while (to_divide.size() > 0) {
        if (DEPTH > nodes[to_divide.front()].depth) {
            int n = nodes.size();
            subdivide(to_divide.front());
            to_divide.push(n); // 0
            to_divide.push(n + 1); // 1
        }
        to_divide.pop();
    }
}

void Quadtree::uniform_divide(int DEPTH) {
    nodes.clear();
    Node2D n;
    n.parent = -1;
    n.depth = 1;
    nodes.push_back(n);
    std::queue<int> to_divide;
    to_divide.push(0);
    while (to_divide.size() > 0) {
        if (DEPTH > nodes[to_divide.front()].depth) {
            int n = nodes.size();
            subdivide(to_divide.front());
            to_divide.push(n); // 0 0
            to_divide.push(n + 1); // 1 0
            to_divide.push(n + 2); // 1 1
            to_divide.push(n + 3); // 0 1
        }
        to_divide.pop();
    }
    //std::cout << *this;
}

void Octree::uniform_divide(int DEPTH) {
    nodes.clear();
    Node3D n;
    n.parent = -1;
    n.depth = 1;
    nodes.push_back(n);
    std::queue<int> to_divide;
    to_divide.push(0);
    while (to_divide.size() > 0) {
        if (DEPTH > nodes[to_divide.front()].depth) {
            int n = nodes.size();
            subdivide(to_divide.front());
            to_divide.push(n); // 0 0 0
            to_divide.push(n + 1); // 1 0 0
            to_divide.push(n + 2); // 1 1 0
            to_divide.push(n + 3); // 0 1 0
            to_divide.push(n + 4); // 0 0 1
            to_divide.push(n + 5); // 1 0 1
            to_divide.push(n + 6); // 1 1 1
            to_divide.push(n + 7); // 0 1 1
        }
        to_divide.pop();
    }
    //std::cout << *this;
}

void BinaryTree::load_points(const std::vector<double> &pts) {
    for (auto &val: pts) add_point(val);
}

void BinaryTree::load_points(const std::vector<double> &pts, const std::vector<double> &wt) {
    for (int i = 0; i < pts.size(); i++) add_point(pts[i], wt[i]);
}

void BinaryTree::load_points(const std::vector<double> &x, const std::vector<double> &y, std::function<double(double, double)> &fun) {
    for (int i = 0; i < x.size(); i++) add_point(fun(x[i], y[i]));
}

void BinaryTree::load_points(const std::vector<double> &x, const std::vector<double> &y, std::function<double(double, double)> &fun, std::function<double(double, double)> &fun_wt) {
    for (int i = 0; i < x.size(); i++) add_point(fun(x[i], y[i]), fun_wt(x[i], y[i]));
}

void BinaryTree::load_points(double x[], double y[], const std::function<double(double, double)> &fun, const std::function<double(double, double)> &fun_wt, int size) {
    for (int i = 0; i < size; i++) add_point(fun(x[i], y[i]), fun_wt(x[i], y[i]));
}

void Quadtree::load_points(double x[], double y[], const std::function<double(double, double)> &fun_x, const std::function<double(double, double)> &fun_y, const std::function<double(double, double)> &fun_wt, int size) {
    for (int i = 0; i < size; i++) add_point(fun_x(x[i], y[i]), fun_y(x[i], y[i]), fun_wt(x[i], y[i]));
}

void BinaryTree::load_points(std::string path, std::string prefix, std::vector<std::string> filename_x) {
    std::vector<double> x;
    int N = filename_x.size();
    std::cout << "Adding points..." << std::endl;
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        for (int j = 0; j < x.size(); j++) {
            add_point(x[j]);
        }
        progress = (double)(i + 1) / N;
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << 100 * progress << "     ";
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
}

void Quadtree::load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y) {
    std::vector<double> x, y;
    std::string del = "";
    for (int i = 0; i < 300; i++) del += "\b";
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    std::cout << "Adding points..." << std::endl;
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        for (int j = 0; j < x.size(); j++) {
            add_point(x[j], y[j]);
        }
        progress = (double)(i + 1) / N;
        std::cout << del << 100 * progress << "          ";
    }
    std::cout << del << "                              " << del;
}

void Quadtree::load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::function<double(double, double, double)> fun1, std::function<double(double, double, double)> fun2, std::function<bool(double, double, double)> cond) {
    std::vector<double> x, y, z;
    std::string del = "";
    for (int i = 0; i < 300; i++) del += "\b";
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_z.size()) N = filename_z.size();
    std::cout << "Adding points..." << std::endl;
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_z[i], z) == -1) {
            std::cout << "Filename " << filename_z[i] << " not found!" << std::endl;
            continue;
        }
        for (int j = 0; j < x.size(); j++) {
            if (cond(x[j], y[j], z[j])) {
                double num1 = fun1(x[j], y[j], z[j]);
                double num2 = fun2(x[j], y[j], z[j]);
                add_point(num1, num2);
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << del << 100 * progress << "          ";
    }
    std::cout << del << "                              " << del;
}

void Quadtree::load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::function<double(double, double, double)> fun1, std::function<double(double, double, double)> fun2, std::function<bool(double, double, double, double, double)> cond, double cond_param1, double cond_param2) {
    std::vector<double> x, y, z;
    std::string del = "";
    for (int i = 0; i < 300; i++) del += "\b";
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_z.size()) N = filename_z.size();
    std::cout << "Adding points..." << std::endl;
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_z[i], z) == -1) {
            std::cout << "Filename " << filename_z[i] << " not found!" << std::endl;
            continue;
        }
        for (int j = 0; j < x.size(); j++) {
            if (cond(x[j], y[j], z[j], cond_param1, cond_param2)) {
                double num1 = fun1(x[j], y[j], z[j]);
                double num2 = fun2(x[j], y[j], z[j]);
                add_point(num1, num2);
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << del << 100 * progress << "          ";
    }
    std::cout << del << "                              " << del;
}

void load_points(std::vector<Quadtree>& q, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::function<double(double, double, double)> fun1, std::function<double(double, double, double)> fun2, std::function<bool(double, double, double, double, double)> cond, std::vector<std::pair<double, double>> conditions) {
    std::vector<double> x, y, z;
    std::string del = "";
    for (int i = 0; i < 300; i++) del += "\b";
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_z.size()) N = filename_z.size();
    std::cout << "Adding points..." << std::endl;
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_z[i], z) == -1) {
            std::cout << "Filename " << filename_z[i] << " not found!" << std::endl;
            continue;
        }
        for (int j = 0; j < x.size(); j++) {
            for (int k = 0; k < q.size(); k++) {
                if (cond(x[j], y[j], z[j], conditions[k].first, conditions[k].second)) {
                    double num1 = fun1(x[j], y[j], z[j]);
                    double num2 = fun2(x[j], y[j], z[j]);
                    q[k].add_point(num1, num2);
                }
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << del << 100 * progress << "          ";
    }
    std::cout << del << "                              " << del;
}

void Octree::load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z) {
    std::vector<double> x, y, z;
    std::string del = "";
    for (int i = 0; i < 300; i++) del += "\b";
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_z.size()) N = filename_z.size();
    std::cout << "Adding points..." << std::endl;
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_z[i], z) == -1) {
            std::cout << "Filename " << filename_z[i] << " not found!" << std::endl;
            continue;
        }
        for (int j = 0; j < x.size(); j++) {
            add_point(x[j], y[j], z[j]);
        }
        progress = (double)(i + 1) / N;
        std::cout << del << 100 * progress << "          ";
    }
    std::cout << del << "                              " << del;
}

void Octree::load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::function<double(double, double, double)> fun1, std::function<double(double, double, double)> fun2, std::function<double(double, double, double)> fun3, std::function<bool(double, double, double)> cond) {
    std::vector<double> x, y, z;
    std::string del = "";
    for (int i = 0; i < 300; i++) del += "\b";
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_z.size()) N = filename_z.size();
    std::cout << "Adding points..." << std::endl;
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_z[i], z) == -1) {
            std::cout << "Filename " << filename_z[i] << " not found!" << std::endl;
            continue;
        }
        for (int j = 0; j < x.size(); j++) {
            if (cond(x[j], y[j], z[j])) {
                double num1 = fun1(x[j], y[j], z[j]);
                double num2 = fun2(x[j], y[j], z[j]);
                double num3 = fun3(x[j], y[j], z[j]);
                add_point(num1, num2, num3);
            }
            else {
                std::cout << del << "Point (" << x[j] << ", " << y[j] << ", " << z[j] << ") removed!                " << std::endl;
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << del << 100 * progress << "          ";
    }
    std::cout << del << "                              " << del;
}

void BinaryTree::load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_weight) {
    std::vector<double> x, w;
    int N = filename_x.size();
    if (N > filename_weight.size()) N = filename_weight.size();
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        int read_weights = ::read_from_binary(path, prefix, filename_weight[i], w);
        if (read_weights == -1) {
            std::cout << "Filename " << filename_weight[i] << " not found! Using the default weights." << std::endl;
            for (int j = 0; j < x.size(); j++) {
                add_point(x[i]);
            }
        }
        else {
            for (int j = 0; j < x.size(); j++) {
                add_point(x[j], w[j]);
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << backspace << 100 * progress << "                                          ";
    }
    std::cout << backspace;
}

void Quadtree::load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_weight) {
    std::vector<double> x, y, w;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_weight.size()) N = filename_weight.size();
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        int read_weights = ::read_from_binary(path, prefix, filename_weight[i], w);
        if (read_weights == -1) {
            std::cout << "Filename " << filename_weight[i] << " not found! Using the default weights." << std::endl;
            for (int j = 0; j < x.size(); j++) {
                add_point(x[i], y[i]);
            }
        }
        else {
            for (int j = 0; j < x.size(); j++) {
                add_point(x[j], y[j], w[j]);
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << backspace << 100 * progress << "                                             ";
    }
    std::cout << backspace;
}

void Octree::load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_z, std::vector<std::string> filename_weight) {
    std::vector<double> x, y, z, w;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_z.size()) N = filename_z.size();
    if (N > filename_weight.size()) N = filename_weight.size();
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_z[i], z) == -1) {
            std::cout << "Filename " << filename_z[i] << " not found!" << std::endl;
            continue;
        }
        int read_weights = ::read_from_binary(path, prefix, filename_weight[i], w);
        if (read_weights == -1) {
            std::cout << "Filename " << filename_weight[i] << " not found! Using the default weights." << std::endl;
            for (int j = 0; j < x.size(); j++) {
                add_point(x[i], y[i], z[j]);
            }
        }
        else {
            for (int j = 0; j < x.size(); j++) {
                add_point(x[j], y[j], z[j], w[j]);
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << backspace << 100 * progress << "                                             ";
    }
    std::cout << backspace;
}

double BinaryTree::count_total(int i) {
    double res = 0;
    std::queue<int> to_check;
    to_check.push(i);
    while (!to_check.empty()) {
        if (nodes[to_check.front()].child[0] >= 0) {
            to_check.push(nodes[to_check.front()].child[0]);
            to_check.push(nodes[to_check.front()].child[1]);
        }
        else {
            res += nodes[to_check.front()].weight;
        }
        to_check.pop();
    }
    return res;
}

double Quadtree::count_total(int i) {
    double res = 0;
    std::queue<int> to_check;
    to_check.push(i);
    while (!to_check.empty()) {
        if (nodes[to_check.front()].child[0] >= 0) {
            for (int i = 0; i < children_per_node; i++) to_check.push(nodes[to_check.front()].child[i]);
        }
        else {
            res += nodes[to_check.front()].weight;
        }
        to_check.pop();
    }
    return res;
}

double Octree::count_total(int i) {
    double res = 0;
    std::queue<int> to_check;
    to_check.push(i);
    while (!to_check.empty()) {
        if (nodes[to_check.front()].child[0] >= 0) {
            for (int i = 0; i < children_per_node; i++) to_check.push(nodes[to_check.front()].child[i]);
        }
        else {
            res += nodes[to_check.front()].weight;
        }
        to_check.pop();
    }
    return res;
}

double BinaryTree::count_total() {
    return count_total(0);
}

double Quadtree::count_total() {
    return count_total(0);
}

double Octree::count_total() {
    return count_total(0);
}

void BinaryTree::merge(const double& merge_weight, int index) {
    Node1D& current_node = nodes[index];
    if (current_node.child[0] >= 0) { // there's nodes below
        for (int i = 0; i < 2; i++) {
            if (nodes[current_node.child[i]].child[0] >= 0) { // if the nodes below contain more nodes below, continue deeper down the tree
                merge(merge_weight, nodes[index].child[i]);
            }
        }
        bool good_to_merge = true;
        for (int i = 0; i < 2; i++) { // is the current node good to merge? Does it only have four children with no deeper structure?
            if (nodes[current_node.child[i]].child[0] >= 0) {
                good_to_merge = false;
                break;
            }
        }
        double weight = nodes[current_node.child[0]].weight + nodes[current_node.child[1]].weight;
        good_to_merge = good_to_merge && (weight < merge_weight);
        if (good_to_merge) { // no deeper structure AND weight is smaller than the threshold
            for (int i = 0; i < 2; i++) {
                nodes[current_node.child[i]].parent = -1; // severe the child from the parent
                current_node.child[i] = -1; // severe the parent from the child
            }
            current_node.weight = weight;
        }
    }
}

void Quadtree::merge(const double& merge_weight, int index) {
    Node2D& current_node = nodes[index];
    if (current_node.child[0] >= 0) { // there's nodes below
        for (int i = 0; i < children_per_node; i++) {
            if (nodes[current_node.child[i]].child[0] >= 0) { // if the nodes below contain more nodes below, continue deeper down the tree
                merge(merge_weight, nodes[index].child[i]);
            }
        }
        bool good_to_merge = true;
        for (int i = 0; i < children_per_node; i++) { // is the current node good to merge? Does it only have four children with no deeper structure?
            if (nodes[current_node.child[i]].child[0] >= 0) {
                good_to_merge = false;
                break;
            }
        }
        double weight = 0;
        for (int i = 0; i < children_per_node; i++) weight += nodes[current_node.child[i]].weight;
        good_to_merge = good_to_merge && (weight < merge_weight);
        if (good_to_merge) { // no deeper structure AND weight is smaller than the threshold
            for (int i = 0; i < children_per_node; i++) {
                nodes[current_node.child[i]].parent = -1; // severe the child from the parent
                current_node.child[i] = -1; // severe the parent from the child
            }
            current_node.weight = weight;
        }
    }
}

void Octree::merge(const double& merge_weight, int index) {
    Node3D& current_node = nodes[index];
    if (current_node.child[0] >= 0) { // there's nodes below
        for (int i = 0; i < children_per_node; i++) {
            if (nodes[current_node.child[i]].child[0] >= 0) { // if the nodes below contain more nodes below, continue deeper down the tree
                merge(merge_weight, nodes[index].child[i]);
            }
        }
        bool good_to_merge = true;
        for (int i = 0; i < children_per_node; i++) { // is the current node good to merge? Does it only have children but no grandchildren?
            if (nodes[current_node.child[i]].child[0] >= 0) {
                good_to_merge = false;
                break;
            }
        }
        double weight = 0;
        for (int i = 0; i < children_per_node; i++) weight += nodes[current_node.child[i]].weight;
        good_to_merge = good_to_merge && (weight < merge_weight);
        if (good_to_merge) { // no deeper structure AND weight is smaller than the threshold
            for (int i = 0; i < children_per_node; i++) {
                nodes[current_node.child[i]].parent = -1; // severe the child from the parent
                current_node.child[i] = -1; // severe the parent from the child
            }
            current_node.weight = weight;
        }
    }
}

void Quadtree::merge(const double& merge_weight, int index, double min_depth) {
    Node2D& current_node = nodes[index];
    if (current_node.child[0] >= 0) { // there's nodes below
        for (int i = 0; i < children_per_node; i++) {
            if (nodes[current_node.child[i]].child[0] >= 0) { // if the nodes below contain more nodes below, continue deeper down the tree
                merge(merge_weight, nodes[index].child[i], min_depth);
            }
        }
        bool good_to_merge = true;
        for (int i = 0; i < children_per_node; i++) { // is the current node good to merge? Does it only have children but no grandchildren?
            if (nodes[current_node.child[i]].child[0] >= 0) {
                good_to_merge = false;
                break;
            }
        }
        double weight = 0;
        for (int i = 0; i < children_per_node; i++) weight += nodes[current_node.child[i]].weight;
        good_to_merge = good_to_merge && (weight < merge_weight) && (current_node.depth >= min_depth);
        if (good_to_merge) { // no deeper structure AND weight is smaller than the threshold
            for (int i = 0; i < children_per_node; i++) {
                nodes[current_node.child[i]].parent = -1; // severe the child from the parent
                current_node.child[i] = -1; // severe the parent from the child
            }
            current_node.weight = weight;
        }
    }
}

void Octree::merge(const double& merge_weight, int index, double min_depth) {
    Node3D& current_node = nodes[index];
    if (current_node.child[0] >= 0) { // there's nodes below
        for (int i = 0; i < children_per_node; i++) {
            if (nodes[current_node.child[i]].child[0] >= 0) { // if the nodes below contain more nodes below, continue deeper down the tree
                merge(merge_weight, nodes[index].child[i], min_depth);
            }
        }
        bool good_to_merge = true;
        for (int i = 0; i < children_per_node; i++) { // is the current node good to merge? Does it only have 8 children but no grandchildren?
            if (nodes[current_node.child[i]].child[0] >= 0) {
                good_to_merge = false;
                break;
            }
        }
        double weight = 0;
        for (int i = 0; i < children_per_node; i++) weight += nodes[current_node.child[i]].weight;
        good_to_merge = good_to_merge && (weight < merge_weight) && (current_node.depth >= min_depth);
        if (good_to_merge) { // no deeper structure AND weight is smaller than the threshold
            for (int i = 0; i < children_per_node; i++) {
                nodes[current_node.child[i]].parent = -1; // severe the child from the parent
                current_node.child[i] = -1; // severe the parent from the child
            }
            current_node.weight = weight;
        }
    }
}

void BinaryTree::cleanup() {
    int index = nodes.size() - 1;
    while (index > 0) {
        if (nodes[index].parent == -1) {
            if (index < nodes.size() - 1) {
                nodes[index] = nodes.back();
                for (int i = 0; i < 2; i++) {
                    if (nodes[nodes[index].parent].child[i] == nodes.size() - 1) {
                        nodes[nodes[index].parent].child[i] = index;
                    }
                    if (nodes[index].child[i] > 0) nodes[nodes[index].child[i]].parent = index;
                }
            }
            nodes.pop_back();
        }
        index--;
    }
}

void Quadtree::cleanup() {
    int index = nodes.size() - 1;
    while (index > 0) { // traverse the nodes back to front
        //std::cout << index << std::endl;
        //if (nodes.back().parent == 153) {
        //    std::cout << "Uh oh..." << std::endl;
        //}
        //if (index == 56) {
        //    std::cout << nodes.back().parent << " " << nodes.back().child[0] << " " << nodes.back().child[1] << std::endl;
        //}
        if (nodes[index].parent == -1) { // whoops, this child is severed from the parent - it needs to be cleaned up
            //std::cout << index << std::endl;
            if (index < nodes.size() - 1) { // if the severed child is the back of the nodes, simply pop. no need to move anything
                nodes[index] = nodes.back(); // if not, swap the back node to current position
                for (int i = 0; i < children_per_node; i++) { // and adjust back's parent's corresponding child's index to current index.
                    if (nodes[nodes[index].parent].child[i] == nodes.size() - 1) { // this is because we don't quite know which child corresponds to the index we just moved
                        nodes[nodes[index].parent].child[i] = index;
                        //break;
                    }
                    if (nodes[index].child[i] > 0) nodes[nodes[index].child[i]].parent = index;
                }
            }
            nodes.pop_back();
        }
        index--;
    }
}

void Octree::cleanup() {
    int index = nodes.size() - 1;
    while (index > 0) { // traverse the nodes back to front
        //std::cout << index << std::endl;
        //if (nodes.back().parent == 153) {
        //    std::cout << "Uh oh..." << std::endl;
        //}
        //if (index == 56) {
        //    std::cout << nodes.back().parent << " " << nodes.back().child[0] << " " << nodes.back().child[1] << std::endl;
        //}
        if (nodes[index].parent == -1) { // whoops, this child is severed from the parent - it needs to be cleaned up
            //std::cout << index << std::endl;
            if (index < nodes.size() - 1) { // if the severed child is the back of the nodes, simply pop. no need to move anything
                nodes[index] = nodes.back(); // if not, swap the back node to current position
                for (int i = 0; i < children_per_node; i++) { // and adjust back's parent's corresponding child's index to current index.
                    if (nodes[nodes[index].parent].child[i] == nodes.size() - 1) { // this is because we don't quite know which child corresponds to the index we just moved
                        nodes[nodes[index].parent].child[i] = index;
                        //break;
                    }
                    if (nodes[index].child[i] > 0) nodes[nodes[index].child[i]].parent = index;
                }
            }
            nodes.pop_back();
        }
        index--;
    }
}

void BinaryTree::merge(double threshold_weight) {
    double total_weight = count_total();
    std::cout << "Merging excess nodes..." << std::endl;
    merge(threshold_weight * total_weight, 0); // keep merging nodes with weights smaller than the threshold weight
    std::cout << "Cleaning up..." << std::endl;
    cleanup(); // clean up the leftovers after merged nodes
}

void Quadtree::merge(double threshold_weight) {
    double total_weight = count_total();
    std::cout << "Merging excess nodes..." << std::endl;
    merge(threshold_weight * total_weight, 0); // keep merging nodes with weights smaller than the threshold weight
    std::cout << "Cleaning up..." << std::endl;
    cleanup(); // clean up the leftovers after merged nodes
}

void Octree::merge(double threshold_weight) {
    double total_weight = count_total();
    std::cout << "Merging excess nodes..." << std::endl;
    merge(threshold_weight * total_weight, 0); // keep merging nodes with weights smaller than the threshold weight
    std::cout << "Cleaning up..." << std::endl;
    cleanup(); // clean up the leftovers after merged nodes
}

void Quadtree::merge(double threshold_weight, double min_depth) {
    double total_weight = count_total();
    std::cout << "Merging excess nodes..." << std::endl;
    merge(threshold_weight * total_weight, 0, min_depth); // keep merging nodes with weights smaller than the threshold weight
    std::cout << "Cleaning up..." << std::endl;
    cleanup(); // clean up the leftovers after merged nodes
}

void Octree::merge(double threshold_weight, double min_depth) {
    double total_weight = count_total();
    std::cout << "Merging excess nodes..." << std::endl;
    merge(threshold_weight * total_weight, 0, min_depth); // keep merging nodes with weights smaller than the threshold weight
    std::cout << "Cleaning up..." << std::endl;
    cleanup(); // clean up the leftovers after merged nodes
}

void BinaryTree::weigh(int node_index, double x_scaled, double current_depth) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        nodes[node_index].weight *= scale(0.5 * (XMIN + XMAX));
    }
    else {
        current_depth /= 2;
        weigh(nodes[node_index].child[0], x_scaled, current_depth);
        weigh(nodes[node_index].child[1], x_scaled + current_depth, current_depth);
    }
}

void Quadtree::weigh(int node_index, double x_scaled, double y_scaled, double current_depth) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double YMIN = ymin + (ymax - ymin) * y_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        double YMAX = ymin + (ymax - ymin) * (y_scaled + current_depth);
        nodes[node_index].weight *= scale_x(0.5 * (XMIN + XMAX)) * scale_y(0.5 * (YMIN + YMAX));
    }
    else {
        current_depth /= 2;
        weigh(nodes[node_index].child[0], x_scaled, y_scaled, current_depth);
        weigh(nodes[node_index].child[1], x_scaled + current_depth, y_scaled, current_depth);
        weigh(nodes[node_index].child[2], x_scaled + current_depth, y_scaled + current_depth, current_depth);
        weigh(nodes[node_index].child[3], x_scaled, y_scaled + current_depth, current_depth);
    }
}

void Octree::weigh(int node_index, double x_scaled, double y_scaled, double z_scaled, double current_depth) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double YMIN = ymin + (ymax - ymin) * y_scaled;
        double ZMIN = zmin + (zmax - zmin) * z_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        double YMAX = ymin + (ymax - ymin) * (y_scaled + current_depth);
        double ZMAX = zmin + (zmax - zmin) * (z_scaled + current_depth);
        nodes[node_index].weight *= scale_x(0.5 * (XMIN + XMAX)) * scale_y(0.5 * (YMIN + YMAX)) * scale_z(0.5 * (ZMIN + ZMAX));
    }
    else {
        current_depth /= 2;
        weigh(nodes[node_index].child[0], x_scaled, y_scaled, z_scaled, current_depth);
        weigh(nodes[node_index].child[1], x_scaled + current_depth, y_scaled, z_scaled, current_depth);
        weigh(nodes[node_index].child[2], x_scaled + current_depth, y_scaled + current_depth, z_scaled, current_depth);
        weigh(nodes[node_index].child[3], x_scaled, y_scaled + current_depth, z_scaled, current_depth);
        weigh(nodes[node_index].child[4], x_scaled, y_scaled, z_scaled + current_depth, current_depth);
        weigh(nodes[node_index].child[5], x_scaled + current_depth, y_scaled, z_scaled + current_depth, current_depth);
        weigh(nodes[node_index].child[6], x_scaled + current_depth, y_scaled + current_depth, z_scaled + current_depth, current_depth);
        weigh(nodes[node_index].child[7], x_scaled, y_scaled + current_depth, z_scaled + current_depth, current_depth);
    }
}

void BinaryTree::write_node(std::ofstream& output, int node_index, const double& total_weight) {
    if (nodes[node_index].child[0] < 0) {
        output << node_index << " " << nodes[node_index].weight / total_weight << std::endl;
    }
    else {
        output << node_index;
        for (int i = 0; i < 2; i++) output << " " << nodes[node_index].child[i];
        output << std::endl;
        for (int i = 0; i < 2; i++) write_node(output, nodes[node_index].child[i], total_weight);
    }
}

void Quadtree::write_node(std::ofstream& output, int node_index, const double& total_weight) {
    if (nodes[node_index].child[0] < 0) {
        output << node_index << " " << nodes[node_index].weight / total_weight << std::endl;
    }
    else {
        output << node_index;
        for (int i = 0; i < children_per_node; i++) output << " " << nodes[node_index].child[i];
        output << std::endl;
        for (int i = 0; i < children_per_node; i++) write_node(output, nodes[node_index].child[i], total_weight);
    }
}

void Octree::write_node(std::ofstream& output, int node_index, const double& total_weight) {
    if (nodes[node_index].child[0] < 0) {
        output << node_index << " " << nodes[node_index].weight / total_weight << std::endl;
    }
    else {
        output << node_index;
        for (int i = 0; i < children_per_node; i++) output << " " << nodes[node_index].child[i];
        output << std::endl;
        for (int i = 0; i < children_per_node; i++) write_node(output, nodes[node_index].child[i], total_weight);
    }
}

void BinaryTree::weigh() {
    std::cout << "Reweighing the probabilities...";
    weigh(0, 0, 1);
}

void Quadtree::weigh() {
    std::cout << "Reweighing the probabilities...";
    weigh(0, 0, 0, 1);
}

void Octree::weigh() {
    std::cout << "Reweighing the probabilities...";
    weigh(0, 0, 0, 0, 1);
}

void BinaryTree::save_structure(std::string path, std::string prefix, std::string filename) {
    std::ofstream output(path + prefix + "_" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the bintree to the " << filename << " ...";
        output << xmin << " " << xmax << " " << std::endl; // bounds
        output << nodes.size() << std::endl; // total number of nodes
        write_node(output, 0, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void BinaryTree::save_structure(std::string path, std::string filename) {
    while (path.back() == '/') path.pop_back();
    std::ofstream output(path + "/" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the bintree to the " << filename << " ...";
        output << xmin << " " << xmax << " " << std::endl; // bounds
        output << nodes.size() << std::endl; // total number of nodes
        write_node(output, 0, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void Quadtree::save_structure(std::string path, std::string prefix, std::string filename) {
    std::ofstream output(path + prefix + "_" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the quadtree to the " << filename << " ...";
        output << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl; // bounds
        output << nodes.size() << std::endl; // total number of nodes
        write_node(output, 0, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void Quadtree::save_structure(std::string path, std::string filename) {
    while (path.back() == '/') path.pop_back();
    std::ofstream output(path + "/" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the quadtree to the " << filename << " ...";
        output << xmin << " " << ymin << " " << xmax << " " << ymax << std::endl; // bounds
        output << nodes.size() << std::endl; // total number of nodes
        write_node(output, 0, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void Octree::save_structure(std::string path, std::string prefix, std::string filename) {
    std::ofstream output(path + prefix + "_" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the quadtree to the " << filename << " ...";
        output << xmin << " " << ymin << " " << zmin << " " << xmax << " " << ymax << " " << zmax << std::endl; // bounds
        output << nodes.size() << std::endl; // total number of nodes
        write_node(output, 0, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void Octree::save_structure(std::string path, std::string filename) {
    while (path.back() == '/') path.pop_back();
    std::ofstream output(path + "/" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the quadtree to the " << filename << " ...";
        output << xmin << " " << ymin << " " << zmin << " " << xmax << " " << ymax << " " << zmax << std::endl; // bounds
        output << nodes.size() << std::endl; // total number of nodes
        write_node(output, 0, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void BinaryTree::write_leaf(std::ofstream& output, int node_index, double x_scaled, double current_depth, const double& total_weight) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        double weight = nodes[node_index].weight / total_weight;
        double dx = XMAX - XMIN;
        double pdf = weight / dx;
        output << XMIN << " " << XMAX << " " << weight << " " << pdf << std::endl;
    }
    else {
        current_depth /= 2;
        write_leaf(output, nodes[node_index].child[0], x_scaled, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[1], x_scaled + current_depth, current_depth, total_weight);
    }
}

void Quadtree::write_leaf(std::ofstream& output, int node_index, double x_scaled, double y_scaled, double current_depth, const double& total_weight) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double YMIN = ymin + (ymax - ymin) * y_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        double YMAX = ymin + (ymax - ymin) * (y_scaled + current_depth);
        double weight = nodes[node_index].weight / total_weight;
        double dxdy = (XMAX - XMIN) * (YMAX - YMIN);
        double pdf = weight / dxdy;
        output << XMIN << " " << YMIN << " " << XMAX << " " << YMAX << " " << weight << " " << pdf << std::endl;
    }
    else {
        current_depth /= 2;
        write_leaf(output, nodes[node_index].child[0], x_scaled, y_scaled, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[1], x_scaled + current_depth, y_scaled, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[2], x_scaled + current_depth, y_scaled + current_depth, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[3], x_scaled, y_scaled + current_depth, current_depth, total_weight);
    }
}

void Octree::write_leaf(std::ofstream& output, int node_index, double x_scaled, double y_scaled, double z_scaled, double current_depth, const double& total_weight) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double YMIN = ymin + (ymax - ymin) * y_scaled;
        double ZMIN = zmin + (zmax - zmin) * z_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        double YMAX = ymin + (ymax - ymin) * (y_scaled + current_depth);
        double ZMAX = zmin + (zmax - zmin) * (z_scaled + current_depth);
        double weight = nodes[node_index].weight / total_weight;
        double dxdydz = (XMAX - XMIN) * (YMAX - YMIN) * (ZMAX - ZMIN);
        double pdf = weight / dxdydz;
        output << XMIN << " " << YMIN << " " << ZMIN << " " << XMAX << " " << YMAX << " " << ZMAX << " " << weight << " " << pdf << std::endl;
    }
    else {
        current_depth /= 2;
        write_leaf(output, nodes[node_index].child[0], x_scaled, y_scaled, z_scaled, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[1], x_scaled + current_depth, y_scaled, z_scaled, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[2], x_scaled + current_depth, y_scaled + current_depth, z_scaled, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[3], x_scaled, y_scaled + current_depth, z_scaled, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[4], x_scaled, y_scaled, z_scaled + current_depth, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[5], x_scaled + current_depth, y_scaled, z_scaled + current_depth, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[6], x_scaled + current_depth, y_scaled + current_depth, z_scaled + current_depth, current_depth, total_weight);
        write_leaf(output, nodes[node_index].child[7], x_scaled, y_scaled + current_depth, z_scaled + current_depth, current_depth, total_weight);
    }
}

void BinaryTree::save_leaves(std::string path, std::string prefix, std::string filename) {
    auto copy = *this; copy.weigh();
    std::ofstream output(path + prefix + "_" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the tiles to file " << filename << " ...";
        copy.write_leaf(output, 0, 0, 1, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void BinaryTree::save_leaves(std::string path, std::string filename) {
    auto copy = *this; copy.weigh();
    while (path.back() == '/') path.pop_back();
    std::ofstream output(path + "/" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the tiles to file " << filename << " ...";
        copy.write_leaf(output, 0, 0, 1, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void Quadtree::save_leaves(std::string path, std::string prefix, std::string filename) {
    auto copy = *this; copy.weigh();
    std::ofstream output(path + prefix + "_" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the tiles to file " << filename << " ...";
        copy.write_leaf(output, 0, 0, 0, 1, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void Quadtree::save_leaves(std::string path, std::string filename) {
    auto copy = *this; copy.weigh();
    while (path.back() == '/') path.pop_back();
    std::ofstream output(path + "/" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the tiles to file " << filename << " ...";
        copy.write_leaf(output, 0, 0, 0, 1, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void Octree::save_leaves(std::string path, std::string prefix, std::string filename) {
    auto copy = *this; copy.weigh();
    std::ofstream output(path + prefix + "_" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the tiles to file " << filename << " ...";
        copy.write_leaf(output, 0, 0, 0, 0, 1, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void Octree::save_leaves(std::string path, std::string filename) {
    auto copy = *this; copy.weigh();
    while (path.back() == '/') path.pop_back();
    std::ofstream output(path + "/" + filename);
    double total_weight = count_total();
    if (output) {
        std::cout << "Writing the tiles to file " << filename << " ...";
        copy.write_leaf(output, 0, 0, 0, 0, 1, total_weight);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

std::pair<double, double> Quadtree::center(int node_index) {
    std::pair<double, double> res;
    res.first = std::pow(0.5, depth);
    res.second = std::pow(0.5, depth);
    while (nodes[node_index].parent >= 0) {
        if (nodes[nodes[node_index].parent].child[1] == node_index) {
            res.first += std::pow(0.5, nodes[node_index].depth);
        }
        else if (nodes[nodes[node_index].parent].child[2] == node_index) {
            res.first += std::pow(0.5, nodes[node_index].depth);
            res.second += std::pow(0.5, nodes[node_index].depth);
        }
        else {
            res.second += std::pow(0.5, nodes[node_index].depth);
        }
        node_index = nodes[node_index].parent;
    }
    res.first = xmin + (xmax - xmin) * res.first;
    res.second = ymin + (ymax - ymin) * res.second;
    return res;
}

void Quadtree::center_and_dxdy(int node_index, double& center_x, double& center_y, double& dxdy) {
    center_x = 0;
    center_y = 0;
    int node_index_orig = node_index;
    while (nodes[node_index].parent >= 0) {
        if (nodes[nodes[node_index].parent].child[1] == node_index) {
            center_x += std::pow(0.5, nodes[node_index].depth);
        }
        else if (nodes[nodes[node_index].parent].child[2] == node_index) {
            center_x += std::pow(0.5, nodes[node_index].depth);
            center_y += std::pow(0.5, nodes[node_index].depth);
        }
        else {
            center_y += std::pow(0.5, nodes[node_index].depth);
        }
        node_index = nodes[node_index].parent;
    }
    double XMIN = xmin + (xmax - xmin) * center_x;
    double YMIN = ymin + (ymax - ymin) * center_y;
    double XMAX = xmin + (xmax - xmin) * (center_x + std::pow(0.5, nodes[node_index_orig].depth));
    double YMAX = ymin + (ymax - ymin) * (center_y + std::pow(0.5, nodes[node_index_orig].depth));
    dxdy = (XMAX - XMIN) * (YMAX - YMIN);
    center_x = 0.5 * (XMIN + XMAX);
    center_y = 0.5 * (YMIN + YMAX);
}

void BinaryTree::traverse_nodes(std::map<double, double>& pdf_weight, int node_index, double x_scaled, double current_depth, const double& total_weight) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        double weight = nodes[node_index].weight / total_weight;
        double dx = XMAX - XMIN;
        double pdf = weight / dx * scale(0.5 * (XMIN + XMAX));
        if (pdf_weight.find(pdf) == pdf_weight.end()) pdf_weight[pdf] = weight;
        else pdf_weight[pdf] += weight;
    }
    else {
        current_depth /= 2;
        traverse_nodes(pdf_weight, nodes[node_index].child[0], x_scaled, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[1], x_scaled + current_depth, current_depth, total_weight);
    }
}

void Quadtree::traverse_nodes(std::map<double, double>& pdf_weight, int node_index, double x_scaled, double y_scaled, double current_depth, const double& total_weight) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double YMIN = ymin + (ymax - ymin) * y_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        double YMAX = ymin + (ymax - ymin) * (y_scaled + current_depth);
        double weight = nodes[node_index].weight / total_weight;
        double dxdy = (XMAX - XMIN) * (YMAX - YMIN);
        double pdf = weight / dxdy * scale_x(0.5 * (XMIN + XMAX)) * scale_y(0.5 * (YMIN + YMAX));
        if (pdf_weight.find(pdf) == pdf_weight.end()) pdf_weight[pdf] = weight;
        else pdf_weight[pdf] += weight;
    }
    else {
        current_depth /= 2;
        traverse_nodes(pdf_weight, nodes[node_index].child[0], x_scaled, y_scaled, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[1], x_scaled + current_depth, y_scaled, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[2], x_scaled + current_depth, y_scaled + current_depth, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[3], x_scaled, y_scaled + current_depth, current_depth, total_weight);
    }
}

void Octree::traverse_nodes(std::map<double, double>& pdf_weight, int node_index, double x_scaled, double y_scaled, double z_scaled, double current_depth, const double& total_weight) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double YMIN = ymin + (ymax - ymin) * y_scaled;
        double ZMIN = zmin + (zmax - zmin) * z_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        double YMAX = ymin + (ymax - ymin) * (y_scaled + current_depth);
        double ZMAX = zmin + (zmax - zmin) * (z_scaled + current_depth);
        double weight = nodes[node_index].weight / total_weight;
        double dxdydz = (XMAX - XMIN) * (YMAX - YMIN) * (ZMAX - ZMIN);
        double pdf = weight / dxdydz * scale_x(0.5 * (XMIN + XMAX)) * scale_y(0.5 * (YMIN + YMAX)) * scale_z(0.5 * (ZMIN + ZMAX));
        if (pdf_weight.find(pdf) == pdf_weight.end()) pdf_weight[pdf] = weight;
        else pdf_weight[pdf] += weight;
    }
    else {
        current_depth /= 2;
        traverse_nodes(pdf_weight, nodes[node_index].child[0], x_scaled, y_scaled, z_scaled, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[1], x_scaled + current_depth, y_scaled, z_scaled, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[2], x_scaled + current_depth, y_scaled + current_depth, z_scaled, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[3], x_scaled, y_scaled + current_depth, z_scaled, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[4], x_scaled, y_scaled, z_scaled + current_depth, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[5], x_scaled + current_depth, y_scaled, z_scaled + current_depth, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[6], x_scaled + current_depth, y_scaled + current_depth, z_scaled + current_depth, current_depth, total_weight);
        traverse_nodes(pdf_weight, nodes[node_index].child[7], x_scaled, y_scaled + current_depth, z_scaled + current_depth, current_depth, total_weight);
    }
}

void BinaryTree::prob_to_pdf(std::string path, std::string prefix, std::string filename) {
    //auto copy = *this; copy.weigh();
    std::map<double, double> data; // (pdf, tile weight)
    double total_weight = count_total();
    traverse_nodes(data, 0, 0, 1, total_weight);
    std::map<double, double> mapping; // (CDF, pdf value)
    double weight_acc = 0;
    mapping[0] = 0;
    for (auto& tiles : data) {
        weight_acc += tiles.second;
        if (tiles.second > 0) mapping[weight_acc] = tiles.first;
        else mapping[weight_acc] += tiles.first;
    }
    std::ofstream output(path + prefix + "_" + filename);
    output << std::setprecision(12);
    if (output) {
        std::cout << "Writing the CDF mapping to file " << filename << " ...";
        for (auto& cdf : mapping) {
            output << cdf.first << " " << cdf.second << std::endl;
        }
    }
    std::cout << " done!" << std::endl;
    output.close();
}

void BinaryTree::prob_to_pdf(std::string path, std::string filename) {
    //auto copy = *this; copy.weigh();
    std::map<double, double> data; // (pdf, tile weight)
    double total_weight = count_total();
    traverse_nodes(data, 0, 0, 1, total_weight);
    std::map<double, double> mapping; // (CDF, pdf value)
    double weight_acc = 0;
    mapping[0] = 0;
    for (auto& tiles : data) {
        weight_acc += tiles.second;
        if (tiles.second > 0) mapping[weight_acc] = tiles.first;
        else mapping[weight_acc] += tiles.first;
    }
    while (path.back() == '/') path.pop_back();
    std::ofstream output(path + "/" + filename);
    output << std::setprecision(12);
    if (output) {
        std::cout << "Writing the CDF mapping to file " << filename << " ...";
        for (auto& cdf : mapping) {
            output << cdf.first << " " << cdf.second << std::endl;
        }
    }
    std::cout << " done!" << std::endl;
    output.close();
}

void Quadtree::prob_to_pdf(std::string path, std::string prefix, std::string filename) {
    //auto copy = *this; copy.weigh();
    std::map<double, double> data; // (pdf, tile weight)
    double total_weight = count_total();
    traverse_nodes(data, 0, 0, 0, 1, total_weight);
    std::map<double, double> mapping; // (CDF, pdf value)
    double weight_acc = 0;
    mapping[0] = 0;
    for (auto& tiles : data) {
        weight_acc += tiles.second;
        if (tiles.second > 0) mapping[weight_acc] = tiles.first;
        else mapping[weight_acc] += tiles.first;
    }
    std::ofstream output(path + prefix + "_" + filename);
    output << std::setprecision(12);
    if (output) {
        std::cout << "Writing the CDF mapping to file " << filename << " ...";
        for (auto& cdf : mapping) {
            output << cdf.first << " " << cdf.second << std::endl;
        }
    }
    std::cout << " done!" << std::endl;
    output.close();
}

void Quadtree::prob_to_pdf(std::string path_out, std::string filename) {
    //auto copy = *this; copy.weigh();
    std::map<double, double> data; // (pdf, tile weight)
    double total_weight = count_total();
    traverse_nodes(data, 0, 0, 0, 1, total_weight);
    std::map<double, double> mapping; // (CDF, pdf value)
    double weight_acc = 0;
    mapping[0] = 0;
    for (auto& tiles : data) {
        weight_acc += tiles.second;
        if (tiles.second > 0) mapping[weight_acc] = tiles.first;
        else mapping[weight_acc] += tiles.first;
    }
    while (path_out.back() == '/') path_out.pop_back();
    std::ofstream output(path_out + "/" + filename);
    output << std::setprecision(12);
    if (output) {
        std::cout << "Writing the CDF mapping to file " << filename << " ...";
        for (auto& cdf : mapping) {
            output << cdf.first << " " << cdf.second << std::endl;
        }
    }
    std::cout << " done!" << std::endl;
    output.close();
}

void Octree::prob_to_pdf(std::string path, std::string prefix, std::string filename) {
    auto copy = *this; copy.weigh();
    std::map<double, double> data; // (pdf, tile weight)
    double total_weight = copy.count_total();
    copy.traverse_nodes(data, 0, 0, 0, 0, 1, total_weight);
    std::map<double, double> mapping; // (CDF, pdf value)
    double weight_acc = 0;
    mapping[0] = 0;
    for (auto& tiles : data) {
        weight_acc += tiles.second;
        if (tiles.second > 0) mapping[weight_acc] = tiles.first;
        else mapping[weight_acc] += tiles.first;
    }
    std::ofstream output(path + prefix + "_" + filename);
    output << std::setprecision(12);
    if (output) {
        std::cout << "Writing the CDF mapping to file " << filename << " ...";
        for (auto& cdf : mapping) {
            output << cdf.first << " " << cdf.second << std::endl;
        }
    }
    std::cout << " done!" << std::endl;
    output.close();
}

void Octree::prob_to_pdf(std::string path, std::string filename) {
    auto copy = *this; copy.weigh();
    std::map<double, double> data; // (pdf, tile weight)
    double total_weight = copy.count_total();
    copy.traverse_nodes(data, 0, 0, 0, 0, 1, total_weight);
    std::map<double, double> mapping; // (CDF, pdf value)
    double weight_acc = 0;
    mapping[0] = 0;
    for (auto& tiles : data) {
        weight_acc += tiles.second;
        if (tiles.second > 0) mapping[weight_acc] = tiles.first;
        else mapping[weight_acc] += tiles.first;
    }
    while (path.back() == '/') path.pop_back();
    std::ofstream output(path + "/" + filename);
    output << std::setprecision(12);
    if (output) {
        std::cout << "Writing the CDF mapping to file " << filename << " ...";
        for (auto& cdf : mapping) {
            output << cdf.first << " " << cdf.second << std::endl;
        }
    }
    std::cout << " done!" << std::endl;
    output.close();
}

void BinaryTree::adjust_bounds(double& XMIN, double& XMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x) {
    std::vector<double> x;
    int N = filename_x.size();

    XMIN = 1e10;
    XMAX = -1e10;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        for (auto& x : x) {
            if (XMIN > x) XMIN = x;
            if (XMAX < x) XMAX = x;
        }
    }
    XMIN = XMIN;
    XMAX = XMAX;
    double dx = XMAX - XMIN;
    XMIN -= dx / 100;
    XMAX += dx / 100;

    int N_bins = 1000;

    std::map<int, double> bins;
    for (int i = 0; i < N_bins; i++) bins[i] = 0;

    double total_weight = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }

        for (auto& x : x) {
            int index = std::floor(N_bins * (x - XMIN) / (XMAX - XMIN));
            bins[index] += 1;
            total_weight += 1;
        }
    }

    std::map<double, std::tuple<double, double, double>> bins_pdf; // <pdf, (xmin, xmax, weight)>

    for (auto& bin : bins) {
        dx = ((bin.first + 1) * (XMAX - XMIN) / N_bins) - (bin.first * (XMAX - XMIN) / N_bins);
        bins_pdf[bin.second / (total_weight * dx)] = std::tuple<double, double, double>(XMIN + bin.first * (XMAX - XMIN) / N_bins, XMIN + (bin.first + 1) * (XMAX - XMIN) / N_bins, bin.second / total_weight);
    }

    total_weight = 0;
    XMIN = 1e10; XMAX = -1e10;
    for (auto bin = bins_pdf.rbegin(); bin != bins_pdf.rend(); bin++) {
        total_weight += std::get<2>(bin->second);
        if (XMIN > std::get<0>(bin->second)) XMIN = std::get<0>(bin->second);
        if (XMAX < std::get<1>(bin->second)) XMAX = std::get<1>(bin->second);
        if (total_weight > threshold_probability) break;
    }
    std::cout << total_weight << std::endl;
}

void Quadtree::adjust_bounds_x(double& XMIN, double& XMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x) {
    std::vector<double> x;
    int N = filename_x.size();

    XMIN = 1e10;
    XMAX = -1e10;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        for (auto& x : x) {
            if (XMIN > x) XMIN = x;
            if (XMAX < x) XMAX = x;
        }
    }
    XMIN = scale_x(XMIN);
    XMAX = scale_x(XMAX);
    double dx = XMAX - XMIN;
    XMIN -= dx / 100;
    XMAX += dx / 100;

    int N_bins = 1000;

    std::map<int, double> bins;
    for (int i = 0; i < N_bins; i++) bins[i] = 0;

    double total_weight = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }

        for (auto& x : x) {
            int index = std::floor(N_bins * (scale_x(x) - XMIN) / (XMAX - XMIN));
            bins[index] += 1;
            total_weight += 1;
        }
    }

    std::map<double, std::tuple<double, double, double>> bins_pdf; // <pdf, (xmin, xmax, weight)>

    for (auto& bin : bins) {
        dx = ((bin.first + 1) * (XMAX - XMIN) / N_bins) - (bin.first * (XMAX - XMIN) / N_bins);
        bins_pdf[bin.second / (total_weight * dx)] = std::tuple<double, double, double>(XMIN + bin.first * (XMAX - XMIN) / N_bins, XMIN + (bin.first + 1) * (XMAX - XMIN) / N_bins, bin.second / total_weight);
    }

    total_weight = 0;
    XMIN = 1e10; XMAX = -1e10;
    for (auto bin = bins_pdf.rbegin(); bin != bins_pdf.rend(); bin++) {
        total_weight += std::get<2>(bin->second);
        if (XMIN > std::get<0>(bin->second)) XMIN = std::get<0>(bin->second);
        if (XMAX < std::get<1>(bin->second)) XMAX = std::get<1>(bin->second);
        if (total_weight > threshold_probability) break;
    }
    std::cout << total_weight << std::endl;
}

void BinaryTree::adjust_bounds(double& XMIN, double& XMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_w) {
    std::vector<double> x, w;
    int N = filename_x.size();

    XMIN = 1e10;
    XMAX = -1e10;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        for (auto& x : x) {
            if (XMIN > x) XMIN = x;
            if (XMAX < x) XMAX = x;
        }
    }
    XMIN = XMIN;
    XMAX = XMAX;
    double dx = XMAX - XMIN;
    XMIN -= dx / 100;
    XMAX += dx / 100;

    int N_bins = 1000;

    std::map<int, double> bins;
    for (int i = 0; i < N_bins; i++) bins[i] = 0;

    double total_weight = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_w[i], w) == -1) {
            std::cout << "Filename " << filename_w[i] << " not found!" << std::endl;
            continue;
        }

        for (int i = 0; i < x.size(); i++) {
            int index = std::floor(N_bins * (x[i] - XMIN) / (XMAX - XMIN));
            bins[index] += w[i];
            total_weight += w[i];
        }
    }

    std::map<double, std::tuple<double, double, double>> bins_pdf; // <pdf, (xmin, xmax, weight)>

    for (auto& bin : bins) {
        dx = ((bin.first + 1) * (XMAX - XMIN) / N_bins) - (bin.first * (XMAX - XMIN) / N_bins);
        bins_pdf[bin.second / (total_weight * dx)] = std::tuple<double, double, double>(XMIN + bin.first * (XMAX - XMIN) / N_bins, XMIN + (bin.first + 1) * (XMAX - XMIN) / N_bins, bin.second / total_weight);
    }

    total_weight = 0;
    XMIN = 1e10; XMAX = -1e10;
    for (auto bin = bins_pdf.rbegin(); bin != bins_pdf.rend(); bin++) {
        total_weight += std::get<2>(bin->second);
        if (XMIN > std::get<0>(bin->second)) XMIN = std::get<0>(bin->second);
        if (XMAX < std::get<1>(bin->second)) XMAX = std::get<1>(bin->second);
        if (total_weight > threshold_probability) break;
    }
    std::cout << total_weight << std::endl;
}

void Quadtree::adjust_bounds_x(double& XMIN, double& XMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_w) {
    std::vector<double> x, w;
    int N = filename_x.size();

    XMIN = 1e10;
    XMAX = -1e10;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        for (auto& x : x) {
            if (XMIN > x) XMIN = x;
            if (XMAX < x) XMAX = x;
        }
    }
    XMIN = scale_x(XMIN);
    XMAX = scale_x(XMAX);
    double dx = XMAX - XMIN;
    XMIN -= dx / 100;
    XMAX += dx / 100;

    int N_bins = 1000;

    std::map<int, double> bins;
    for (int i = 0; i < N_bins; i++) bins[i] = 0;

    double total_weight = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_w[i], w) == -1) {
            std::cout << "Filename " << filename_w[i] << " not found!" << std::endl;
            continue;
        }

        for (int i = 0; i < x.size(); i++) {
            int index = std::floor(N_bins * (scale_x(x[i]) - XMIN) / (XMAX - XMIN));
            bins[index] += w[i];
            total_weight += w[i];
        }
    }

    std::map<double, std::tuple<double, double, double>> bins_pdf; // <pdf, (xmin, xmax, weight)>

    for (auto& bin : bins) {
        dx = ((bin.first + 1) * (XMAX - XMIN) / N_bins) - (bin.first * (XMAX - XMIN) / N_bins);
        bins_pdf[bin.second / (total_weight * dx)] = std::tuple<double, double, double>(XMIN + bin.first * (XMAX - XMIN) / N_bins, XMIN + (bin.first + 1) * (XMAX - XMIN) / N_bins, bin.second / total_weight);
    }

    total_weight = 0;
    XMIN = 1e10; XMAX = -1e10;
    for (auto bin = bins_pdf.rbegin(); bin != bins_pdf.rend(); bin++) {
        total_weight += std::get<2>(bin->second);
        if (XMIN > std::get<0>(bin->second)) XMIN = std::get<0>(bin->second);
        if (XMAX < std::get<1>(bin->second)) XMAX = std::get<1>(bin->second);
        if (total_weight > threshold_probability) break;
    }
    std::cout << total_weight << std::endl;
}

void Quadtree::adjust_bounds_y(double& YMIN, double& YMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_y) {
    std::vector<double> y;
    int N = filename_y.size();

    YMIN = 1e10;
    YMAX = -1e10;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        for (auto& y : y) {
            if (YMIN > y) YMIN = y;
            if (YMAX < y) YMAX = y;
        }
    }
    YMIN = scale_y(YMIN);
    YMAX = scale_y(YMAX);
    double dy = YMAX - YMIN;
    YMIN -= dy / 100;
    YMAX += dy / 100;

    int N_bins = 1000;

    std::map<int, double> bins;
    for (int i = 0; i < N_bins; i++) bins[i] = 0;

    double total_weight = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }

        for (auto& y : y) {
            int index = std::floor(N_bins * (scale_y(y) - YMIN) / (YMAX - YMIN));
            bins[index] += 1;
            total_weight += 1;
        }
    }

    std::map<double, std::tuple<double, double, double>> bins_pdf; // <pdf, (xmin, xmax, weight)>

    for (auto& bin : bins) {
        dy = ((bin.first + 1) * (YMAX - YMIN) / N_bins) - (bin.first * (YMAX - YMIN) / N_bins);
        bins_pdf[bin.second / (total_weight * dy)] = std::tuple<double, double, double>(YMIN + bin.first * (YMAX - YMIN) / N_bins, YMIN + (bin.first + 1) * (YMAX - YMIN) / N_bins, bin.second / total_weight);
    }

    total_weight = 0;
    YMIN = 1e10; YMAX = -1e10;
    for (auto bin = bins_pdf.rbegin(); bin != bins_pdf.rend(); bin++) {
        total_weight += std::get<2>(bin->second);
        if (YMIN > std::get<0>(bin->second)) YMIN = std::get<0>(bin->second);
        if (YMAX < std::get<1>(bin->second)) YMAX = std::get<1>(bin->second);
        if (total_weight > threshold_probability) break;
    }
    std::cout << total_weight << std::endl;
}

void Quadtree::adjust_bounds_y(double& YMIN, double& YMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_y, std::vector<std::string> filename_w) {
    std::vector<double> y, w;
    int N = filename_y.size();

    YMIN = 1e10;
    YMAX = -1e10;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        for (auto& y : y) {
            if (YMIN > y) YMIN = y;
            if (YMAX < y) YMAX = y;
        }
    }
    YMIN = scale_y(YMIN);
    YMAX = scale_y(YMAX);
    double dy = YMAX - YMIN;
    YMIN -= dy / 100;
    YMAX += dy / 100;

    int N_bins = 1000;

    std::map<int, double> bins;
    for (int i = 0; i < N_bins; i++) bins[i] = 0;

    double total_weight = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_w[i], w) == -1) {
            std::cout << "Filename " << filename_w[i] << " not found!" << std::endl;
            continue;
        }

        for (int i = 0; i < y.size(); i++) {
            int index = std::floor(N_bins * (scale_y(y[i]) - YMIN) / (YMAX - YMIN));
            bins[index] += w[i];
            total_weight += w[i];
        }
    }

    std::map<double, std::tuple<double, double, double>> bins_pdf; // <pdf, (xmin, xmax, weight)>

    for (auto& bin : bins) {
        dy = ((bin.first + 1) * (YMAX - YMIN) / N_bins) - (bin.first * (YMAX - YMIN) / N_bins);
        bins_pdf[bin.second / (total_weight * dy)] = std::tuple<double, double, double>(YMIN + bin.first * (YMAX - YMIN) / N_bins, YMIN + (bin.first + 1) * (YMAX - YMIN) / N_bins, bin.second / total_weight);
    }

    total_weight = 0;
    YMIN = 1e10; YMAX = -1e10;
    for (auto bin = bins_pdf.rbegin(); bin != bins_pdf.rend(); bin++) {
        total_weight += std::get<2>(bin->second);
        if (YMIN > std::get<0>(bin->second)) YMIN = std::get<0>(bin->second);
        if (YMAX < std::get<1>(bin->second)) YMAX = std::get<1>(bin->second);
        if (total_weight > threshold_probability) break;
    }
    std::cout << total_weight << std::endl;
}

void Quadtree::adjust_bounds(double& XMIN, double& XMAX, double& YMIN, double& YMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y) {
    adjust_bounds_x(XMIN, XMAX, threshold_probability, path, prefix, filename_x);
    adjust_bounds_y(YMIN, YMAX, threshold_probability, path, prefix, filename_y);
}

void Quadtree::adjust_bounds(double& XMIN, double& XMAX, double& YMIN, double& YMAX, double threshold_probability, std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_w) {
    adjust_bounds_x(XMIN, XMAX, threshold_probability, path, prefix, filename_x, filename_w);
    adjust_bounds_y(YMIN, YMAX, threshold_probability, path, prefix, filename_y, filename_w);
}

void BinaryTree::fix_depth(int index, int depth) {
    nodes[index].depth = depth;
    if (nodes[index].child[0] > 0) {
        fix_depth(nodes[index].child[0], depth + 1);
        fix_depth(nodes[index].child[1], depth + 1);
    }
}

void Quadtree::fix_depth(int index, int depth) {
    nodes[index].depth = depth;
    if (nodes[index].child[0] > 0) {
        fix_depth(nodes[index].child[0], depth + 1);
        fix_depth(nodes[index].child[1], depth + 1);
        fix_depth(nodes[index].child[2], depth + 1);
        fix_depth(nodes[index].child[3], depth + 1);
    }
}

void Octree::fix_depth(int index, int depth) {
    nodes[index].depth = depth;
    if (nodes[index].child[0] > 0) {
        for (int i = 0; i < children_per_node; i++) fix_depth(nodes[index].child[i], depth + 1);
    }
}

void BinaryTree::fix_depth() {
    fix_depth(0, 0);
}

void Quadtree::fix_depth() {
    fix_depth(0, 0);
}

void Octree::fix_depth() {
    fix_depth(0, 0);
}

void BinaryTree::load_tree(std::string path, std::string prefix, std::string filename, std::string scale) {
    std::fstream input(path + prefix + "_" + filename);
    std::string line;
    std::getline(input, line);
    std::stringstream iss(line);
    iss >> xmin >> xmax;
    std::getline(input, line);
    nodes.resize(std::stoi(line));
    while (std::getline(input, line)) {
        std::stringstream iss(line);
        double a, b, c;
        iss >> a >> b;
        if (iss.rdbuf()->in_avail()) {
            iss >> c;
            nodes[(int)a].child[0] = b;
            nodes[(int)b].parent = a;

            nodes[(int)a].child[1] = c;
            nodes[(int)c].parent = a;
        }
        else {
            nodes[(int)a].weight = b;
        }
    }
    set_scales(scale);
    fix_depth();
}

BinaryTree load_tree(std::string path, std::string prefix, std::string filename, std::string scale) {
    BinaryTree res;
    res.load_tree(path, prefix, filename, scale);
    return res;
}

void Quadtree::load_tree(std::string path, std::string prefix, std::string filename, std::string scale_X, std::string scale_Y) {
    std::fstream input(path + prefix + "_" + filename);
    std::string line;
    std::getline(input, line);
    std::stringstream iss(line);
    iss >> xmin >> ymin >> xmax >> ymax;
    std::getline(input, line);
    nodes.resize(std::stoi(line));
    while (std::getline(input, line)) {
        std::stringstream iss(line);
        double a, b, c, d, e;
        iss >> a >> b;
        if (iss.rdbuf()->in_avail()) {
            iss >> c >> d >> e;
            nodes[(int)a].child[0] = b;
            nodes[(int)b].parent = a;

            nodes[(int)a].child[1] = c;
            nodes[(int)c].parent = a;

            nodes[(int)a].child[2] = d;
            nodes[(int)d].parent = a;

            nodes[(int)a].child[3] = e;
            nodes[(int)e].parent = a;
        }
        else {
            nodes[(int)a].weight = b;
        }
    }
    set_scales(scale_X, scale_Y);
    fix_depth();
}

void Octree::load_tree(std::string path, std::string prefix, std::string filename, std::string scale_X, std::string scale_Y, std::string scale_Z) {
    std::fstream input(path + prefix + "_" + filename);
    std::string line;
    std::getline(input, line);
    std::stringstream iss(line);
    iss >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax;
    std::getline(input, line);
    nodes.resize(std::stoi(line));
    while (std::getline(input, line)) {
        std::stringstream iss(line);
        double a, b, c, d, e, f, g, h, i;
        iss >> a >> b;
        if (iss.rdbuf()->in_avail()) {
            iss >> c >> d >> e >> f >> g >> h >> i;
            nodes[(int)a].child[0] = b;
            nodes[(int)b].parent = a;

            nodes[(int)a].child[1] = c;
            nodes[(int)c].parent = a;

            nodes[(int)a].child[2] = d;
            nodes[(int)d].parent = a;

            nodes[(int)a].child[3] = e;
            nodes[(int)e].parent = a;

            nodes[(int)a].child[4] = f;
            nodes[(int)e].parent = a;

            nodes[(int)a].child[5] = g;
            nodes[(int)e].parent = a;

            nodes[(int)a].child[6] = h;
            nodes[(int)e].parent = a;

            nodes[(int)a].child[7] = i;
            nodes[(int)e].parent = a;
        }
        else {
            nodes[(int)a].weight = b;
        }
    }
    set_scales(scale_X, scale_Y, scale_Z);
    fix_depth();
}

void Octree::slice_check_x(int N3D_index, Quadtree& q, std::map<int, int>& N3D_to_n2D, double x0_scaled, double x_scaled, double y_scaled, double z_scaled, double depth) {
    const auto& N = nodes[N3D_index];
    if (N.child[0] > 0) {
        int n = q.nodes.size();
        q.subdivide(N3D_to_n2D[N3D_index]);
        if (x0_scaled < 0.5) {
            N3D_to_n2D[N.child[0]] = n;
            N3D_to_n2D[N.child[3]] = n + 1;
            N3D_to_n2D[N.child[7]] = n + 2;
            N3D_to_n2D[N.child[4]] = n + 3;
            x0_scaled = 2 * x0_scaled;
            depth /= 2;
            slice_check_x(N.child[0], q, N3D_to_n2D, x0_scaled, x_scaled, y_scaled, z_scaled, depth);
            slice_check_x(N.child[3], q, N3D_to_n2D, x0_scaled, x_scaled, y_scaled + depth, z_scaled, depth);
            slice_check_x(N.child[4], q, N3D_to_n2D, x0_scaled, x_scaled, y_scaled + depth, z_scaled + depth, depth);
            slice_check_x(N.child[7], q, N3D_to_n2D, x0_scaled, x_scaled, y_scaled, z_scaled + depth, depth);
        }
        else {
            N3D_to_n2D[N.child[1]] = n;
            N3D_to_n2D[N.child[2]] = n + 1;
            N3D_to_n2D[N.child[6]] = n + 2;
            N3D_to_n2D[N.child[5]] = n + 3;
            x0_scaled = 2 * (x0_scaled - 0.5);
            depth /= 2;
            slice_check_x(N.child[1], q, N3D_to_n2D, x0_scaled, x_scaled + depth, y_scaled, z_scaled, depth);
            slice_check_x(N.child[2], q, N3D_to_n2D, x0_scaled, x_scaled + depth, y_scaled + depth, z_scaled, depth);
            slice_check_x(N.child[5], q, N3D_to_n2D, x0_scaled, x_scaled + depth, y_scaled + depth, z_scaled + depth, depth);
            slice_check_x(N.child[6], q, N3D_to_n2D, x0_scaled, x_scaled + depth, y_scaled, z_scaled + depth, depth);
        }
    }
    else {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double YMIN = ymin + (ymax - ymin) * y_scaled;
        double ZMIN = zmin + (zmax - zmin) * z_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + depth);
        double YMAX = ymin + (ymax - ymin) * (y_scaled + depth);
        double ZMAX = zmin + (zmax - zmin) * (z_scaled + depth);
        double dxdy = (YMAX - YMIN) * (ZMAX - ZMIN);
        double dxdydz = (XMAX - XMIN) * dxdy;
        double pdf = nodes[N3D_index].weight / dxdydz;
        q.nodes[N3D_to_n2D[N3D_index]].weight = pdf * dxdy;
    }
}

Quadtree Octree::slice(int coordinate, double value) {
    Quadtree res;
    if (coordinate == 0) {
        res.set_bounds(ymin, ymax, zmin, zmax);
        res.set_scales(scale_x, scale_y);
        double x_scaled = std::min(1., std::max(0., (value - xmin) / (xmax - xmin)));
        Node2D n;
        n.depth = 1;
        n.parent = -1;
        res.nodes.push_back(n);
        std::map<int, int> N3D_to_n2D;
        N3D_to_n2D[0] = 0;
        slice_check_x(0, res, N3D_to_n2D, x_scaled, 0, 0, 0, 1);
    }
    return res;
}

int find_index(double x, std::vector<std::pair<std::pair<double, double>, double>>& pdf) {
    int index1 = 0;
    int index2 = pdf.size() - 1;
    if (pdf[index1].first.first >= x) return index1;
    else if (pdf[index2].first.second <= x) return index2;
    else {
        while (index1 < index2 - 1) {
            int index_mid = (index1 + index2) / 2;
            if ((pdf[index_mid].first.first <= x) && (pdf[index_mid].first.second > x)) return index_mid;
            else if (pdf[index_mid].first.first > x) index2 = index_mid;
            else if (pdf[index_mid].first.second < x) index1 = index_mid;
        }
        if ((pdf[index1].first.first <= x) && (pdf[index1].first.second > x)) return index1;
        if ((pdf[index2].first.first <= x) && (pdf[index2].first.second > x)) return index2;
    }
}

double Quadtree::collect_weight_X(double XMIN, double XMAX, std::vector<std::pair<std::pair<double, double>, double>>& pdf_X) {
    int index1 = find_index(XMIN, pdf_X);
    int index2 = find_index(XMAX, pdf_X);
    double weight = 0;
    if (index1 == index2) {
        weight = pdf_X[index1].second * (XMAX - XMIN) / (pdf_X[index1].first.second - pdf_X[index1].first.first);
    }
    else if (index2 == index1 + 1) {
        double x = XMIN;
        double x1 = pdf_X[index1].first.first;
        double x2 = pdf_X[index1].first.second;
        weight += pdf_X[index1].second * (x2 - x) / (x2 - x1);

        x = XMAX;
        x1 = pdf_X[index2].first.first;
        x2 = pdf_X[index2].first.second;
        weight += pdf_X[index2].second * (x - x1) / (x2 - x1);
    }
    else {
        double x = XMIN;
        double x1 = pdf_X[index1].first.first;
        double x2 = pdf_X[index1].first.second;
        weight += pdf_X[index1].second * (x2 - x) / (x2 - x1);

        x = XMAX;
        x1 = pdf_X[index2].first.first;
        x2 = pdf_X[index2].first.second;
        weight += pdf_X[index2].second * (x - x1) / (x2 - x1);
        for (int i = index1 + 1; i < index2; i++) {
            weight += pdf_X[i].second;
        }
    }
    return weight;
}

double Quadtree::collect_weight_Y(double YMIN, double YMAX, std::vector<std::pair<std::pair<double, double>, double>>& pdf_Y) {
    int index1 = find_index(YMIN, pdf_Y);
    int index2 = find_index(YMAX, pdf_Y);
    double weight = 0;
    if (index1 == index2) {
        weight = pdf_Y[index1].second * (YMAX - YMIN) / (pdf_Y[index1].first.second - pdf_Y[index1].first.first);
    }
    else if (index2 == index1 + 1) {
        double y = YMIN;
        double y1 = pdf_Y[index1].first.first;
        double y2 = pdf_Y[index1].first.second;
        weight += pdf_Y[index1].second * (y2 - y) / (y2 - y1);

        y = YMAX;
        y1 = pdf_Y[index2].first.first;
        y2 = pdf_Y[index2].first.second;
        weight += pdf_Y[index2].second * (y - y1) / (y2 - y1);
    }
    else {
        double y = YMIN;
        double y1 = pdf_Y[index1].first.first;
        double y2 = pdf_Y[index1].first.second;
        weight += pdf_Y[index1].second * (y2 - y) / (y2 - y1);

        y = YMAX;
        y1 = pdf_Y[index2].first.first;
        y2 = pdf_Y[index2].first.second;
        weight += pdf_Y[index2].second * (y - y1) / (y2 - y1);
        for (int i = index1 + 1; i < index2; i++) {
            weight += pdf_Y[i].second;
        }
    }
    return weight;
}

void Quadtree::fill_marginalized(int node_index, double current_depth, double x_scaled, double y_scaled, std::vector<std::pair<std::pair<double, double>, double>>& pdf_rho, std::vector<std::pair<std::pair<double, double>, double>>& pdf_v) {
    if (nodes[node_index].child[0] < 0) {
        double XMIN = xmin + (xmax - xmin) * x_scaled;
        double XMAX = xmin + (xmax - xmin) * (x_scaled + current_depth);
        double YMIN = ymin + (ymax - ymin) * y_scaled;
        double YMAX = ymin + (ymax - ymin) * (y_scaled + current_depth);
        //std::cout << node_index << std::endl;
        nodes[node_index].weight = collect_weight_X(XMIN, XMAX, pdf_rho) * collect_weight_Y(YMIN, YMAX, pdf_v);
    }
    else {
        current_depth /= 2;
        fill_marginalized(nodes[node_index].child[0], current_depth, x_scaled, y_scaled, pdf_rho, pdf_v);
        fill_marginalized(nodes[node_index].child[1], current_depth, x_scaled + current_depth, y_scaled, pdf_rho, pdf_v);
        fill_marginalized(nodes[node_index].child[2], current_depth, x_scaled + current_depth, y_scaled + current_depth, pdf_rho, pdf_v);
        fill_marginalized(nodes[node_index].child[3], current_depth, x_scaled, y_scaled + current_depth, pdf_rho, pdf_v);
    }
}

void Quadtree::fill_marginalized(std::vector<std::pair<std::pair<double, double>, double>>& pdf_rho, std::vector<std::pair<std::pair<double, double>, double>>& pdf_v, std::string scale_X, std::string scale_Y) {
    set_scales(scale_X, scale_Y);
    fill_marginalized(0, 1, 0, 0, pdf_rho, pdf_v);
}

void write_leaf(std::ofstream& output, Quadtree &q1, Quadtree& q2, int node_index, double x_scaled, double y_scaled, double current_depth, const double& total_weight1, const double& total_weight2) {
    if (q1.nodes[node_index].child[0] < 0) {
        double XMIN = q1.xmin + (q1.xmax - q1.xmin) * x_scaled;
        double YMIN = q1.ymin + (q1.ymax - q1.ymin) * y_scaled;
        double XMAX = q1.xmin + (q1.xmax - q1.xmin) * (x_scaled + current_depth);
        double YMAX = q1.ymin + (q1.ymax - q1.ymin) * (y_scaled + current_depth);
        double weight1 = q1.nodes[node_index].weight / total_weight1;
        double weight2 = q2.nodes[node_index].weight / total_weight2;
        double dxdy = (std::pow(10, XMAX) - std::pow(10, XMIN)) * (YMAX - YMIN);
        double pdf1 = weight1 / dxdy;
        double pdf2 = weight2 / dxdy;
        output << XMIN << " " << YMIN << " " << XMAX << " " << YMAX << " " << pdf1 - pdf2 << std::endl;
    }
    else {
        current_depth /= 2;
        write_leaf(output, q1, q2, q1.nodes[node_index].child[0], x_scaled, y_scaled, current_depth, total_weight1, total_weight2);
        write_leaf(output, q1, q2, q1.nodes[node_index].child[1], x_scaled + current_depth, y_scaled, current_depth, total_weight1, total_weight2);
        write_leaf(output, q1, q2, q1.nodes[node_index].child[2], x_scaled + current_depth, y_scaled + current_depth, current_depth, total_weight1, total_weight2);
        write_leaf(output, q1, q2, q1.nodes[node_index].child[3], x_scaled, y_scaled + current_depth, current_depth, total_weight1, total_weight2);
    }
}

void save_difference(std::string path, std::string prefix, std::string filename, Quadtree& q1, Quadtree& q2) {
    std::ofstream output(path + prefix + "_" + filename);
    double total_weight1 = q1.count_total();
    double total_weight2 = q2.count_total();
    if (output) {
        std::cout << "Writing the tiles to file " << filename << " ...";
        write_leaf(output, q1, q2, 0, 0, 0, 1, total_weight1, total_weight2);
    }
    output.close();
    std::cout << " done!" << std::endl;
}

void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, Quadtree &q1, BinaryTree &X, BinaryTree &Y) {
    std::vector<double> x, y;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    std::cout << "Adding points..." << std::endl;
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        for (int j = 0; j < x.size(); j++) {
            q1.add_point(x[j], y[j]);
            X.add_point(x[j]);
            Y.add_point(y[j]);
        }
        progress = (double)(i + 1) / N;
        std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b" << 100 * progress << "     ";
    }
    std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
}

void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_weight, Quadtree& Q, BinaryTree& X, BinaryTree& Y) {
    std::vector<double> x, y, w;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_weight.size()) N = filename_weight.size();
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << backspace << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << backspace << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        if ((filename_weight[i] != filename_x[i]) && (filename_weight[i] != filename_y[i])) {
            int read_weights = ::read_from_binary(path, prefix, filename_weight[i], w);
            if (read_weights == -1) {
                std::cout << backspace << "Filename " << filename_weight[i] << " not found! Using the default weights." << std::endl;
                for (int j = 0; j < x.size(); j++) {
                    Q.add_point(x[i], y[i]);
                    X.add_point(x[i]);
                    Y.add_point(y[i]);
                }
            }
            else {
                for (int j = 0; j < x.size(); j++) {
                    Q.add_point(x[j], y[j], w[j]);
                    X.add_point(x[j], w[j]);
                    Y.add_point(y[j], w[j]);
                }
            }
        }
        else if (filename_weight[i] == filename_x[i]) {
            for (int j = 0; j < x.size(); j++) {
                Q.add_point(x[j], y[j], x[j]);
                X.add_point(x[j], x[j]);
                Y.add_point(y[j], x[j]);
            }
        }
        else if (filename_weight[i] == filename_y[i]) {
            for (int j = 0; j < x.size(); j++) {
                Q.add_point(x[j], y[j], y[j]);
                X.add_point(x[j], y[j]);
                Y.add_point(y[j], y[j]);
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << backspace << 100 * progress << "                                             ";
    }
    std::cout << backspace;
}

void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_weight1, std::vector<std::string> filename_weight2, Quadtree& Q1, Quadtree& Q2, BinaryTree& X1, BinaryTree& X2, BinaryTree& Y1, BinaryTree& Y2) {
    // weight1 is volume, weight2 is mass
    std::vector<double> x, y, w2;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_weight1.size()) N = filename_weight1.size();
    if (N > filename_weight2.size()) N = filename_weight2.size();
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << backspace << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << backspace << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        int read_weights2 = ::read_from_binary(path, prefix, filename_weight2[i], w2);
        if (read_weights2 == -1) {
            std::cout << backspace << "Filename " << filename_weight2[i] << " not found! Using the default weights." << std::endl;
            for (int j = 0; j < x.size(); j++) {
                Q1.add_point(x[i], y[i]);
                X1.add_point(x[i]);
                Y1.add_point(y[i]);
                Q2.add_point(x[i], y[i]);
                X2.add_point(x[i]);
                Y2.add_point(y[i]);
            }
        }
        else {
            for (int j = 0; j < x.size(); j++) {
                Q1.add_point(x[j], y[j]);
                X1.add_point(x[j]);
                Y1.add_point(y[j]);
                Q2.add_point(x[j], y[j], w2[j]);
                X2.add_point(x[j], w2[j]);
                Y2.add_point(y[j], w2[j]);
            }
        }
            
        progress = (double)(i + 1) / N;
        std::cout << backspace << 100 * progress << "                                             ";
    }
    std::cout << backspace;
}

void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::string> filename_weight1, std::vector<std::string> filename_weight2, std::vector<std::string> filename_weight3, Quadtree& Q1, Quadtree& Q2, Quadtree& Q3, BinaryTree& X1, BinaryTree& X2, BinaryTree& X3, BinaryTree& Y1, BinaryTree& Y2, BinaryTree& Y3) {
    // weight1 is always volume
    std::vector<double> x, y, w2, w3;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    if (N > filename_weight1.size()) N = filename_weight1.size();
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    double progress = 0;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << backspace << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << backspace << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        int read_weights2 = ::read_from_binary(path, prefix, filename_weight2[i], w2);
        int read_weights3 = ::read_from_binary(path, prefix, filename_weight3[i], w3);
        if ((read_weights2 == -1) && (read_weights3 == -1)) {
            std::cout << backspace << "Filename " << filename_weight2[i] << " not found! Using the default weights." << std::endl;
            std::cout << backspace << "Filename " << filename_weight3[i] << " not found! Using the default weights." << std::endl;
            for (int j = 0; j < x.size(); j++) {
                Q1.add_point(x[j], y[j]);
                X1.add_point(x[j]);
                Y1.add_point(y[j]);
                Q2.add_point(x[j], y[j]);
                X2.add_point(x[j]);
                Y2.add_point(y[j]);
                Q3.add_point(x[j], y[j]);
                X3.add_point(x[j]);
                Y3.add_point(y[j]);
            }
        }
        else if ((read_weights2 == 0) && (read_weights3 == 0)) {
            for (int j = 0; j < x.size(); j++) {
                Q1.add_point(x[j], y[j]);
                X1.add_point(x[j]);
                Y1.add_point(y[j]);
                Q2.add_point(x[j], y[j], w2[j]);
                X2.add_point(x[j], w2[j]);
                Y2.add_point(y[j], w2[j]);
                Q3.add_point(x[j], y[j], w3[j]);
                X3.add_point(x[j], w3[j]);
                Y3.add_point(y[j], w3[j]);
            }
        }
        else if (read_weights2 == 0) {
            for (int j = 0; j < x.size(); j++) {
                Q1.add_point(x[j], y[j]);
                X1.add_point(x[j]);
                Y1.add_point(y[j]);
                Q2.add_point(x[j], y[j], w2[j]);
                X2.add_point(x[j], w2[j]);
                Y2.add_point(y[j], w2[j]);
                Q3.add_point(x[j], y[j]);
                X3.add_point(x[j]);
                Y3.add_point(y[j]);
            }
        }
        else if (read_weights3 == 0) {
            for (int j = 0; j < x.size(); j++) {
                Q1.add_point(x[j], y[j]);
                X1.add_point(x[j]);
                Y1.add_point(y[j]);
                Q2.add_point(x[j], y[j]);
                X2.add_point(x[j]);
                Y2.add_point(y[j]);
                Q3.add_point(x[j], y[j], w3[j]);
                X3.add_point(x[j], w3[j]);
                Y3.add_point(y[j], w3[j]);
            }
        }
        else {
            for (int j = 0; j < x.size(); j++) {
                Q1.add_point(x[j], y[j]);
                X1.add_point(x[j]);
                Y1.add_point(y[j]);
                Q2.add_point(x[j], y[j], w2[j]);
                X2.add_point(x[j], w2[j]);
                Y2.add_point(y[j], w2[j]);
                Q3.add_point(x[j], y[j], w3[j]);
                X3.add_point(x[j], w3[j]);
                Y3.add_point(y[j], w3[j]);
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << backspace << 100 * progress << "                                             ";
    }
    std::cout << backspace;
}
/*
void load_points(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector< std::function<double(double, double)>> weight_fun, std::vector<std::function<double(double, double)>> field_1D, std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> field_2D, std::vector<BinaryTree>& data_1D, std::vector<Quadtree>& data_2D) {
    std::vector<double> x, y;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    double progress = 0;
    std::cout << backspace << "Loading points...                  " << std::endl;
    for (int i = 0; i < N; i++) {
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << backspace << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << backspace << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        for (int j = 0; j < x.size(); j++) {
            if (weight_fun.size() > 0) {
                for (int w = 0; w < weight_fun.size(); w++) {
                    for (int k = 0; k < field_1D.size(); k++) {
                        data_1D[w * field_1D.size() + k].add_point(field_1D[k](x[j], y[j]), weight_fun[w](x[j], y[j]));
                    }
                    for (int k = 0; k < field_2D.size(); k++) {
                        data_2D[w * field_2D.size() + k].add_point(field_2D[k].first(x[j], y[j]), field_2D[k].second(x[j], y[j]), weight_fun[w](x[j], y[j]));
                    }
                }
            }
            else {
                for (int k = 0; k < field_1D.size(); k++) {
                    data_1D[k].add_point(field_1D[k](x[j], y[j]));
                }
                for (int k = 0; k < field_2D.size(); k++) {
                    data_2D[k].add_point(field_2D[k].first(x[j], y[j]), field_2D[k].second(x[j], y[j]));
                }
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << backspace << 100 * progress << "                                             ";
    }
    std::cout << backspace;
}

void load_points(std::string path, std::vector<int> frames, std::string fieldname_x, std::string fieldname_y, std::vector< std::function<double(double, double)>> weight_fun, std::vector<std::function<double(double, double)>> field_1D, std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> field_2D, std::vector<BinaryTree>& data_1D, std::vector<Quadtree>& data_2D) {
    std::vector<double> x, y;
    for (int i = 0; i < frames.size(); i++) {
        int frame = frames[i];
        std::cout << "Loading points from frame " << frame << std::endl;
        auto cpu_list = get_cpu_list(path, frame);
        for (int c = 0; c < 1024; c++) {
            std::cout << "Loading points from frame " << frame << ", ";
            std::cout << "cpu" << string_pad_left(c, 0, 4) << ", ";
            std::string filename = cpu_filename(path, frame, c);
            for (auto & groupname:cpu_list[c]) {
                std::cout << groupname << std::endl;
                x = read_from_hdf5(filename, groupname, fieldname_x);
                y = read_from_hdf5(filename, groupname, fieldname_y);
                for (int j = 0; j < x.size(); j++) {
                    if (weight_fun.size() > 0) {
                        for (int w = 0; w < weight_fun.size(); w++) {
                            for (int k = 0; k < field_1D.size(); k++) {
                                data_1D[w * field_1D.size() + k].add_point(field_1D[k](x[j], y[j]), weight_fun[w](x[j], y[j]));
                            }
                            for (int k = 0; k < field_2D.size(); k++) {
                                data_2D[w * field_2D.size() + k].add_point(field_2D[k].first(x[j], y[j]), field_2D[k].second(x[j], y[j]), weight_fun[w](x[j], y[j]));
                            }
                        }
                    }
                    else {
                        for (int k = 0; k < field_1D.size(); k++) {
                            data_1D[k].add_point(field_1D[k](x[j], y[j]));
                        }
                        for (int k = 0; k < field_2D.size(); k++) {
                            data_2D[k].add_point(field_2D[k].first(x[j], y[j]), field_2D[k].second(x[j], y[j]));
                        }
                    }
                }
            }
        }
    }
}*/

void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::string scale_X, std::string scale_Y, Quadtree& Q, BinaryTree& X, BinaryTree& Y) {
    std::string backspace = "";
    for (int i = 0; i < 100; i++) backspace = backspace + "\b";
    Q.set_scales(scale_X, scale_Y);
    X.set_scales(scale_X);
    Y.set_scales(scale_Y);
    std::vector<double> x, y;
    double XMIN = 1e10;
    double YMIN = 1e10;
    double XMAX = -1e10;
    double YMAX = -1e10;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    for (int i = 0; i < N; i++) {
        std::cout << backspace << "Reading " << filename_x[i] << "...             ";
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << backspace << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << backspace << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        bool adjusted = false;
        for (int j = 0; j < x.size(); j++) {
            if (XMIN > x[j]) {
                XMIN = x[j];
                adjusted = true;
            }
            if (YMIN > y[j]) {
                YMIN = y[j];
                adjusted = true;
            }
            if (XMAX < x[j]) {
                XMAX = x[j];
                adjusted = true;
            }
            if (YMAX < y[j]) {
                YMAX = y[j];
                adjusted = true;
            }
        }
        if (adjusted) std::cout << backspace << "Bounds adjusted: (" << XMIN << ", " << XMAX << "), (" << YMIN << ", " << YMAX << ")                 ";
    }
    std::cout << std::endl;
    Q.set_bounds(XMIN, XMAX, YMIN, YMAX);
    X.set_bounds(XMIN, XMAX);
    Y.set_bounds(YMIN, YMAX);
    Q.nodes.clear();
    X.nodes.clear();
    Y.nodes.clear();
    Node2D n2; Node1D n1;
    n2.parent = -1;
    n2.depth = 1;
    n1.parent = -1;
    n1.depth = 1;
    Q.nodes.push_back(n2);
    X.nodes.push_back(n1);
    Y.nodes.push_back(n1);
}

void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::string scale_X, std::string scale_Y, Quadtree& Q1, Quadtree& Q2, BinaryTree& X1, BinaryTree& X2, BinaryTree& Y1, BinaryTree& Y2) {
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    Q1.set_scales(scale_X, scale_Y);
    X1.set_scales(scale_X);
    Y1.set_scales(scale_Y);
    Q2.set_scales(scale_X, scale_Y);
    X2.set_scales(scale_X);
    Y2.set_scales(scale_Y);
    std::vector<double> x, y;
    double XMIN = 1e10;
    double YMIN = 1e10;
    double XMAX = -1e10;
    double YMAX = -1e10;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    for (int i = 0; i < N; i++) {
        //std::cout << backspace << "Reading " << filename_x[i] << "...             ";
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << backspace << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << backspace << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        bool adjusted = false;
        for (int j = 0; j < x.size(); j++) {
            if (XMIN > x[j]) {
                XMIN = x[j];
                adjusted = true;
            }
            if (YMIN > y[j]) {
                YMIN = y[j];
                adjusted = true;
            }
            if (XMAX < x[j]) {
                XMAX = x[j];
                adjusted = true;
            }
            if (YMAX < y[j]) {
                YMAX = y[j];
                adjusted = true;
            }
        }
        std::cout << backspace << "Bounds adjusted: (" << XMIN << ", " << XMAX << "), (" << YMIN << ", " << YMAX << ") (" << std::round(10000 * (double)i / N) / 100 << "%)      ";
    }
    std::cout << backspace << "Bounds adjusted: (" << XMIN << ", " << XMAX << "), (" << YMIN << ", " << YMAX << ")                      " << std::endl;
    Q1.set_bounds(XMIN, XMAX, YMIN, YMAX);
    X1.set_bounds(XMIN, XMAX);
    Y1.set_bounds(YMIN, YMAX);
    Q2.set_bounds(XMIN, XMAX, YMIN, YMAX);
    X2.set_bounds(XMIN, XMAX);
    Y2.set_bounds(YMIN, YMAX);
    Q1.nodes.clear();
    X1.nodes.clear();
    Y1.nodes.clear();
    Q2.nodes.clear();
    X2.nodes.clear();
    Y2.nodes.clear();
    Node2D n2; Node1D n1;
    n2.parent = -1;
    n2.depth = 1;
    n1.parent = -1;
    n1.depth = 1;
    Q1.nodes.push_back(n2);
    X1.nodes.push_back(n1);
    Y1.nodes.push_back(n1);
    Q2.nodes.push_back(n2);
    X2.nodes.push_back(n1);
    Y2.nodes.push_back(n1);
}

void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::string scale_X, std::string scale_Y, Quadtree& Q1, Quadtree& Q2, Quadtree& Q3, BinaryTree& X1, BinaryTree& X2, BinaryTree& X3, BinaryTree& Y1, BinaryTree& Y2, BinaryTree& Y3) {
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    Q1.set_scales(scale_X, scale_Y);
    X1.set_scales(scale_X);
    Y1.set_scales(scale_Y);
    Q2.set_scales(scale_X, scale_Y);
    X2.set_scales(scale_X);
    Y2.set_scales(scale_Y);
    Q3.set_scales(scale_X, scale_Y);
    X3.set_scales(scale_X);
    Y3.set_scales(scale_Y);
    std::vector<double> x, y;
    double XMIN = 1e10;
    double YMIN = 1e10;
    double XMAX = -1e10;
    double YMAX = -1e10;
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    for (int i = 0; i < N; i++) {
        //std::cout << backspace << "Reading " << filename_x[i] << "...             ";
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        bool adjusted = false;
        for (int j = 0; j < x.size(); j++) {
            if (XMIN > x[j]) {
                XMIN = x[j];
                adjusted = true;
            }
            if (YMIN > y[j]) {
                YMIN = y[j];
                adjusted = true;
            }
            if (XMAX < x[j]) {
                XMAX = x[j];
                adjusted = true;
            }
            if (YMAX < y[j]) {
                YMAX = y[j];
                adjusted = true;
            }
        }
        std::cout << backspace << "Bounds adjusted: (" << XMIN << ", " << XMAX << "), (" << YMIN << ", " << YMAX << ") (" << std::round(10000 * (double)i / N) / 100 << "%)      ";
    }
    std::cout << backspace << "Bounds adjusted: (" << XMIN << ", " << XMAX << "), (" << YMIN << ", " << YMAX << ")                      " << std::endl;
    Q1.set_bounds(XMIN, XMAX, YMIN, YMAX);
    X1.set_bounds(XMIN, XMAX);
    Y1.set_bounds(YMIN, YMAX);
    Q2.set_bounds(XMIN, XMAX, YMIN, YMAX);
    X2.set_bounds(XMIN, XMAX);
    Y2.set_bounds(YMIN, YMAX);
    Q3.set_bounds(XMIN, XMAX, YMIN, YMAX);
    X3.set_bounds(XMIN, XMAX);
    Y3.set_bounds(YMIN, YMAX);
    Q1.nodes.clear();
    X1.nodes.clear();
    Y1.nodes.clear();
    Q2.nodes.clear();
    X2.nodes.clear();
    Y2.nodes.clear();
    Q3.nodes.clear();
    X3.nodes.clear();
    Y3.nodes.clear();
    Node2D n2; Node1D n1;
    n2.parent = -1;
    n2.depth = 1;
    n1.parent = -1;
    n1.depth = 1;
    Q1.nodes.push_back(n2);
    X1.nodes.push_back(n1);
    Y1.nodes.push_back(n1);
    Q2.nodes.push_back(n2);
    X2.nodes.push_back(n1);
    Y2.nodes.push_back(n1);
    Q3.nodes.push_back(n2);
    X3.nodes.push_back(n1);
    Y3.nodes.push_back(n1);
}

void initialize(std::string path, std::string prefix, std::vector<std::string> filename_x, std::vector<std::string> filename_y, std::vector<std::function<double(double, double)>> field_1D, std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> field_2D, std::vector<std::string> scale_1D, std::vector<std::pair<std::string, std::string>> scale_2D, std::vector<BinaryTree>& data_1D, std::vector<Quadtree>& data_2D) {
    std::string backspace = "";
    for (int i = 0; i < 300; i++) backspace = backspace + "\b";
    std::vector<double> x, y;
    std::vector<double> XMIN, XMAX, Z_XMIN, Z_XMAX, Z_YMIN, Z_YMAX;
    for (int i = 0; i < field_1D.size(); i++) {
        XMIN.push_back(1e10);
        XMAX.push_back(-1e10);
    }
    for (int i = 0; i < field_2D.size(); i++) {
        Z_XMIN.push_back(1e10);
        Z_XMAX.push_back(-1e10);
        Z_YMIN.push_back(1e10);
        Z_YMAX.push_back(-1e10);
    }
    int N = filename_x.size();
    if (N > filename_y.size()) N = filename_y.size();
    double progress = 0;
    std::cout << backspace << "Adjusting bounds...                " << std::endl;
    for (int i = 0; i < N; i++) {
        //std::cout << backspace << "Reading " << filename_x[i] << "...             ";
        if (::read_from_binary(path, prefix, filename_x[i], x) == -1) {
            std::cout << "Filename " << filename_x[i] << " not found!" << std::endl;
            continue;
        }
        if (::read_from_binary(path, prefix, filename_y[i], y) == -1) {
            std::cout << "Filename " << filename_y[i] << " not found!" << std::endl;
            continue;
        }
        bool adjusted = false;
        for (int j = 0; j < x.size(); j++) {
            for (int k = 0; k < field_1D.size(); k++) {
                if (XMIN[k] > field_1D[k](x[j], y[j])) {
                    XMIN[k] = field_1D[k](x[j], y[j]);
                    adjusted = true;
                }
                if (XMAX[k] < field_1D[k](x[j], y[j])) {
                    XMAX[k] = field_1D[k](x[j], y[j]);
                    adjusted = true;
                }
            }
            for (int k = 0; k < field_2D.size(); k++) {
                if (Z_XMIN[k] > field_2D[k].first(x[j], y[j])) {
                    Z_XMIN[k] = field_2D[k].first(x[j], y[j]);
                    adjusted = true;
                }
                if (Z_XMAX[k] < field_2D[k].first(x[j], y[j])) {
                    Z_XMAX[k] = field_2D[k].first(x[j], y[j]);
                    adjusted = true;
                }
                if (Z_YMIN[k] > field_2D[k].second(x[j], y[j])) {
                    Z_YMIN[k] = field_2D[k].second(x[j], y[j]);
                    adjusted = true;
                }
                if (Z_YMAX[k] < field_2D[k].second(x[j], y[j])) {
                    Z_YMAX[k] = field_2D[k].second(x[j], y[j]);
                    adjusted = true;
                }
            }
        }
        progress = (double)(i + 1) / N;
        std::cout << backspace << 100 * progress << "                    ";
        //std::cout << backspace << "Bounds adjusted: (" << Q1.transform_x(XMIN) << ", " << Q1.transform_x(XMAX) << "), (" << Q1.transform_y(YMIN) << ", " << Q1.transform_y(YMAX) << ") (" << std::round(10000 * (double)i / N) / 100 << "%)      ";
    }
    int weight_x_size = 0;
    int weight_z_size = 0;
    if (field_1D.size() > 0) weight_x_size = data_1D.size() / field_1D.size();
    if (field_2D.size() > 0)  weight_z_size = data_2D.size() / field_2D.size();
    Node2D n2; Node1D n1;
    n2.parent = -1;
    n2.depth = 1;
    n1.parent = -1;
    n1.depth = 1;
    for (int k = 0; k < field_1D.size(); k++) {
        for (int w = 0; w < weight_x_size; w++) {
            data_1D[w * field_1D.size() + k].set_scales(scale_1D[k]);
            data_1D[w * field_1D.size() + k].set_bounds(XMIN[k], XMAX[k]);
            data_1D[w * field_1D.size() + k].nodes.push_back(n1);
        }
    }
    for (int k = 0; k < field_2D.size(); k++) {
        for (int w = 0; w < weight_z_size; w++) {
            data_2D[w * field_2D.size() + k].set_scales(scale_2D[k].first, scale_2D[k].second);
            data_2D[w * field_2D.size() + k].set_bounds(Z_XMIN[k], Z_XMAX[k], Z_YMIN[k], Z_YMAX[k]);
            data_2D[w * field_2D.size() + k].nodes.push_back(n2);
        }
    }
}
/*
void initialize(std::string path_to_sim, std::vector<int> frames, std::string fieldname_x, std::string fieldname_y, std::vector<std::function<double(double, double)>> field_1D, std::vector<std::pair<std::function<double(double, double)>, std::function<double(double, double)>>> field_2D, std::vector<std::string> scale_1D, std::vector<std::pair<std::string, std::string>> scale_2D, std::vector<BinaryTree>& data_1D, std::vector<Quadtree>& data_2D) {
    std::vector<double> x, y;
    std::vector<double> XMIN, XMAX, Z_XMIN, Z_XMAX, Z_YMIN, Z_YMAX;
    for (int i = 0; i < field_1D.size(); i++) {
        XMIN.push_back(1e10);
        XMAX.push_back(-1e10);
    }
    for (int i = 0; i < field_2D.size(); i++) {
        Z_XMIN.push_back(1e10);
        Z_XMAX.push_back(-1e10);
        Z_YMIN.push_back(1e10);
        Z_YMAX.push_back(-1e10);
    }
    int N = frames.size();
    for (int i = 0; i < N; i++) {
	    int frame = frames[i];
        std::cout << "Adjusting bounds from frame " << frame << std::endl;
        auto cpu_files = get_cpu_list(path_to_sim, frame);
        for (int c = 0; c < 1024; c++) {
            std::cout << "Adjusting bounds from frame " << frame << ", cpu" << string_pad_left(c, 0, 4) << ", ";
            std::string filename = path_to_sim + "/DD" + string_pad_left(frame, 0, 4) + "/data" + string_pad_left(frame, 0, 4) + ".cpu" + string_pad_left(c, 0, 4);
            for (auto & groupname : cpu_files[c]) {
                std::cout << groupname;
                x = read_from_hdf5(filename, groupname, fieldname_x);
                y = read_from_hdf5(filename, groupname, fieldname_y);
                bool adjusted = false;
                for (int j = 0; j < x.size(); j++) {
                    for (int k = 0; k < field_1D.size(); k++) {
                        if (XMIN[k] > field_1D[k](x[j], y[j])) {
                            XMIN[k] = field_1D[k](x[j], y[j]);
                            adjusted = true;
                        }
                        if (XMAX[k] < field_1D[k](x[j], y[j])) {
                            XMAX[k] = field_1D[k](x[j], y[j]);
                            adjusted = true;
                        }
                    }
                    for (int k = 0; k < field_2D.size(); k++) {
                        if (Z_XMIN[k] > field_2D[k].first(x[j], y[j])) {
                            Z_XMIN[k] = field_2D[k].first(x[j], y[j]);
                            adjusted = true;
                        }
                        if (Z_XMAX[k] < field_2D[k].first(x[j], y[j])) {
                            Z_XMAX[k] = field_2D[k].first(x[j], y[j]);
                            adjusted = true;
                        }
                        if (Z_YMIN[k] > field_2D[k].second(x[j], y[j])) {
                            Z_YMIN[k] = field_2D[k].second(x[j], y[j]);
                            adjusted = true;
                        }
                        if (Z_YMAX[k] < field_2D[k].second(x[j], y[j])) {
                            Z_YMAX[k] = field_2D[k].second(x[j], y[j]);
                            adjusted = true;
                        }
                    }
                }
                if (adjusted) std::cout << ", adjusted!" << std::endl;
                else std::cout << std::endl;
            }
        }
    }
    int weight_x_size = 0;
    int weight_z_size = 0;
    if (field_1D.size() > 0) weight_x_size = data_1D.size() / field_1D.size();
    if (field_2D.size() > 0)  weight_z_size = data_2D.size() / field_2D.size();
    Node2D n2; Node1D n1;
    n2.parent = -1;
    n2.depth = 1;
    n1.parent = -1;
    n1.depth = 1;
    for (int k = 0; k < field_1D.size(); k++) {
        for (int w = 0; w < weight_x_size; w++) {
            data_1D[w * field_1D.size() + k].set_scales(scale_1D[k]);
            data_1D[w * field_1D.size() + k].set_bounds(XMIN[k], XMAX[k]);
            data_1D[w * field_1D.size() + k].nodes.push_back(n1);
        }
    }
    for (int k = 0; k < field_2D.size(); k++) {
        for (int w = 0; w < weight_z_size; w++) {
            data_2D[w * field_2D.size() + k].set_scales(scale_2D[k].first, scale_2D[k].second);
            data_2D[w * field_2D.size() + k].set_bounds(Z_XMIN[k], Z_XMAX[k], Z_YMIN[k], Z_YMAX[k]);
            data_2D[w * field_2D.size() + k].nodes.push_back(n2);
        }
    }
}
*/