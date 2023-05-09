//
// Created by Cameron Aidan McEleney on 08/05/2023.
//

#include "NewSolver.h"

#include <iostream>
#include <vector>
#include <functional>
#include <tuple>
#include <fstream>

using namespace std;

const int num_spins = 100;
const double gyro = 1; // Gyromagnetic ratio
const double h0 = 0.1; // External field strength
const double hEx = 43.5;  // Exchange field strength
const double hp = 3e-3; // Pumping field strength

double hX(const vector<double>& mx, int position) {
    double mx_left = (position > 0) ? mx[position - 1] : 0;
    double mx_right = (position < num_spins - 1) ? mx[position + 1] : 0;
    return hEx * mx_left + hEx * mx_right + hp;
}

double hY(const vector<double>& my, int position) {
    double my_left = (position > 0) ? my[position - 1] : 0;
    double my_right = (position < num_spins - 1) ? my[position + 1] : 0;
    return hEx * my_left + hEx * my_right;
}

double hZ(const vector<double>& mz, int position) {
    double mz_left = (position > 0) ? mz[position - 1] : 0;
    double mz_right = (position < num_spins - 1) ? mz[position + 1] : 0;
    return hEx * mz_left + hEx * mz_right + h0;
}

double mx(const vector<double>& mx, const vector<double>& my, const vector<double>& mz, int position) {
    return -1.0 * gyro * (my[position] * hZ(mz, position) - mz[position] * hY(my, position));
}

double my(const vector<double>& mx, const vector<double>& my, const vector<double>& mz, int position) {
    return gyro * (mx[position] * hZ(mz, position) - mz[position] * hX(mx, position));
}

double mz(const vector<double>& mx, const vector<double>& my, const vector<double>& mz, int position) {
    return -1.0 * gyro * (mx[position] * hY(my, position) - my[position] * hX(mx, position));
}


tuple<double, double, double> RK2(double t, double h, int position,
    const vector<double>& mx0, const vector<double>& my0, const vector<double>& mz0,
    function<double(const vector<double>&, const vector<double>&, const vector<double>&, int)> fx,
    function<double(const vector<double>&, const vector<double>&, const vector<double>&, int)> fy,
    function<double(const vector<double>&, const vector<double>&, const vector<double>&, int)> fz) {

    double k1_x = h * fx(mx0, my0, mz0, position);
    double k1_y = h * fy(mx0, my0, mz0, position);
    double k1_z = h * fz(mx0, my0, mz0, position);

    vector<double> mx1 = mx0, my1 = my0, mz1 = mz0;
    mx1[position] += k1_x / 2;
    my1[position] += k1_y / 2;
    mz1[position] += k1_z / 2;

    double k2_x = h * fx(mx1, my1, mz1, position);
    double k2_y = h * fy(mx1, my1, mz1, position);
    double k2_z = h * fz(mx1, my1, mz1, position);

    double mx_total = mx0[position] + k2_x;
    double my_total = my0[position] + k2_y;
    double mz_total = mz0[position] + k2_z;

    return make_tuple(mx_total, my_total, mz_total);
}


tuple<vector<double>, vector<double>, vector<double>>
solveODE(vector<double> mx_init, vector<double> my_init, vector<double> mz_init, double t, double h,
         function<double(const vector<double>&, const vector<double>&, const vector<double>&, int)> fx,
         function<double(const vector<double>&, const vector<double>&, const vector<double>&, int)> fy,
         function<double(const vector<double>&, const vector<double>&, const vector<double>&, int)> fz) {

    for (int spin = 0; spin < num_spins; spin++) {
        tie(mx_init[spin + 1], my_init[spin + 1], mz_init[spin + 1]) = RK2(t, h, spin, mx_init, my_init, mz_init,
                                                                        fx, fy, fz);
    }

    return make_tuple(mx_init, my_init, mz_init);
}

void save_to_csv(const string& file_path,
                 const vector<double>& xs,
                 const vector<double>& ys,
                 const vector<double>& zs) {
    ofstream file(file_path);

    if (file.is_open()) {
        file << "Spin,mx,my,mz\n";
        for (int spin = 0; spin < num_spins; ++spin) {
            file << spin + 1 << ",";
            file << xs[spin] << ",";
            file << ys[spin] << ",";
            file << zs[spin] << "\n";
        }
        file.close();
    } else {
        cerr << "Unable to open file";
    }
}

void NewSolver::solver_main() {
    double t0 = 0;
    double tmax = 0.7e-9;
    double stepsize = 1e-15;
    double steps = tmax / stepsize;

    function<double(const vector<double>&, const vector<double>&, const vector<double>&, int)> mx_func = mx;
    function<double(const vector<double>&, const vector<double>&, const vector<double>&, int)> my_func = my;
    function<double(const vector<double>&, const vector<double>&, const vector<double>&, int)> mz_func = mz;

    // Initialize magnetization vectors for all spins
    vector<double> mx_init(num_spins, 0.01);
    vector<double> my_init(num_spins, 0.01);
    vector<double> mz_init(num_spins, 0.98);

    vector<double> xs(num_spins);
    vector<double> ys(num_spins);
    vector<double> zs(num_spins);

    for (int i = 0; i < steps; i++) {
        double t = t0 + i * stepsize;
        auto [spin_xs, spin_ys, spin_zs] = solveODE(mx_init, my_init, mz_init, t, stepsize,
                                                     mx_func, my_func, mz_func);

        xs = spin_xs;
        ys = spin_ys;
        zs = spin_zs;
    }

    string file_path = "/Users/cameronmceleney/CLionProjects/Data/2023-05-08/magnetization.csv"; // Change this to your desired file path
    save_to_csv(file_path, xs, ys, zs);

}

