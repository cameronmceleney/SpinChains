//
// Created by Cameron Aidan McEleney on 08/05/2023.
//

#include "../include/NewSolver.h"

#include <iostream>
#include <vector>
#include <functional>
#include <tuple>
#include <fstream>
#include <cmath>

using namespace std;

const int num_spins = 4000;
const double gyro = 28.8e9 * 2 * M_PI; // Gyromagnetic ratio
const double h0 = 0.1; // External field strength
const double hEx = 43.5;  // Exchange field strength
const double hp = 3e-3; // Pumping field strength

const bool useLLG = true;

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
    double mz_right = (position < num_spins - 1) ? mz[position + 1]: 0;
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

    string file_path = "D:/Data/2023-05-09/test.csv"; // Change this to your desired file path
    save_to_csv(file_path, xs, ys, zs);

}




double hX_old(int spin, double mx0LHS, double mx0RHS, double t0) {
    
    double hX0; // The effective field (H_eff) component acting upon each spin
    bool hasStaticDrive = false; // Whether the driving field is static
    double drivingFreq = 42.5e9 * 2 * M_PI;
    
    if (spin >= 1 && spin <= 200) {
        // The pulse of input energy will be restricted to being along the x-direction, and it will only be generated within the driving region
        if (hasStaticDrive)
            hX0 = hEx * mx0LHS + hEx * mx0RHS + hp;
        else if (!hasStaticDrive)
            hX0 = hEx * mx0LHS + hEx * mx0RHS + hp * cos(drivingFreq * t0);
    } else
        // All spins along x which are not within the driving region
        hX0 = hEx * mx0LHS + hEx * mx0RHS;
    
    return hX0;
}

double hY_old(double my0LHS, double my0RHS) {
    double hY0 = hEx * my0LHS + hEx * my0RHS;
    
    return hY0;
}

double hZ_old(double mz0LHS, double mz0RHS) {
    double hZ0 = hEx * mz0LHS + hEx * mz0RHS + h0;
    
    return hZ0;
}

double mx_old(double mx0MID, double my0MID, double mz0MID, double gilbert, double hX0, double hY0, double hZ0) {
    
    double mxK; // The magnetic moment component along the x-direction for the first stage of RK2
    
    if (useLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mxK = gyro * (- (gilbert * hY0 * mx0MID * my0MID) + hY0 * mz0MID - hZ0 * (my0MID + gilbert * mx0MID * mz0MID) + gilbert * hX0 * (pow(my0MID,2) + pow(mz0MID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mxK = -1.0 * gyro * (my0MID * hZ0 - mz0MID * hY0);
    }
    
    return mxK;
}

double my_old(double mx0MID, double my0MID, double mz0MID, double gilbert, double hX0, double hY0, double hZ0) {
    
    double myK;
    
    if (useLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        myK = gyro * (-(hX0 * mz0MID) + hZ0 * (mx0MID - gilbert * my0MID * mz0MID) + gilbert * (hY0 * pow(mx0MID,2) - hX0 * mx0MID * my0MID + hY0 * pow(mz0MID,2)));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        myK =        gyro * (mx0MID * hZ0 - mz0MID * hX0);
    }
    
    return myK;
}

double mz_old(double mx0MID, double my0MID, double mz0MID, double gilbert, double hX0, double hY0, double hZ0) {
    
    double mzK; 
    
    if (useLLG) {
        // The magnetic moment components' coupled equations (obtained from magDynamics equation) with the parameters for the first stage of RK2.
        mzK = gyro * (hX0 * my0MID + gilbert * hZ0 * (pow(mx0MID,2) + pow(my0MID,2)) - gilbert*hX0*mx0MID*mz0MID - hY0 * (mx0MID + gilbert * my0MID * mz0MID));
    } else {
        // The magnetic moment components' coupled equations (obtained from the torque equation) with the parameters for the first stage of RK2.
        mzK = -1.0 * gyro * (mx0MID * hY0 - my0MID * hX0);
    }
    
    return mzK;
}


void NewSolver::rework_old_solver() {
    
    double totalTime = 0;
    double stepsize = 1e-15;
    double stepsizeHalf = stepsize / 2;
    double gilbert = 1e-4;

    vector<double> mx0{0};
    vector<double> my0{0};
    vector<double> mz0{0};
    double _mxInit = 0.0, _myInit = 0.0, _mzInit = 1.0;
    vector<double> mxInitCond(num_spins, _mxInit), myInitCond(num_spins, _myInit), mzInitCond(num_spins, _mzInit);
    mx0.insert(mx0.end(), mxInitCond.begin(), mxInitCond.end());
    my0.insert(my0.end(), myInitCond.begin(), myInitCond.end());
    mz0.insert(mz0.end(), mzInitCond.begin(), mzInitCond.end());
    mx0.push_back(0);
    my0.push_back(0);
    mz0.push_back(0);
    
    for (int iteration = 0; iteration <= 7e5; iteration++) {

        double t0 = totalTime, t0HalfStep = totalTime + stepsizeHalf;

        // The estimate of the slope for the x/y/z-axis magnetic moment component at the midpoint; mx1 = mx0 + (h * k1 / 2) etc
        std::vector<double> mx1(num_spins + 2, 0), my1(num_spins + 2, 0), mz1(num_spins + 2, 0);

        // Excludes the 0th and last spins as they will always be zero-valued (end, pinned spins)
        for (int spin = 1; spin <= num_spins; spin++) {
            // RK2 Stage 1. Takes initial conditions as inputs.

            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the first stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx0MID = mx0[spin], mx0LHS = mx0[spinLHS], mx0RHS = mx0[spinRHS];
            double my0MID = my0[spin], my0LHS = my0[spinLHS], my0RHS = my0[spinRHS];
            double mz0MID = mz0[spin], mz0LHS = mz0[spinLHS], mz0RHS = mz0[spinRHS];

            double hX0 = hX_old(spin, mx0LHS, mx0RHS, t0);
            double hY0 = hY_old(my0LHS, my0RHS);
            double hZ0 = hZ_old(mz0LHS, mz0RHS);
            
            double mxK1 = mx_old(mx0MID, my0MID, mz0MID, gilbert, hX0, hY0, hZ0);
            double myK1 = my_old(mx0MID, my0MID, mz0MID, gilbert, hX0, hY0, hZ0);
            double mzK1 = mz_old(mx0MID, my0MID, mz0MID, gilbert, hX0, hY0, hZ0);
                    
            // Find (m0 + k1/2) for each spin, which is used in the next stage.
            mx1[spin] = mx0MID + stepsizeHalf * mxK1;
            my1[spin] = my0MID + stepsizeHalf * myK1;
            mz1[spin] = mz0MID + stepsizeHalf * mzK1;
        }
        // The estimations of the m-components values for the next iteration.
        std::vector<double> mx2(num_spins + 2,0), my2(num_spins + 2,0), mz2(num_spins + 2,0);

        for (int spin = 1; spin <= num_spins; spin++) {

            int spinLHS = spin - 1, spinRHS = spin + 1;

            // The m-components for the first stage for the: current spin site (MID); site to the left (LHS); site to the right (RHS)
            double mx1MID = mx1[spin], mx1LHS = mx1[spinLHS], mx1RHS = mx1[spinRHS];
            double my1MID = my1[spin], my1LHS = my1[spinLHS], my1RHS = my1[spinRHS];
            double mz1MID = mz1[spin], mz1LHS = mz1[spinLHS], mz1RHS = mz1[spinRHS];

            double hX1 = hX_old(spin, mx1LHS, mx1RHS, t0);
            double hY1 = hY_old(my1LHS, my1RHS);
            double hZ1 = hZ_old(mz1LHS, mz1RHS);
            
            double mxK2 = mx_old(mx1MID, my1MID, mz1MID, gilbert, hX1, hY1, hZ1);
            double myK2 = my_old(mx1MID, my1MID, mz1MID, gilbert, hX1, hY1, hZ1);
            double mzK2 = mz_old(mx1MID, my1MID, mz1MID, gilbert, hX1, hY1, hZ1);

            mx2[spin] = mx0[spin] + stepsize * mxK2;
            my2[spin] = my0[spin] + stepsize * myK2;
            mz2[spin] = mz0[spin] + stepsize * mzK2;
        }
        // Everything below here is part of the class method, but not the internal RK2 stage loops.

        /**
         * Removes (possibly) large arrays as they can lead to memory overloads later in main.cpp. Failing to clear
         * these between loop iterations sometimes led to incorrect values cropping up.
         */
        mx0.clear();
        my0.clear();
        mz0.clear();
        mx1.clear();
        my1.clear();
        mz1.clear();

        //Sets the final value of the current iteration of the loop to be the starting value of the next loop.
        mx0 = mx2;
        my0 = my2;
        mz0 = mz2;
        

        totalTime += stepsize;
    }

    string file_path = "D:/Data/2023-05-09/test.csv"; // Change this to your desired file path
    save_to_csv(file_path, mx0, my0, mz0);
}

