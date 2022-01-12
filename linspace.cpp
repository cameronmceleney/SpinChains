#include "linspace.h"

void linspace::setStart(double start) {
    // Start point of the interval
    start_in = start;
}

void linspace::setEnd(double end) {
    // Stop point of the interval
    stop_in = end;
}

void linspace::setNum(int num) {
    // 'Number' of evenly spaced samples to be returned over the given interval
    num_in = num;
}

std::vector<double> linspace::generate_array(void)
{
    std::vector<double> linspace_array;

    if (num_in == 0) {
        // if no range in inputted then empty vector is returned
        return linspace_array;
    }

    if (num_in == 1)
    {
        // if range is one then start value is the only element of the array returned
        linspace_array.push_back(start_in);
        return linspace_array;
    }

    // All other cases are dealt with here. Output will be of length 'num'

    // 'Delta' finds the stepsize between values
    double delta = (stop_in - start_in) / (num_in - 1);

    for(int i=0; i < num_in-1; ++i)
    {
        // Adds elements to list. Method is very similar to Euler's numerical method
        linspace_array.push_back(start_in + delta * i);
    }
    // Includes endpoint
    linspace_array.push_back(stop_in);
    return linspace_array;
}
