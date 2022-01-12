#include "linspace.h"

void linspace::set_values(double start, double stop, int num, bool endpoint){

    // Start point of the interval.
    start_in = start;

    // Stop point of the interval.
    stop_in = stop;

    // 'Number' of evenly spaced samples to be returned over the given interval.
    num_in = num;

    //If true, stop is the last sample. Otherwise, it is not included. Default is true.
    endpoint_in = endpoint;
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

    // All other cases are dealt with here. Output will be of length 'num' if endpoint = True

    // 'Delta' finds the stepsize between values
    double delta = (stop_in - start_in) / (num_in - 1);

    for(int i=0; i < num_in-1; ++i)
    {
        // Adds elements to list. Method is very similar to Euler's numerical method
        linspace_array.push_back(start_in + delta * i);
    }
    // If true, stop is the last sample; otherwise it's excluded
    if (endpoint_in == true){
        linspace_array.push_back(stop_in);
    }

    return linspace_array;
}
