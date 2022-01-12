#include "linspace.h"

void linspace::setStart(double start) {
    start_in = start;
}

void linspace::setEnd(double end) {
    end_in = end;
}

void linspace::setNum(int num) {
    num_in = num;
}

std::vector<double> linspace::generate_array(void)
{
    std::vector<double> linspace_array;


    if (num_in == 0) {
        return linspace_array;
    }

    if (num_in == 1)
    {
        linspace_array.push_back(start_in);
        return linspace_array;
    }

    double delta = (end_in - start_in) / (num_in - 1);

    for(int i=0; i < num_in-1; ++i)
    {
        linspace_array.push_back(start_in + delta * i);
    }
    linspace_array.push_back(end_in); // I want to ensure that start and end
    // are exactly the same as the input
    return linspace_array;
}
