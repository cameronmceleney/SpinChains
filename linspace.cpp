#include "linspace.h"
linspace::linspace(double start, double end, int num)
{
    start_in = start;
    end_in = end;
    num_in = num;
}

std::vector<double> linspace::linspace(double start, double end, int num)
{
    std::vector<double> linspace;


    if (num_in == 0) {
        return linspace;
    }

    if (num_in == 1)
    {
        linspace.push_back(start_in);
        return linspace;
    }

    double delta = (end_in - start_in) / (num_in - 1);

    for(int i=0; i < num_in-1; ++i)
    {
        linspace.push_back(start_in + delta * i);
    }
    linspace.push_back(end_in); // I want to ensure that start and end
    // are exactly the same as the input
    return linspace;
}
