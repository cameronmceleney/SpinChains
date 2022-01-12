
#ifndef SPINCHAINS_LINSPACE_H
#define SPINCHAINS_LINSPACE_H

#include <vector>

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

    std::vector<double> linspace;

    auto start = static_cast<double>(start_in);
    auto end = static_cast<double>(end_in);
    auto num = static_cast<double>(num_in);

    if (num == 0) {
        return linspace;
    }

    if (num == 1)
    {
        linspace.push_back(start);
        return linspace;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspace.push_back(start + delta * i);
    }
    linspace.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspace;
}
#endif //SPINCHAINS_LINSPACE_H
