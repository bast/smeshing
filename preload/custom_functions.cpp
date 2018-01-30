// provides std::max
#include <algorithm>

double h_function(const double coefficients[])
{
    // in this example this is the distance across strait
    // using view angle of 90 degs
    return 0.5*coefficients[0];
}

double g_function(const double distance)
{
    // this is to make sure we do not end up with zero distance
    // and then try to divide by zero later
    double result = std::max(0.5, distance);

    return result;
}
