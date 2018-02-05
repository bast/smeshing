// provides std::max
#include <algorithm>

// The resolution is expressed as distance using the same (arbitrary)
// units as used by the boundary and polygons - this means that
// larger resolution number means that points are farther apart.

// Resolution in point r is defined as min(f(r, p)), where the miminum
// is taken over all boundary points p for a particular boundary point p,
// f is given as f(r, p) = g(d(r, p)) + h(c_p).
// d(r, p) is the distance r to p and function h(c_p) depends on
// coefficients c_p of a boundary point p. The number of coefficients
// per point and their meaning can be freely specified and interpreted.

// Below you are asked to specify functions g and h.
// You have two restrictions:
// 1) You have to respect is that g should not decrease for an increasing d.
//    In other words, for an increasing distance the resolution should not
//    decrease.
// 2) The sum g + h should never become zero since the code will divide by
//    the distance.

// This function only depends on the distance to a boundary point but not
// on coefficients at the boundary point.
double g_function(const double distance)
{
    // this is to make sure we do not end up with zero distance
    // and then try to divide by zero later
    double result = std::max(0.5, distance);

    return result;
}

// The code will give you all coefficients for a point in h_function
// and then you can use and combine them freely.
double h_function(const double coefficients[])
{
    // in this example we simply return the first coefficient
    return coefficients[0];
}
