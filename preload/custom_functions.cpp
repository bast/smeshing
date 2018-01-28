double h_function(const double coefficients[])
{
    return 0.1*coefficients[0];
}

double g_function(const double distance)
{
    double scale_factor = 0.99;
    return scale_factor * distance;
}
