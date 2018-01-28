double h_function(const double coefficients[])
{
    return coefficients[0] + 0.5*coefficients[1];
}

double g_function(const double distance)
{
    double scale_factor = 0.99;
    return scale_factor * distance;
}
