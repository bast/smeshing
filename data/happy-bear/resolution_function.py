def resolution_function(d_gp, d_pp, q):
    '''
    d_gp: distance from grid point
          to nearest polygon point
          (same units as polygon data)
    d_pp: distance from nearest polygon point
          to nearest polygon point across strait
          taking into account the view vector
          (same units as polygon data)
    q: dictionary of quantities at nearest polygon point
       keys are defined in the configuration file
       (in other words, you define what these mean)
    '''
    s = 0.9

    result = s*d_gp + d_pp/q['number of points across strait']

    if result < 1.0:
        result = 1.0

    return result
