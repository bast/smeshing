#ifndef SMESHING_H_INCLUDED
#define SMESHING_H_INCLUDED

#ifndef SMESHING_API
#include "smeshing_export.h"
#define SMESHING_API smeshing_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

SMESHING_API
void get_resolution(const int num_points,
                    const double x[],
                    const double y[],
                          double resolutions[],
                    const bool use_tanh,
                    const int num_reference_points,
                    const double reference_x[],
                    const double reference_y[],
                    const double nearest_distance_at_coastline_point[],
                    const int flanders_indices[]);

#ifdef __cplusplus
}
#endif

#endif /* SMESHING_H_INCLUDED */
