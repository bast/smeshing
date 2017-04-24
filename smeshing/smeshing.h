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
double get_resolution(const double x,
                      const double y,
                      const bool use_tanh,
                      const int num_points,
                      const double points_x[],
                      const double points_y[],
                      const int flanders_indices[]);

#ifdef __cplusplus
}
#endif

#endif /* SMESHING_H_INCLUDED */
