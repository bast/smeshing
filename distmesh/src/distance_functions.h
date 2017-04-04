#ifndef DISTANCE_FUNCTIONS_H_INCLUDED
#define DISTANCE_FUNCTIONS_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

double dsegment(const double x0,
                const double y0,
                const double p1x,
                const double p1y,
                const double p2x,
                const double p2y);

double vdsegment(const int num_points,
                 const double ps_x[],
                 const double ps_y[],
                 const int num_vertices,
                 const double vs_x[],
                 const double vs_y[],
                       double distances[]);

#ifdef __cplusplus
}
#endif

#endif /* DISTANCE_FUNCTIONS_H_INCLUDED */
