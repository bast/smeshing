#ifndef DISTANCE_FUNCTIONS_H_INCLUDED
#define DISTANCE_FUNCTIONS_H_INCLUDED

double dellipse(double x0, double y0, double a, double b);

double dellipsoid(double x0, double y0, double z0, double a, double b, double c);

double dsegment(double x0, double y0, double p1x, double p1y, double p2x, double p2y);

void dhseqr_(char*,char*,int*,int*,int*,double*,int*,double*,
             double*,void*,int*,double*,void*,int*);

#endif /* DISTANCE_FUNCTIONS_H_INCLUDED */
