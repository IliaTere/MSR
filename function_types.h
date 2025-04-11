#ifndef FUNCTION_TYPES_H
#define FUNCTION_TYPES_H

#include "common_types.h"

double F_IJ(int nx, int ny, double hx, double hy, double a, double с, int i, int j, double (*f)(double, double));
void fill_B(int nx, int ny, double hx, double hy, double a, double c, double* B, double (*f)(double, double), int p, int k);
double f_0(double, double);
double f_1(double x, double);
double f_2(double, double y);
double f_3(double x, double y);
double f_4(double x, double y);
double f_5(double x, double y);
double f_6(double x, double y);
double f_7(double x, double y);

class Functions {
public:
    double (*f)(double, double);
    void select_f(int func_id);
};

double r1(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k);
double r2(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k);
double r3(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k);
double r4(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k);

#endif // FUNCTION_TYPES_H 