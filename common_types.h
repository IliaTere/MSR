#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <stdio.h>

enum class Status {
    success,
    error
};

struct Args{
    double a;
    double b;
    double c;
    double d;
    double eps;
    int* I;
    double* A;
    double* B;
    double* x;
    double* r;
    double* u;
    double* v;
    int nx;
    int ny;
    int maxit;
    int p;
    int k;
    double (*f)(double, double);
    int its = 0;
    double t1 = 0;
    double t2 = 0;
    double res_1 = 0;
    double res_2 = 0;
    double res_3 = 0;
    double res_4 = 0;
    Status status = Status::success;
};

#endif // COMMON_TYPES_H 