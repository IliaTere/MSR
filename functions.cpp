#include "all_includes.h"

double f_0(double, double) {
    return 1;
}

double f_1(double x, double) {
    return x;
}

double f_2(double, double y) {
    return y;
}

double f_3(double x, double y) {
    return x + y;
}

double f_4(double x, double y) {
    return sqrt(x*x + y*y);
}

double f_5(double x, double y) {
    return x*x + y*y;
}

double f_6(double x, double y) {
    return exp(x*x - y*y);
}

double f_7(double x, double y) {
    return 1. / (25*(x*x + y*y) + 1);
}

void Functions::select_f(int func_id) {
    if (func_id == 0) {
        f = f_0;
    } else if (func_id == 1) {
        f = f_1;
    } else if (func_id == 2) {
        f = f_2;
    } else if (func_id == 3) {
        f = f_3;
    } else if (func_id == 4) {
        f = f_4;
    } else if (func_id == 5) {
        f = f_5;
    } else if (func_id == 6) {
        f = f_6;
    } else if (func_id == 7) {
        f = f_7;
    }
}

void displayVector(int vectorSize, double* dataArray) {
    for (int idx = 0; idx < vectorSize; ++idx) {
        printf("%.6lf ", dataArray[idx]);
    }
    printf("\n");
}
