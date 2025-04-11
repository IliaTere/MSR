#include "all_includes.h"

static double* results = nullptr;

int init_reduce_sum(int p) {
    results = new double[p];
    if (results == nullptr) {
        return -1;
    }

    return 0;
}

double reduce_sum_det(int p, int k, double s) {
    double sum = 0; int l;
    results[k] = s;
    reduce_sum<int>(p);
    for (l = 0; l < p; ++l) {
        sum += results[l];
    }

    reduce_sum<int>(p);
    return sum;
}

void free_results() {
    delete[] results;
}
