#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H

#include <stdio.h>
#include <pthread.h>
#include <stddef.h>
#include "math.h"
#include "parallel_utils.h"
#include "common_types.h"
#include <string>

void matrix_mult_vector_msr(int n, double* A, int* I, double* x, double* y, int p, int k);
int minimal_errors_msr_matrix(int n, double* A, int* I, double* b, double* x,
    double* r, double* u, double* v, double eps, int maxit, int p, int k);

int minimal_errors_msr_matrix_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, 
    double eps, int maxit, int maxsteps, int p, int k);

void ij2l(int nx, int, int i, int j, int& l);
void l2ij(int nx, int, int& i, int& j, int l);
int get_len_msr(int nx, int ny);
int get_off_diag(int nx, int ny, int i, int j, int* I_ij = nullptr);
int get_len_msr_off_diag(int nx, int ny);
int allocate_msr_matrix(int nx, int ny, double** p_A, int** p_I);
void fill_I(int nx, int ny, int* I);
void fill_A_ij(int nx, int ny, double hx, double hy, int i, int j, double* A_diag, double* A_off_diag);
void fill_A(int nx, int ny, double hx, double hy, int* I, double* A, int p, int k);
int check_symm(int nx, int ny, int* I, double* A, double eps, int p, int k);

void apply_preconditioner_msr_matrix(int n, double* A, int* I, double* v1, double* v2, int flag, int p, int k);
void solve_rsystem(int n, int* I, double* U, double* b, double* x, double w, int p, int k);
void solve_lsystem(int n, int* I, double* U, double* b, double* x, double w, int p, int k);
bool step(int n, double* A, int* I, double* x, double* r, double* u, double* v, double prec, int p, int k);

void displayVector(int vectorSize, double* dataArray);
bool isNumber(std::string& str);

#endif // MATRIX_OPERATIONS_H 