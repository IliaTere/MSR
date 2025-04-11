#include <stdio.h>
#include <pthread.h>
#include <stddef.h>
#include "math.h"

void matrix_mult_vector_msr(int n, double* A, int* I, double* x, double* y, int p, int k);
int minimal_errors_msr_matrix(int n, double* A, int* I, double* b, double* x,
    double* r, double* u, double* v, double eps, int maxit, int p, int k);

int minimal_errors_msr_matrix_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, 
    double eps, int maxit, int maxsteps, int p, int k);

void thread_rows(int n, int p, int k, int& i1, int& i2);
double scalar_product(int n, double* x, double* y, int p, int k);
void mult_sub_vector(int n, double* x, double* y, double tau, int p, int k);
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

void apply_preconditioner_msr_matrix(int n, double* A, int* I, double* v1, double* v2, int flag, int p, int k);
void solve_rsystem(int n, int* I, double* U, double* b, double* x, double w, int p, int k);
void solve_lsystem(int n, int* I, double* U, double* b, double* x, double w, int p, int k);
bool step(int n, double* A, int* I, double* x, double* r, double* u, double* v, double prec, int p, int k);
void* solution(void* ptr);
double get_cpu_time();

int init_reduce_sum(int p);
double reduce_sum_det(int p, int k, double s);
void free_results();

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

template<class T>
void sum(T* r, T* a, int n) {
    for (int i = 0; i < n; ++i) {
        r[i] += a[i];
    }    
}

template<class T>
void max(T* r, T* a, int n) {
    for (int i = 0; i < n; ++i) {
        r[i] = std::max(r[i], a[i]);
    }  
}

template<class T>
void reduce_sum(int p, T* a = nullptr, int n = 0, void (*func)(T*, T*, int) = &sum) {
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static T* r = nullptr;
    int i;

    if (p <= 1) {
        return;
    }

    pthread_mutex_lock(&m);

    if (r == nullptr) {
        r = a;
    } else {
        func(r, a, n);
    }

    t_in++;
    if (t_in >= p) {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
        
    } else {
        while (t_in < p) {
            pthread_cond_wait(&c_in, &m);
        }
    }

    if (r != a) {
        for (i = 0; i < n; ++i) {
            a[i] = r[i];
        }
    } 

    t_out++;
    if (t_out >= p) {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);        
    } else {
        while (t_out < p) {
            pthread_cond_wait(&c_out, &m);
        }
    }

    pthread_mutex_unlock(&m);
}
