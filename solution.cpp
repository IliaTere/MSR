#include <sys/sysinfo.h>
#include "all_includes.h"

void* solution(void* ptr) {
    Args* args = (Args*)ptr;
    double a = args->a; double b = args->b; double c = args->c; double d = args->d;
    double eps = args->eps; int* I = args->I; double* A = args->A; double* B = args->B;
    double* x = args->x; double* r = args->r; double* u = args->u; double* v = args->v;
    int maxit = args->maxit; int nx = args->nx; int ny = args->ny;
    int p = args->p; int k = args->k; double (*f)(double, double) = args->f;

    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (k % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();

    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    double hx = (b - a) / nx;
    double hy = (d - c) / ny;
    int N = (nx + 1) * (ny + 1);
    
    fill_A(nx, ny, hx, hy, I, A, p, k);
    fill_B(nx, ny, hx, hy, a, c, B, f, p, k); 

    int maxsteps = 300; // гиперпараметр
    args->t1 = get_cpu_time();
    // Использование метода минимальных ошибок с блочным симметричным предобуславливателем Зейделя
    int its = minimal_errors_msr_matrix_full(N, A, I, B, x, r, u, v, eps, maxit, maxsteps, p, k);  
    args->t1 = get_cpu_time() - args->t1;
    args->its = its;

    args->t2 = get_cpu_time();
    double res_1 = r1(nx, ny, a, c, hx, hy, x, f, p, k);
    double res_2 = r2(nx, ny, a, c, hx, hy, x, f, p, k);
    double res_3 = r3(nx, ny, a, c, hx, hy, x, f, p, k);
    double res_4 = r4(nx, ny, a, c, hx, hy, x, f, p, k);
    args->t2 = get_cpu_time() - args->t2;

    args->res_1 = res_1;
    args->res_2 = res_2;
    args->res_3 = res_3;
    args->res_4 = res_4;

    reduce_sum<int>(p);
    return nullptr;
}
