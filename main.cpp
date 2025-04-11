#include <string>
#include <cstring>
#include "all_includes.h"
#include <fenv.h>
#include <iostream>
#include <stdexcept>


int main(int argc, char* argv[]) {
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    
    
    if (argc != 11) {
        std::cerr << "Error: Expected 10 command-line arguments." << std::endl;
        std::cerr << "Usage: " << argv[0] << " a b c d nx ny k epsilon max_iterations threads" << std::endl;
        return 1;
    }

    double a, b, c, d, eps;
    int nx, ny, k, max_its, p;
    
    try {
        a = std::stod(argv[1]);
        b = std::stod(argv[2]);
        c = std::stod(argv[3]);
        d = std::stod(argv[4]);
        nx = std::stoi(argv[5]);
        ny = std::stoi(argv[6]);
        k = std::stoi(argv[7]);
        eps = std::stod(argv[8]);
        max_its = std::stoi(argv[9]);
        p = std::stoi(argv[10]);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Error: Invalid argument format. All parameters must be valid numbers." << std::endl;
        std::cerr << "Usage: " << argv[0] << " a b c d nx ny k epsilon max_iterations threads" << std::endl;
        return 1;
    } catch (const std::out_of_range& e) {
        std::cerr << "Error: Number out of range." << std::endl;
        return 1;
    }
    
    int* I = nullptr;
    double* A = nullptr;
    if (allocate_msr_matrix(nx, ny, &A, &I)) { 
        std::cerr << "Error: Failed to allocate MSR matrix." << std::endl;
        return 2; 
    }
    
    init_reduce_sum(p);
    
    int n = (nx + 1) * (ny + 1);
    
    double* B = new double[n];
    double* x = new double[n];
    double* r = new double[n];
    double* u = new double[n];
    double* v = new double[n];

    fill_I(nx, ny, I);

    memset(x, 0, n * sizeof(double));

    Functions func;
    func.select_f(k);
    double (*f)(double, double) = func.f;

    Args* args = new Args[p];
    pthread_t* threads = new pthread_t[p];
        
    for (int i = 1; i < p; ++i) {
        args[i].a = a;
        args[i].b = b;
        args[i].c = c;
        args[i].d = d;
        args[i].eps = eps;
        args[i].I = I;
        args[i].A = A;
        args[i].B = B;
        args[i].x = x;
        args[i].r = r;
        args[i].u = u;
        args[i].v = v;
        args[i].nx = nx;
        args[i].ny = ny;
        args[i].maxit = max_its;
        args[i].p = p;
        args[i].k = i;
        args[i].f = f;

        pthread_create(&threads[i], nullptr, &::solution, &args[i]); 
    }

    args[0].a = a;
    args[0].b = b;
    args[0].c = c;
    args[0].d = d;
    args[0].eps = eps;
    args[0].I = I;
    args[0].A = A;
    args[0].B = B;
    args[0].x = x;
    args[0].r = r;
    args[0].u = u;
    args[0].v = v;
    args[0].nx = nx;
    args[0].ny = ny;
    args[0].maxit = max_its;
    args[0].p = p;
    args[0].k = 0;
    args[0].f = f;
    
    ::solution(&args[0]);

    for (int i = 1; i < p; ++i) {
        pthread_join(threads[i], nullptr);
    }    

    int its = args[0].its;
    double r1 = args[0].res_1;
    double r2 = args[0].res_2;
    double r3 = args[0].res_3;
    double r4 = args[0].res_4;
    double t1 = args[0].t1;
    double t2 = args[0].t2;

    const int task = 6;

    printf(
        "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\n"
        "      It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
        argv[0], task, 
        r1, r2, r3, r4, 
        t1, t2, 
        its, eps, k, 
        nx, ny, p);


    free_results();
    delete[] I;
    delete[] A;
    delete[] B;
    delete[] x;
    delete[] r;
    delete[] u;
    delete[] v;
    delete[] args;
    delete[] threads;

    return 0;
}
