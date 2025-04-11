#ifndef PARALLEL_UTILS_H
#define PARALLEL_UTILS_H

#include <pthread.h>
#include <stddef.h>
#include <algorithm>

void thread_rows(int n, int p, int k, int& i1, int& i2);
double scalar_product(int n, double* x, double* y, int p, int k);
void mult_sub_vector(int n, double* x, double* y, double tau, int p, int k);
void* solution(void* ptr);
double get_cpu_time();

int init_reduce_sum(int p);
double reduce_sum_det(int p, int k, double s);
void free_results();

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

#endif // PARALLEL_UTILS_H 