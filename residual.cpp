#include "all_includes.h"
#include <algorithm>

double r1(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    const int gridSize = (nx+1) * (ny+1);
    int startIdx, endIdx;
    int rowIdx, colIdx;
    
    thread_rows(gridSize, p, k, startIdx, endIdx);
    
    double maxError = -1;
    
    for (int idx = startIdx; idx < endIdx; ++idx) {
        l2ij(nx, ny, rowIdx, colIdx, idx);

        if (rowIdx == nx || colIdx == ny) {
            continue;
        }

        const double x1 = a + (rowIdx + 2.0/3.0) * hx;
        const double y1 = c + (colIdx + 1.0/3.0) * hy;
        const double x2 = a + (rowIdx + 1.0/3.0) * hx;
        const double y2 = c + (colIdx + 2.0/3.0) * hy;
        
        const double node1 = x[idx];
        const double node2 = x[idx + 1];
        const double node3 = x[idx + 1 + nx + 1];
        const double node4 = x[idx + nx + 1];
        
        const double exactVal1 = f(x1, y1);
        const double approxVal1 = (node1 + node2 + node3) / 3.0;
        const double error1 = fabs(exactVal1 - approxVal1);
        
        const double exactVal2 = f(x2, y2);
        const double approxVal2 = (node1 + node4 + node3) / 3.0;
        const double error2 = fabs(exactVal2 - approxVal2);
        
        const double localMax = std::max(error1, error2);
        maxError = std::max(maxError, localMax);
    }

    reduce_sum(p, &maxError, 1, &max);
    
    return maxError;
}

double r2(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    const int gridSize = (nx+1) * (ny+1);
    int startIdx, endIdx;
    int rowIdx, colIdx;
    
    thread_rows(gridSize, p, k, startIdx, endIdx);
    
    double errorSum = 0.0;
    
    for (int idx = startIdx; idx < endIdx; ++idx) {
        l2ij(nx, ny, rowIdx, colIdx, idx);
        
        if (rowIdx == nx || colIdx == ny) {
            continue;
        }

        const double px1 = a + (rowIdx + 2.0/3.0) * hx;
        const double py1 = c + (colIdx + 1.0/3.0) * hy;
        const double px2 = a + (rowIdx + 1.0/3.0) * hx;
        const double py2 = c + (colIdx + 2.0/3.0) * hy;
        
        const double valAtNode = x[idx];
        const double valAtRightNode = x[idx + 1];
        const double valAtDiagNode = x[idx + 1 + nx + 1];
        const double valAtTopNode = x[idx + nx + 1];
        
        const double triangleError1 = fabs(f(px1, py1) - (valAtNode + valAtRightNode + valAtDiagNode) / 3.0);
        const double triangleError2 = fabs(f(px2, py2) - (valAtNode + valAtTopNode + valAtDiagNode) / 3.0);
        
        errorSum += triangleError1 + triangleError2;
    }

    // Combine results from all threads
    double totalError = reduce_sum_det(p, k, errorSum);
    return (hx * hy * totalError) / 2.0;
}

double r3(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    const int gridSize = (nx+1) * (ny+1);
    int startIdx, endIdx;
    int rowIdx, colIdx;
    
    thread_rows(gridSize, p, k, startIdx, endIdx);
    
    double nodeMaxError = -1.0;
    
    for (int nodeIdx = startIdx; nodeIdx < endIdx; ++nodeIdx) {
        l2ij(nx, ny, rowIdx, colIdx, nodeIdx);
        
        const double exactValue = f(a + rowIdx*hx, c + colIdx*hy);
        const double approxValue = x[nodeIdx];
        
        const double nodeError = fabs(exactValue - approxValue);
        nodeMaxError = std::max(nodeMaxError, nodeError);
    }

    reduce_sum(p, &nodeMaxError, 1, &max);
    
    return nodeMaxError;
}

double r4(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    const int gridSize = (nx+1) * (ny+1);
    int startIdx, endIdx;
    int rowIdx, colIdx;
    
    thread_rows(gridSize, p, k, startIdx, endIdx);
    double errorAccumulator = 0.0;
    
    for (int nodeIdx = startIdx; nodeIdx < endIdx; ++nodeIdx) {
        l2ij(nx, ny, rowIdx, colIdx, nodeIdx);
        
        const double physX = a + rowIdx * hx;
        const double physY = c + colIdx * hy;
        
        errorAccumulator += fabs(f(physX, physY) - x[nodeIdx]);
    }

    double totalError = reduce_sum_det(p, k, errorAccumulator);
    
    return hx * hy * totalError;
}
