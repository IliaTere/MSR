#include "inc.h"
#include <algorithm>

// Renamed function: r1 → computeMaxTriangleError
double r1(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    const int gridSize = (nx+1) * (ny+1);
    int startIdx, endIdx;
    int rowIdx, colIdx;
    
    // Get thread's portion of work
    thread_rows(gridSize, p, k, startIdx, endIdx);
    
    // Initial error value
    double maxError = -1;
    
    // Process assigned elements
    for (int idx = startIdx; idx < endIdx; ++idx) {
        // Convert linear index to 2D coordinates
        l2ij(nx, ny, rowIdx, colIdx, idx);
        
        // Skip boundary elements
        if (rowIdx == nx || colIdx == ny) {
            continue;
        }

        // Calculate triangle interpolation points
        const double x1 = a + (rowIdx + 2.0/3.0) * hx;
        const double y1 = c + (colIdx + 1.0/3.0) * hy;
        const double x2 = a + (rowIdx + 1.0/3.0) * hx;
        const double y2 = c + (colIdx + 2.0/3.0) * hy;
        
        // Get node values for interpolation
        const double node1 = x[idx];
        const double node2 = x[idx + 1];
        const double node3 = x[idx + 1 + nx + 1];
        const double node4 = x[idx + nx + 1];
        
        // Calculate interpolation and error for first triangle
        const double exactVal1 = f(x1, y1);
        const double approxVal1 = (node1 + node2 + node3) / 3.0;
        const double error1 = fabs(exactVal1 - approxVal1);
        
        // Calculate interpolation and error for second triangle
        const double exactVal2 = f(x2, y2);
        const double approxVal2 = (node1 + node4 + node3) / 3.0;
        const double error2 = fabs(exactVal2 - approxVal2);
        
        // Update maximum error
        const double localMax = std::max(error1, error2);
        maxError = std::max(maxError, localMax);
    }

    // Combine results from all threads
    reduce_sum(p, &maxError, 1, &max);
    
    return maxError;
}

// Renamed function: r2 → calculateIntegratedTriangleError
double r2(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    const int gridSize = (nx+1) * (ny+1);
    int startIdx, endIdx;
    int rowIdx, colIdx;
    
    // Get thread's portion of work
    thread_rows(gridSize, p, k, startIdx, endIdx);
    
    // Accumulate error
    double errorSum = 0.0;
    
    // Process assigned elements
    for (int idx = startIdx; idx < endIdx; ++idx) {
        // Convert linear index to 2D coordinates
        l2ij(nx, ny, rowIdx, colIdx, idx);
        
        // Skip boundary elements
        if (rowIdx == nx || colIdx == ny) {
            continue;
        }

        // Calculate triangle sample points
        const double px1 = a + (rowIdx + 2.0/3.0) * hx;
        const double py1 = c + (colIdx + 1.0/3.0) * hy;
        const double px2 = a + (rowIdx + 1.0/3.0) * hx;
        const double py2 = c + (colIdx + 2.0/3.0) * hy;
        
        // Get nodal values
        const double valAtNode = x[idx];
        const double valAtRightNode = x[idx + 1];
        const double valAtDiagNode = x[idx + 1 + nx + 1];
        const double valAtTopNode = x[idx + nx + 1];
        
        // Calculate errors at sample points
        const double triangleError1 = fabs(f(px1, py1) - (valAtNode + valAtRightNode + valAtDiagNode) / 3.0);
        const double triangleError2 = fabs(f(px2, py2) - (valAtNode + valAtTopNode + valAtDiagNode) / 3.0);
        
        // Add to total error
        errorSum += triangleError1 + triangleError2;
    }

    // Combine results from all threads
    double totalError = reduce_sum_det(p, k, errorSum);
    
    // Scale by area element
    return (hx * hy * totalError) / 2.0;
}

// Renamed function: r3 → findMaxNodeError
double r3(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    const int gridSize = (nx+1) * (ny+1);
    int startIdx, endIdx;
    int rowIdx, colIdx;
    
    // Get thread's work range
    thread_rows(gridSize, p, k, startIdx, endIdx);
    
    // Track maximum error
    double nodeMaxError = -1.0;
    
    // Process each node in this thread's range
    for (int nodeIdx = startIdx; nodeIdx < endIdx; ++nodeIdx) {
        // Convert to grid coordinates
        l2ij(nx, ny, rowIdx, colIdx, nodeIdx);
        
        // Calculate exact and approximate values
        const double exactValue = f(a + rowIdx*hx, c + colIdx*hy);
        const double approxValue = x[nodeIdx];
        
        // Update maximum error
        const double nodeError = fabs(exactValue - approxValue);
        nodeMaxError = std::max(nodeMaxError, nodeError);
    }

    // Combine results across threads
    reduce_sum(p, &nodeMaxError, 1, &max);
    
    return nodeMaxError;
}

// Renamed function: r4 → computeIntegratedNodeError
double r4(int nx, int ny, double a, double c, double hx, double hy, double* x, double (*f)(double, double), int p, int k) {
    const int gridSize = (nx+1) * (ny+1);
    int startIdx, endIdx;
    int rowIdx, colIdx;
    
    // Divide work among threads
    thread_rows(gridSize, p, k, startIdx, endIdx);
    
    // Initialize error accumulator
    double errorAccumulator = 0.0;
    
    // Process nodes in thread's range
    for (int nodeIdx = startIdx; nodeIdx < endIdx; ++nodeIdx) {
        // Get 2D coordinates
        l2ij(nx, ny, rowIdx, colIdx, nodeIdx);
        
        // Calculate physical coordinates
        const double physX = a + rowIdx * hx;
        const double physY = c + colIdx * hy;
        
        // Add absolute error at this node
        errorAccumulator += fabs(f(physX, physY) - x[nodeIdx]);
    }

    // Combine partial sums from all threads
    double totalError = reduce_sum_det(p, k, errorAccumulator);
    
    // Scale by area element for integration
    return hx * hy * totalError;
}
