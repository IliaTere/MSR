#include <string>
#include <cstring>
#include "inc.h"
#include <fenv.h>
#include <iostream>
#include <stdexcept>

// Display vector values for debugging
void displayVector(int vectorSize, double* dataArray) {
    for (int idx = 0; idx < vectorSize; ++idx) {
        printf("%.6lf ", dataArray[idx]);
    }
    printf("\n");
}

int main(int argc, char* argv[]) {
    // Enable floating-point exceptions
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    
    
    // Check for correct number of arguments
    if (argc != 11) {
        std::cerr << "Error: Expected 10 command-line arguments." << std::endl;
        std::cerr << "Usage: " << argv[0] << " a b c d nx ny k epsilon max_iterations threads" << std::endl;
        return 1;
    }
    
    // Parse command line arguments directly
    double domainStartX = std::stod(argv[1]);
    double domainEndX = std::stod(argv[2]);
    double domainStartY = std::stod(argv[3]);
    double domainEndY = std::stod(argv[4]);
    int gridPointsX = std::stoi(argv[5]);
    int gridPointsY = std::stoi(argv[6]);
    int functionIdentifier = std::stoi(argv[7]);
    double solverTolerance = std::stod(argv[8]);
    int maxSolverIterations = std::stoi(argv[9]);
    int numThreads = std::stoi(argv[10]);
    
    // Allocate and initialize MSR matrix structure
    int* rowIndices = nullptr;
    double* matrixValues = nullptr;
    if (allocate_msr_matrix(gridPointsX, gridPointsY, &matrixValues, &rowIndices)) { 
        std::cerr << "Error: Failed to allocate MSR matrix." << std::endl;
        return 2; 
    }
    
    // Initialize reduction operations for parallel execution
    init_reduce_sum(numThreads);
    
    // Calculate total grid size
    int totalGridPoints = (gridPointsX + 1) * (gridPointsY + 1);
    
    // Allocate arrays for the solution process
    double* rightHandSide = new double[totalGridPoints];
    double* solutionVector = new double[totalGridPoints];
    double* residual = new double[totalGridPoints];
    double* auxVector1 = new double[totalGridPoints];
    double* auxVector2 = new double[totalGridPoints];

    // Initialize matrix structure
    fill_I(gridPointsX, gridPointsY, rowIndices);

    // Initialize solution vector to zeros
    memset(solutionVector, 0, totalGridPoints * sizeof(double));

    // Set up function to approximate
    Functions functionProvider;
    functionProvider.select_f(functionIdentifier);
    double (*targetFunction)(double, double) = functionProvider.f;

    // Initialize thread arguments and thread IDs
    Args* threadArgs = new Args[numThreads];
    pthread_t* threadIds = new pthread_t[numThreads];
        
    // Create worker threads (skipping the first one for main thread)
    for (int threadIdx = 1; threadIdx < numThreads; ++threadIdx) {
        // Fill thread arguments
        threadArgs[threadIdx].a = domainStartX;
        threadArgs[threadIdx].b = domainEndX;
        threadArgs[threadIdx].c = domainStartY;
        threadArgs[threadIdx].d = domainEndY;
        threadArgs[threadIdx].eps = solverTolerance;
        threadArgs[threadIdx].I = rowIndices;
        threadArgs[threadIdx].A = matrixValues;
        threadArgs[threadIdx].B = rightHandSide;
        threadArgs[threadIdx].x = solutionVector;
        threadArgs[threadIdx].r = residual;
        threadArgs[threadIdx].u = auxVector1;
        threadArgs[threadIdx].v = auxVector2;
        threadArgs[threadIdx].nx = gridPointsX;
        threadArgs[threadIdx].ny = gridPointsY;
        threadArgs[threadIdx].maxit = maxSolverIterations;
        threadArgs[threadIdx].p = numThreads;
        threadArgs[threadIdx].k = threadIdx;
        threadArgs[threadIdx].f = targetFunction;

        // Create the thread - using the solution function from inc.h
        pthread_create(&threadIds[threadIdx], nullptr, &::solution, &threadArgs[threadIdx]); 
    }

    // Configure main thread arguments
    threadArgs[0].a = domainStartX;
    threadArgs[0].b = domainEndX;
    threadArgs[0].c = domainStartY;
    threadArgs[0].d = domainEndY;
    threadArgs[0].eps = solverTolerance;
    threadArgs[0].I = rowIndices;
    threadArgs[0].A = matrixValues;
    threadArgs[0].B = rightHandSide;
    threadArgs[0].x = solutionVector;
    threadArgs[0].r = residual;
    threadArgs[0].u = auxVector1;
    threadArgs[0].v = auxVector2;
    threadArgs[0].nx = gridPointsX;
    threadArgs[0].ny = gridPointsY;
    threadArgs[0].maxit = maxSolverIterations;
    threadArgs[0].p = numThreads;
    threadArgs[0].k = 0;
    threadArgs[0].f = targetFunction;
    
    // Execute computation in the main thread
    ::solution(&threadArgs[0]);

    // Wait for all worker threads to complete
    for (int threadIdx = 1; threadIdx < numThreads; ++threadIdx) {
        pthread_join(threadIds[threadIdx], nullptr);
    }    

    // Extract results from the main thread's arguments
    int iterationCount = threadArgs[0].its;
    double residualNorm1 = threadArgs[0].res_1;
    double residualNorm2 = threadArgs[0].res_2;
    double residualNorm3 = threadArgs[0].res_3;
    double residualNorm4 = threadArgs[0].res_4;
    double setupTime = threadArgs[0].t1;
    double solveTime = threadArgs[0].t2;

    const int task = 6;

    printf(
        "%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\n"
        "      It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
        argv[0], task, 
        residualNorm1, residualNorm2, residualNorm3, residualNorm4, 
        setupTime, solveTime, 
        iterationCount, solverTolerance, functionIdentifier, 
        gridPointsX, gridPointsY, numThreads);

    // Uncomment for debugging
    // displayVector(totalGridPoints, solutionVector);
    // displayVector(totalGridPoints, rightHandSide);    

    // Clean up resources
    free_results();
    delete[] rowIndices;
    delete[] matrixValues;
    delete[] rightHandSide;
    delete[] solutionVector;
    delete[] residual;
    delete[] auxVector1;
    delete[] auxVector2;
    delete[] threadArgs;
    delete[] threadIds;

    return 0;
}
