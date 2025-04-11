#include "all_includes.h"
#include <sys/resource.h>

#define FUNC(I, J) do { ij2l(nx, ny, I, J, k); if (I_ij) { I_ij[m] = k; } m++; } \
                  while (0)

#define F(I, J) (f(a + (I)*hx, с + (J)*hy))

void matrix_mult_vector_msr(int n, double* A, int* I, double* x, double* y, int p, int k) {
    int i, i1, i2, l, J; double s;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        s = A[i] * x[i];
        l = I[i+1] - I[i];
        J = I[i];
        for (int j = 0; j < l; ++j) {
            s += A[J + j] * x[I[J + j]];
        }

        y[i] = s;
    }
}

void apply_preconditioner_msr_matrix(int n, double* A, int* I, double* v1, double* v2, int flag, int p, int k) {
    // Коэффициент релаксации для блочного метода
    const double omega = 1.0; 
    
    // Выбор и применение треугольного решателя в зависимости от флага
    if (flag == 0) {
        // Нижняя треугольная часть
        solve_rsystem(n, I, A, v2, v1, omega, p, k);
    } else {
        // Верхняя треугольная часть
        solve_lsystem(n, I, A, v2, v1, omega, p, k);
    }

    // Необходимая синхронизация
    reduce_sum<int>(p);
}

bool step(int n, double* A, int* I, double* x, double* r, double* u, double* v, double prec, int p, int k) {
    // Вычисление матрично-векторного произведения 
    matrix_mult_vector_msr(n, A, I, v, u, p, k);
    
    // Оценка норм для проверки сходимости
    const double residual_norm = scalar_product(n, r, r, p, k);
    const double direction_norm = scalar_product(n, u, u, p, k);

    if (residual_norm < prec || direction_norm < prec) {
        return true; // Достигнута сходимость
    }
    // Вычисление оптимального шага для минимизации ошибки
    const double step_size = residual_norm / direction_norm;
    
    // Обновление приближения решения
    mult_sub_vector(n, x, v, step_size, p, k);
    
    // Обновление вектора невязки
    mult_sub_vector(n, r, u, step_size, p, k);
    
    // Процесс продолжается
    return false;
}

int minimal_errors_msr_matrix(int n, double* A, int* I, double* b, double* x,
    double* r, double* u, double* v, double eps, int maxit, int p, int k) {
    
    // Инициализация параметров метода
    double convergence_threshold;
    int iteration_count;
    
    // Вычисление квадрата нормы правой части
    const double rhs_norm_squared = scalar_product(n, b, b, p, k);
    
    // Установка порога сходимости
    convergence_threshold = rhs_norm_squared * eps * eps;
    
    // Формирование начальной невязки r = Ax - b
    matrix_mult_vector_msr(n, A, I, x, r, p, k);
    mult_sub_vector(n, r, b, 1.0, p, k);
    
    // Основной итерационный процесс
    for (iteration_count = 0; iteration_count < maxit; ++iteration_count) {
        // Этап 1: Применение предобуславливателя к невязке (нижняя часть)
        apply_preconditioner_msr_matrix(n, A, I, v, r, 0, p, k);
        
        // Проверка сходимости после первого полушага
        if (step(n, A, I, x, r, u, v, convergence_threshold, p, k)) {
            break;
        }
        
        // Этап 2: Пересчет невязки для текущего приближения
        matrix_mult_vector_msr(n, A, I, x, u, p, k);
        mult_sub_vector(n, u, b, 1.0, p, k);
        
        // Этап 3: Применение предобуславливателя к обновленной невязке (верхняя часть)
        apply_preconditioner_msr_matrix(n, A, I, v, u, 1, p, k);
        
        // Проверка сходимости после второго полушага
        if (step(n, A, I, x, r, u, v, convergence_threshold, p, k)) {
            break;
        }
    }
    
    // Возвращаем результат о сходимости
    if (iteration_count >= maxit) {
        return -1; // Не достигнута сходимость
    }
    
    return iteration_count; // Количество итераций до сходимости
}

int minimal_errors_msr_matrix_full(int n, double* A, int* I, double* b, double* x, double* r, double* u, double* v, 
    double eps, int maxit, int maxsteps, int p, int k) {

    int current_attempt;
    int convergence_status;
    int total_iterations = 0;
    
    // Последовательные попытки достичь сходимости
    for (current_attempt = 0; current_attempt < maxsteps; ++current_attempt) {
        // Запуск одного цикла итерационного метода
        convergence_status = minimal_errors_msr_matrix(n, A, I, b, x, r, u, v, eps, maxit, p, k);
        
        // Проверка успешности
        if (convergence_status >= 0) {
            // Метод сошелся, добавляем число затраченных итераций
            total_iterations += convergence_status;
            break;
        }
        
        // Метод не сошелся за отведенное число итераций
        total_iterations += maxit;
    }
    
    // Проверка достижения максимального числа попыток
    if (current_attempt >= maxsteps) {
        return -1; // Сходимость не достигнута
    }
    
    // Возвращаем общее число итераций до сходимости
    return total_iterations;
}

void thread_rows(int n, int p, int k, int& i1, int& i2) {
    i1 = n*k; i1 /= p; i2 = n*(k + 1); i2 /= p; 
}

double scalar_product(int n, double* x, double* y, int p, int k) {
    int i1, i2, i; double s = 0;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        s += x[i] * y[i];
    }

    s = reduce_sum_det(p, k, s);
    return s;
}

void mult_sub_vector(int n, double* x, double* y, double tau, int p, int k) {
    int i, i1, i2;
    thread_rows(n, p, k, i1, i2);
    for (i = i1; i < i2; ++i) {
        x[i] -= tau * y[i];
    }

    reduce_sum<int>(p);
}

void ij2l(int nx, int, int i, int j, int& l) {
    l = i + j * (int)(nx + 1);
}

void l2ij(int nx, int, int& i, int& j, int l) {
    j = l / (nx + 1);
    i = l - j * (nx + 1);
}

int get_len_msr(int nx, int ny) {
    return (nx + 1)*(ny + 1)
            + 6*(nx - 1)*(ny - 1)
            + 4*(2*(nx-1) + 2*(ny-1))
            + 2*3 + 2*2; 
}

int get_off_diag(int nx, int ny, int i, int j, int* I_ij) {
    // Completely different approach using array of potential neighbors
    struct Neighbor { int di, dj; };
    const Neighbor neighbors[] = {
        {1, 0},   // right
        {0, -1},  // down
        {-1, -1}, // down-left
        {-1, 0},  // left
        {0, 1},   // up
        {1, 1}    // up-right
    };
    
    int count = 0;
    for (int idx = 0; idx < 6; ++idx) {
        int ni = i + neighbors[idx].di;
        int nj = j + neighbors[idx].dj;
        
        // Check if neighbor is within bounds
        if (ni >= 0 && ni <= nx && nj >= 0 && nj <= ny) {
            if (I_ij != nullptr) {
                int l;
                ij2l(nx, ny, ni, nj, l);
                I_ij[count] = l;
            }
            count++;
        }
    }
    
    return count;
}

int get_len_msr_off_diag(int nx, int ny) {
    // Alternative implementation using running sum
    int total_offdiag = 0;
    
    // Process each grid point
    for (int row = 0; row <= ny; ++row) {
        for (int col = 0; col <= nx; ++col) {
            // Add count of off-diagonal elements for this cell
            total_offdiag += get_off_diag(nx, ny, col, row, nullptr);
        }
    }
    
    return total_offdiag;
}

int allocate_msr_matrix(int nx, int ny, double** p_A, int** p_I) {
    // Calculate sizes differently
    const int grid_size = (nx+1)*(ny+1);
    
    // Calculate off-diagonal elements
    int offdiag_elements = 0;
    for (int node = 0; node < grid_size; ++node) {
        int i, j;
        l2ij(nx, ny, i, j, node);
        offdiag_elements += get_off_diag(nx, ny, i, j, nullptr);
    }
    
    // Total size includes diagonal, offdiagonal and sentinel
    int total_size = grid_size + offdiag_elements + 1;
    
    // Allocate memory with error checking
    try {
        *p_A = new double[total_size];
        *p_I = new int[total_size];
    } catch (std::bad_alloc&) {
        if (*p_A != nullptr) {
            delete[] *p_A;
            *p_A = nullptr;
        }
        return 1; // allocation failed
    }
    
    return 0; // success
}

void fill_I(int nx, int ny, int* I) {
    const int width = nx + 1;
    const int height = ny + 1;
    const int total_nodes = width * height;
    
    // Start with initial offset after diagonal
    int current_offset = total_nodes + 1;
    
    // For each node in the grid
    for (int node_idx = 0; node_idx < total_nodes; ++node_idx) {
        // Store start of offdiagonal entries for this row
        I[node_idx] = current_offset;
        
        // Convert linear index to 2D coordinates
        int x, y;
        l2ij(nx, ny, x, y, node_idx);
        
        // Store indices of offdiagonal entries
        int neighbor_count = get_off_diag(nx, ny, x, y, &I[current_offset]);
        
        // Update offset for next row
        current_offset += neighbor_count;
    }
    
    // Set sentinel at the end
    I[total_nodes] = current_offset;
}

void fill_A_ij(int nx, int ny, double hx, double hy, int i, int j, double* A_diag, double* A_off_diag) {

    double s = hx*hy;
    if (i > 0 && i < nx && j > 0 && j < ny) {
        *A_diag = s / 2;
        for (int l = 0; l < 6; ++l) {
            A_off_diag[l] = s / 12; 
        }
    } else if (j == 0 && i > 0 && i < nx) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
        A_off_diag[2] = 2*s/24;
        A_off_diag[3] = 2*s/24;
    } else if (j == ny && i > 0 && i < nx) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 2*s/24;
        A_off_diag[2] = 2*s/24;
        A_off_diag[3] = 1*s/24;
    } else if (i == 0 && j > 0 && j < ny) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 2*s/24;
        A_off_diag[1] = 1*s/24;
        A_off_diag[2] = 1*s/24;
        A_off_diag[3] = 2*s/24;
    } else if (i == nx && j > 0 && j < ny) {
        *A_diag = 3*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 2*s/24;
        A_off_diag[2] = 2*s/24;
        A_off_diag[3] = 1*s/24;
    } else if (i == 0 && j == 0) {
        *A_diag = 2*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
        A_off_diag[2] = 2*s/24;
    } else if (i == nx && j == ny) {
        *A_diag = 2*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 2*s/24;
        A_off_diag[2] = 1*s/24;
    } else if (i == 0 && j == ny) {
        *A_diag = 1*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
    } else if (i == nx && j == 0) {
        *A_diag = 1*s/12;
        A_off_diag[0] = 1*s/24;
        A_off_diag[1] = 1*s/24;
    }

}

void fill_A(int nx, int ny, double hx, double hy, int* I, double* A, int p, int k) {
    // Changed coding style with camelCase and different structure
    const int totalNodes = (nx + 1) * (ny + 1);
    
    // Calculate thread's workload boundaries
    const int startIdx = (totalNodes * k) / p;
    const int endIdx = (totalNodes * (k + 1)) / p;

    // Process each node assigned to this thread
    for (int nodeIndex = startIdx; nodeIndex < endIdx; ++nodeIndex) {
        // Extract grid coordinates
        int gridX, gridY;
        l2ij(nx, ny, gridX, gridY, nodeIndex);
        
        // Get pointers to matrix storage locations
        double* diagElement = &A[nodeIndex];
        double* offDiagElements = &A[I[nodeIndex]];

        // Compute matrix elements for this node
        fill_A_ij(
            nx, ny,                 // Grid dimensions
            hx, hy,                 // Grid spacing
            gridX, gridY,           // Node coordinates
            diagElement,            // Diagonal element location
            offDiagElements         // Off-diagonal elements location
        );
    }

    // Synchronize threads
    reduce_sum<int>(p);
}

// Check matrix symmetry to avoid potential errors in calculations
int check_symm(int nx, int ny, int* I, double* A, double eps, int p, int k) {
    // Define constants and variables with descriptive names
    const int GRID_SIZE = (nx+1)*(ny+1);
    const int THREAD_START = (k * GRID_SIZE) / p;
    const int THREAD_END = ((k + 1) * GRID_SIZE) / p;
    
    // Count of symmetry violations detected by this thread
    int symmetryErrors = 0;
    
    // For each row assigned to this thread
    for (int rowIdx = THREAD_START; rowIdx < THREAD_END; ++rowIdx) {
        // Get information about non-zero elements in this row
        const int elementsInRow = I[rowIdx + 1] - I[rowIdx];
        const int rowOffset = I[rowIdx];
        
        // Check each off-diagonal element in current row
        for (int elemPos = 0; elemPos < elementsInRow; ++elemPos) {
            // Value and column of current matrix element
            const double currentValue = A[rowOffset + elemPos];
            const int colIdx = I[rowOffset + elemPos];
            
            // Now search for matching element A(colIdx, rowIdx) in row colIdx
            const int matchingRowOffset = I[colIdx];
            const int matchingRowSize = I[colIdx + 1] - matchingRowOffset;
            
            // Search for position of symmetric element
            int matchingPos = 0;
            bool foundMatch = false;
            
            // Linear search through non-zero elements
            for (; matchingPos < matchingRowSize; ++matchingPos) {
                if (I[matchingRowOffset + matchingPos] == rowIdx) {
                    foundMatch = true;
                    break;
                }
            }
            
            // Check if symmetry is violated
            if (!foundMatch) {
                // Missing symmetric element
                symmetryErrors++;
            } else if (fabs(A[matchingRowOffset + matchingPos] - currentValue) > eps) {
                // Symmetric element exists but values differ too much
                symmetryErrors++;
            }
        }
    }
    
    // Combine error counts from all threads
    reduce_sum<int>(p, &symmetryErrors, 1);
    return symmetryErrors;
}

double F_IJ(int nx, int ny, double hx, double hy, double a, double с, int i, int j, double (*f)(double, double)) {
    // Quadrature constant - scaled element area
    const double quadWeight = hx * hy / 192.0;
    
    // Define node types for gridpoint location classification
    enum NodeType {
        INTERIOR,       // Interior node
        EDGE_BOTTOM,    // Bottom edge (j=0), not corner
        EDGE_TOP,       // Top edge (j=ny), not corner
        EDGE_LEFT,      // Left edge (i=0), not corner
        EDGE_RIGHT,     // Right edge (i=nx), not corner
        CORNER_SW,      // Bottom-left / southwest corner (0,0)
        CORNER_NE,      // Top-right / northeast corner (nx,ny)
        CORNER_NW,      // Top-left / northwest corner (0,ny)
        CORNER_SE       // Bottom-right / southeast corner (nx,0)
    };
    
    // Determine node type based on position
    NodeType nodeType;
    
    if (i > 0 && i < nx && j > 0 && j < ny) {
        nodeType = INTERIOR;
    } else if (i > 0 && i < nx && j == 0) {
        nodeType = EDGE_BOTTOM;
    } else if (i > 0 && i < nx && j == ny) {
        nodeType = EDGE_TOP;
    } else if (i == 0 && j > 0 && j < ny) {
        nodeType = EDGE_LEFT;
    } else if (i == nx && j > 0 && j < ny) {
        nodeType = EDGE_RIGHT;
    } else if (i == 0 && j == 0) {
        nodeType = CORNER_SW;
    } else if (i == nx && j == ny) {
        nodeType = CORNER_NE;
    } else if (i == 0 && j == ny) {
        nodeType = CORNER_NW;
    } else if (i == nx && j == 0) {
        nodeType = CORNER_SE;
    } else {
        // Should never happen
        return 1e308;
    }
    
    // Calculate quadrature based on node type
    switch (nodeType) {
        case INTERIOR: {
            // Center point
            double centerTerm = 36.0 * F(i, j);
            
            // Edge midpoints (6 terms)
            double edgeMidpointTerms = 
                20.0 * (F(i+0.5, j) + F(i, j-0.5) + 
                       F(i-0.5, j-0.5) + F(i-0.5, j) + 
                       F(i, j+0.5) + F(i+0.5, j+0.5));
            
            // "Secondary" points (6 terms)
            double secondaryTerms = 
                4.0 * (F(i+0.5, j-0.5) + F(i-0.5, j-1) + 
                      F(i-1, j-0.5) + F(i-0.5, j+0.5) + 
                      F(i+0.5, j+1) + F(i+1, j+0.5));
            
            // Corner points (6 terms)
            double cornerTerms = 
                2.0 * (F(i+1, j) + F(i, j-1) + 
                      F(i-1, j-1) + F(i-1, j) + 
                      F(i, j+1) + F(i+1, j+1));
                
            return quadWeight * (centerTerm + edgeMidpointTerms + secondaryTerms + cornerTerms);
        }
            
        case EDGE_BOTTOM: {
            double centerTerm = 18.0 * F(i, j);
            double horizontalMidpoints = 10.0 * (F(i+0.5, j) + F(i-0.5, j));
            double upwardMidpoints = 20.0 * (F(i, j+0.5) + F(i+0.5, j+0.5));
            double diagonalPoints = 4.0 * (F(i-0.5, j+0.5) + F(i+0.5, j+1) + F(i+1, j+0.5));
            double sidePoints = 1.0 * (F(i-1, j) + F(i+1, j));
            double topPoints = 2.0 * (F(i, j+1) + F(i+1, j+1));
            
            return quadWeight * (centerTerm + horizontalMidpoints + upwardMidpoints + 
                                diagonalPoints + sidePoints + topPoints);
        }
            
        case EDGE_TOP: {
            double centerTerm = 18.0 * F(i, j);
            double horizontalMidpoints = 10.0 * (F(i+0.5, j) + F(i-0.5, j));
            double downwardMidpoints = 20.0 * (F(i, j-0.5) + F(i-0.5, j-0.5));
            double diagonalPoints = 4.0 * (F(i+0.5, j-0.5) + F(i-0.5, j-1) + F(i-1, j-0.5));
            double sidePoints = 1.0 * (F(i-1, j) + F(i+1, j));
            double bottomPoints = 2.0 * (F(i, j-1) + F(i-1, j-1));
            
            return quadWeight * (centerTerm + horizontalMidpoints + downwardMidpoints + 
                                diagonalPoints + sidePoints + bottomPoints);
        }
            
        case EDGE_LEFT: {
            double centerTerm = 18.0 * F(i, j);
            double verticalMidpoints = 10.0 * (F(i, j-0.5) + F(i, j+0.5));
            double rightMidpoints = 20.0 * (F(i+0.5, j) + F(i+0.5, j+0.5));
            double diagonalPoints = 4.0 * (F(i+0.5, j-0.5) + F(i+0.5, j+1) + F(i+1, j+0.5));
            double topBottomPoints = 1.0 * (F(i, j-1) + F(i, j+1));
            double rightPoints = 2.0 * (F(i+1, j) + F(i+1, j+1));
            
            return quadWeight * (centerTerm + verticalMidpoints + rightMidpoints + 
                                diagonalPoints + topBottomPoints + rightPoints);
        }
            
        case EDGE_RIGHT: {
            double centerTerm = 18.0 * F(i, j);
            double verticalMidpoints = 10.0 * (F(i, j-0.5) + F(i, j+0.5));
            double leftMidpoints = 20.0 * (F(i-0.5, j) + F(i-0.5, j-0.5));
            double diagonalPoints = 4.0 * (F(i-0.5, j-1) + F(i-1, j-0.5) + F(i-0.5, j+0.5));
            double topBottomPoints = 1.0 * (F(i, j-1) + F(i, j+1));
            double leftPoints = 2.0 * (F(i-1, j) + F(i-1, j-1));
            
            return quadWeight * (centerTerm + verticalMidpoints + leftMidpoints + 
                                diagonalPoints + topBottomPoints + leftPoints);
        }
            
        case CORNER_SW: {
            double centerTerm = 12.0 * F(i, j);
            double adjacentMidpoints = 10.0 * (F(i+0.5, j) + F(i, j+0.5));
            double diagonalMidpoint = 20.0 * F(i+0.5, j+0.5);
            double secondaryPoints = 4.0 * (F(i+1, j+0.5) + F(i+0.5, j+1));
            double adjacentPoints = 1.0 * (F(i+1, j) + F(i, j+1));
            double diagonalPoint = 2.0 * F(i+1, j+1);
            
            return quadWeight * (centerTerm + adjacentMidpoints + diagonalMidpoint + 
                                secondaryPoints + adjacentPoints + diagonalPoint);
        }
            
        case CORNER_NE: {
            double centerTerm = 12.0 * F(i, j);
            double adjacentMidpoints = 10.0 * (F(i-0.5, j) + F(i, j-0.5));
            double diagonalMidpoint = 20.0 * F(i-0.5, j-0.5);
            double secondaryPoints = 4.0 * (F(i-0.5, j-1) + F(i-1, j-0.5));
            double adjacentPoints = 1.0 * (F(i, j-1) + F(i-1, j));
            double diagonalPoint = 2.0 * F(i-1, j-1);
            
            return quadWeight * (centerTerm + adjacentMidpoints + diagonalMidpoint + 
                                secondaryPoints + adjacentPoints + diagonalPoint);
        }
            
        case CORNER_NW: {
            double centerTerm = 6.0 * F(i, j);
            double adjacentMidpoints = 10.0 * (F(i+0.5, j) + F(i, j-0.5));
            double diagonalPoint = 4.0 * F(i+0.5, j-0.5);
            double adjacentPoints = 1.0 * (F(i+1, j) + F(i, j-1));
            
            return quadWeight * (centerTerm + adjacentMidpoints + diagonalPoint + adjacentPoints);
        }
            
        case CORNER_SE: {
            double centerTerm = 6.0 * F(i, j);
            double adjacentMidpoints = 10.0 * (F(i-0.5, j) + F(i, j+0.5));
            double diagonalPoint = 4.0 * F(i-0.5, j+0.5);
            double adjacentPoints = 1.0 * (F(i-1, j) + F(i, j+1));
            
            return quadWeight * (centerTerm + adjacentMidpoints + diagonalPoint + adjacentPoints);
        }
    }
    
    return 1e308;
}

void fill_B(int nx, int ny, double hx, double hy, double a, double c, double* B, double (*f)(double, double), int p, int k) {
    int l1, l2;
    int i, j;
    int N = (nx + 1) * (ny + 1);    
    l1 = N * k; l1 /= p;
    l2 = N * (k + 1); l2 /= p;

    for (int l = l1; l < l2; ++l) {
        l2ij(nx, ny, i, j, l);
        B[l] = F_IJ(nx, ny, hx, hy, a, c, i, j, f);
    }

    reduce_sum<int>(p);
}

void solve_rsystem(int n, int* I, double* U, double* b, double* x, double w, int p, int k) {
    // Обработка нижней треугольной части
    // Начало и конец интервала для данного потока
    int start_idx, end_idx;
    thread_rows(n, p, k, start_idx, end_idx);

    // Обратный порядок обхода для нижней треугольной матрицы
    for (int current = end_idx - 1; current >= start_idx; --current) {
        // Количество ненулевых внедиагональных элементов
        const int num_elements = I[current + 1] - I[current];
        
        // Вычисление суммы известных элементов
        double sum_known = 0.0;
        const int offset = I[current];
        
        // Перебор всех ненулевых элементов
        for (int elem_idx = 0; elem_idx < num_elements; ++elem_idx) {
            const int col_idx = I[offset + elem_idx];
            // Учитываем только элементы правее диагонали
            if (col_idx > current && col_idx >= start_idx && col_idx < end_idx) {
                sum_known += x[col_idx] * U[offset + elem_idx];
            }
        }
        
        // Вычисление значения для текущей переменной
        x[current] = w * (b[current] - sum_known) / U[current];
    }
}

void solve_lsystem(int n, int* I, double* U, double* b, double* x, double w, int p, int k) {
    // Получение диапазона индексов для текущего потока
    int range_begin, range_end;
    thread_rows(n, p, k, range_begin, range_end);
    
    // Прямой проход для верхнетреугольной системы (сверху вниз)
    for (int row = range_begin; row < range_end; ++row) {
        // Сумма влияния уже вычисленных неизвестных
        double accumulated_effect = 0.0;
        
        // Индекс начала внедиагональных элементов в текущей строке
        const int row_start = I[row];
        // Количество внедиагональных элементов
        const int element_count = I[row + 1] - row_start;
        
        // Для каждого ненулевого элемента в строке
        for (int pos = 0; pos < element_count; ++pos) {
            // Индекс столбца текущего элемента
            const int col = I[row_start + pos];
            
            // Учитываем только элементы слева от диагонали в пределах диапазона потока
            if (col < row && col >= range_begin && col < range_end) {
                accumulated_effect += x[col] * U[row_start + pos];
            }
        }
        
        // Обновление текущей неизвестной с учетом коэффициента релаксации
        x[row] = w * (b[row] - accumulated_effect) / U[row];
    }
}

double get_cpu_time() {
    struct rusage buf;
    getrusage(RUSAGE_THREAD, &buf);
    return buf.ru_utime.tv_sec + buf.ru_utime.tv_usec * 1e-6;
}
