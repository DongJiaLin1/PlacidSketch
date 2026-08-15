/**
 * @file Param.h
 * @brief XSketch algorithm parameters and regression functions
 */

#ifndef _PARAM_H_
#define _PARAM_H_

#include <bits/stdc++.h>
#include "Common/Matrix.h"

// Time window size
const uint32_t window_size = 10000;

// Memory size in MB
const int memory_p = 300;

// Memory ratio for Stage1
const double stage_ratio_p = 0.8;

// Variance threshold for slope coefficient
const double var_thres_p = 0.1;

// Number of windows for linear regression
const int P = 7;

// Polynomial degree (0 = constant model)
const int K = 0;

// Number of windows in Stage1
const int S = 4;

const int S_p = S;

// Mean squared error threshold
const double error_thres_p = 0.5;

// Number of slots per bucket
const int bucket_size_p = 4;

// Potential threshold for replacement
const double potential_thres_p = 0.25;

// Coefficient matrix for linear regression (P windows, K+1 coefficients)
double calcu_matrix[K + 1][P] = {};

// Coefficient matrix for Stage1 (S windows, K+1 coefficients)
double calcu_matrix_try[K + 1][S_p] = {};

/**
 * @brief Report slot for storing flow information
 */
template<typename ID_TYPE>
class Report_Slot {
public:
    ID_TYPE id;              // Flow ID
    uint32_t start_window;   // Start window
    uint32_t end_window;     // End window

    Report_Slot() {}

    Report_Slot(ID_TYPE _id, uint32_t st, uint32_t et): 
        id(_id), start_window(st), end_window(et) {}

    ~Report_Slot() {}

    bool operator < (const Report_Slot &r) {
        if (end_window - start_window == r.end_window - r.start_window) {
            if (start_window == r.start_window)
                return id < r.id;
            return start_window < r.start_window;
        }
        return end_window - start_window < r.end_window - r.start_window;
    }
};

/**
 * @brief Initialize coefficient matrices for linear regression
 */
void init_matrix() {
    double X[P][K + 1] = {};
    double XX[K + 1][K + 1] = {};
    double temp[(K + 1) * (K + 1)];

    // Build design matrix X[i][j] = i^j
    for (int i = 0; i <= P - 1; ++i) {
        for (int j = 0; j <= K; ++j) {
            X[i][j] = pow(i, j);
        }
    }

    // Compute X'X
    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= K; ++j) {
            for (int k = 0; k <= P - 1; ++k) {
                XX[i][j] += X[k][i] * X[k][j];
            }
        }
    }

    // Copy to temp array
    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= K; ++j) {
            temp[i * (K + 1) + j] = XX[i][j];
        }
    }

    // Compute matrix inverse
    inverse(K + 1, temp);

    // Update XX to inverse matrix
    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= K; ++j) {
            XX[i][j] = temp[i * (K + 1) + j];
        }
    }

    // Compute (X'X)^{-1}X'
    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= P - 1; ++j) {
            for (int k = 0; k <= K; ++k) {
                calcu_matrix[i][j] += XX[i][k] * X[j][k];
            }
        }
    }

    // Matrix computation for Stage1
    double Z[S_p][K + 1] = {};

    for (int i = 0; i <= S_p - 1; ++i) {
        for (int j = 0; j <= K; ++j) {
            X[i][j] = pow(i, j);
        }
    }

    memset(XX, 0, sizeof(XX));
    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= K; ++j) {
            for (int k = 0; k <= S_p - 1; ++k) {
                XX[i][j] += X[k][i] * X[k][j];
            }
        }
    }

    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= K; ++j) {
            temp[i * (K + 1) + j] = XX[i][j];
        }
    }

    inverse(K + 1, temp);

    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= K; ++j) {
            XX[i][j] = temp[i * (K + 1) + j];
        }
    }

    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= S_p - 1; ++j) {
            for (int k = 0; k <= K; ++k) {
                calcu_matrix_try[i][j] += XX[i][k] * Z[j][k];
            }
        }
    }
}

/**
 * @brief Linear regression for P windows
 */
void linear_regressing(double* y, double* b) {
    memset(b, 0, (K + 1) * sizeof(double));
    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= P - 1; ++j) {
            b[i] += calcu_matrix[i][j] * y[j];
        }
    }
}

/**
 * @brief Linear regression for S windows
 */
void linear_regressing_try(double* y, double* b) {
    memset(b, 0, (K + 1) * sizeof(double));
    for (int i = 0; i <= K; ++i) {
        for (int j = 0; j <= S_p - 1; ++j) {
            b[i] += calcu_matrix[i][j] * y[j];
        }
    }
}

#endif
