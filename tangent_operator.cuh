#ifndef TANGENT_OPERATOR_CUH
#define TANGENT_OPERATOR_CUH

#include <cuda_runtime.h>

// Material properties structure
struct MaterialProperties {
    double E;        // Young's modulus
    double nu;       // Poisson's ratio
    double H;        // Hardening modulus
    double beta;     // Back stress parameter
    double sigma_y0; // Initial yield stress
};

// CUDA kernel function declaration
__global__ void computeTangentOperator(
    const double* strain,
    const double* stress,
    const double* alpha,
    const double* R,
    double* tangent_operator,
    double E,
    double nu,
    double H,
    double beta,
    double sigma_y0,
    int num_elements
);

// Host function declaration
void launchTangentOperatorComputation(
    double* h_tangent_operator,
    const double* h_strain,
    const double* h_stress,
    const double* h_alpha,
    const double* h_R,
    const MaterialProperties& props,
    int num_elements
);

#endif // TANGENT_OPERATOR_CUH 