#include "tangent_operator.cuh"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdio.h>
#include <math.h>

// Constants
#define BLOCK_SIZE 256
#define NUM_COMPONENTS 6

// Function prototypes
__device__ void computeElasticTangent(double* C_alg, const MaterialProperties& props);
__device__ void computeDeviatoricStress(const double* stress, double* dev_stress);
__device__ double computeEquivalentStress(const double* XI);
__device__ void computePlasticTangent(double* C_alg, const double* XI, double a, const double* nhat, double delta_lambda, const MaterialProperties& props);

// Device function: Compute elastic prediction
__device__ void computeElasticPrediction(
    double* sig_tr,
    const double* current_stress,
    const double* strain_increment,
    const MaterialProperties& props
) {
    double Ce[36];
    computeElasticTangent(Ce, props);
    
    // Convert engineering shear strain to tensor shear strain
    double strain_tensor[6];
    for (int i = 0; i < 3; i++) {
        strain_tensor[i] = strain_increment[i];
    }
    for (int i = 3; i < 6; i++) {
        strain_tensor[i] = strain_increment[i] / 2.0;  // Convert engineering shear strain to tensor shear strain
    }
    
    // Matrix multiplication: sig_tr = current_stress + Ce * strain_tensor
    for (int i = 0; i < 6; i++) {
        sig_tr[i] = current_stress[i];
        for (int j = 0; j < 6; j++) {
            sig_tr[i] += Ce[i + j*6] * strain_tensor[j];
        }
    }
}

// Device function: Compute deviatoric stress
__device__ void computeDeviatoricStress(const double* stress, double* dev_stress) {
    // Compute mean stress
    double mean_stress = (stress[0] + stress[1] + stress[2]) / 3.0;
    
    // Compute deviatoric stress (using vectorized operations)
    #pragma unroll
    for (int i = 0; i < 3; i++) {
        dev_stress[i] = stress[i] - mean_stress;
    }
    #pragma unroll
    for (int i = 3; i < 6; i++) {
        dev_stress[i] = stress[i];
    }
}

// Device function: Compute relative stress
__device__ void computeRelativeStress(
    double* XI,
    const double* s_tr,
    const double* alpha
) {
    for (int i = 0; i < 6; i++) {
        XI[i] = s_tr[i] - alpha[i];
    }
}

// Device function: Compute equivalent stress
__device__ double computeEquivalentStress(const double* XI) {
    double a = 0.0;
    // Use vectorized operations
    #pragma unroll
    for (int i = 0; i < 3; i++) {
        a += XI[i] * XI[i];
    }
    #pragma unroll
    for (int i = 3; i < 6; i++) {
        a += 2.0 * XI[i] * XI[i];
    }
    return sqrt(a);
}

// Device function: Compute elastic tangent operator
__device__ void computeElasticTangent(
    double* C_alg,
    const MaterialProperties& props
) {
    double mu = props.E / (2.0 * (1.0 + props.nu));
    double lambda = props.E * props.nu / ((1.0 + props.nu) * (1.0 - 2.0 * props.nu));
    
    // Initialize tangent operator
    for (int i = 0; i < 36; i++) {
        C_alg[i] = 0.0;
    }
    
    // Fill elastic tangent operator (column-major storage)
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            C_alg[i + j*6] = lambda;
        }
        C_alg[j + j*6] += 2.0 * mu;
    }
    
    for (int i = 3; i < 6; i++) {
        C_alg[i + i*6] = mu;
    }
}

// Device function: Compute plastic tangent operator
__device__ void computePlasticTangent(
    double* C_alg,
    const double* XI,
    double a,
    const double* nhat,
    double delta_lambda,
    const MaterialProperties& props
) {
    double mu = props.E / (2.0 * (1.0 + props.nu));
    
    // Compute elastic tangent operator
    double Ce[36];
    computeElasticTangent(Ce, props);
    
    // Compute nhat*nhat' term (column-major storage)
    double nhat_nhat[36];
    for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 6; i++) {
            nhat_nhat[i + j*6] = nhat[i] * nhat[j];
        }
    }
    
    // Compute Cep
    double Cep[36];
    double factor = 2.0 * mu / (1.0 + props.H/(2.0 * mu));
    for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 6; i++) {
            Cep[i + j*6] = Ce[i + j*6] - factor * nhat_nhat[i + j*6];
        }
    }
    
    // Compute psi
    double psi = 2.0 * mu * delta_lambda / a;
    
    // Correctly compute Ia matrix
    double Ia[36];
    for (int i = 0; i < 36; i++) {
        Ia[i] = 0.0;
    }
    for (int i = 0; i < 3; i++) {
        Ia[i + i*6] = 1.0;  // First 3 diagonal elements
    }
    for (int i = 3; i < 6; i++) {
        Ia[i + i*6] = 0.5;  // Last 3 diagonal elements
    }
    
    // ID*ID' term remains unchanged
    double ID_ID[36];
    for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 6; i++) {
            ID_ID[i + j*6] = (i < 3 && j < 3) ? 1.0 : 0.0;
        }
    }
    
    // Compute final tangent operator
    for (int j = 0; j < 6; j++) {
        for (int i = 0; i < 6; i++) {
            double temp = Ia[i + j*6] - (1.0/3.0) * ID_ID[i + j*6] - nhat_nhat[i + j*6];
            C_alg[i + j*6] = Cep[i + j*6] - 2.0 * mu * psi * temp;
        }
    }
}

// Device function: Compute stress norm
__device__ double computeStressNorm(const double* stress) {
    double norm = 0.0;
    for (int i = 0; i < 6; i++) {
        norm += stress[i] * stress[i];
    }
    return sqrt(norm);
}

// Compute elastic stiffness matrix
__device__ void computeElasticStiffness(
    double* C,
    const MaterialProperties& props
) {
    double E = props.E;
    double nu = props.nu;
    double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));

    // Initialize stiffness matrix
    for (int i = 0; i < 36; i++) {
        C[i] = 0.0;
    }

    // Fill stiffness matrix (column-major storage)
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            C[i + j*6] = lambda;
        }
        C[j + j*6] += 2 * mu;
    }

    for (int i = 3; i < 6; i++) {
        C[i + i*6] = mu;
    }
}

// Compute plastic flow direction
__device__ void computePlasticFlow(
    double* n,
    const double* stress,
    const double* alpha,
    const double R,
    const MaterialProperties& props
) {
    double sigma_dev[6];
    double sigma_eff[6];
    double norm = 0.0;
    
    // Compute deviatoric stress
    computeDeviatoricStress(stress, sigma_dev);
    
    // Compute relative stress
    for (int i = 0; i < 6; i++) {
        sigma_eff[i] = sigma_dev[i] - alpha[i];
    }
    
    // Compute norm
    norm = computeEquivalentStress(sigma_eff);
    
    // Compute flow direction
    if (norm > 1e-10) {
        for (int i = 0; i < 6; i++) {
            n[i] = sigma_eff[i] / norm;
        }
    } else {
        for (int i = 0; i < 6; i++) {
            n[i] = 0.0;
        }
    }
}

// CUDA kernel: Compute tangent operator
__global__ void computeTangentOperator(
    const double* strain_increment,
    const double* current_stress,
    const double* alpha,
    const double* R,
    double* tangent_operator,
    double E,
    double nu,
    double H,
    double beta,
    double sigma_y0,
    int num_elements
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_elements) return;

    // Set material properties
    MaterialProperties props;
    props.E = E;
    props.nu = nu;
    props.H = H;
    props.beta = beta;
    props.sigma_y0 = sigma_y0;

    // Get current element data
    const double* curr_strain_inc = &strain_increment[idx * 6];
    const double* curr_stress = &current_stress[idx * 6];
    const double* curr_alpha = &alpha[idx * 6];
    double curr_R = R[idx];
    double* curr_tangent = &tangent_operator[idx * 36];

    // Debug output for first element
    if (idx == 0) {
        printf("=== CUDA Debug (Element 0) ===\n");
        printf("strain_increment (before conversion): ");
        for(int i = 0; i < 6; i++) printf("%.6e ", curr_strain_inc[i]);
        printf("\n");
    }

    // Compute elastic prediction stress
    double sig_tr[6];
    computeElasticPrediction(sig_tr, curr_stress, curr_strain_inc, props);

    // Debug output for first element
    if (idx == 0) {
        printf("sig_tr: ");
        for(int i = 0; i < 6; i++) printf("%.6e ", sig_tr[i]);
        printf("\n");
    }

    // Compute deviatoric stress
    double s_tr[6];
    computeDeviatoricStress(sig_tr, s_tr);

    // Compute relative stress
    double XI[6];
    computeRelativeStress(XI, s_tr, curr_alpha);

    // Compute equivalent stress
    double a = computeEquivalentStress(XI);

    // Debug output for first element - add a and R output
    if (idx == 0) {
        printf("Element 0: a = %.6e, R = %.6e\n", a, curr_R);
        printf("Yield condition (a <= R): %d\n", a <= curr_R);
    }

    // Compute tangent operator
    double mu = props.E / (2.0 * (1.0 + props.nu));
    
    if (a <= curr_R) {
        // Elastic case
        computeElasticTangent(curr_tangent, props);
        if (idx == 0) {
            printf("Elastic case: Using elastic tangent operator\n");
        }
    } else {
        // Plastic case
        double nhat[6];
        for (int i = 0; i < 6; i++) {
            nhat[i] = XI[i] / a;
        }
        
        // Compute plastic multiplier
        double delta_lambda = (a - curr_R) / (2.0 * mu + props.H);
        
        // Compute tangent operator
        computePlasticTangent(curr_tangent, XI, a, nhat, delta_lambda, props);
        
        // Debug output for first element
        if (idx == 0) {
            printf("Plastic case:\n");
            printf("delta_lambda = %.6e\n", delta_lambda);
            printf("nhat = [%.6e, %.6e, %.6e, %.6e, %.6e, %.6e]\n",
                   nhat[0], nhat[1], nhat[2], nhat[3], nhat[4], nhat[5]);
        }
    }
}

// Host function: Allocate memory and launch kernel function
void launchTangentOperatorComputation(
    double* h_tangent_operator,
    const double* h_strain,
    const double* h_stress,
    const double* h_alpha,
    const double* h_R,
    const MaterialProperties& props,
    int num_elements
) {
    // Initialize all pointers to NULL
    double *d_tangent_operator = NULL;
    double *d_strain = NULL;
    double *d_stress = NULL;
    double *d_alpha = NULL;
    double *d_R = NULL;
    cudaError_t error;
    int num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;
    
    // Allocate device memory
    error = cudaMalloc(&d_tangent_operator, num_elements * 36 * sizeof(double));
    if (error != cudaSuccess) {
        printf("CUDA memory allocation failed (tangent_operator): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    error = cudaMalloc(&d_strain, num_elements * 6 * sizeof(double));
    if (error != cudaSuccess) {
        printf("CUDA memory allocation failed (strain): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    error = cudaMalloc(&d_stress, num_elements * 6 * sizeof(double));
    if (error != cudaSuccess) {
        printf("CUDA memory allocation failed (stress): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    error = cudaMalloc(&d_alpha, num_elements * 6 * sizeof(double));
    if (error != cudaSuccess) {
        printf("CUDA memory allocation failed (alpha): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    error = cudaMalloc(&d_R, num_elements * sizeof(double));
    if (error != cudaSuccess) {
        printf("CUDA memory allocation failed (R): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    // Copy data to device
    error = cudaMemcpy(d_strain, h_strain, num_elements * 6 * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess) {
        printf("CUDA memory copy failed (strain): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    error = cudaMemcpy(d_stress, h_stress, num_elements * 6 * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess) {
        printf("CUDA memory copy failed (stress): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    error = cudaMemcpy(d_alpha, h_alpha, num_elements * 6 * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess) {
        printf("CUDA memory copy failed (alpha): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    error = cudaMemcpy(d_R, h_R, num_elements * sizeof(double), cudaMemcpyHostToDevice);
    if (error != cudaSuccess) {
        printf("CUDA memory copy failed (R): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    // Launch kernel function
    computeTangentOperator<<<num_blocks, BLOCK_SIZE>>>(
        d_strain,
        d_stress,
        d_alpha,
        d_R,
        d_tangent_operator,
        props.E,
        props.nu,
        props.H,
        props.beta,
        props.sigma_y0,
        num_elements
    );
    
    // Check kernel execution errors
    error = cudaGetLastError();
    if (error != cudaSuccess) {
        printf("CUDA kernel function execution error: %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    // Wait for GPU computation to complete
    error = cudaDeviceSynchronize();
    if (error != cudaSuccess) {
        printf("CUDA synchronization error: %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    // Copy results back to host
    error = cudaMemcpy(h_tangent_operator, d_tangent_operator, num_elements * 36 * sizeof(double), cudaMemcpyDeviceToHost);
    if (error != cudaSuccess) {
        printf("CUDA memory copy failed (results): %s\n", cudaGetErrorString(error));
        goto cleanup;
    }
    
    printf("CUDA computation completed successfully!\n");

cleanup:
    // Safely free device memory
    if (d_tangent_operator != NULL) {
        cudaFree(d_tangent_operator);
        d_tangent_operator = NULL;
    }
    if (d_strain != NULL) {
        cudaFree(d_strain);
        d_strain = NULL;
    }
    if (d_stress != NULL) {
        cudaFree(d_stress);
        d_stress = NULL;
    }
    if (d_alpha != NULL) {
        cudaFree(d_alpha);
        d_alpha = NULL;
    }
    if (d_R != NULL) {
        cudaFree(d_R);
        d_R = NULL;
    }
    
    // Final error check
    cudaError_t lastError = cudaGetLastError();
    if (lastError != cudaSuccess) {
        printf("CUDA Error in launchTangentOperatorComputation: %s\n", 
               cudaGetErrorString(lastError));
    }
} 