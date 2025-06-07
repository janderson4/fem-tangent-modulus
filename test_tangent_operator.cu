#include <stdio.h>
#include <cuda_runtime.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <unistd.h>  // For getcwd

// Include tangent operator header
#include "tangent_operator.cuh"

// Function to read input data
bool readInputData(const char* filename, double* strain, double* stress, double* alpha, double* R, int num_elements) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error: Cannot open input file %s\n", filename);
        return false;
    }

    printf("Successfully opened %s, reading %d elements...\n", filename, num_elements);

    char line[512];
    int elements_read = 0;
    
    while (elements_read < num_elements && fgets(line, sizeof(line), file) != NULL) {
        // Skip empty lines
        if (strlen(line) <= 1) continue;
        
        int elem_index;
        
        // Try to read stress line
        if (sscanf(line, "Stress %d: %lf %lf %lf %lf %lf %lf",
            &elem_index, &stress[elements_read*6], &stress[elements_read*6+1], 
            &stress[elements_read*6+2], &stress[elements_read*6+3], 
            &stress[elements_read*6+4], &stress[elements_read*6+5]) == 7) {
            
            // Verify element index
            if (elem_index != elements_read + 1) {
                printf("Warning: Expected element %d, but read element %d\n", elements_read+1, elem_index);
            }
            
            // Read alpha line
            if (fgets(line, sizeof(line), file) == NULL ||
                sscanf(line, "Alpha %d: %lf %lf %lf %lf %lf %lf",
                &elem_index, &alpha[elements_read*6], &alpha[elements_read*6+1], 
                &alpha[elements_read*6+2], &alpha[elements_read*6+3], 
                &alpha[elements_read*6+4], &alpha[elements_read*6+5]) != 7) {
                printf("Error reading alpha data for element %d\n", elements_read+1);
                printf("Alpha line content: %s", line);
                fclose(file);
                return false;
            }
            
            // Read R line
            if (fgets(line, sizeof(line), file) == NULL ||
                sscanf(line, "R %d: %lf", &elem_index, &R[elements_read]) != 2) {
                printf("Error reading R data for element %d\n", elements_read+1);
                printf("R line content: %s", line);
                fclose(file);
                return false;
            }
            
            elements_read++;
            
            // Progress indicator
            if (elements_read % 100 == 0) {
                printf("Read %d elements so far...\n", elements_read);
            }
        }
    }

    fclose(file);
    
    if (elements_read != num_elements) {
        printf("Error: Expected %d elements, but only read %d elements\n", num_elements, elements_read);
        return false;
    }
    
    printf("Successfully read all %d elements from %s\n", num_elements, filename);
    return true;
}

// Save results
bool saveResults(
    const char* filename,
    const double* tangent_operator,
    int num_elements
) {
    printf("\nStarting to save results to file: %s\n", filename);
    FILE* file = fopen(filename, "w");
    if (!file) {
        printf("Error: Cannot create output file %s\n", filename);
        return false;
    }

    for (int i = 0; i < num_elements; i++) {
        fprintf(file, "Element %d:\n", i + 1);
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 6; k++) {
                fprintf(file, "%e ", tangent_operator[i*36 + j*6 + k]);
            }
            fprintf(file, "\n");
        }
        fprintf(file, "\n");
    }

    fclose(file);
    printf("Results saved to file\n");
    return true;
}

// Verify results
bool verifyResults(
    const double* tangent_operator,
    int num_elements
) {
    printf("\nStarting result verification...\n");
    
    // Check symmetry
    for (int i = 0; i < num_elements; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 6; k++) {
                double diff = fabs(tangent_operator[i*36 + j*6 + k] - 
                                 tangent_operator[i*36 + k*6 + j]);
                if (diff > 1e-10) {
                    printf("Warning: Tangent operator for element %d is not symmetric (position[%d,%d])\n", i, j, k);
                    return false;
                }
            }
        }
    }
    
    printf("Result verification passed\n");
    return true;
}

// Compute elastic strain
void computeElasticStrain(
    const double* stress,
    double* strain,
    double E,
    double nu,
    int num_elements
) {
    double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));
    
    for (int i = 0; i < num_elements; i++) {
        // Compute volumetric stress
        double p = (stress[i*6] + stress[i*6+1] + stress[i*6+2]) / 3.0;
        
        // Compute deviatoric stress
        double s[6];
        s[0] = stress[i*6] - p;
        s[1] = stress[i*6+1] - p;
        s[2] = stress[i*6+2] - p;
        s[3] = stress[i*6+3];
        s[4] = stress[i*6+4];
        s[5] = stress[i*6+5];
        
        // Compute strain
        strain[i*6] = s[0] / (2 * mu) + p / (3 * lambda);
        strain[i*6+1] = s[1] / (2 * mu) + p / (3 * lambda);
        strain[i*6+2] = s[2] / (2 * mu) + p / (3 * lambda);
        strain[i*6+3] = s[3] / (2 * mu);
        strain[i*6+4] = s[4] / (2 * mu);
        strain[i*6+5] = s[5] / (2 * mu);
    }
}

int main() {
    printf("Program execution started...\n");
    
    // Declare all variables at function start
    MaterialProperties props;
    
    // Scan input file to determine actual data count
    printf("Scanning input file to determine total data count...\n");
    FILE* scan_file = fopen("cuda_input.txt", "r");
    if (!scan_file) {
        printf("Error: Cannot open cuda_input.txt for scanning\n");
        return 1;
    }
    
    int num_elements = 0;
    char line[512];
    while (fgets(line, sizeof(line), scan_file)) {
        if (strstr(line, "Stress ") != NULL && strstr(line, ":") != NULL) {
            num_elements++;
        }
    }
    fclose(scan_file);
    
    printf("Found %d integration points in input file\n", num_elements);
    printf("Will process all %d integration points\n", num_elements);
    double *strain = NULL;
    double *stress = NULL;
    double *alpha = NULL;
    double *R = NULL;
    double *tangent_operator = NULL;
         FILE *output_file = NULL;  // Add output_file variable declaration
     const char* output_filename = "cuda_output.txt";  // Move to beginning declaration
    
    // Material parameter variable declarations
    double E = 200e3;        // 200 GPa
    double nu = 0.3;         // Poisson's ratio
    double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));
    
    // Directory path variable declarations
    char *cwd = NULL;
    char *cwd2 = NULL;
    
    
    
    // Print current working directory
    cwd = getcwd(NULL, 0);
    if (cwd) {
        printf("Current working directory: %s\n", cwd);
        free(cwd);
        cwd = NULL;
    } else {
        printf("Failed to get current directory\n");
    }
    
    // Set material parameters - exactly match MATLAB test case
    props.E = 200e3;        // 200 GPa - match MATLAB
    props.nu = 0.3;         // Poisson's ratio - match MATLAB
    props.H = 20e3;         // 20 GPa - match MATLAB
    props.beta = 0.0;       // Pure kinematic hardening - match MATLAB
    props.sigma_y0 = 300.0; // Not used, replaced by R
    
    printf("Material parameters set:\n");
    printf("E = %e\n", props.E);
    printf("nu = %e\n", props.nu);
    printf("H = %e\n", props.H);
    printf("beta = %e\n", props.beta);
    printf("sigma_y0 = %e\n", props.sigma_y0);

    printf("\nSetting number of elements: %d\n", num_elements);

    // Allocate memory
    strain = (double*)malloc(num_elements * 6 * sizeof(double));
    stress = (double*)malloc(num_elements * 6 * sizeof(double));
    alpha = (double*)malloc(num_elements * 6 * sizeof(double));
    R = (double*)malloc(num_elements * sizeof(double));
    tangent_operator = (double*)malloc(num_elements * 36 * sizeof(double));

    if (!strain || !stress || !alpha || !R || !tangent_operator) {
        printf("Error: Memory allocation failed\n");
        goto cleanup;
    }

    // Read input data
    printf("\n=== STEP 1: Reading input data from cuda_input.txt ===\n");
    if (!readInputData("cuda_input.txt", strain, stress, alpha, R, num_elements)) {
        printf("CRITICAL ERROR: Failed to read cuda_input.txt file!\n");
        printf("This file is required for proper comparison with MATLAB results.\n");
        printf("Please check if the file exists and has correct format.\n");
        cwd = getcwd(NULL, 0);
        if (cwd) {
            printf("Current directory: %s\n", cwd);
            free(cwd);
            cwd = NULL;
        }
        goto cleanup;
    }
    printf("✓ Step 1: Successfully read input data\n");
    
    printf("Successfully read input data from cuda_input.txt\n");
    
    // Calculate strain from stress data (like compare_results.m does)
    printf("Calculating strain from stress data...\n");
    
    for (int i = 0; i < num_elements; i++) {
        // Calculate volumetric stress
        double p = (stress[i*6] + stress[i*6+1] + stress[i*6+2]) / 3.0;
        
        // Calculate deviatoric stress and strain
        for (int j = 0; j < 3; j++) {
            double s = stress[i*6 + j] - p;
            strain[i*6 + j] = s / (2 * mu) + p / (3 * lambda);
        }
        for (int j = 3; j < 6; j++) {
            strain[i*6 + j] = stress[i*6 + j] / (2 * mu);
        }
    }
    printf("Strain calculation completed\n");

    // All elements use data from cuda_input.txt for accurate comparison
    printf("\n=== All elements using real data from cuda_input.txt ===\n");
    printf("First element data preview:\n");
    printf("Strain: ");
    for(int j = 0; j < 6; j++) printf("%.3e ", strain[j]);
    printf("\n");
    printf("Stress: ");
    for(int j = 0; j < 6; j++) printf("%.3e ", stress[j]);
    printf("\n");
    printf("Alpha: ");
    for(int j = 0; j < 6; j++) printf("%.3e ", alpha[j]);
    printf("\n");
    printf("R: %.3e\n", R[0]);

    // Start computation
    printf("\n=== STEP 2: Starting tangent operator computation ===\n");
    launchTangentOperatorComputation(
        tangent_operator,
        strain,
        stress,
        alpha,
        R,
        props,
        num_elements
    );
    printf("✓ Step 2: Computation completed\n");

    // Print first element result
    printf("\n=== Element 0 Tangent Operator Result ===\n");
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            printf("%.6e ", tangent_operator[i*6 + j]);
        }
        printf("\n");
    }

    // Verify results
    if (!verifyResults(tangent_operator, num_elements)) {
        printf("Warning: Result verification failed\n");
    }

    // Save results to file
    printf("\n=== STEP 3: Saving results to file ===\n");
    cwd2 = getcwd(NULL, 0);
    if (cwd2) {
        printf("Current directory before writing: %s\n", cwd2);
        free(cwd2);
        cwd2 = NULL;
    }
    printf("Attempting to create file: %s\n", output_filename);
    
    output_file = fopen(output_filename, "w");
    if (!output_file) {
        printf("Error: Failed to open output file '%s'\n", output_filename);
        perror("fopen error");
        goto cleanup;
    }
    
    printf("File opened successfully, writing results for all %d elements...\n", num_elements);
    
    // Write results for all elements in the format expected by compare_results.m
    for (int elem = 0; elem < num_elements; elem++) {
        fprintf(output_file, "Element %d:\n", elem + 1);
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                fprintf(output_file, "%e ", tangent_operator[elem*36 + i*6 + j]);
            }
            fprintf(output_file, "\n");
        }
        fprintf(output_file, "\n");
    }
    
            fflush(output_file);  // Ensure data is written to disk
        fclose(output_file);
        output_file = NULL;  // Set pointer to NULL to prevent double close
    printf("Complete results successfully saved to: %s\n", output_filename);
    printf("Total elements written: %d\n", num_elements);
    printf("✓ Step 3: File writing completed successfully\n");

    

    printf("\nProgram execution completed\n");

cleanup:
    // Free memory
    if (strain) free(strain);
    if (stress) free(stress);
    if (alpha) free(alpha);
    if (R) free(R);
    if (tangent_operator) free(tangent_operator);
    if (output_file) {  // Only close if file pointer is not NULL
        fclose(output_file);
        output_file = NULL;
    }
    // cwd and cwd2 are freed immediately after use, no need to free again here

    return 0;
} 