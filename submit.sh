#!/bin/bash
#SBATCH -p gpu-turing
#SBATCH --gres gpu:1
#SBATCH --output=tangent_operator_%j.out
#SBATCH --error=tangent_operator_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G

### ---------------------------------------
### BEGINNING OF EXECUTION
### ---------------------------------------

echo "Starting at $(date '+%Y-%m-%d %H:%M:%S')"
echo

# Show available modules
echo "Available modules:"
module avail

# Load CUDA 12.4
echo "Loading CUDA 12.4..."
module load cuda/12.4

# Show loaded modules
echo "Loaded modules:"
module list

# Show CUDA compiler path
echo "CUDA compiler path:"
which nvcc

# Clean and compile
make clean
make

# Run program
./test_tangent_operator

echo "Ending at $(date '+%Y-%m-%d %H:%M:%S')"
echo "Computation completed, results saved to cuda_results.txt" 