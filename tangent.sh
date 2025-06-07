#!/bin/bash
#SBATCH -p gpu-turing
#SBATCH --gres gpu:1

### ---------------------------------------
### BEGINNING OF EXECUTION
### ---------------------------------------

echo "Starting at `date`"
echo
make

echo
echo Output from tangent operator
echo ----------------
./test_tangent_operator

echo
echo Profiling tangent operator
echo ----------------
nsys profile --trace=cuda,nvtx,osrt \
             --output=tangent_profile \
             --force-overwrite=true \
             ./test_tangent_operator 