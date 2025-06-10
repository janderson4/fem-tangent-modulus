NVCC = nvcc
# Basic optimization flags
NVCC_FLAGS = -O3 -arch=sm_60 -diag-suppress 177
# Debug and profiling flags (add lineinfo for better Nsight analysis)
PROFILE_FLAGS = -lineinfo -g

TARGET = test_tangent_operator
TARGET_PROFILE = test_tangent_operator_profile
SRCS = test_tangent_operator.cu tangent_operator.cu
OBJS = $(SRCS:.cu=.o)
OBJS_PROFILE = $(SRCS:.cu=_profile.o)

all: $(TARGET)

profile: $(TARGET_PROFILE)

# Regular optimized build
$(TARGET): $(OBJS)
	$(NVCC) $(NVCC_FLAGS) -o $@ $^

# Profile build with debug info for Nsight
$(TARGET_PROFILE): $(OBJS_PROFILE)
	$(NVCC) $(NVCC_FLAGS) $(PROFILE_FLAGS) -o $@ $^

%.o: %.cu
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@

%_profile.o: %.cu
	$(NVCC) $(NVCC_FLAGS) $(PROFILE_FLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(TARGET_PROFILE) $(OBJS) $(OBJS_PROFILE)

run: $(TARGET)
	./$(TARGET)

# Profile runs
profile-compute: $(TARGET_PROFILE)
	mkdir -p ../nsight_results
	export TMPDIR=/tmp/nsight-compute-lock-$$USER; \
	mkdir -p "$$TMPDIR"; \
	ncu -o ../nsight_results/tangent_compute_report \
	    -f \
	    --set full \
	    --kernel-name "computeTangentOperator" \
	    --launch-skip-before-match 0 \
	    ./$(TARGET_PROFILE)

profile-systems: $(TARGET_PROFILE)
	mkdir -p ../nsight_results
	nsys profile --trace=cuda,nvtx,osrt \
	             --output=../nsight_results/tangent_timeline \
	             --force-overwrite=true \
	             ./$(TARGET_PROFILE)

.PHONY: all clean run profile profile-compute profile-systems
