# User defined parameters
KERNELS 	  			 = omp3
SUBPROJECT				 = cg
COMPILER    			 = INTEL
MPI								 = no
SILO 							 = no
DECOMP						 = TILES
ARCH_COMPILER_CC   = mpicc
ARCH_COMPILER_CPP  = mpic++
OPTIONS		  			 = -DENABLE_PROFILING 

# Compiler-specific flags
CFLAGS_INTEL			 = -qopenmp -no-prec-div -std=gnu99 -DINTEL \
										 -Wall -qopt-report=5 #-xhost
CFLAGS_INTEL_KNL	 = -O3 -qopenmp -no-prec-div -std=gnu99 -DINTEL \
										 -xMIC-AVX512 -Wall -qopt-report=5
CFLAGS_GCC				 = -std=gnu99 -fopenmp -march=native -Wall #-std=gnu99
CFLAGS_GCC_KNL   	 = -O3 -fopenmp -std=gnu99 \
										 -mavx512f -mavx512cd -mavx512er -mavx512pf #-fopt-info-vec-all
CFLAGS_GCC_POWER   = -O3 -mcpu=power8 -mtune=power8 -fopenmp -std=gnu99
CFLAGS_CRAY				 = -lrt -hlist=a
CFLAGS_XL					 = -O3 -qsmp=omp
CFLAGS_XL_OMP4		 = -qsmp -qoffload
CFLAGS_CLANG_OMP4  = -O3 -Wall -fopenmp-targets=nvptx64-nvidia-cuda -fopenmp-nonaliased-maps \
										 -fopenmp=libomp --cuda-path=$(CUDAROOT) -DCLANG
CFLAGS_PGI				 = -O3 -fast -mp

ifeq ($(KERNELS), cuda)
  CHECK_CUDA_ROOT = yes
endif
ifeq ($(COMPILER), CLANG_OMP4)
  CHECK_CUDA_ROOT = yes
endif

ifeq ($(CHECK_CUDA_ROOT), yes)
ifeq ("${CUDAROOT}", "")
$(error "$$CUDAROOT is not set, please set this to the root of your CUDA install.")
endif
endif

ifeq ($(DEBUG), yes)
  OPTIONS += -O0 -g -DDEBUG 
endif

ifeq ($(MPI), yes)
  OPTIONS += -DMPI
endif

ifeq ($(DECOMP), TILES)
OPTIONS += -DTILES
endif
ifeq ($(DECOMP), ROWS)
OPTIONS += -DROWS
endif
ifeq ($(DECOMP), COLS)
OPTIONS += -DCOLS
endif

# Default compiler
ARCH_LINKER    		= $(ARCH_COMPILER_CC)
ARCH_FLAGS     		= $(CFLAGS_$(COMPILER))
ARCH_LDFLAGS   		= $(ARCH_FLAGS) -lm 
ARCH_BUILD_DIR 		= ../obj/nodes/
ARCH_DIR       		= ..

ifeq ($(SILO), YES)
 OPTIONS += -DENABLE_SILO
 ARCH_LDFLAGS += -L/Applications/VisIt.app//Contents/Resources/2.10.2/darwin-x86_64/lib -lsiloh5
 ARCH_FLAGS += -I/Applications/VisIt.app//Contents/Resources/2.10.2/darwin-x86_64/include/silo/include/ 
endif

ifeq ($(KERNELS), cuda)
include Makefile.cuda
endif

# Get specialised kernels
SRC  			 = $(wildcard *.c)
SRC  			+= $(wildcard $(KERNELS)/$(SUBPROJECT)/*.c)
SRC  			+= $(wildcard $(ARCH_DIR)/$(KERNELS)/*.c)
SRC 			+= $(subst main.c,, $(wildcard $(ARCH_DIR)/*.c))
SRC_CLEAN  = $(subst $(ARCH_DIR)/,,$(SRC))
OBJS 			+= $(patsubst %.c, $(ARCH_BUILD_DIR)/%.o, $(SRC_CLEAN))

nodes: make_build_dir $(OBJS) Makefile
	$(ARCH_LINKER) $(OBJS) $(ARCH_LDFLAGS) -o nodes.$(KERNELS)

# Rule to make controlling code
$(ARCH_BUILD_DIR)/%.o: %.c Makefile 
	$(ARCH_COMPILER_CC) $(ARCH_FLAGS) $(OPTIONS) -c $< -o $@

$(ARCH_BUILD_DIR)/%.o: $(ARCH_DIR)/%.c Makefile 
	$(ARCH_COMPILER_CC) $(ARCH_FLAGS) $(OPTIONS) -c $< -o $@

make_build_dir:
	@mkdir -p $(ARCH_BUILD_DIR)/
	@mkdir -p $(ARCH_BUILD_DIR)/$(KERNELS)/$(SUBPROJECT)

clean:
	rm -rf $(ARCH_BUILD_DIR)/* nodes.$(KERNELS) *.vtk *.bov *.dat \
		*.optrpt *.cub *.ptx *.silo

