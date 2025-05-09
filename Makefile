# conda environment name
CONDA_ENV_NAME ?= lrgsgenv
# paths
PATH_BUILD = build/
PATH_DOCS  = docs/
PATH_DAT   = data/
PATH_SRC   = src/
PATH_TEST  = test/
PATH_TOOLS = tools/
#
PATH_LRGSGLIB = $(PATH_SRC)lrgsglib/
PATH_SH       = $(PATH_TOOLS)bash/
#
PATH_CCORE  = $(PATH_LRGSGLIB)Ccore/
PATH_GTPTCH = $(PATH_LRGSGLIB)gt_patches/
#
PATH_CCORE_BIN  = $(PATH_CCORE)bin/
PATH_GTPTCH_CPP = $(PATH_GTPTCH)cpp/
PATH_SFMT       = $(PATH_CCORE)SFMT/
PATH_STATSYS    = $(PATH_CCORE)statsys/
#
PATH_RBIM = $(PATH_STATSYS)RBIsingM/
PATH_SRW  = $(PATH_STATSYS)signedRw/
PATH_VM   = $(PATH_STATSYS)voterM/
#
PATH_RBIM_BASE  = $(PATH_RBIM)base/
PATH_RBIM_SIMC  = $(PATH_RBIM)simulatorC/
PATH_RBIM_STORE = $(PATH_RBIM)storer/
PATH_SRW_LATT   = $(PATH_SRW)Lattices/
#
DIRS_TO_MAKE = $(PATH_DAT) $(PATH_CCORE_BIN)
# conda paths
CONDA_PREFIX = $(shell conda info --root)/envs/lrgsgenv
CONDA_BIN    = $(CONDA_PREFIX)/bin
# conda activation
ACTIVATE_D   = $(CONDA_PREFIX)/etc/conda/activate.d
DEACTIVATE_D = $(CONDA_PREFIX)/etc/conda/deactivate.d
#
CONFIG_SCRIPT_GEN        = $(PATH_SH)generate_config.sh
CONFIG_SCRIPT_PATH       = $(PATH_SH)config_env.sh
UNCONFIG_SCRIPT_PATH     = $(PATH_SH)unconfig_env.sh
CUSTOM_ACTIVATE_SCRIPT   = $(ACTIVATE_D)/custom_env_setup.sh
CUSTOM_DEACTIVATE_SCRIPT = $(DEACTIVATE_D)/custom_env_cleanup.sh
# python includes and libraries
PYTHON_INC = $(shell python3 -m pybind11 --includes)
PYTHON_LIB = $(shell python3-config --ldflags)
export PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig#:$(PKG_CONFIG_PATH)
# C filenames
FN_RBIMSIM0 = IsingSimulator0
FN_RBIMSIM1 = IsingSimulator1
FN_RBIMSIM2 = IsingSimulator2
FN_RBIMSIM3 = IsingSimulator3
FN_RBIMSIM4 = IsingSimulator4
FN_RBIMSIM5 = IsingSimulator5
FN_VMSIM0   = voter_model
FN_LRGSGLIB = LRGSG_utils sfmtrng 
SRC_RBIM    = LRGSG_rbim
SRC_VM      = LRGSG_vm
SFMTSRC     = SFMT
# 
RBIMSIM0.c = $(addsuffix .c, $(FN_RBIMSIM0))
RBIMSIM1.c = $(addsuffix .c, $(FN_RBIMSIM1))
RBIMSIM2.c = $(addsuffix .c, $(FN_RBIMSIM2))
RBIMSIM3.c = $(addsuffix .c, $(FN_RBIMSIM3))
RBIMSIM4.c = $(addsuffix .c, $(FN_RBIMSIM4))
RBIMSIM5.c = $(addsuffix .c, $(FN_RBIMSIM5))
VMSIM0.c   = $(addsuffix .c, $(FN_VMSIM0))
#
SRCCFILES.c     = $(addsuffix .c, $(FN_LRGSGLIB))
SRCCFILESRBIM.c = $(addsuffix .c, $(SRC_RBIM))
SRCCFILESVM.c   = $(addsuffix .c, $(SRC_VM))
SFMTFILES.c     = $(addsuffix .c, $(SFMTSRC))
#
PATH_SRCC_FILES = $(addprefix $(PATH_CCORE), $(SRCCFILES.c))
PATH_SRCC_RBIM  = $(addprefix $(PATH_RBIM_SIMC), $(SRCCFILESRBIM.c))
PATH_SFMT_FILES = $(addprefix $(PATH_SFMT), $(SFMTFILES.c))
PATH_SRCC_VM    = $(addprefix $(PATH_VM), $(SRCCFILESVM.c))
#
PATHSRS0.c := $(addprefix $(PATH_RBIM_SIMC), $(RBIMSIM0.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_RBIM) $(PATH_SFMT_FILES)
PATHSRS1.c := $(addprefix $(PATH_RBIM_SIMC), $(RBIMSIM1.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_RBIM) $(PATH_SFMT_FILES)
PATHSRS2.c := $(addprefix $(PATH_RBIM_SIMC), $(RBIMSIM2.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_RBIM) $(PATH_SFMT_FILES)
PATHSRS3.c := $(addprefix $(PATH_RBIM_SIMC), $(RBIMSIM3.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_RBIM) $(PATH_SFMT_FILES)
PATHSRS4.c := $(addprefix $(PATH_RBIM_SIMC), $(RBIMSIM4.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_RBIM) $(PATH_SFMT_FILES)
PATHSRS5.c := $(addprefix $(PATH_RBIM_SIMC), $(RBIMSIM5.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_RBIM) $(PATH_SFMT_FILES)
PATHSRSVM0.c := $(addprefix $(PATH_VM), $(VMSIM0.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_VM) $(PATH_SFMT_FILES)
#
# FILES.o = ${FILES.c:.c=.o}
#
PROGRAMN0 = $(addprefix $(PATH_CCORE_BIN), $(FN_RBIMSIM0))
PROGRAMN1 = $(addprefix $(PATH_CCORE_BIN), $(FN_RBIMSIM1))
PROGRAMN2 = $(addprefix $(PATH_CCORE_BIN), $(FN_RBIMSIM2))
PROGRAMN3 = $(addprefix $(PATH_CCORE_BIN), $(FN_RBIMSIM3))
PROGRAMN4 = $(addprefix $(PATH_CCORE_BIN), $(FN_RBIMSIM4))
PROGRAMN5 = $(addprefix $(PATH_CCORE_BIN), $(FN_RBIMSIM5))
PROGRAMNVM0 = $(addprefix $(PATH_CCORE_BIN), $(FN_VMSIM0))
#
GCC := $(CONDA_PREFIX)/bin/gcc
CPP := $(CONDA_PREFIX)/bin/g++
CXXFLAGS      = -O3 -Wall -shared -std=c++14 -fPIC
GFLAGS        = -g
OFLAGS        = -O3
DSFMTFLAG     = -DSFMT_MEXP=19937
LMFLAG        = -lm
WFLAG1        = -Wall
WFLAG2        = -Wextra
OPENMPFLAG    = -fopenmp
WFLAGS        = ${WFLAG1} ${WFLAG2}
INC_PATH1     = -I${PATH_CCORE}
INC_PATH_SFMT = -I${PATH_SFMT}
INC_PATHS     = ${INC_PATH1} ${INC_PATH_SFMT}
ALLFLAGS      = ${GFLAGS} ${OFLAGS} ${WFLAGS} ${DSFMTFLAG} ${INC_PATHS} ${OPENMPFLAG}

configure: setup chmod_scripts echo_paths create_dirs
c_make: ${PROGRAMN0} ${PROGRAMN1} ${PROGRAMN2} ${PROGRAMN3} ${PROGRAMN4} ${PROGRAMN5}
all: configure c_make #sub_make chmod_scripts

generate_config_script:
	@echo "Generating config script..."
	@chmod +x $(CONFIG_SCRIPT_GEN)
	@bash $(CONFIG_SCRIPT_GEN)

make_rootf:
	@echo "Creating .isrootf file..."
	@echo "YES" > .isrootf

setup_conda_activate:
	@echo "Creating activate.d directory..."
	@mkdir -p $(ACTIVATE_D)
	@echo "Creating custom activation script..."
	@echo 'source $(CONFIG_SCRIPT_PATH)' > $(CUSTOM_ACTIVATE_SCRIPT)
	@echo "Making custom activation script executable..."
	@chmod +x $(CUSTOM_ACTIVATE_SCRIPT)
	@echo "Setup conda environment activate.d complete."

# Setup target to create deactivation script
setup_conda_deactivate:
	@echo "Creating deactivate.d directory..."
	@mkdir -p $(DEACTIVATE_D)
	@echo "Creating custom deactivation script..."
	@echo 'source $(UNCONFIG_SCRIPT_PATH)' > $(CUSTOM_DEACTIVATE_SCRIPT)
	@echo "Making custom deactivation script executable..."
	@chmod +x $(CUSTOM_DEACTIVATE_SCRIPT)
	@echo "Setup conda environment deactivate.d complete."

setup: make_rootf generate_config_script setup_conda_activate setup_conda_deactivate create_dirs
	@# Find the gcc compiler binary
	$(eval GCC_BIN := $(shell find $(CONDA_BIN) -name 'x86_64-conda_cos*-linux-gnu-cc' | head -n 1))
	@# Find the g++ compiler binary
	$(eval GPP_BIN := $(shell find $(CONDA_BIN) -name 'x86_64-conda_cos*-linux-gnu-cpp' | head -n 1))

	@# Remove existing gcc symlink if it exists and create a new one
	@if [ -L $(CONDA_BIN)/gcc ]; then \
	    rm $(CONDA_BIN)/gcc; \
	fi
	@ln -s $(GCC_BIN) $(CONDA_BIN)/gcc;

	@# Remove existing g++ symlink if it exists and create a new one
	@if [ -L $(CONDA_BIN)/g++ ]; then \
	    rm $(CONDA_BIN)/g++; \
	fi
	@ln -s $(GPP_BIN) $(CONDA_BIN)/g++;

	@echo "Using GCC at $(GCC_BIN)"
	@echo "Using G++ at $(GPP_BIN)"
	@if [ -z "$(GCC_BIN)" ] || [ -z "$(GPP_BIN)" ]; then \
		echo "Error: Unable to find required compilers in $(CONDA_BIN)"; exit 1; \
	fi

echo_paths:
	@echo "LRGSG_ROOT = $$LRGSG_ROOT"
	@echo "LRGSG_BUILD = $$LRGSG_BUILD"
	@echo "LRGSG_DATA = $$LRGSG_DATA"
	@echo "LRGSG_IPYNB = $$LRGSG_IPYNB"
	@echo "LRGSG_SRC = $$LRGSG_SRC"
	@echo "LRGSG_LIB = $$LRGSG_LIB"
	@echo "LRGSG_LIB_CCORE = $$LRGSG_LIB_CCORE"
	@echo "LRGSG_LIB_CBIN = $$LRGSG_LIB_CBIN"
	@echo "LRGSG_LIB_CONFIG = $$LRGSG_LIB_CONFIG"
	@echo "LRGSG_LIB_GT_PATCHES = $$LRGSG_LIB_GT_PATCHES"
	@echo "LRGSG_LIB_NX_PATCHES = $$LRGSG_LIB_NX_PATCHES"
	@echo "LRGSG_LIB_STOCPROC = $$LRGSG_LIB_STOCPROC"
	@echo "LRGSG_TEST = $$LRGSG_TEST"
	@echo "LRGSG_TOOLS = $$LRGSG_TOOLS"

${PROGRAMN0}: ${PATHSRS0.c}
	${GCC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN1}: ${PATHSRS1.c}
	${GCC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN2}: ${PATHSRS2.c}
	${GCC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN3}: ${PATHSRS3.c}
	${GCC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN4}: ${PATHSRS4.c}
	${GCC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN5}: ${PATHSRS5.c}
	${GCC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

chmod_scripts:
	find $(PATH_SH) -type f -name '*.sh' -exec chmod +x {} \;

# Target that creates necessary directories
create_dirs:
	@mkdir -p $(DIRS_TO_MAKE)

sub_make:
	$(MAKE) -C $(PATH_GTPTCH_CPP)
	$(MAKE) -C $(PATH_SRW_LATT)
	$(MAKE) -C $(PATH_RBIM_BASE)
	$(MAKE) -C $(PATH_RBIM_STORE)

DEBRIS = a.out *~ 
RM_FR  = rm -fr

clean_programs:
	${RM_FR} ${PROGRAMN0}
	${RM_FR} ${PROGRAMN1}
	${RM_FR} ${PROGRAMN2}
	${RM_FR} ${PROGRAMN3}
	${RM_FR} ${PROGRAMN4}
	${RM_FR} ${PROGRAMN5}

clean_directories:
	$(MAKE) -C $(PATH_GTPTCH_CPP) clean
	$(MAKE) -C $(PATH_SRW_LATT) clean
	$(MAKE) -C $(PATH_RBIM_BASE) clean
	$(MAKE) -C $(PATH_RBIM_STORE) clean

clean_files:
	rm -f $(RW_TARGET)
	rm -f *.o main

clean_activations:
	@rm -rf $(ACTIVATE_D) $(DEACTIVATE_D)

clean: clean_programs clean_directories clean_files clean_activations

# Prevent duplication of 'all' and 'clean' by removing previous definitions
.PHONY: all clean setup chmod_scripts create_dirs make_rootf sub_make