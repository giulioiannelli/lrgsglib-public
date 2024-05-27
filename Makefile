# exports
CONDA_PREFIX = $(shell conda info --root)/envs/lrgsgenv
CONDA_BIN := $(CONDA_PREFIX)/bin
export PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig:$(PKG_CONFIG_PATH)
# directories
MAKEDIR := data/
# paths
PATH_SH = tools/bash
PATH_LRGSG = src/lrgsglib/
PATH_CCORE  = $(PATH_LRGSG)Ccore/
PATH_statsys = $(PATH_CCORE)statsys/
PATH_statsys_LATTICES = $(PATH_statsys)signedRw/Lattices/
PATH_statsys_rbim_base = $(PATH_statsys)RBIsingM/base/

PATH_GTPTCH = $(PATH_LRGSG)gt_patches/cpp/
PATH_SFMT = $(PATH_CCORE)SFMT/
#
RBIMSIMULATOR0_PN = IsingSimulator0
RBIMSIMULATOR1_PN = IsingSimulator1
VMSIMULATOR_PN = voter_model
LRGSGSRCC = LRGSG_utils sfmtrng 
ISINGLRGSRC = ${LRGSGSRCC} LRGSG_rbim
VOTERLRGSRC = ${LRGSGSRCC} LRGSG_vm
SFMTSRC = SFMT

CXXFLAGS = -O3 -Wall -shared -std=c++14 -fPIC

# Python includes and libraries
PYTHON_INC = $(shell python3 -m pybind11 --includes)
PYTHON_LIB = $(shell python3-config --ldflags)

# Targets
RW_TARGET = $(PATH_statsys)random_walk$(shell python3-config --extension-suffix)
RW_SOURCES = $(PATH_statsys)random_walk.cpp
#
#
RBIMSIM0.c = $(addsuffix .c, $(RBIMSIMULATOR0_PN))
RBIMSIM1.c = $(addsuffix .c, $(RBIMSIMULATOR1_PN))
VMSIM.c = $(addsuffix .c, $(VMSIMULATOR_PN))
#
SRCCFILES.c = $(addsuffix .c, $(ISINGLRGSRC))
SRCCFILESVM.c = $(addsuffix .c, $(VOTERLRGSRC))
SFMTFILES.c = $(addsuffix .c, $(SFMTSRC))
#
SRCCRS0.c = ${RBIMSIM0.c} ${SRCCFILES.c}
SRCCRS1.c = ${RBIMSIM1.c} ${SRCCFILES.c}
SRCCVM1.c = ${VMSIM.c} ${SRCCFILESVM.c}
#
PATHSRS0.c := $(addprefix $(PATH_CCORE), $(SRCCRS0.c)) $(addprefix $(PATH_SFMT), $(SFMTFILES.c))
PATHSRS1.c := $(addprefix $(PATH_CCORE), $(SRCCRS1.c)) $(addprefix $(PATH_SFMT), $(SFMTFILES.c))
PATHSRS2.c := $(addprefix $(PATH_CCORE), $(SRCCVM1.c)) $(addprefix $(PATH_SFMT), $(SFMTFILES.c))
#
# FILES.o = ${FILES.c:.c=.o}
#
PROGRAMN0  = $(addprefix $(PATH_CCORE), $(RBIMSIMULATOR0_PN))
PROGRAMN1 = $(addprefix $(PATH_CCORE), $(RBIMSIMULATOR1_PN))
PROGRAMN2 = $(addprefix $(PATH_CCORE), $(VMSIMULATOR_PN))
#
GCC := $(CONDA_PREFIX)/bin/gcc
CPP := $(CONDA_PREFIX)/bin/g++
GFLAGS    = -g
OFLAGS    = -O3
DSFMTFLAG = -DSFMT_MEXP=19937
LMFLAG    = -lm
WFLAG1    = -Wall
WFLAG2    = -Wextra
WFLAGS    = ${WFLAG1} ${WFLAG2}
INC_PATH1 = -Isrc/LRSG_package/Ccore
# INC_PATH2 = -Isrc
# INC_PATH3 = -Idep/SFMT
INC_PATHS = ${INC_PATH1} #${INC_PATH2} ${INC_PATH3}
ALLFLAGS  = ${GFLAGS} ${OFLAGS} ${WFLAGS} ${INC_PATHS} ${DSFMTFLAG} ${WFLAGS} ${INC_PATHS}

setup:
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


all: setup ${PROGRAMN0} ${PROGRAMN1} ${PROGRAMN2} chmod_scripts create_dirs make_rootf sub_make

${PROGRAMN0}: ${PATHSRS0.c}
	${GCC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN1}: ${PATHSRS1.c}
	${GCC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN2}: ${PATHSRS2.c}
	${GCC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

# $(RW_TARGET): $(RW_SOURCES)
# 	$(CPP) $(CXXFLAGS) $(PYTHON_INC) $(RW_SOURCES) -o $(RW_TARGET) $(PYTHON_LIB)

chmod_scripts:
	find $(PATH_SH) -type f -name '*.sh' -exec chmod +x {} \;

make_rootf:
	echo "YES" > .isrootf

# Target that creates necessary directories
create_dirs:
	@mkdir -p $(MAKEDIR)

sub_make:
	$(MAKE) -C $(PATH_GTPTCH)
	$(MAKE) -C $(PATH_statsys_LATTICES)
	$(MAKE) -C $(PATH_statsys_rbim_base)

DEBRIS = a.out *~ 
RM_FR  = rm -fr

clean:
	${RM_FR} ${FILES.o} ${FILES.o} ${DEBRIS}
	$(MAKE) -C $(PATH_GTPTCH) clean
	$(MAKE) -C $(PATH_statsys_LATTICES) clean
	$(MAKE) -C $(PATH_statsys_rbim_base) clean
	rm -f $(RW_TARGET)


# Prevent duplication of 'all' and 'clean' by removing previous definitions
.PHONY: all clean setup chmod_scripts create_dirs make_rootf sub_make