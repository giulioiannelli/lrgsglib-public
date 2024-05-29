# exports
CONDA_PREFIX = $(shell conda info --root)/envs/lrgsgenv
CONDA_BIN := $(CONDA_PREFIX)/bin
# Python includes and libraries
PYTHON_INC = $(shell python3 -m pybind11 --includes)
PYTHON_LIB = $(shell python3-config --ldflags)
export PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig:$(PKG_CONFIG_PATH)
# paths
PATH_BUILD = build/
PATH_DOCS =	docs/
PATH_DAT = data/
PATH_SRC = src/
PATH_TEST = test/
PATH_TOOLS = tools/

PATH_LRGSGLIB = $(PATH_SRC)lrgsglib/
PATH_SH = $(PATH_TOOLS)bash/
#
PATH_CCORE  = $(PATH_LRGSGLIB)Ccore/
PATH_GTPTCH = $(PATH_LRGSGLIB)gt_patches/
#
PATH_CCORE_BIN  = $(PATH_CCORE)bin/
PATH_GTPTCH_CPP = $(PATH_GTPTCH)cpp/
PATH_SFMT = $(PATH_CCORE)SFMT/
PATH_STATSYS = $(PATH_CCORE)statsys/
#
PATH_RBIM = $(PATH_STATSYS)RBIsingM/
PATH_SRW = $(PATH_STATSYS)signedRw/
PATH_VM = $(PATH_STATSYS)voterM/
#
PATH_RBIM_BASE = $(PATH_RBIM)base/
PATH_RBIM_SIMC = $(PATH_RBIM)simulatorC/
PATH_RBIM_STORE =$(PATH_RBIM)storer/
PATH_SRW_LATT = $(PATH_SRW)Lattices/
#
# directories
#
DIRS_TO_MAKE := $(PATH_DAT)
#
RBIMSIMULATOR0_PN = IsingSimulator0
RBIMSIMULATOR1_PN = IsingSimulator1
VMSIMULATOR_PN = voter_model
SRC_LRGSGLIB = LRGSG_utils sfmtrng 
SRC_RBIM = LRGSG_rbim
SRC_VM = LRGSG_vm
SFMTSRC = SFMT
#
CXXFLAGS = -O3 -Wall -shared -std=c++14 -fPIC
#
#
RBIMSIM0.c = $(addsuffix .c, $(RBIMSIMULATOR0_PN))
RBIMSIM1.c = $(addsuffix .c, $(RBIMSIMULATOR1_PN))
VMSIM.c = $(addsuffix .c, $(VMSIMULATOR_PN))
#
SRCCFILES.c = $(addsuffix .c, $(SRC_LRGSGLIB))
SRCCFILESRBIM.c = $(addsuffix .c, $(SRC_RBIM))
SRCCFILESVM.c = $(addsuffix .c, $(SRC_VM))
SFMTFILES.c = $(addsuffix .c, $(SFMTSRC))
#
SRCCRS0.c = ${RBIMSIM0.c}
SRCCRS1.c = ${RBIMSIM1.c}
SRCCVM1.c = ${VMSIM.c}
#
PATH_SRCC_FILES = $(addprefix $(PATH_CCORE), $(SRCCFILES.c))
PATH_SRCC_RBIM = $(addprefix $(PATH_RBIM_SIMC), $(SRCCFILESRBIM.c))
PATH_SFMT_FILES = $(addprefix $(PATH_SFMT), $(SFMTFILES.c))
PATH_SRCC_VM = $(addprefix $(PATH_VM), $(SRCCFILESVM.c))
PATHSRS0.c := $(addprefix $(PATH_RBIM_SIMC), $(RBIMSIM0.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_RBIM) $(PATH_SFMT_FILES)
PATHSRS1.c := $(addprefix $(PATH_RBIM_SIMC), $(RBIMSIM1.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_RBIM) $(PATH_SFMT_FILES)
PATHSRS2.c := $(addprefix $(PATH_VM), $(VMSIM.c)) $(PATH_SRCC_FILES) $(PATH_SRCC_VM) $(PATH_SFMT_FILES)
#
# FILES.o = ${FILES.c:.c=.o}
#
PROGRAMN0  = $(addprefix $(PATH_CCORE_BIN), $(RBIMSIMULATOR0_PN))
PROGRAMN1 = $(addprefix $(PATH_CCORE_BIN), $(RBIMSIMULATOR1_PN))
PROGRAMN2 = $(addprefix $(PATH_CCORE_BIN), $(VMSIMULATOR_PN))
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
INC_PATH1 = -I${PATH_CCORE}
INC_PATHS = ${INC_PATH1}
ALLFLAGS  = ${GFLAGS} ${OFLAGS} ${WFLAGS} ${DSFMTFLAG} ${INC_PATHS} 

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

chmod_scripts:
	find $(PATH_SH) -type f -name '*.sh' -exec chmod +x {} \;

make_rootf:
	echo "YES" > .isrootf

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

clean:
	${RM_FR} ${FILES.o} ${FILES.o} ${DEBRIS}
	$(MAKE) -C $(PATH_GTPTCH_CPP) clean
	$(MAKE) -C $(PATH_SRW_LATT) clean
	$(MAKE) -C $(PATH_RBIM_BASE) clean
	$(MAKE) -C $(PATH_RBIM_STORE) clean
	rm -f $(RW_TARGET)


# Prevent duplication of 'all' and 'clean' by removing previous definitions
.PHONY: all clean setup chmod_scripts create_dirs make_rootf sub_make