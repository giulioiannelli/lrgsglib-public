# exports
export PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig:$(PKG_CONFIG_PATH)
# directories
MAKEDIR := data/
# paths
PATH_SH = tools/bash
PATH_LRGSG = src/LRGSG_package/
PATH_CCORE  = $(PATH_LRGSG)/Ccore/
PATH_GTPTCH = $(PATH_LRGSG)/gt_patches/cpp/
PATH_SFMT = $(PATH_CCORE)/SFMT/
#
RBIMSIMULATOR0_PN = IsingSimulator0
RBIMSIMULATOR1_PN = IsingSimulator1
VMSIMULATOR_PN = voter_model
LRGSGSRCC = LRGSG_utils sfmtrng 
ISINGLRGSRC = ${LRGSGSRCC} LRGSG_rbim
VOTERLRGSRC = ${LRGSGSRCC} LRGSG_vm
SFMTSRC = SFMT
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
CC        = gcc
GFLAGS    = -g
OFLAGS    = -O3
DSFMTFLAG = -DSFMT_MEXP=19937
LMFLAG    = -lm
WFLAG1    = -Wall
WFLAG2    = -Wextra
WFLAGS    = ${WFLAG1} ${WFLAG2}
INC_PATH1 = -Isrc/LRSG_package/
# INC_PATH2 = -Isrc
# INC_PATH3 = -Idep/SFMT
INC_PATHS = ${INC_PATH1} #${INC_PATH2} ${INC_PATH3}
ALLFLAGS  = ${GFLAGS} ${OFLAGS} ${WFLAGS} ${INC_PATHS} ${DSFMTFLAG} ${WFLAGS} ${INC_PATHS}



all: ${PROGRAMN0} ${PROGRAMN1} ${PROGRAMN2} chmod_scripts create_dirs make_rootf

${PROGRAMN0}: ${PATHSRS0.c}
	${CC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN1}: ${PATHSRS1.c}
	${CC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN2}: ${PATHSRS2.c}
	${CC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

chmod_scripts:
	find $(PATH_SH) -type f -name '*.sh' -exec chmod +x {} \;

make_rootf:
	echo "YES" > .isrootf

# Target that creates necessary directories
create_dirs:
	@mkdir -p $(MAKEDIR)

DEBRIS = a.out *~ 
RM_FR  = rm -fr

clean:
	${RM_FR} ${FILES.o} ${FILES.o} ${DEBRIS}


all:
	$(MAKE) -C $(PATH_GTPTCH)

clean:
	$(MAKE) -C $(PATH_GTPTCH) clean
