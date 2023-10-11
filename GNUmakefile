# Minimal makefile for hnn-main.c

PATH_SRC  = src/LRGSG_package/
PATH_SFMT = src/LRGSG_package/SFMT/
SRCCFILES.c = IsingSimulator.c LRGSG_utils.c sfmtrng.c
SFMTFILES.c = SFMT.c

SUBDIRS := dep/ising-model/

FILES.c = ${SRCCFILES.c} ${SFMTFILES.c}
FILES.o = ${FILES.c:.c=.o}

PATHS.c := $(addprefix $(PATH_SRC), $(SRCCFILES.c)) $(addprefix $(PATH_SFMT), $(SFMTFILES.c))

# create directories
DIRS := data/l2d_sq_ising/graphs/ src/LRGSG_package/dump 
$(shell mkdir -p $(DIRS))

PROGRAMN  = $(addprefix $(PATH_SRC), IsingSimulator)
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



all: ${PROGRAMN}
${PROGRAMN}: ${PATHS.c}
	${CC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}


DEBRIS = a.out *~ 
RM_FR  = rm -fr

clean:
	${RM_FR} ${FILES.o} ${FILES.o} ${DEBRIS}




