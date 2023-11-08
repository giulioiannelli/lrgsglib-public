# Minimal makefile for hnn-main.c
#
# create directories
DIRS := src/LRGSG_package/dump 
$(shell mkdir -p $(DIRS))
#
PATH_SRC  = src/LRGSG_package/
PATH_SFMT = src/LRGSG_package/SFMT/
#
RBIMSIMULATOR0_PN = IsingSimulator0
RBIMSIMULATOR1_PN = IsingSimulator1
LRGSGSRCC = LRGSG_utils LRGSG_rbim sfmtrng
SFMTSRCC = SFMT
#
#
RBIMSIM0.c = $(addsuffix .c, $(RBIMSIMULATOR0_PN))
RBIMSIM1.c = $(addsuffix .c, $(RBIMSIMULATOR1_PN))
#
SRCCFILES.c = $(addsuffix .c, $(LRGSGSRCC))
SFMTFILES.c = $(addsuffix .c, $(SFMTSRCC))
#
SRCCRS0.c = ${RBIMSIM0.c} ${SRCCFILES.c}
SRCCRS1.c = ${RBIMSIM1.c} ${SRCCFILES.c}
#
PATHSRS0.c := $(addprefix $(PATH_SRC), $(SRCCRS0.c)) $(addprefix $(PATH_SFMT), $(SFMTFILES.c))
PATHSRS1.c := $(addprefix $(PATH_SRC), $(SRCCRS1.c)) $(addprefix $(PATH_SFMT), $(SFMTFILES.c))
#
# FILES.o = ${FILES.c:.c=.o}
#
PROGRAMN0  = $(addprefix $(PATH_SRC), $(RBIMSIMULATOR0_PN))
PROGRAMN1 = $(addprefix $(PATH_SRC), $(RBIMSIMULATOR1_PN))
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



all: ${PROGRAMN0} ${PROGRAMN1}

${PROGRAMN0}: ${PATHSRS0.c}
	${CC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

${PROGRAMN1}: ${PATHSRS1.c}
	${CC} ${ALLFLAGS} -o $@ $^ ${LMFLAG}

DEBRIS = a.out *~ 
RM_FR  = rm -fr

clean:
	${RM_FR} ${FILES.o} ${FILES.o} ${DEBRIS}




