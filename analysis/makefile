COMPILER = g++
INCLUDE := `root-config --libs --cflags`

SRC_DIR := ${PROJECT_DIR}/src
SRC := ${shell find ${SRC_DIR} -name "*.cpp"}

VPATH := ${SRC_DIR} ${SRC_DIR}/histograms
INCLUDE += ${addprefix -I, ${VPATH}}

.PHONY: all clean

all: anal_PtRapidity anal_Centrality anal_EmissionTime

% : %.cpp ${SRC}
	${COMPILER} $^ -o $@.exe ${INCLUDE} 

clean:
	rm *.exe
