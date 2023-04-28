COMPILER = g++
INCLUDE := `root-config --libs --cflags`

SRC_DIR := ${PROJECT_DIR}/src 
SRC := ${shell find ${SRC_DIR} -name "*.cpp"}


VPATH := ${SRC_DIR}
INCLUDE += -I${SRC_DIR}

.PHONY: all clean

all: anal

anal : anal.cpp Histograms.cpp ${SRC}
	${COMPILER} $^ -o $@.exe ${INCLUDE} 

clean:
	rm *.exe
