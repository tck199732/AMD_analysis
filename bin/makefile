GCC := g++
INCLUDE := `root-config --cflags --libs`


SRC_DIR := ${PROJECT_DIR}/src
SUB_DIR := ${SRC_DIR}/e15190

VPATH := ${SRC_DIR} ${SUB_DIR}
INCLUDE += -I${SRC_DIR} -I${SRC_DIR}/e15190
SRC := ${shell find ${SRC_DIR} -name "*.cpp"}

.PHONY: all clean

all : amd2root filter_e15190

amd2root : amd2root.cpp 
	${GCC} -o $@.exe ${PROJECT_DIR}/bin/$< ${INCLUDE}

filter_e15190 : filter_e15190.cpp ${SRC} 
	${GCC} -o $@.exe ${PROJECT_DIR}/bin/$^ ${INCLUDE}

clean:
	rm -f ${PROJECT_DIR}/bin/*.o ${PROJECT_DIR}/bin/*.exe