GCC := g++
INC := `root-config --cflags --libs`
VPATH := ${PROJECT_DIR}/src
INC += -I${VPATH}

.PHONY: all clean

all : amd2root filter_e15190

amd2root : amd2root.cpp 
	${GCC} -o ${PROJECT_DIR}/bin/$@.exe ${PROJECT_DIR}/bin/$< ${INC}

filter_e15190 : filter_e15190.cpp 
	${GCC} -o ${PROJECT_DIR}/bin/$@.exe ${PROJECT_DIR}/bin/$< ${INC}

clean:
	rm -f ${PROJECT_DIR}/bin/*.o ${PROJECT_DIR}/bin/*.exe