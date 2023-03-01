GCC := g++
INC := `root-config --cflags --libs`
VPATH := ${PROJECT_DIR}/src

.PHONY: all clean

all : amd2root filter_e15190

amd2root : amd2root.cpp
	${GCC} -o ${PROJECT_DIR}/bin/amd2root.exe ${PROJECT_DIR}/bin/amd2root ${INC}

filter_e15190 : filter_e15190.cpp
	${GCC} -o ${PROJECT_DIR}/bin/filter_e15190.exe ${PROJECT_DIR}/bin/filter_e15190 ${INC}

clean:
	rm -f ${PROJECT_DIR}/bin/*.o ${PROJECT_DIR}/bin/*.exe