COMPILER = g++
INCLUDE := `root-config --libs --cflags`
VPATH := ${PROJECT_DIR}/src
INCLUDE += -I${VPATH}

.PHONY: all clean

all: anal

% : %.cpp
	${COMPILER} $< -o $@.exe -I${INCLUDE} 

clean:
	rm *.exe
