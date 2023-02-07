COMPILER = g++
INCLUDE := ${PROJECT_DIR}/src
CFLAGS := -fPIC -O2

.PHONY: clean

all: anal anal_raw

anal:
	${COMPILER} anal.cpp -o anal.exe -I${INCLUDE} `root-config --libs --cflags`
anal_raw : 
	${COMPILER} anal_raw.cpp -o anal_raw.exe -I${INCLUDE} `root-config --libs --cflags`
clean:
	rm *.exe
