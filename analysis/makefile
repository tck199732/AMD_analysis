COMPILER = g++
INCLUDE := ${PROJECT_DIR}/src
CFLAGS := -fPIC -O2
# $(info    ${INCLUDE})

.PHONY: clean

all: anal

anal:
	${COMPILER} anal.cpp -o anal.exe -I${INCLUDE} `root-config --libs --cflags`

clean:
	rm *.exe
