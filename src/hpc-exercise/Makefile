PROGRAM=hpc_exercise
OBJS=utils/mat.o utils/mat_util.o main.o

CXX = g++
CXXFLAGS = -std=c++0x -fopenmp -Wno-unused-result -march=native #-mtune=native #-O3
#for Mac
#CXX = clang++
#CXXFLAGS = -std=c++0x -Xpreprocessor -fopenmp -lomp -Wno-unused-result -Ofast -march=native -mtune=native

.SUFFIXES: .cpp .o

.PHONY: all
all: depend $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(PROGRAM) $^

.c.o:
	$(CXX) $(CXXFLAGS) -c $<

.PHONY: clean
clean:
	$(RM) $(PROGRAM) $(OBJS) depend.inc

.PHONY: depend
depend: $(OBJS:.o=.cpp)
	-@ $(RM) depend.inc
	-@ for i in $^; do cpp -MM $$i | sed "s/\ [_a-zA-Z0-9][_a-zA-Z0-9]*\.cpp//g" >> depend.inc; done

-include depend.inc
