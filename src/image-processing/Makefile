PROGRAM=image_processing
OBJS=utils/image.o utils/image_util.o main.o

CXX = g++
CXXFLAGS = -std=c++0x -Wno-unused-result -pipe -fopenmp -march=native -O3
#for O1
#CXXFLAGS = -std=c++0x -Wno-unused-result -pipe -fopenmp -march=native -O1

#for O2 without auto vectorization
#CXXFLAGS = -std=c++0x -Wno-unused-result -pipe -fopenmp -march=native -O2

#for O2 with auto vectorizion
#CXXFLAGS = -std=c++0x -Wno-unused-result -pipe -fopenmp -march=native -O2 -ftree-vectorize

#for asm
#CXXFLAGS = -std=c++0x -Wno-unused-result -fopenmp -g -O2 -march=native -ftree-vectorize -save-temps=obj #-pipe
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
	$(RM) $(PROGRAM) $(OBJS) depend.inc *.s *.ii

.PHONY: depend
depend: $(OBJS:.o=.cpp)
	-@ $(RM) depend.inc
	-@ for i in $^; do cpp -MM $$i | sed "s/\ [_a-zA-Z0-9][_a-zA-Z0-9]*\.cpp//g" >> depend.inc; done

-include depend.inc
