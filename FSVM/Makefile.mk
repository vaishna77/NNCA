CC			=g++
KERNEL		= -DUSE_Matern # -DUSE_Matern, -DUSE_Gaussian
DIM		= -DUSE_DIM2 # -DUSE_DIM2, -DUSE_DIM4
MATVEC = -DUSE_AFMMnD  # -DUSE_directMatVec, -DUSE_AFMMnD
CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17 -I$(EIGEN_PATH)
LDFLAGS		=-fopenmp -std=c++17
SOURCES		=./testFSVM.cpp
OBJECTS		=$(SOURCES:.cpp=.o)
EXECUTABLE	=./testFSVM

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $(KERNEL) $(DIM) $(MATVEC) $< -o $@

clean:
	rm a.out testFSVM *.o
