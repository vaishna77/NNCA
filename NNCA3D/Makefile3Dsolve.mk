CC			=g++
CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17 -I$(EIGEN_PATH)
LDFLAGS		=-fopenmp -std=c++17
SOURCES		=./testFMM3Dsolve.cpp
OBJECTS		=$(SOURCES:.cpp=.o)
EXECUTABLE	=./testFMM3Dsolve

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $< -o $@

clean:
	rm a.out testFMM3Dsolve *.o
