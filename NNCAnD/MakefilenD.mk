CC			=/opt/homebrew/bin/g++-13
CFLAGS		=-c -Wall -O4 -fopenmp -std=c++17
LDFLAGS		=-fopenmp -std=c++17 -I$(EIGEN_PATH)
SOURCES		=./testFMMnD.cpp
OBJECTS		=$(SOURCES:.cpp=.o)
EXECUTABLE	=./testFMMnD

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
		$(CC) $(CFLAGS) $< -o $@

clean:
	rm a.out testFMMnD *.o
