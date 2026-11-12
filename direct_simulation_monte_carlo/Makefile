#Defining the compiler.
CC=g++

#Compilation time flags.
CFLAGS=-Wall -Ofast -funroll-loops -march=native -g

#Libraries to link into the executable.
LIBS=-lm

dsmc: particle.h cellList.h hardSphere.h hardSphere.cpp main.cpp 
	$(CC) $(CFLAGS) main.cpp -o dsmc $(LIBS)
