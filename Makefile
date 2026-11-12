#Defining the compiler.
CC=g++

#Compilation time flags.
CFLAGS=-Wall --std=c++11 -funroll-loops -march=native -g -pg
## CFLAGS=-Wall -std=c++11 -Ofast -funroll-loops -march=native -g -pg 
#Libraries to link into the executable.
LIBS=-lm

edmd: particle.h hardSphere.h hardSphere.cpp main.cpp 
	$(CC) $(CFLAGS) main.cpp -o edmd $(LIBS)
