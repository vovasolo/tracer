CXX = g++
CXXFLAGS = -O2 -Wall -std=c++11 -march=native
INCLUDES = -I/usr/include/eigen3 
ROOTINCLUDES = -I$(shell root-config --incdir)
ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

example: main.o tracer.o
	$(CXX) -o example main.o tracer.o $(ROOTLIBS) -lGeom -lGeomPainter -lGeomBuilder

main.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(ROOTINCLUDES) main.cpp
	
tracer.o:
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(ROOTINCLUDES) tracer.cpp
	
clean:
	rm -f *.o

