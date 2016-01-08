VPATH=../src
INCPATH=../include

CXX=g++
CXXFLAGS=-Wall -g -std=c++0x -O4

LIBLAPACK=-llapack -lblas

OBJS=lapack.o vector.o matrix.o sqmat.o davidson.o main.o

all: main

clean:
	rm -r $(OBJS) main

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -I$(INCPATH) -o $@ -c $<

lapack.o: $(INCPATH)/lapack.hpp
vector.o: $(INCPATH)/vector.hpp
matrix.o: $(INCPATH)/matrix.hpp
sqmat.o: $(INCPATH)/sqmat.hpp
davidson.o: $(INCPATH)/operator.hpp $(INCPATH)/davidson.hpp
main.o: vector.o matrix.o sqmat.o davidson.o

main: $(OBJS)
	$(CXX) $(CXXFLAGS) -o main $(OBJS) $(LIBLAPACK)
