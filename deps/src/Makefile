CXX=g++

EIGENROOT=/shared/eigen/3.3.4
BOOSTROOT=/shared/boost/1.67.0
LIBINTROOT=/export/home/david/programs/install/libint

JULIAROOT=/export/home/david/programs/install/julia-1.5.0
JLCXXROOT=/export/home/david/programs/install/cxxwrap-julia

CXXFLAGS=-std=c++1z -O0 -g -I${EIGENROOT} -I${BOOSTROOT}/include -I'${JULIAROOT}/include/julia' -fPIC -L'${JULIAROOT}/lib' -Wl,--export-dynamic -Wl,-rpath,'${JULIAROOT}/lib' -Wl,-rpath,'${JULIAROOT}/lib/julia' -ljulia -I${JLCXXROOT}/include -L${JLCXXROOT}/lib -lcxxwrap_julia -L${LIBINTROOT}/lib -lint2 -shared -fPIC

EXE=jeri

all: $(EXE)

jeri: jeri.cpp
	$(CXX) $(CXXFLAGS) -c jeri.cpp 
	$(CXX) $(CXXFLAGS) jeri.o -o libjeri.so

.PHONY: clean
clean:
	rm -f $(EXE)
