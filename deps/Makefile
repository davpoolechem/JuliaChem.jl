#CXXFLAGS=-std=c++1z -O3 -ffast-math -I${EIGEN_ROOT} -I${BOOST_ROOT}/include -I${JLCXX_ROOT}/include -I${JULIA_ROOT}/include/julia -fPIC 
CXXFLAGS=-std=c++1z -O3 -ffast-math -fPIC 

#LDFLAGS=-${CXXFLAGS} -L${JULIA_ROOT}/lib -Wl,--export-dynamic -Wl,-rpath,'${JULIA_ROOT}/lib' -Wl,-rpath,'${JULIA_ROOT}/lib/julia' -ljulia -L${JLCXX_ROOT}/lib -lcxxwrap_julia -L${LIBINT_ROOT}/lib -lint2 -shared -fPIC
LDFLAGS=-${CXXFLAGS} -L${JULIA_ROOT}/lib -Wl,--export-dynamic -Wl,-rpath,'${JULIA_ROOT}/lib' -Wl,-rpath,'${JULIA_ROOT}/lib/julia' -ljulia -lcxxwrap_julia -lint2 -shared -fPIC

EXE=jeri

all: $(EXE)

jeri: src/jeri.cpp
	$(CXX) $(CXXFLAGS) -c src/jeri.cpp 
	mv jeri.o src/
	$(CXX) $(LDFLAGS) src/jeri.o -o libjeri.so

.PHONY: clean
clean:
	rm libjeri.so
	rm src/jeri.o
