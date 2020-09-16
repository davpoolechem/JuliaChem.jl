CXXFLAGS=-std=c++1z -O3 -ffast-math -fPIC 
LDFLAGS=-${CXXFLAGS} -L${JULIA_ROOT}/lib -Wl,--export-dynamic -Wl,-rpath,'${JULIA_ROOT}/lib' -Wl,-rpath,'${JULIA_ROOT}/lib/julia' -ljulia -L${JLCXX_ROOT}/lib -lcxxwrap_julia -L${LIBINT_ROOT}/lib -lint2 -shared -fPIC

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
