CXXFLAGS=-std=c++1z -O3 -ffast-math -fPIC 
LDFLAGS=-${CXXFLAGS} -Wl,--export-dynamic -Wl,--no-as-needed -Wl,-rpath,'${JULIA_ROOT}/lib' -Wl,-rpath,'${JULIA_ROOT}/lib/julia' -ljulia -Wl,-rpath,${JLCXX_ROOT}/lib -lcxxwrap_julia -Wl,-rpath,${LIBINT_ROOT}/lib -lint2 -shared -fPIC

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
