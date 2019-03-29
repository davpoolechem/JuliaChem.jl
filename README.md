# JuliChem
A electronic structure theory program written in Julia, designed to combine
the ease of use of interactive Julia with the performance of statically
compiled Julia to offer users the best of both worlds for quantum chemistry
computations.

# How to run
JuliChem can be run in one of three ways:
1. Interactively, though the REPL. This can be done by loading the Main.JuliChem
module into the REPL and running

julia> JuliChem.julia_main(["path/to/data/file.dat"])

where file.dat is the data file, discussed in more detail below.

2. Through the shell as a statically compiled binary. The binary can be created
by simply running the build-julichem script, which builds the binary to
builddir/julichem. The binary can then be executed by running

builddir/julichem path/to/data/file.dat

in the shell, where file.dat is the data file.

3. Through the shell dynamically. This can be done by running

julia julichem_shell.jl path/to/data/file.dat

in the shell, where file.dat is the data file. Note that running JuliChem
through the shell dynamically yields significantly lower performance compared
to running JuliChem through the shell via statically compiled binary.

# Input information
A few steps must be taken to create a proper input for JuliChem:
1. JuliChem uses a .jl script to input calculation flag and geometry
information. Such a script must be created for the calculation to be performed.
Examples of .jl input files can be seen in the examples/ directory.

2. A JuliChem run pulls information from the .jl input file via input.jl,
acquiring information from the file specified by the "input_file" variable.
So, the value of input_file must be set to "path/to/desired/input.jl".

3. Finally, as of now, JuliChem also requires the manual entry of a few pieces
of data: the nuclear energy, and the one- and two- electron integral lists. These
must be included in a separate .dat data file, and the .dat file must be used as
the input argument for execution of the JuliChem program. Example data files
can be see in the examples/ directory.
