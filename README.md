# Welcome to JuliChem!
JuliChem is an electronic structure theory program written in Julia, designed to combine
Julia's dynamic and interactive nature with its high-performance capabilities to offer users the best of both worlds for quantum chemistry computations.

# How to run
JuliChem can be run in one of two ways:
1. Interactively, though the REPL. This can be done by loading the Main.JuliChem
module (located in src/core/julichem.jl) into the REPL and running

julia> JuliChem.julia_main()

2. Through the shell dynamically. This can be done by running

JULIA_NUM_THREADS=x julia julichem_shell.jl

in the shell, where x is the desired number of threads to execute the program
with.

# Input information
A few steps must be taken to create a proper input environment for JuliChem:
1. JuliChem uses a .jl script to input calculation flag and geometry
information. As of now, this .jl script also contains a few other pieces
of data: the nuclear energy, and the one- and two- electron integral lists.
Such a script must be created for the calculation to be performed.
Examples of .jl input files can be seen in the examples/ directory.

2. A JuliChem run pulls information from the .jl input file via InputFile.jl,
acquiring information from the file specified by the "input_file" variable.
So, the value of input_file must be set to "path/to/desired/input.jl".

3. JuliChem is highly modularized, allowing users to construct complex
calculation routines from simple building-block calculations. This is done
via filling out the script() function in InputScript.jl. Information about
filling out this script function, as well as the modules that can be used
to construct scripts, can be found in the documentation.

4. Once an input file is selected and an input script is written, JuliChem
can be executed through one of the two methods listed above.

# Documentation
JuliChem uses the Documenter.jl package to allow for the generation of its
documentation. This is done by going into the docs/ directory and running

julia make_xxx.jl

in the shell, where make_xxx can refer to one of two make scripts:
1. make_user.jl returns information about the different user modules that exist
within JuliChem, their role in a calculation script, and the function used to utilize
the module in the script. This is the recommended make script for users of JuliChem.

2. make_dev.jl returns information about all modules in JuliChem, and all functions
that exist for each module. This effectively returns documentation about every
component of the code entirely, and is thus the recommended make script for JuliChem
developers.
