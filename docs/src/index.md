# Table of Contents
```@contents
```

# Basics

## Dependencies

The use of JuliaChem requires an MPI library to be downloaded. The core wave
function computation modules use hybrid MPI/threaded parallelism in their schemes,
so this necessitates an MPI library to take advantage of them. As MPI is
required to run JuliaChem, the MPI.jl package must also be installed via the Julia
package manager.

Aside from this, JuliaChem has no other external dependencies.

## Running JuliaChem

The first step to running JuliaChem is generating an input file to process. The
input file is itself a .jl Julia file. It contains initializations of
variables that serve as "flags" that control specific parts of the calculation.
These initializations occur within specific functions associated with a group
of flags. Additionally, the input file contains the molecular coordinates of
the system, as well as the atomic number of the atom associated with each set
of coordinates. Optionally, the input file can also contain information to be
read in, including one- and two- electron integrals. Example input files can be
found in the example_inputs/ directory; and more information regarding the
specific groups of flags that can be set, can be found in the "Flags" section
of the documentation.

Once an input file is created, the calculation can be performed on that input file.
JuliaChem directly executes the functions defined in the input file to read the
input file information, rather than parsing it. Thus, the input file must be
defined for the calculation in a special manner, which is done via the
JCInputFile.assign function. JCInputFile.assign can be called either via
the shell (using julia -e) or via the REPL with the command:

>import JCInputFile; JCInputFile.assign(path/to/input/file.jl)

This marks the selected input file as the file to use for JuliaChem computations.
The input file can be changed by simply rerunning the command above.

With the input file selected, JuliaChem calculations can now be performed on
that input file. JuliaChem consists of a large amount of modules, each of
which have a specific role in a calculation. Full JuliaChem calculations
are performed by executing each desired module in the desired sequence. By
implementing the code this way, smaller module calculations can be chained
together arbitrarily to execute more complex calculations, in a building-block
fashion. Details about specific modules and how to execute them can be seen
in the "Modules" section of the documentation.

Module sequences can be executed in multiple ways. It is possible to run the
desired sequence of modules in order as individual commands in the REPL. This
method allows for individual analysis and control of the outputs of each
module calculation. It is also possible, and recommended for general calculations,
to execute a sequence of modules via a handcrafted .jl script. The desired
calculation can then be called each time the .jl script is executed, either via
REPL or via shell. Example script files can be found in the example_scripts
directory. *Note that the call to JCInputFile.assign cannot occur within the
script file; these two must occur separately!*

Finally, the MPI requirement must be account for. In both running on the REPL and
via shell, the MPI library must be initialized before module calculations start,
and finalized once the module calculations are complete. If running module calculations
in sequence through the REPL, Julia must be started on multiple processes via mpirun
to do so:

>mpirun -np <nprocs> julia

When running a JuliaChem script via shell, this must also be executed via mpirun:

>mpirun -np <nprocs> julia

If one wished to take full advantage of the hybrid parallelism present within
JuliaChem, the JULIA_NUM_THREADS environmental variable can also be defined:

>JULIA_NUM_THREADS=<nthreads> mpirun -np <nprocs> julia

For ease of use to the user, the above command exists in a shell script provided
to the user. The juliachem.sh enables automatic marking of the input file
and execution of the above statement through the following shell command:

>./juliachem.sh </path/to/script/file.jl> <path/to/input/file.jl> <nprocs> <nthreads>

# Modules

The different modules available for use in JuliaChem can be seen in the
Table of Contents above; clicking on a module will take you to its'
documentation section.

A given module Foo can be "executed" in the script by:
1. Adding "using Foo" to the module list in InputScript.jl
2. Adding Foo.run to the script() function in InputScript.jl.

In this manner, the different modules can be combined, building-block style, to
create a more expansive computation within a single calculation. Different
modules' run functions may require certain inputs and output certain outputs;
these are discussed in the module's individual documentation sections.

## Input

```@meta
CurrentModule = Input
```

```@docs
Input
run()
```
## Molecule

```@meta
CurrentModule = Molecule
```

```@docs
Molecule
run(coord::Array{Float64,2})
```

## RHF

```@meta
CurrentModule = RHF
```

```@docs
RHF
run(flags::Flags)
```

## Properties

```@meta
CurrentModule = Properties
```

```@docs
Properties
run(scf::Data,flags::Flags)
```

# Flags

The input file contains different flags that control various aspects of the
calculation. These flags can be divided into certain subsections, which can
be seen in the table of contents at the beginning of the manual. Clicking on
a section will take you to that section's available flags.

## Basis Set Flags

```@meta
CurrentModule = InputStructs
```

```@docs
Basis_Flags
```

## Hartree-Fock Flags

```@meta
CurrentModule = InputStructs
```

```@docs
HF_Flags
```
