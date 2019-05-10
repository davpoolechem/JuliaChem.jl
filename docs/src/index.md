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

JuliaChem uses the JSON format to output data from calculations, and future support
is planned for parsing data from JSON files for checkpointing and restarting calculations.
Thus, the JSON.jl package must also be downloaded and installed to use JuliaChem.

Support for SnoopCompile.jl exists within JuliaChem, so SnoopCompile can be taken advantage
of to reduce package loading times when running JuliaChem in the shell. However, this is
an optional dependency and is not required. Future support for compilation via PackageCompiler.jl
is planned; however, this is also optional and is not yet implemented yet.

Documenter.jl is required to view the JuliaChem documentation; but if you are reading
this, then you have already fulfilled that dependency.

## Running JuliaChem

The first step to running JuliaChem is generating an input file to process. The
input file is a .json file. It contains variables that serve as "flags"
that control specific parts of the calculation. Different sections of the json
input file correspond to different groups of flags.
Additionally, the input file contains the molecular coordinates of
the system, as well as the atomic number of the atom associated with each set
of coordinates. Optionally, the input file can also contain information to be
read in, including one- and two- electron integrals. Example input files can be
found in the example_inputs/ directory; and more information regarding the
specific groups of flags that can be set, can be found in the "Flags" section
of the documentation.

With the input file created, JuliaChem calculations can now be performed on
that input file. JuliaChem consists of a large amount of modules, each of
which have a specific role in a calculation. Full JuliaChem calculations
are performed by executing each desired module in a defined sequence. By
implementing the code this way, smaller module-based calculations can be chained
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
directory.

Finally, the MPI requirement must be accounted for. While Julia can be run interactively
over MPI, doing so comes with a few quirks. To run the Julia REPL over MPI, type the
following command into the shell:

```
mpirun -np <nprocs> --xterm <list of ranks, one per proc> julia
```

This will open a new Julia REPL for each proc. To fully initialize MPI, the
MPI.Init() function must then be executed in each Julia REPL. This is
easily handled by putting the line

```
import MPI; MPI.Init()
```

into your Julia startup file.

Deviating from these instructions leads to the aforementioned quirks:

-Not including the --xterm flag will cause the shell command to hang, as there is
no output for the Julia REPL to go to.
-Including only a subset of the processes in the --xterm flag list will cause
the REPL to open for each rank in the subset, but the REPLs will immediately
crash. This is the case even if MPI.Init() is included in the startup file.
I'll admit I have no idea why this is the case.

So while running Julia via the REPL over MPI required a bit of knowledge
to work, it can indeed be done.

To run a Julia script in Julia via the shell, type the following command into
the shell:
```
mpirun -np <nprocs> julia <path/to/script.jl>
```

If one wishes to take full advantage of the hybrid parallelism present within
JuliaChem, the JULIA_NUM_THREADS environmental variable can also be defined:

```
JULIA_NUM_THREADS=<nthreads> mpirun -np <nprocs> julia <path/to/script.jl>
```

For ease of use to the user, the above command exists in a shell script provided
to the user. The juliachem.sh shell script enables automatic marking of the input
file and execution of the above statement through the following shell command:

```
./juliachem.sh <path/to/script/file.jl> <path/to/input/file.jl> <nprocs> <nthreads>
```

# Modules

The different modules available for use in JuliaChem can be seen in the
Table of Contents above; clicking on a module will take you to its'
documentation section.

A given module Foo can be "executed" in the script by:
1. Adding "using Foo" to the module list in YourInputScript.jl
2. Adding Foo.run to the script() function in YourInputScript.jl.

In this manner, the different modules can be combined, building-block style, to
create a more expansive computation within a single calculation. Different
modules' run functions may require certain inputs and output certain outputs;
these are discussed in the module's individual documentation sections.

## Input

```@meta
CurrentModule = JCInput
```

```@docs
JCInput
run(args::String)
```
## Molecule

```@meta
CurrentModule = JCMolecule
```

```@docs
JCMolecule
run(input_info::Dict{String,Dict{String,Any}})
```

## RHF

```@meta
CurrentModule = JCRHF
```

```@docs
JCRHF
run(input_info::Dict{String,Dict{String,Any}}, basis::Basis)
```

## Properties

```@meta
CurrentModule = JCProperties
```

```@docs
JCProperties
run(scf::Data,input_info::Dict{String,Dict{String,Any}})
```

# Flags

The input file contains different flags that control various aspects of the
calculation. These flags can be divided into certain subsections, which can
be seen in the table of contents at the beginning of the manual. Clicking on
a section will take you to that section's available flags.

## Calculation Control Flags

```@meta
CurrentModule = JCStructs
```

```@docs
Ctrl_Flags
```

## Basis Set Flags

```@meta
CurrentModule = JCStructs
```

```@docs
Basis_Flags
```

## Hartree-Fock Flags

```@meta
CurrentModule = JCStructs
```

```@docs
SCF_Flags
```
