# Table of Contents
```@contents
```

# Basics

## Dependencies

The use of JuliaChem requires an MPI library to be downloaded. The core wave
function computation modules use MPI parallelism in their schemes,
so this necessitates an MPI library to take advantage of them.  

JuliaChem also uses the SIMINT integral library to perform electron repulsion 
integral calculations. This requires installation of both SIMINT itself; and
the CMake software package, which is used to build the interface between 
JuliaChem and SIMINT.

## Building JuliaChem

The first step to building JuliaChem is to download it. To do this, you
can download it from GitHub using the following Julia commands:

1. import Pkg

2. Pkg.add(PackageSpec(url="https://github.com/davpoolechem/JuliaChem.jl"))

This is download the JuliaChem package to your computer.

The next step is building the interface to SIMINT,
currently called Julia Electron Repulsion Integrals (JERI). This is done as 
follows:

1. Install SIMINT.

2. Define the environmental variable SIMINT, which is the directory of your
SIMINT installation.

3. Go to the deps/ folder and run "julia build.jl" 

4. Profit! If done correctly, this should compile correctly and lead to a
libjeri shared library in the deps folder. This library contains the
interface functions between JuliaChem and SIMINT.

5. Finally, define an environmental variable named JERIPATH containing the 
directory to libjeri.so, and add JERIPATH to your LD_LIBRARY_PATH environment
variable. 

Now, you can use JuliaChem in any script you wish simply by importing the 
JuliaChem module:

using JuliaChem

## Running JuliaChem

The first step to running JuliaChem is generating an input file to process. The
input file is a .json file. It contains variables that serve as "flags"
that control specific parts of the calculation. Different sections of the json
input file correspond to different groups of flags. Examples can be seen in
the example_inputs folder.

Currently, the nuclear repulsion energies and one-electron integrals
must be read in; however, work in underway to remove this requirement in
the future.

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
module calculation. However, running through the REPL prevents parallel
runs.

It is also possible, and recommended for general calculations,
to execute a sequence of modules via a handcrafted .jl script. The desired
calculation can then be called each time the .jl script is executed, either via
REPL or via shell. Example script files can be found in the example_scripts
directory.

To run a Julia script in Julia via the shell, type the following command into
the shell:
```
mpirun -np <nprocs> julia <path/to/script.jl>
```

# Modules

The different modules available for use in JuliaChem can be seen in the
Table of Contents above; clicking on a module will take you to its'
documentation section.

A given module Foo can be "executed" in the script by:
1. Adding "using Foo" to the module list in YourInputScript.jl
2. Adding Foo.run to the script() function in YourInputScript.jl.

In this manner, the different modules can be combined, building-block style, to
create a more expansive computation within a single script. Different
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
## Basis 

```@meta
CurrentModule = JCBasis
```

```@docs
JCBasis
run(molecule, model)

## RHF

```@meta
CurrentModule = JCRHF
```

```@docs
JCRHF
run(basis, molecule, keywords)
```

# Flags

The input file contains different flags that control various aspects of the
calculation. These flags can be divided into certain subsections, which can
be seen in the table of contents at the beginning of the manual. Clicking on
a section will take you to that section's available flags.
