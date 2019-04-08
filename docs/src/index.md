# Table of Contents
```@contents
```

# Modules

The different modules available for use in JuliChem can be seen in the
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
