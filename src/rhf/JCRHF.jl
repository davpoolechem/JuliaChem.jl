"""
     module JCRHF
The module required for computation of the wave function using the *Restricted
Hartree-Fock* (RHF) method in a Self-Consistent Field (SCF) calculation. This
module will be used often, as the RHF wave function is often the zeroth-order
wave function for closed-shell systems.
"""
module JCRHF

Base.include(@__MODULE__,"RHFSCF.jl")

using JCStructs

using MPI
using JSON

"""
    run(flags::Flags)

Execute the JuliaChem RHF algorithm.

One input variable is required:
1. flags = The calculation flags from the input file.

One variable is output:
1. scf = Data saved from the SCF calculation.

Thus, proper use of the RHF.run() function would look like this:

```
scf = RHF.run(flags)
```
"""
function run(input_info::Dict{String,Dict{String,Any}})
    comm=MPI.COMM_WORLD

    if (MPI.Comm_rank(comm) == 0)
        println("--------------------------------------------------------------------------------------")
        println("                       ========================================          ")
        println("                         RESTRICTED CLOSED-SHELL HARTREE-FOCK            ")
        println("                       ========================================          ")
        println("")
    end

    #set up rhf flags
    basis_info::Dict{String,Any} = input_info["Basis Flags"]
    basis_flags::Basis_Flags = Basis_Flags(basis_info["norb"], basis_info["nocc"])

    scf_info::Dict{String,Any} = input_info["SCF Flags"]
    scf_flags::SCF_Flags = SCF_Flags(scf_info["niter"], scf_info["dele"],
                                        scf_info["rmsd"], scf_info["direct"])

    rhf_flags::RHF_Flags = RHF_Flags(basis_flags,scf_flags)

    #set up values to read in if not doing direct
    read_in::Dict{String,Any} = Dict([])

    merge!(read_in, input_info["Enuc"])
    merge!(read_in, input_info["Overlap"])
    merge!(read_in, input_info["Kinetic Energy"])
    merge!(read_in, input_info["Nuclear Attraction"])
    merge!(read_in, input_info["Two-Electron"])

    #read_in::Dict{String,Dict{String,Any}} = Dict([
    #                                                ("Enuc",input_info["Enuc"]),
    #                                                ("Ovr",input_info["Overlap"]),
    #                                                ("Kei",input_info["Kinetic Energy"]),
    #                                                ("Nai",input_info["Nuclear Attraction"]),
    #                                                ("Tei",input_info["Two-Electron"])
    #                                             ])
    #GC.enable(false)
    if (scf_flags.DIRECT == false)
        scf = rhf_energy(rhf_flags, read_in)
    end
    #GC.enable(true)
    #GC.gc()

    if (MPI.Comm_rank(comm) == 0)
        println("                       ========================================          ")
        println("                             END RESTRICTED CLOSED-SHELL                 ")
        println("                                     HARTREE-FOCK                        ")
        println("                       ========================================          ")
    end

    json_output = open("test.json","w")
        output_fock = Dict([("Structure","Fock"),("Data",scf.Fock)])
        output_density = Dict([("Structure","Density"),("Data",scf.Density)])
        output_coeff = Dict([("Structure","Coeff"),("Data",scf.Coeff)])
        if (MPI.Comm_rank(comm) == 0)
            write(json_output,JSON.json(output_fock))
            write(json_output,JSON.json(output_density))
            write(json_output,JSON.json(output_coeff))
        end

        if (typeof(scf) == RHFRestartData)
            output_hcore = Dict([("Structure","Hcore"),("Data",scf.H)])
            output_ortho = Dict([("Structure","Ortho"),("Data",scf.Ortho)])
            output_iter = Dict([("Structure","Iteration"),("Data",scf.iter)])
            if (MPI.Comm_rank(comm) == 0)
                write(json_output,JSON.json(output_hcore))
                write(json_output,JSON.json(output_ortho))
                write(json_output,JSON.json(output_iter))
            end
        end
    close(json_output)

    return scf
end
export run

#=
function run(flags::RHF_Flags, restart::RHFRestartData)
    comm=MPI.COMM_WORLD

    if (MPI.Comm_rank(comm) == 0)
        println("--------------------------------------------------------------------------------------")
        println("                       ========================================          ")
        println("                         RESTRICTED CLOSED-SHELL HARTREE-FOCK            ")
        println("                       ========================================          ")
        println("")
    end

    #GC.enable(false)
    scf = rhf_energy(flags,restart)
    #GC.enable(true)
    #GC.gc()

    if (MPI.Comm_rank(comm) == 0)
        println("                       ========================================          ")
        println("                             END RESTRICTED CLOSED-SHELL                 ")
        println("                                     HARTREE-FOCK                        ")
        println("                       ========================================          ")
    end

    output_fock = Dict([("Structure","Fock"),("Data",scf.Fock)])
    output_density = Dict([("Structure","Density"),("Data",scf.Density)])
    output_coeff = Dict([("Structure","Coeff"),("Data",scf.Coeff)])
    if (MPI.Comm_rank(comm) == 0)
        json_output = open("test.json","w")
            write(json_output,JSON.json(output_fock))
            write(json_output,JSON.json(output_density))
            write(json_output,JSON.json(output_coeff))
        close(json_output)
    end

    return scf
end
export run
=#

end
