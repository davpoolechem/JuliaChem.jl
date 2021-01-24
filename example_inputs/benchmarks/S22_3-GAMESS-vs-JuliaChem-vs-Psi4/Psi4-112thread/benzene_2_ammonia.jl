#== import modules here ==#
using PyCall
psi4 = pyimport("psi4")

#== build molecule ==#
psi4.set_memory("500 MB")

input_file = "/gpfs/jlse-fs0/users/davpoolechem/projects/Julia/Julia-Sandbox/QC-Programs/Psi4/S22_3/benzene_2_ammonia.xyz"

io = open(input_file, "r")
input_file_string = read(io, String)
close(io)

mol = psi4.geometry(input_file_string)
mol.reset_point_group("c1")

options = Dict( 
  "REFERENCE" => "RHF", 
  "BASIS" => "6-311++G(2d,2p)", 
  "SCF_TYPE" => "DIRECT", 
  "GUESS" => "CORE", 
  "DF_SCF_GUESS" => false, 
  "PUREAM" => false, 
  "PROPERTIES" => [],
  "PRINT" => 0,
  "E_CONVERGENCE" => 1.0,
  "D_CONVERGENCE" => 2e-5
)

psi4.set_options(options)
psi4.set_num_threads(112)

psi4.energy("scf/6-311++G(2d,2p)")
