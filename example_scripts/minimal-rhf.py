#=============================#
#== put needed modules here ==#
#=============================#
import sys
import os

import julia
from julia import JuliaChem
from julia import Base 

#================================#
#== JuliaChem execution script ==#
#================================#
dir_path = os.path.dirname(os.path.realpath(__file__))
exec(open(dir_path+"/minimal-rhf-repl.py").read())

JuliaChem.initialize()
scf = minimal_rhf(sys.argv[1])
JuliaChem.finalize()
