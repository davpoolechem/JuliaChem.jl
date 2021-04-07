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
exec(open(dir_path+"/rhf-properties-xyz-repl.py").read())

JuliaChem.initialize()
energy, properties = rhf_properties_xyz(sys.argv[1])
JuliaChem.finalize()
