#/bin/bash

source ~/.bash_profile

input_dir="example_inputs/S22_3/"
#basis_sets=("${input_dir}6-311++G_2d_2p" "${input_dir}6-311++G_d_p-J")
basis_sets=("${input_dir}6-311++G_2d_2p")

for dir in ${basis_sets[@]}; do
  for file in "$dir"/*.json; do
    name=$(echo "$file" | cut -f 1 -d '.')
    echo $file
    compute-julia-pbs-benchmark $name skylake_8180 1 1 112
  done
done

#compute-julia-pbs-benchmark example_inputs/S22_3/6-311++G_2d_2p/uracil_trimer skylake_8180 1 1 112
