#/bin/bash

source ~/.bash_profile

input_dir="example_inputs/S22_3/"
#basis_sets=("${input_dir}6-311++G_2d_2p" "${input_dir}6-311++G_d_p-J")
basis_sets=("${input_dir}6-311++G_2d_2p")

for dir in ${basis_sets[@]}; do
  for file in "$dir"/*.json; do
    name=$(echo "$file" | cut -f 1 -d '.')
    echo $file
    compute-juliachem-benchmark $name haswell 3 1 36
  done
done
