#/bin/bash

source ~/.bash_profile

basis_sets=("6-311++G_2d_2p" "6-311++G_d_p-J")

for dir in ${basis_sets[@]}; do
  for file in "$dir"/*; do
    name=$(echo "$file" | cut -f 1 -d '.')
    echo $name
    compute-julia-pbs-benchmark $name skylake_8180 1 1 36
  done
done
