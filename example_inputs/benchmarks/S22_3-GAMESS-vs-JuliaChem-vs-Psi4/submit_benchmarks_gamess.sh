#/bin/bash

source ~/.bash_profile

input_dir="S22_3/"
#basis_sets=("${input_dir}6-311++G_2d_2p" "${input_dir}6-311++G_d_p-J")
basis_sets=("${input_dir}6-311++G_2d_2p")

for dir in ${basis_sets[@]}; do
  for file in "$dir"/*.inp; do
    name=$(echo "$file" | cut -f 1 -d '.')
    echo $name
    compute-pbs $name skylake_8180 1 1 112 
  done
done
    
#compute-pbs S22_3/6-311++G_2d_2p/ammonia_trimer skylake_8180 1 1 112 
