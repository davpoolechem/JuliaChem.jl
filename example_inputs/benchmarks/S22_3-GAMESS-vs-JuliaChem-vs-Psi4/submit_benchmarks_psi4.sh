#/bin/bash

source ~/.bash_profile

#basis_sets=("6-311++G_2d_2p" "6-311++G_d_p-J")
basis_sets=("6-311++G_2d_2p_MKL")

function compute-psi4-pbs() { 
  echo "#!/bin/bash" > $1.sh                                                    
  echo "#" >> $1.sh                                                             
  echo "source /etc/profile.d/modules.sh" >> $1.sh                              
  echo "source ~/.bashrc" >> $1.sh                                              
  echo "#" >> $1.sh                                                             
  echo "export OMP_NUM_THREADS=$5" >> $1.sh                                   
  echo "export MKL_NUM_THREADS=1" >> $1.sh                                   
  echo "#" >> $1.sh                                                             
  echo "julia --check-bounds=no --math-mode=fast --optimize=3 --inline=yes --compiled-modules=yes $1.jl" >> $1.sh
                                                                                
  chmod 755 $1.sh                                                               
  qsub -O $1 -A CHM135 -q $2 -t 360 -n $3 $1.sh                                 
}     

for dir in ${basis_sets[@]}; do
  for file in "$dir"/*; do
    name=$(echo "$file" | cut -f 1 -d '.')
    echo $name
    compute-psi4-pbs $name skylake_8180 1 1 112
  done
done
    
#compute-psi4-pbs 6-311++G_2d_2p/ammonia_trimer skylake_8180 1 1 112
