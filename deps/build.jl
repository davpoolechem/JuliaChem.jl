using Eigen_jll
using boost_jll
using libcxxwrap_julia_jll
using libint_jll

function build_jeri()
  eigen_root = Eigen_jll.artifact_dir
  boost_root = boost_jll.artifact_dir
  jlcxx_root = libcxxwrap_julia_jll.artifact_dir
  libint_root = "/export/home/david/programs/install/libint" 

  if ispath("build") rm("build", recursive=true) end
  mkdir("build")
  cd("build")
  run(`cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++  
    -DEIGEN_PATH=$eigen_root -DBOOST_PATH=$boost_root 
    -DJLCXX_PATH=$jlcxx_root -DLIBINT_PATH=$libint_root ../`)
  run(`make`)
  run(`make install`)
end

build_jeri()
