using Eigen_jll
using boost_jll
using libcxxwrap_julia_jll

function build_jeri()
  eigen_root = Eigen_jll.artifact_dir
  boost_root = boost_jll.artifact_dir
  jlcxx_root = libcxxwrap_julia_jll.artifact_dir
  libint_root = ENV["LIBINT_ROOT"] 

  if !ispath("build") mkdir("build") end
  cd("build")
  run(`cmake -DEIGEN_INPUT=$eigen_root -DBOOST_INPUT=$boost_root 
    -DJLCXX_INPUT=$jlcxx_root -DLIBINT_INPUT=$libint_root ../`)
  run(`make`)
  run(`make install`)
end

build_jeri()
