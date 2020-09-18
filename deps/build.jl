using Eigen_jll
using boost_jll
using libcxxwrap_julia_jll

function build_jeri()
  eigen_root = Eigen_jll.LIBPATH
  boost_root = boost_jll.LIBPATH
  jlcxx_root = libcxxwrap_julia_jll.LIBPATH
  libint_root = ENV["LIBINT_ROOT"] 

  if !ispath("build") mkdir("build") end
  cd("build")
  run(`cmake -DEIGEN_INPUT=$eigen_root -DBOOST_INPUT=$boost_root 
    -DJLCXX_INPUT=$jlcxx_root -DLIBINT_INPUT=$libint_root ../`)
  run(`make`)
  run(`make install`)
end

build_jeri()
