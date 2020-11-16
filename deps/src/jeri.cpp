#include "jeri-core.hpp"
#include "jeri-oei.hpp"
#include "jeri-prop.hpp"
#include "jeri-tei.hpp"

JLCXX_MODULE define_jeri(jlcxx::Module& mod) {
  //-- initialize/finalize functions --//
  mod.method("initialize", &initialize);
  mod.method("finalize", &finalize);

  //-- atom information --//
  mod.add_type<libint2::Atom>("Atom")
    .method("create_atom", &create_atom);
  jlcxx::stl::apply_stl<libint2::Atom>(mod);

  //-- shell information --//
  mod.add_type<libint2::Shell>("Shell")
    .method("create_shell", &create_shell);
  jlcxx::stl::apply_stl<libint2::Shell>(mod);

  //-- basis set information --//
  mod.add_type<libint2::BasisSet>("BasisSet")
    .constructor<const std::vector<libint2::Atom>&, 
      const std::vector<std::vector<libint2::Shell> >& >()
    .constructor<const libint2::BasisSet&>();
  mod.method("copy_basis", &copy_basis);
  
  //-- shell pair information --//
  mod.add_type<libint2::ShellPair>("ShellPair");
  jlcxx::stl::apply_stl<libint2::ShellPair>(mod);
  mod.method("precompute_shell_pair_data", &precompute_shell_pair_data);
 
  //-- oei engine information --//
  //mod.add_type<libint2::Engine>("LibIntEngine");
  mod.add_type<OEIEngine>("OEIEngine")
    //.constructor<const std::vector<libint2::Atom>&, 
    //  const std::vector<std::vector<libint2::Shell> >& >()
    .constructor<const std::vector<libint2::Atom>&,
      const libint2::BasisSet&, julia_int>() 
    //.method("basis", &OEIEngine::basis)
    .method("compute_overlap_block", &OEIEngine::compute_overlap_block)
    .method("compute_overlap_grad_block", &OEIEngine::compute_overlap_grad_block)
    .method("compute_kinetic_block", &OEIEngine::compute_kinetic_block)
    .method("compute_kinetic_grad_block", &OEIEngine::compute_kinetic_grad_block)
    .method("compute_nuc_attr_block", &OEIEngine::compute_nuc_attr_block)
    .method("compute_nuc_attr_grad_block", &OEIEngine::compute_nuc_attr_grad_block);

  //-- prop engine information --//
  mod.add_type<PropEngine>("PropEngine")
    .constructor<const std::vector<libint2::Atom>&,
      const libint2::BasisSet&>() 
    .method("compute_dipole_block", &PropEngine::compute_dipole_block);

  //-- tei engine information --//
  mod.add_type<TEIEngine>("TEIEngine")
    .constructor<const libint2::BasisSet&, 
      const std::vector<libint2::ShellPair>& >() 
    .method("compute_eri_block", &TEIEngine::compute_eri_block);
} 

