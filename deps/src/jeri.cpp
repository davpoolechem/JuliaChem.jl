#include "jeri-core.hpp"
#include "jeri-oei.hpp"

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
  mod.add_type<libint2::BasisSet>("BasisSet");

  //-- oei engine information --//
  //mod.add_type<libint2::Engine>("LibIntEngine");
  mod.add_type<OEIEngine>("OEIEngine")
    .constructor<const std::vector<libint2::Atom>&, 
      const std::vector<std::vector<libint2::Shell> >& >()
    .method("basis", &OEIEngine::basis)
    .method("compute_overlap_block", &OEIEngine::compute_overlap_block)
    .method("compute_kinetic_block", &OEIEngine::compute_kinetic_block)
    .method("compute_nuc_attr_block", &OEIEngine::compute_nuc_attr_block);
} 
