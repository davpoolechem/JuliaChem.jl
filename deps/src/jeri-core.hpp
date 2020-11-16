#ifndef JERI_CORE_H
#define JERI_CORE_H

#include <libint2.hpp>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>

#include <vector>

typedef int64_t julia_int;

//-----------------------------------//
//-- initialize/finalize functions --//
//-----------------------------------//
void initialize() {
  libint2::initialize(); 
}

void finalize() {
  libint2::finalize();
}

//--------------------------------//
//-- Map libint2::Atom to Julia --// 
//--------------------------------//
template<> struct jlcxx::IsMirroredType<libint2::Atom> : std::false_type { };
template<> struct jlcxx::IsMirroredType<std::vector<libint2::Atom> > : 
  std::false_type { };

libint2::Atom create_atom(julia_int t_atomic_number, double coords[3]) {
  return libint2::Atom{ static_cast<int>(t_atomic_number), coords[0], coords[1], coords[2] };
}

//---------------------------------//
//-- Map libint2::Shell to Julia --// 
//---------------------------------//
template<> struct jlcxx::IsMirroredType<libint2::Shell> : std::false_type { };
template<> struct jlcxx::IsMirroredType<std::vector<libint2::Shell> > : 
  std::false_type { };

libint2::Shell create_shell(julia_int ang_mom, 
  const jlcxx::ArrayRef<double> t_exps,
  const jlcxx::ArrayRef<double> t_coeffs, double t_atom_center[3]) {

  libint2::svector<double> exps(t_exps.size());
  for (int iexp = 0; iexp != t_exps.size(); ++iexp)
    exps[iexp] = t_exps[iexp];

  libint2::svector<double> coeffs(t_coeffs.size());
  for (int icoeff = 0; icoeff != t_coeffs.size(); ++icoeff)
    coeffs[icoeff] = t_coeffs[icoeff];

  libint2::Shell::Contraction shell_contract{ static_cast<int>(ang_mom), 
    false, { coeffs } };
  libint2::svector<libint2::Shell::Contraction> shell_contract_vec(1);
  shell_contract_vec[0] = std::move(shell_contract);

  std::array<double, 3> atom_center = { t_atom_center[0], t_atom_center[1],
    t_atom_center[2] };

  return libint2::Shell(exps, shell_contract_vec, atom_center); 
}

//------------------------------------//
//-- Map libint2::BasisSet to Julia --// 
//------------------------------------//
template<> struct jlcxx::IsMirroredType<libint2::BasisSet> : std::false_type { };

libint2::BasisSet copy_basis(const libint2::BasisSet& basis_set) {
  libint2::BasisSet copy = basis_set;
  return copy;
}

//-------------------------------------//
//-- Map libint2::ShellPair to Julia --// 
//-------------------------------------//
template<> struct jlcxx::IsMirroredType<libint2::ShellPair> : 
  std::false_type { };
template<> struct jlcxx::IsMirroredType<std::vector<libint2::ShellPair> > : 
  std::false_type { };

void precompute_shell_pair_data(std::vector<libint2::ShellPair>& shpdata,
  const libint2::BasisSet& basis_set) {

  int nshells = basis_set.size();

  for (int ash = 0; ash != nshells; ++ash) {
    for (int bsh = 0; bsh <= ash; ++bsh) {
      shpdata.emplace_back(libint2::ShellPair(
        basis_set[ash], basis_set[bsh],
        std::log(std::numeric_limits<double>::epsilon()/1e10)));
    }
  }
}

#endif /* JERI_CORE_H */
