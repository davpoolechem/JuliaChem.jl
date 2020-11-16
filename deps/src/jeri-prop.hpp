#ifndef JERI_PROP_H
#define JERI_PROP_H

#include <libint2.hpp>
#include <jlcxx/jlcxx.hpp>
//#include <jlcxx/stl.hpp>

#include <iostream>
#include <vector>

typedef int64_t julia_int;

//------------------------------------------------------------------------//
//-- C++ JERI engine: small wrapper allowing LibInt to be used in Julia --//
//------------------------------------------------------------------------//
class PropEngine {
  const libint2::BasisSet* m_basis_set;
  //const std::vector<libint2::Atom>* m_atoms;
  //std::vector<julia_int> m_shell2atom; 

  libint2::Engine m_dipole_eng;

public:
  //-- ctors and dtors --//
  PropEngine(const std::vector<libint2::Atom>& t_atoms, 
    const libint2::BasisSet& t_basis_set) 
    : m_basis_set(&t_basis_set),
      //m_atoms(&t_atoms), 
      //m_shell2atom(t_basis_set.shell2atom(t_atoms)),
      m_dipole_eng(libint2::Operator::emultipole1, 
        m_basis_set->max_nprim(),
        m_basis_set->max_l(),
        0)
  {
    //m_dipole_eng.set_params(libint2::make_point_charges(t_atoms));
  }

  ~PropEngine() { };

  //-- setters and getters --//

  //-- member functions --//
  void compute_dipole_block(jlcxx::ArrayRef<double> dipole_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    m_dipole_eng.compute((*m_basis_set)[ash-1], (*m_basis_set)[bsh-1]);
    for (int xyz = 0; xyz != 3; ++xyz) {
      //if (m_coulomb_eng.results()[0] != nullptr) {
      for (int idx = 0; idx != absize; ++idx) {
        dipole_block[absize*xyz + idx] = m_dipole_eng.results()[xyz+1][idx];
      }
    }
  }
};
 
#endif /* JERI_PROP_H */
