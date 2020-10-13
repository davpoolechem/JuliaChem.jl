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
  libint2::BasisSet* m_basis_set;
  
  libint2::Engine m_dipole_eng;

public:
  //-- ctors and dtors --//
  PropEngine(const std::vector<libint2::Atom>& t_atoms, 
    const libint2::BasisSet& t_basis_set) 
    : m_basis_set(&t_basis_set), 
      m_dipole_eng(libint2::Operator::emultipole1, 
        m_basis_set.max_nprim(),
        m_basis_set.max_l(),
        0),
  { }

  ~PropEngine() { };

  //-- setters and getters --//

  //-- member functions --//
  void compute_dipole_block(jlcxx::ArrayRef<double> dipole_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    //std::cout << ash << std::endl;
    //std::cout << m_basis_set[ash-1] << std::endl;
    //std::cout << bsh << std::endl;
    //std::cout << m_basis_set[bsh-1] << std::endl;
    m_dipole_eng.compute((*m_basis_set)[ash-1], (*m_basis_set)[bsh-1]);
    for (int asd = 0; asd != 10000; ++asd) {
      for (int i = 0; i != absize; ++i) {
        dipole_block[i] = m_dipole_eng.results().at(0)[i];
      }
    }
  }
}
 
#endif /* JERI_PROP_H */