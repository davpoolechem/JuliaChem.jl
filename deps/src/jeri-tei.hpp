#ifndef JERI_TEI_H
#define JERI_TEI_H

#include <libint2.hpp>
#include <jlcxx/jlcxx.hpp>
//#include <jlcxx/stl.hpp>

#include <iostream>
#include <vector>

typedef int64_t julia_int;

//------------------------------------------------------------------------//
//-- C++ JERI engine: small wrapper allowing LibInt to be used in Julia --//
//------------------------------------------------------------------------//
class TEIEngine {
  libint2::BasisSet m_basis_set;
  
  libint2::Engine m_coulomb_eng;

public:
 //-- ctors and dtors --//
  TEIEngine(const libint2::BasisSet& t_basis_set)
    : m_basis_set(t_basis_set),
      m_coulomb_eng(libint2::Operator::coulomb,
        m_basis_set.max_nprim(),
        m_basis_set.max_l(),
        0)
  { }

  ~TEIEngine() { };

  //-- member functions --//
  void compute_eri_block(jlcxx::ArrayRef<double> eri_block, 
    julia_int ash, julia_int bsh, julia_int csh, julia_int dsh, 
    julia_int absize, julia_int cdsize) 
  {
    //if (ash == 160 && bsh == 160) {
    //  if (csh == 159 && dsh == 153) {
    //    std::cout << m_basis_set[ash-1] << std::endl;
    //    std::cout << m_basis_set[bsh-1] << std::endl;
    //    std::cout << m_basis_set[csh-1] << std::endl;
    //    std::cout << m_basis_set[dsh-1] << std::endl;
    //  }
    //}
 
    m_coulomb_eng.compute2<libint2::Operator::coulomb, 
      libint2::BraKet::xx_xx, 0>(m_basis_set[ash-1], m_basis_set[bsh-1],
      m_basis_set[csh-1], m_basis_set[dsh-1]);
    if (m_coulomb_eng.results()[0] != nullptr) {
      for (int i = 0; i != absize*cdsize; ++i) {
        eri_block[i] = m_coulomb_eng.results()[0][i];
      }
    }
  }
};

#endif /* JERI_TEI_H */
