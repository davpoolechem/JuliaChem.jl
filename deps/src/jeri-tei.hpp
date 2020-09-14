#ifndef JERI_TEI_H
#define JERI_TEI_H

#include <libint2.hpp>
#include <jlcxx/jlcxx.hpp>
//#include <jlcxx/stl.hpp>

//#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <limits>
#include <vector>

typedef int64_t julia_int;

//------------------------------------------------------------------------//
//-- C++ JERI engine: small wrapper allowing LibInt to be used in Julia --//
//------------------------------------------------------------------------//
class TEIEngine {
  const libint2::BasisSet* m_basis_set;
  const libint2::ShellPair* m_shellpair_data;
  
  libint2::Engine m_coulomb_eng;

public:
 //-- ctors and dtors --//
  TEIEngine(const libint2::BasisSet& t_basis_set, 
    const std::vector<libint2::ShellPair>& t_shellpair_data)
    : m_basis_set(&t_basis_set),
      m_shellpair_data(t_shellpair_data.data()),
      m_coulomb_eng(libint2::Operator::coulomb,
        m_basis_set->max_nprim(),
        m_basis_set->max_l(),
        0)
  {
    //-- no screening done in engine --// 
    m_coulomb_eng.set_precision(0.0); 
  }
    
  ~TEIEngine() { };

  //-- member functions --//
  void compute_eri_block(jlcxx::ArrayRef<double> eri_block, 
    julia_int ash, julia_int bsh, julia_int csh, julia_int dsh, 
    julia_int bra_idx, julia_int ket_idx,
    julia_int absize, julia_int cdsize) 
  {
    //if (ash == 40 && bsh == 26) {
    //  if (csh == 8 && dsh == 8) {
    //    std::cout << m_basis_set[ash-1] << std::endl;
    //    std::cout << m_basis_set[bsh-1] << std::endl;
    //    std::cout << m_basis_set[csh-1] << std::endl;
    //    std::cout << m_basis_set[dsh-1] << std::endl;
    //  }
   // }

    //assert(ash >= bsh);
    //assert(csh >= dsh);
    //assert(ab_idx >= cd_idx);

    //std::cout << ash-1 << "," << bsh-1 << ";" << ab_idx << std::endl;
    //std::cout << csh-1 << "," << dsh-1 << ";" << cd_idx << std::endl << std::endl;

    m_coulomb_eng.compute2<libint2::Operator::coulomb, 
      libint2::BraKet::xx_xx, 0>((*m_basis_set)[ash-1], (*m_basis_set)[bsh-1],
      (*m_basis_set)[csh-1], (*m_basis_set)[dsh-1],
      &m_shellpair_data[bra_idx-1], &m_shellpair_data[ket_idx-1]);
      
    //assert(m_coulomb_eng.results()[0] != nullptr); 
    if (m_coulomb_eng.results()[0] != nullptr) {
      for (int i = 0; i != absize*cdsize; ++i) {
        eri_block[i] = m_coulomb_eng.results()[0][i];
      }
    } else {
      for (int i = 0; i != absize*cdsize; ++i) {
        eri_block[i] = 0.0; 
      }
    }
  }
};

#endif /* JERI_TEI_H */
