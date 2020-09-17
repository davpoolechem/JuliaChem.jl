#ifndef JERI_OEI_H
#define JERI_OEI_H

#include <libint2.hpp>
#include <jlcxx/jlcxx.hpp>
//#include <jlcxx/stl.hpp>

#include <iostream>
#include <vector>

typedef int64_t julia_int;

//------------------------------------------------------------------------//
//-- C++ JERI engine: small wrapper allowing LibInt to be used in Julia --//
//------------------------------------------------------------------------//
class OEIEngine {
  libint2::BasisSet m_basis_set;
  
  libint2::Engine m_overlap_eng;
  libint2::Engine m_kinetic_eng;
  libint2::Engine m_nuc_attr_eng;

public:
  //-- ctors and dtors --//
  OEIEngine(const std::vector<libint2::Atom>& t_atoms, 
    const libint2::BasisSet& t_basis_set) 
    : m_basis_set(t_basis_set), 
      m_overlap_eng(libint2::Operator::overlap, 
        m_basis_set.max_nprim(),
        m_basis_set.max_l(),
        0),
      m_kinetic_eng(libint2::Operator::kinetic, 
        m_basis_set.max_nprim(),
        m_basis_set.max_l(),
        0),
      m_nuc_attr_eng(libint2::Operator::nuclear, 
        m_basis_set.max_nprim(),
        m_basis_set.max_l(),
        0)
  { 
    m_nuc_attr_eng.set_params(libint2::make_point_charges(t_atoms)); 
  }

  /*
  OEIEngine(const std::vector<libint2::Atom>& t_atoms, 
    const std::vector<std::vector<libint2::Shell> >& t_shells) 
    : m_basis_set(t_atoms, t_shells, "", true)
  { 
    std::cout << "ANALYZE ELEMENT BASES" << std::endl;
    for (libint2::Atom atom : t_atoms) {
      auto Z = atom.atomic_number;
      std::cout << "ATOMIC NUMBER: " << Z << std::endl;
      std::cout << "SIZE: " << t_shells[Z].size() << std::endl;
      
      for (libint2::Shell shell : t_shells[Z]) {
        std::cout << shell << std::endl;
      }
    }
    std::cout << "ANALYZE ATOMS" << std::endl;
    for ( auto atom : t_atoms ) {
      std::cout << "{" << std::endl;
      std::cout << "  " << atom.atomic_number << std::endl;
      std::cout << "  " << atom.x << ", " << atom.y << ", " << atom.z << std::endl;
      std::cout << "}" << std::endl << std::endl;
      //std::cout << "SHELL CONTRACT: " << shell.ncontr() << std::endl;
    }
    std::cout << "ANALYZE BASIS SET" << std::endl;
    for ( auto shell : m_basis_set ) {
      std::cout << shell << std::endl;
      //std::cout << "SHELL CONTRACT: " << shell.ncontr() << std::endl;
    }
    
    //std::cout << m_basis_set.max_nprim() << std::endl;
    //std::cout << m_basis_set.max_l() << std::endl;

    m_overlap_eng = libint2::Engine(libint2::Operator::overlap, 
      m_basis_set.max_nprim(), m_basis_set.max_l(), 0);
    //m_overlap_eng.set_params(precision);

    m_kinetic_eng = libint2::Engine(libint2::Operator::kinetic, 
      m_basis_set.max_nprim(), m_basis_set.max_l(), 0);
    //m_kinetic_eng.set_params();

    m_nuc_attr_eng = libint2::Engine(libint2::Operator::nuclear, 
      m_basis_set.max_nprim(), m_basis_set.max_l(), 0);
    m_nuc_attr_eng.set_params(libint2::make_point_charges(t_atoms));
  }
  */

  ~OEIEngine() { };

  //-- setters and getters --//

  //-- member functions --//
  void compute_overlap_block(jlcxx::ArrayRef<double> S_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    //std::cout << ash << std::endl;
    //std::cout << m_basis_set[ash-1] << std::endl;
    //std::cout << bsh << std::endl;
    //std::cout << m_basis_set[bsh-1] << std::endl;
    m_overlap_eng.compute(m_basis_set[ash-1], m_basis_set[bsh-1]);
    for (int i = 0; i != absize; ++i) {
      S_block[i] = m_overlap_eng.results()[0][i];
    }
  }
 
  void compute_kinetic_block(jlcxx::ArrayRef<double> T_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    m_kinetic_eng.compute(m_basis_set[ash-1], m_basis_set[bsh-1]);
    for (int i = 0; i != absize; ++i) {
      T_block[i] = m_kinetic_eng.results()[0][i];
    }
  }

  void compute_nuc_attr_block(jlcxx::ArrayRef<double> V_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    m_nuc_attr_eng.compute(m_basis_set[ash-1], m_basis_set[bsh-1]);
    for (int i = 0; i != absize; ++i) {
      V_block[i] = m_nuc_attr_eng.results()[0][i];
    }
  }
};

#endif /* JERI_OEI_H */