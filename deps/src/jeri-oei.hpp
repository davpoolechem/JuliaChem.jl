#ifndef JERI_OEI_H
#define JERI_OEI_H

#include <libint2.hpp>
#include <jlcxx/jlcxx.hpp>
//#include <jlcxx/stl.hpp>

#include <iostream>
#include <vector>

typedef int64_t julia_int;

//----------------------------------------------------------------//
//-- small wrapper class to enable Libint one-electron integral --// 
//--                  computation in Julia                      --//
//----------------------------------------------------------------//
class OEIEngine {
  const libint2::BasisSet* m_basis_set;
  const std::vector<libint2::Atom>* m_atoms;
  std::vector<julia_int> m_shell2atom; 

  libint2::Engine m_overlap_eng;
  libint2::Engine m_kinetic_eng;
  libint2::Engine m_nuc_attr_eng;
  
  julia_int m_deriv_order;

public:
  //-- ctors and dtors --//
  OEIEngine(const std::vector<libint2::Atom>& t_atoms, 
    const libint2::BasisSet& t_basis_set, julia_int
    t_deriv_order) 
    : m_basis_set(&t_basis_set),
      m_atoms(&t_atoms), 
      m_shell2atom(t_basis_set.shell2atom(t_atoms)),
      m_deriv_order(t_deriv_order),
      m_overlap_eng(libint2::Operator::overlap, 
        m_basis_set->max_nprim(),
        m_basis_set->max_l(),
        t_deriv_order),
      m_kinetic_eng(libint2::Operator::kinetic, 
        m_basis_set->max_nprim(),
        m_basis_set->max_l(),
        t_deriv_order),
      m_nuc_attr_eng(libint2::Operator::nuclear, 
        m_basis_set->max_nprim(),
        m_basis_set->max_l(),
        t_deriv_order)
  { 
    m_nuc_attr_eng.set_params(libint2::make_point_charges(t_atoms)); 
  }

  ~OEIEngine() { };

  //-- setters and getters --//

  //-- member functions --//
  void compute_overlap_block(jlcxx::ArrayRef<double> S_block, julia_int ash, 
    julia_int bsh, julia_int absize)
  {
    assert(m_deriv_order == 0);

    m_overlap_eng.compute((*m_basis_set)[ash-1], (*m_basis_set)[bsh-1]);
    for (int i = 0; i != absize; ++i) {
      S_block[i] = m_overlap_eng.results()[0][i];
    }
  }
  
  void compute_overlap_grad_block(jlcxx::ArrayRef<double> S_grad_block, 
    julia_int ash, julia_int bsh, julia_int absize)
  {
    constexpr auto nopers = libint2::operator_traits<libint2::Operator::overlap>::nopers; 
    const auto nresults =
      nopers * libint2::num_geometrical_derivatives(m_atoms->size(), m_deriv_order);        
    const auto nderivcenters_shset =
      2 + ((libint2::Operator::overlap == libint2::Operator::nuclear) ? m_atoms->size() : 0);    

    auto atom1 = m_shell2atom[ash-1]; 
    auto atom2 = m_shell2atom[bsh-1]; 

    std::cout << "Shells: " << ash-1 << ", " << bsh-1 << std::endl;
    std::cout << "Atoms: " << atom1 << ", " << atom2 << std::endl; 
    std::cout << "nopers: " << nopers << ", " << nresults << std::endl;
    std::cout << "shset: " << nderivcenters_shset << std::endl; 

    m_overlap_eng.compute((*m_basis_set)[ash-1], (*m_basis_set)[bsh-1]);
    
    std::size_t shellset_idx = 0; 
    for (auto c = 0; c != nderivcenters_shset; ++c) {
      auto atom = (c == 0) ? atom1 : ((c == 1) ? atom2 : c - 2);
      auto op_start = 3 * atom * nopers;
      auto op_fence = op_start + nopers;
      for (auto xyz = 0; xyz != 3;
        ++xyz, op_start += nopers, op_fence += nopers) {
        for (unsigned int op = op_start; op != op_fence;
          ++op, ++shellset_idx) {
      
          std::cout << c << ", " << xyz << " => " 
            << op << ", " << shellset_idx << std::endl;
            
          m_overlap_eng.compute((*m_basis_set)[ash-1], (*m_basis_set)[bsh-1]);
          for (int idx = 0; idx != absize; ++idx) {
            S_grad_block[absize*op + idx] += m_overlap_eng.results().at(shellset_idx)[idx];
          }   
        }
      }
    }
  }

  void compute_kinetic_block(jlcxx::ArrayRef<double> T_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    m_kinetic_eng.compute((*m_basis_set)[ash-1], (*m_basis_set)[bsh-1]);
    for (int i = 0; i != absize; ++i) {
      T_block[i] = m_kinetic_eng.results()[0][i];
    }
  }

  void compute_nuc_attr_block(jlcxx::ArrayRef<double> V_block, julia_int ash, 
    julia_int bsh, julia_int absize) 
  {
    m_nuc_attr_eng.compute((*m_basis_set)[ash-1], (*m_basis_set)[bsh-1]);
    for (int i = 0; i != absize; ++i) {
      V_block[i] = m_nuc_attr_eng.results()[0][i];
    }
  }
};

#endif /* JERI_OEI_H */
