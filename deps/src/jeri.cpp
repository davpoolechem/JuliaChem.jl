#include <libint2.hpp>
#include <jlcxx/jlcxx.hpp>
#include <jlcxx/stl.hpp>

#include <unistd.h>
#include <iostream>
#include <vector>

typedef int64_t julia_int;

//--------------------------------//
//-- Map libint2::Atom to Julia --// 
//--------------------------------//
template<> struct jlcxx::IsMirroredType<libint2::Atom> : std::false_type { };
template<> struct jlcxx::IsMirroredType<std::vector<libint2::Atom> > : 
  std::false_type { };

libint2::Atom create_atom(julia_int t_atomic_number, double coords[3]) {
  return libint2::Atom{ t_atomic_number, coords[0], coords[1], coords[2] };
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

  libint2::Shell::Contraction shell_contract{ ang_mom, false, { coeffs } };
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

//------------------------------------------------------------------------//
//-- C++ JERI engine: small wrapper allowing LibInt to be used in Julia --//
//------------------------------------------------------------------------//
class OEIEngine {
  libint2::BasisSet m_basis_set;
  
  libint2::Engine m_kinetic_eng;
  libint2::Engine m_overlap_eng;
  libint2::Engine m_nuc_attr_eng;

public:
  //-- ctors and dtors --//
  OEIEngine() { initialize(); };
  OEIEngine(const std::vector<libint2::Atom> t_atoms, 
    const std::vector<std::vector<libint2::Shell> > t_shells) 
  { 
    initialize();

  /*  
    std::cout << "ANALYZE ELEMENT BASES" << std::endl;
    for (libint2::Atom atom : t_atoms) {
      auto Z = atom.atomic_number;
      std::cout << "ATOMIC NUMBER: " << Z << std::endl;
      std::cout << "SIZE: " << t_shells[Z].size() << std::endl;
      
      for (libint2::Shell shell : t_shells[Z]) {
        std::cout << shell << std::endl;
      }
    }
    */
    /*
    std::cout << "ANALYZE ATOMS" << std::endl;
    for ( auto atom : t_atoms ) {
      std::cout << "{" << std::endl;
      std::cout << "  " << atom.atomic_number << std::endl;
      std::cout << "  " << atom.x << ", " << atom.y << ", " << atom.z << std::endl;
      std::cout << "}" << std::endl << std::endl;
      //std::cout << "SHELL CONTRACT: " << shell.ncontr() << std::endl;
    }
    */
    m_basis_set = libint2::BasisSet(t_atoms, t_shells, "", true);
    /*
    std::cout << "ANALYZE BASIS SET" << std::endl;
    for ( auto shell : m_basis_set ) {
      std::cout << shell << std::endl;
      //std::cout << "SHELL CONTRACT: " << shell.ncontr() << std::endl;
    }
    */
    
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

  ~OEIEngine() { finalize(); };

  //-- setters and getters --//
  libint2::BasisSet basis() { return m_basis_set; }

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

private:
  //-- private member functions --//
  void initialize() {
    libint2::initialize(); 
  }

  void finalize() {
    libint2::finalize();
  }
};

JLCXX_MODULE define_jeri(jlcxx::Module& mod) {
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

  //-- engine information --//
  //mod.add_type<libint2::Engine>("LibIntEngine");
  mod.add_type<OEIEngine>("OEIEngine")
    .constructor<const std::vector<libint2::Atom>&, 
      const std::vector<std::vector<libint2::Shell> >& >()
    .method("basis", &OEIEngine::basis)
    .method("compute_overlap_block", &OEIEngine::compute_overlap_block)
    .method("compute_kinetic_block", &OEIEngine::compute_kinetic_block)
    .method("compute_nuc_attr_block", &OEIEngine::compute_nuc_attr_block);
} 
